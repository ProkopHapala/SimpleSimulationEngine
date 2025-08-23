import os
import sys
import argparse
import numpy as np

# Ensure repo 'python' root is on sys.path to import OpenCL utilities
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PY_ROOT  = os.path.dirname(_THIS_DIR)               # .../python
if _PY_ROOT not in sys.path:
    sys.path.insert(0, _PY_ROOT)

from pyMolecular.OCL.OpenCLBase import OpenCLBase
from Radiosity3D import rasterize_face_to_elems, unit_cube


def three_quads():
    """
    Three parallel unit quads along X: at x=0 (face 0), x=0.5 (face 1, occluder), x=1 (face 2).
    Returns (V,F) with separate vertices per quad for simplicity.
    """
    V = np.array([
        # x = 0
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 1.0, 1.0],
        [0.0, 0.0, 1.0],
        # x = 0.5 (occluder)
        [0.5, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        # x = 1.0
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [1.0, 1.0, 1.0],
        [1.0, 0.0, 1.0],
    ], dtype=np.float64)
    F = [
        [0, 1, 2, 3],    # face 0
        [4, 5, 6, 7],    # face 1 (occluder)
        [8, 9,10,11],    # face 2
    ]
    return V, F


class TriangleOcclusionOCL(OpenCLBase):
    """
    GPU occlusion between element centers using triangle obstacles.
    One work-item per element i; each tests rays to all j against all triangles.
    """
    def __init__(self, nloc=64, device_index=0):
        super().__init__(nloc=nloc, device_index=device_index)
        self.kernel_name = 'occlusion_matrix'
        # Build program from shared kernel library to keep all kernels in one place
        kpath = os.path.abspath(os.path.join(_PY_ROOT, '..', 'cpp', 'common_resources', 'cl', 'radiosity.cl'))
        ok = self.load_program(kernel_path=kpath, bPrint=False, bMakeHeaders=True)
        assert ok and (self.prg is not None), f"Failed to build kernel from {kpath}"

    @staticmethod
    def triangulate_faces(vertices: np.ndarray, faces: list) -> (np.ndarray, np.ndarray):
        """
        Fan-triangulate polygon faces.
        Returns (tris_xyz, tri_face_ids)
        - tris_xyz: (ntris,3,3) float64
        - tri_face_ids: (ntris,) int32
        """
        tris = []
        tids = []
        for fi, f in enumerate(faces):
            if len(f) < 3: continue
            v0 = vertices[f[0]]
            for k in range(1, len(f)-1):
                v1 = vertices[f[k]]
                v2 = vertices[f[k+1]]
                tris.append([v0, v1, v2])
                tids.append(fi)
        if len(tris) == 0:
            return np.zeros((0,3,3), np.float64), np.zeros((0,), np.int32)
        return np.asarray(tris, dtype=np.float64), np.asarray(tids, dtype=np.int32)

    def build_buffers(self, points_xyz: np.ndarray, point_face_ids: np.ndarray,
                      tris_xyz: np.ndarray, tri_face_ids: np.ndarray):
        """
        Prepare device buffers matching kernel header names.
        points: (n,4) float32 xyz + w=face_id
        tris  : (ntris*3,4) float32 sequence of triangle vertices (A,B,C) with w=face_id
        occ   : (n*n,) float32
        """
        assert points_xyz.ndim == 2 and points_xyz.shape[1] == 3
        assert tris_xyz.ndim == 3 and tris_xyz.shape[1:] == (3,3)
        n = points_xyz.shape[0]
        ntris = tris_xyz.shape[0]

        pts = np.zeros((n,4), dtype=np.float32)
        pts[:,:3] = points_xyz.astype(np.float32, copy=False)
        pts[:,  3] = point_face_ids.astype(np.float32, copy=False)

        tri_seq = np.zeros((ntris*3,4), dtype=np.float32)
        tri_seq[0::3, :3] = tris_xyz[:,0,:].astype(np.float32, copy=False)
        tri_seq[1::3, :3] = tris_xyz[:,1,:].astype(np.float32, copy=False)
        tri_seq[2::3, :3] = tris_xyz[:,2,:].astype(np.float32, copy=False)
        tri_seq[:,3]      = tri_face_ids.repeat(3).astype(np.float32, copy=False)

        occ = np.zeros((n*n,), dtype=np.float32)

        # Allocate device buffers with names EXACTLY matching kernel header
        sizes = {
            'points': pts.nbytes,
            'tris':   tri_seq.nbytes,
            'occ':    occ.nbytes,
        }
        self.try_make_buffers(sizes, suffix="")
        self.toGPU('points', pts)
        self.toGPU('tris',   tri_seq)
        self.toGPU('occ',    occ)  # init to zeros

        # Kernel scalar params
        self.kernel_params = {
            'npoints': np.int32(n),
            'ntris':   np.int32(ntris),
        }
        return n

    def run(self, npoints: int):
        args = self.generate_kernel_args(self.kernel_name, bPrint=False)
        g = (self.roundUpGlobalSize(npoints),)
        self.prg.occlusion_matrix(self.queue, g, (self.nloc,), *args)
        # Download occ
        occ = np.zeros((npoints*npoints,), dtype=np.float32)
        self.fromGPU('occ', occ)
        return occ.reshape((npoints, npoints))


def build_elements_from_faces(V: np.ndarray, F: list, hex_size: float, area_keep_ratio: float, verbosity: int=0):
    # Rasterize faces into elements using existing helper
    elems = []
    for fi, face in enumerate(F):
        pts = V[np.array(face, dtype=int)]
        efs = rasterize_face_to_elems(pts, face_id=fi, hex_size=hex_size, area_keep_ratio=area_keep_ratio, verbosity=verbosity)
        elems.extend(efs)
    assert len(elems) > 0, 'No elements produced; check hex size vs geometry scale'
    centers = np.array([e['center'] for e in elems], dtype=np.float64)
    face_ids = np.array([e['face_id'] for e in elems], dtype=np.int32)
    elem_polys3d = [e.get('polys3d', None) for e in elems]
    return centers, face_ids, elem_polys3d


def visualize_occlusion_lines(vertices: np.ndarray, faces: list, centers: np.ndarray, occ: np.ndarray,
                              elem_polys3d=None, mesh_alpha: float = 0.2, elem_alpha: float = 0.9,
                              title: str = 'Occlusion (visible pairs only)'):
    try:
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"[warn] 3D plotting unavailable ({e}); skipping visualization")
        return

    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection='3d')

    # Draw faces (semi-transparent)
    polys = [[vertices[idx] for idx in f] for f in faces]
    mesh = Poly3DCollection(polys, alpha=float(mesh_alpha), facecolor='lightgray', edgecolor='k', linewidths=0.5)
    ax.add_collection3d(mesh)

    # Draw centers
    ax.scatter(centers[:,0], centers[:,1], centers[:,2], c='tab:blue', s=18, depthshade=True)

    # Optionally draw element tiles
    if (elem_polys3d is not None):
        poly_list = []
        for polyset in elem_polys3d:
            if polyset is None: continue
            for poly in polyset:
                if poly is None or len(poly) < 3: continue
                poly_list.append(poly)
        if len(poly_list):
            emesh = Poly3DCollection(poly_list, edgecolor='k', linewidths=0.3)
            emesh.set_facecolor((1.0, 0.6, 0.0, float(elem_alpha)))  # orange with alpha
            ax.add_collection3d(emesh)

    # Draw non-occluded lines (require both directions visible)
    n = centers.shape[0]
    for i in range(n):
        Pi = centers[i]
        for j in range(i+1, n):
            if (occ[i,j] == 0.0) and (occ[j,i] == 0.0):
                Pj = centers[j]
                ax.plot([Pi[0], Pj[0]], [Pi[1], Pj[1]], [Pi[2], Pj[2]], color='tab:green', linewidth=0.6)

    # axes limits and equal aspect
    all_pts = vertices if vertices is not None and len(vertices)>0 else centers
    mins = all_pts.min(axis=0)
    maxs = all_pts.max(axis=0)
    span = float(np.max(maxs - mins))
    mid  = 0.5*(mins+maxs)
    ax.set_xlim(mid[0]-span/2, mid[0]+span/2)
    ax.set_ylim(mid[1]-span/2, mid[1]+span/2)
    ax.set_zlim(mid[2]-span/2, mid[2]+span/2)
    if hasattr(ax, 'set_box_aspect'):
        ax.set_box_aspect((1,1,1))

    ax.set_title(title)
    ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
    plt.tight_layout()
    plt.show()


def main():
    ap = argparse.ArgumentParser(description='GPU occlusion between element COGs using triangle obstacles (PyOpenCL)')
    ap.add_argument('--hex', type=float, default=0.35, help='Hexagon size (grid step in face UV units)')
    ap.add_argument('--size', type=float, default=0.3, help='Alias for --hex (hexagon size)')
    ap.add_argument('--areaFactor', type=float, default=0.5, help='Raster keep ratio for partial hex (A_clip/A_full)')
    ap.add_argument('--verbosity', type=int, default=1, help='Verbosity level')
    ap.add_argument('--nloc', type=int, default=32, help='Preferred local size (work-group size)')
    ap.add_argument('--geom', type=str, default='three_quads', choices=['three_quads','cube2faces'], help='Which test geometry to use')
    ap.add_argument('--noPlot', action='store_true', help='Disable 3D plotting')
    ap.add_argument('--elements', action='store_true', help='Draw element tiling polygons')
    ap.add_argument('--meshAlpha', type=float, default=0.1, help='Opacity (alpha) for base faces mesh [0..1]')
    ap.add_argument('--elemAlpha', type=float, default=0.5, help='Opacity (alpha) for element tiles [0..1]')
    args = ap.parse_args()

    hex_size = args.size if (args.size is not None) else args.hex

    # Geometry selection
    if args.geom == 'three_quads':
        V, F = three_quads()
    else:
        # Fallback: reuse the two faces from unit_cube() defined in Radiosity3D
        V, F = unit_cube()

    # Rasterize -> element centers (+ polygons for optional plotting)
    centers, face_ids, elem_polys3d = build_elements_from_faces(V, F, hex_size, args.areaFactor, verbosity=args.verbosity)
    n = centers.shape[0]
    if args.verbosity:
        print(f'[info] elements: {n}')

    # Triangulate base faces for obstacles
    tris_xyz, tri_face_ids = TriangleOcclusionOCL.triangulate_faces(V, F)
    ntris = tris_xyz.shape[0]
    if args.verbosity:
        print(f'[info] obstacle triangles: {ntris}')

    # Run GPU occlusion
    ocl = TriangleOcclusionOCL(nloc=args.nloc)
    npoints = ocl.build_buffers(centers, face_ids, tris_xyz, tri_face_ids)
    occ = ocl.run(npoints)

    # Basic stats
    n_blocked = int(occ.sum())
    frac_blocked = n_blocked / float(n*n) if n>0 else 0.0
    print(f'[result] blocked pairs (including both directions): {n_blocked} / {n*n}  ({frac_blocked*100:.3f}%)')
    # Ensure diagonal is zero
    diag_sum = float(np.trace(occ))
    print(f'[sanity] diag sum = {diag_sum:.1f} (expected 0)')

    if not args.noPlot:
        vis_polys = (elem_polys3d if args.elements else None)
        visualize_occlusion_lines(V, F, centers, occ, elem_polys3d=vis_polys,
                                  mesh_alpha=args.meshAlpha, elem_alpha=args.elemAlpha,
                                  title=f'Geometry: {args.geom}  (visible pairs only)')


if __name__ == '__main__':
    main()
