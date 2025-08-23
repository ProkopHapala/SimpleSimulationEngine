import numpy as np
import argparse
import matplotlib.pyplot as plt
from typing import Optional
from matplotlib import colors as mcolors

# Reuse modules from this repo
from PolygonRasterization import HexRasterizerNP

# -----------------------------
# Utilities
# -----------------------------

def _norm(v: np.ndarray, eps: float = 1e-12) -> float:
    return float(np.linalg.norm(v) + 0.0)

def _normalize(v: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    n = np.linalg.norm(v)
    assert n > eps, f"Zero-length vector in normalize (|v|={n})"
    return v / n

# -----------------------------
# Geometry: per-face local frame and projection
# -----------------------------

def build_face_frame(points3: np.ndarray, tol: float = 1e-10, verbosity: int = 0):
    """
    Given a set of 3D vertices (>=3, roughly planar), build an orthonormal basis (U,V,N) and origin O.
    - N = normalized cross of two non-collinear edges
    - U = longest edge projected onto plane and normalized
    - V = N x U
    Returns (O, U, V, N)
    """
    assert points3.shape[0] >= 3 and points3.shape[1] == 3
    O = points3[0]

    # find two non-collinear edges
    edges = [points3[i] - O for i in range(1, points3.shape[0])]
    # pick longest edge as e0
    lens = [np.linalg.norm(e) for e in edges]
    i0 = int(np.argmax(lens))
    e0 = edges[i0]
    # find e1 not collinear with e0
    e1 = None
    for e in edges:
        c = np.linalg.norm(np.cross(e0, e))
        if c > tol:
            e1 = e
            break
    assert e1 is not None, "Degenerate face: all points collinear"

    N = _normalize(np.cross(e0, e1))

    # U: project e0 onto plane to avoid any tiny normal component
    e0p = e0 - N * np.dot(N, e0)
    if np.linalg.norm(e0p) < tol:
        # fall back to e1
        e0p = e1 - N * np.dot(N, e1)
    U = _normalize(e0p)
    V = np.cross(N, U)

    # Planarity sanity check: distance of all points from plane
    d = points3 - O
    off = np.abs(np.dot(d, N))
    max_off = float(np.max(off))
    if verbosity >= 2:
        print(f"[planarity] max distance from plane: {max_off:.3e}")
    assert max_off < 1e-6, f"Face not planar within tolerance, max_off={max_off}"
    return O, U, V, N


def project_to_uv(points3: np.ndarray, O: np.ndarray, U: np.ndarray, V: np.ndarray) -> np.ndarray:
    d = points3 - O  # (m,3)
    u = d @ U  # (m,)
    v = d @ V  # (m,)
    return np.stack([u, v], axis=1)


def uv_to_xyz(com_uv: np.ndarray, O: np.ndarray, U: np.ndarray, V: np.ndarray) -> np.ndarray:
    return O + U * com_uv[0] + V * com_uv[1]

# -----------------------------
# Face rasterization -> 3D elements
# -----------------------------

def rasterize_face_to_elems(points3: np.ndarray, face_id: int, hex_size: float, area_keep_ratio: float, verbosity: int = 0):
    """
    Rasterize a single planar face into hex elements using the 2D rasterizer on the UV plane.
    Returns list of dicts with keys: center(3,), normal(3,), area(float), face_id(int)
    """
    O, U, V, N = build_face_frame(points3, verbosity=verbosity)
    poly_uv = project_to_uv(points3, O, U, V)

    r = HexRasterizerNP(poly_uv, float(hex_size), verbosity=verbosity)
    elems2d = r.rasterize(area_keep_ratio=area_keep_ratio)

    out = []
    for k, el in elems2d.items():  # el has fields: area, com (2D), geoms (list)
        center3 = uv_to_xyz(el.com, O, U, V)
        # Map each 2D geom polygon to 3D using the face frame
        polys3d = []
        for geom in el.geoms:
            if geom is None or len(geom) == 0: continue
            poly3 = np.array([uv_to_xyz(p, O, U, V) for p in geom], dtype=np.float64)
            polys3d.append(poly3)
        out.append({
            'center': center3.astype(np.float64, copy=False),
            'normal': N.astype(np.float64, copy=False),
            'area': float(el.area),
            'face_id': int(face_id),
            'polys3d': polys3d
        })
    return out

# -----------------------------
# Assembly for radiosity solver from Radiosity.py (3D-aware)
# -----------------------------

def assemble_elements(elems: list, eps_front: float, eps_back: float, face_heat_W_per_m2: np.ndarray):
    """
    Build the elements dict expected by Radiosity.py solvers.
    - elems: list of {'center':(3,), 'normal':(3,), 'area':float, 'face_id':int}
    - face_heat_W_per_m2: (n_faces,) power density per face
    Returns dict with keys: center(n,3), normal(n,3), area(n,), eps_front, eps_back, rho_front, rho_back, heat_in(n,), length(n,)
    """
    n = len(elems)
    C = np.array([e['center'] for e in elems], dtype=np.float64)
    N = np.array([e['normal'] for e in elems], dtype=np.float64)
    A = np.array([e['area'] for e in elems], dtype=np.float64)
    face_ids = np.array([e['face_id'] for e in elems], dtype=np.int32)

    eps_f = float(eps_front) * np.ones(n, dtype=np.float64)
    eps_b = float(eps_back) * np.ones(n, dtype=np.float64)
    rho_f = 1.0 - eps_f
    rho_b = 1.0 - eps_b

    Pdens = face_heat_W_per_m2[face_ids.astype(int)]  # per-element density by face
    heat_in = Pdens * A  # total power per element

    P = [e.get('polys3d', None) for e in elems]
    elements = {
        'center': C,
        'normal': N,
        'area': A,
        'eps_front': eps_f,
        'eps_back': eps_b,
        'rho_front': rho_f,
        'rho_back': rho_b,
        'heat_in': heat_in,
        # compatibility placeholder (not used in 3D functions below)
        'length': np.ones(n, dtype=np.float64),
        # optional visualization payload
        'polys3d': P,
    }
    return elements

# -----------------------------
# Radiosity core (3D-aware)
# -----------------------------

def compute_view_factors_3d(elements: dict):
    """
    Geometry-only F using midpoint inverse-square kernel with absolute cosines in 3D.
    F[i,j] ~ (|n_i·rhat_ij| |n_j·rhat_ji| / (pi r^2)) * A_j
    """
    C = elements['center']  # (n,3)
    N = elements['normal']  # (n,3)
    A = elements['area']    # (n,)
    n = C.shape[0]

    R = C[np.newaxis, :, :] - C[:, np.newaxis, :]        # (n,n,3) r_ij = c_j - c_i
    d2 = np.einsum('ijk,ijk->ij', R, R) + 1e-12
    d = np.sqrt(d2)
    Rhat = R / d[:, :, None]

    cos_i = np.abs(np.einsum('ik,ijk->ij', N, Rhat))      # |n_i · rhat_ij|
    cos_j = np.abs(np.einsum('jk,ijk->ij', N, -Rhat))     # |n_j · rhat_ji|

    F = (cos_i * cos_j / (np.pi * d2)) * A[np.newaxis, :]
    np.fill_diagonal(F, 0.0)
    F = np.clip(F, 0.0, 1.0)
    return F


def solve_radiosity_system_3d(elements: dict):
    """
    3D version following Radiosity.py logic with directional kernels and receiver-side reflectivity.
    Solves (I - W) B = heat_in / area
    Returns B (n,)
    """
    C = elements['center']
    N = elements['normal']
    A = elements['area']
    n = C.shape[0]
    I = np.eye(n)

    R = C[np.newaxis, :, :] - C[:, np.newaxis, :]
    d2 = np.einsum('ijk,ijk->ij', R, R) + 1e-12
    d = np.sqrt(d2)
    Rhat = R / d[:, :, None]

    cos_i = np.einsum('ik,ijk->ij', N, Rhat)      # signed cos on receiver i
    cos_j = np.einsum('jk,ijk->ij', N, -Rhat)     # signed cos on emitter j

    cip = np.clip(cos_i, 0.0, None)
    cin = np.clip(-cos_i, 0.0, None)
    cjp = np.clip(cos_j, 0.0, None)
    cjn = np.clip(-cos_j, 0.0, None)

    K_ff = (cip * cjp / (np.pi * d2)) * A[np.newaxis, :]
    K_fb = (cip * cjn / (np.pi * d2)) * A[np.newaxis, :]
    K_bf = (cin * cjp / (np.pi * d2)) * A[np.newaxis, :]
    K_bb = (cin * cjn / (np.pi * d2)) * A[np.newaxis, :]

    for K in (K_ff, K_fb, K_bf, K_bb):
        np.fill_diagonal(K, 0.0)

    F_front_i = K_ff + K_fb
    F_back_i  = K_bf + K_bb
    W = elements['rho_front'][:, None] * F_front_i + elements['rho_back'][:, None] * F_back_i

    M = I - W
    P_density = elements['heat_in'] / (elements['area'] + 1e-12)
    B = np.linalg.solve(M, P_density)
    return B


def calculate_temperatures_3d(B: np.ndarray, elements: dict):
    """
    Thin-sheet balance with two-side eps following Radiosity.py, in 3D.
    Returns T4 (n,)
    """
    C = elements['center']
    N = elements['normal']
    A = elements['area']

    R = C[np.newaxis, :, :] - C[:, np.newaxis, :]
    d2 = np.einsum('ijk,ijk->ij', R, R) + 1e-12
    d = np.sqrt(d2)
    Rhat = R / d[:, :, None]

    cos_i = np.einsum('ik,ijk->ij', N, Rhat)
    cos_j = np.einsum('jk,ijk->ij', N, -Rhat)

    cip = np.clip(cos_i, 0.0, None)
    cin = np.clip(-cos_i, 0.0, None)
    cjp = np.clip(cos_j, 0.0, None)
    cjn = np.clip(-cos_j, 0.0, None)

    K_ff = (cip * cjp / (np.pi * d2)) * A[np.newaxis, :]
    K_fb = (cip * cjn / (np.pi * d2)) * A[np.newaxis, :]
    K_bf = (cin * cjp / (np.pi * d2)) * A[np.newaxis, :]
    K_bb = (cin * cjn / (np.pi * d2)) * A[np.newaxis, :]
    for K in (K_ff, K_fb, K_bf, K_bb):
        np.fill_diagonal(K, 0.0)

    F_front_i = K_ff + K_fb
    F_back_i  = K_bf + K_bb

    H_front = (F_front_i * B[None, :]).sum(axis=1)
    H_back  = (F_back_i  * B[None, :]).sum(axis=1)

    num = (elements['heat_in'] / (elements['area'] + 1e-12)) + elements['eps_front'] * H_front + elements['eps_back'] * H_back
    den = (elements['eps_front'] + elements['eps_back'] + 1e-12)
    T4 = num / den
    return T4

# -----------------------------
# Visualization (3D)
# -----------------------------

def visualize_3d(vertices: np.ndarray, faces: list, centers: np.ndarray, T4: np.ndarray, normals: Optional[np.ndarray] = None,
                 draw_mesh: bool = True, draw_normals: bool = False, labels: str = 'none', elem_polys3d: Optional[list] = None):
    # Import 3D toolkit lazily to avoid issues when running headless
    try:
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        have_3d = True
    except Exception as e:
        print(f"[warn] 3D plotting unavailable ({e}); skipping mesh/normal drawing")
        have_3d = False

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d' if have_3d else None)

    if have_3d and draw_mesh and (vertices is not None) and (faces is not None):
        polys = [[vertices[idx] for idx in f] for f in faces]
        mesh = Poly3DCollection(polys, alpha=0.15, facecolor='lightgray', edgecolor='k', linewidths=0.5)
        ax.add_collection3d(mesh)

    # Scatter elements colored by T4 (or T)
    if labels == 'temp':
        vals = np.maximum(T4, 0.0) ** 0.25
        cbar_label = 'Temperature (a.u.)'
    else:
        vals = T4
        cbar_label = 'Proportional to T^4 (a.u.)'

    if have_3d:
        sc = ax.scatter(centers[:, 0], centers[:, 1], centers[:, 2], c=vals, cmap='inferno', s=16, depthshade=True)
    else:
        sc = ax.scatter(centers[:, 0], centers[:, 1], c=vals, cmap='inferno', s=16)
    cbar = fig.colorbar(sc, ax=ax, shrink=0.7)
    cbar.set_label(cbar_label)

    # Draw per-element polygons colored by vals if available
    if have_3d and (elem_polys3d is not None):
        poly_list = []
        face_cols = []
        norm = mcolors.Normalize(vmin=float(np.min(vals)), vmax=float(np.max(vals)))
        cmap = plt.cm.inferno
        for i, polyset in enumerate(elem_polys3d):
            if polyset is None: continue
            for poly in polyset:
                if poly is None or len(poly) < 3: continue
                poly_list.append(poly)
                face_cols.append(cmap(norm(float(vals[i]))))
        if len(poly_list):
            emesh = Poly3DCollection(poly_list, edgecolor='k', linewidths=0.3)
            emesh.set_facecolor(face_cols)
            emesh.set_alpha(1.0)
            ax.add_collection3d(emesh)

    if have_3d and draw_normals and normals is not None:
        # draw a decimated subset to keep plot readable
        step = max(1, centers.shape[0] // 200)
        P = centers[::step]
        Vv = normals[::step]
        scale = np.mean(np.linalg.norm(centers.max(axis=0) - centers.min(axis=0))) * 0.05
        ax.quiver(P[:, 0], P[:, 1], P[:, 2], Vv[:, 0], Vv[:, 1], Vv[:, 2], length=scale, color='cyan', linewidth=0.5)

    # axes limits
    all_pts = vertices if vertices is not None else centers
    mins = all_pts.min(axis=0)
    maxs = all_pts.max(axis=0)
    span = max(maxs - mins)
    mid = 0.5 * (mins + maxs)
    ax.set_xlim(mid[0] - span/2, mid[0] + span/2)
    ax.set_ylim(mid[1] - span/2, mid[1] + span/2)
    if have_3d:
        ax.set_zlim(mid[2] - span/2, mid[2] + span/2)
        # Equal aspect for all three axes if supported (Matplotlib >=3.3)
        if hasattr(ax, 'set_box_aspect'):
            ax.set_box_aspect((1, 1, 1))

    ax.set_title('Radiosity 3D: element centers colored by temperature' + (" (2D fallback)" if not have_3d else ""))
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    if have_3d:
        ax.set_zlabel('Z')
    plt.tight_layout()
    # if have_3d:
    #     plt.show()
    # else:
    #     # Save 2D fallback to file to avoid GUI requirements
    #     out_png = 'OUT-radiosity3D.png'
    #     plt.savefig(out_png, dpi=120)
    #     plt.close(fig)
    #     print(f"[info] Saved 2D fallback plot to {out_png}")

# -----------------------------
# Minimal test geometry: unit cube
# -----------------------------

def unit_cube():
    V = np.array([
        [0, 0, 0],
        [1, 0, 0],
        [1, 1, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [1, 1, 1],
        [0, 1, 1],
    ], dtype=np.float64)
    F = [
        [0, 1, 2, 3],  # bottom z=0
        #[4, 5, 6, 7],  # top z=1
        #[0, 1, 5, 4],  # y=0
        #[1, 2, 6, 5],  # x=1
        [2, 3, 7, 6],  # y=1
        #[3, 0, 4, 7],  # x=0
    ]
    return V, F

# -----------------------------
# Main
# -----------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='3D Radiosity via 2D hex rasterization per face')
    parser.add_argument('--hex', type=float, default=0.35, help='Hexagon size (grid step in face UV units)')
    parser.add_argument('--size', type=float, default=None, help='Alias for --hex (hexagon size)')
    parser.add_argument('--areaFactor', type=float, default=0.7, help='Raster keep ratio for partial hex (A_clip/A_full)')
    parser.add_argument('--epsFront', type=float, default=0.8, help='Front emissivity')
    parser.add_argument('--epsBack', type=float, default=0.8, help='Back emissivity')
    parser.add_argument('--hotFace', type=int, default=0, help='Index of hot face in test geometry')
    parser.add_argument('--P_hot', type=float, default=100.0, help='Heating power density [W/m^2] on hot face')
    parser.add_argument('--verbosity', type=int, default=1, help='Verbosity level')
    parser.add_argument('--showMatrix', action='store_true', help='Show view factor matrix (may be large)')
    parser.add_argument('--labels', type=str, default='none', choices=['none','T4','temp'], help='Colorbar label mode')
    parser.add_argument('--plotNormals', action='store_true', help='Draw a subset of normals')
    parser.add_argument('--noPlot', action='store_true', help='Disable plotting for headless runs')
    args = parser.parse_args()

    # Geometry: built-in cube for now
    V, F = unit_cube()

    # Heating per face
    nF = len(F)
    face_heat = np.zeros(nF, dtype=np.float64)
    if (0 <= args.hotFace < nF):
        face_heat[args.hotFace] = args.P_hot

    # Resolve hex size (prefer --size if provided for consistency with PolygonRasterization.py)
    hex_size = args.size if (args.size is not None) else args.hex

    # Rasterize each face -> elements
    all_elems = []
    for fi, face in enumerate(F):
        pts = V[np.array(face, dtype=int)]
        elems = rasterize_face_to_elems(pts, face_id=fi, hex_size=hex_size, area_keep_ratio=args.areaFactor, verbosity=args.verbosity)
        if args.verbosity >= 1:
            print(f"face {fi}: +{len(elems)} elements")
        all_elems.extend(elems)

    assert len(all_elems) > 0, "No elements produced; check hex size vs geometry scale"

    elements = assemble_elements(all_elems, eps_front=args.epsFront, eps_back=args.epsBack, face_heat_W_per_m2=face_heat)

    if args.verbosity >= 1:
        A_tot = elements['area'].sum()
        print(f"Assembled {len(all_elems)} elements, total area ~ {A_tot:.3f}")

    # Solve radiosity (3D-aware)
    Fgeom = compute_view_factors_3d(elements)
    B = solve_radiosity_system_3d(elements)
    T4 = calculate_temperatures_3d(B, elements)

    if args.showMatrix and (not args.noPlot):
        # local simple matrix view
        fig, ax = plt.subplots(figsize=(6, 5))
        im = ax.imshow(Fgeom, interpolation='nearest', cmap='viridis', origin='lower')
        ax.set_title('View Factor Matrix (linear, 3D)')
        fig.colorbar(im, ax=ax)
        plt.show()

    # Visualize
    if not args.noPlot:
        visualize_3d(vertices=V, faces=F, centers=elements['center'], T4=T4, normals=elements['normal'],
                     draw_mesh=True, draw_normals=args.plotNormals, labels=args.labels, elem_polys3d=elements.get('polys3d'))
        plt.show()



