import numpy as np
from dataclasses import dataclass
from typing import List, Tuple, Set

@dataclass
class TrussEdge:
    i: int
    j: int
    kind: int = 1

class Truss:
    def __init__(self):
        self.points = np.zeros((0, 3))  # Nx3 array of point coordinates
        self.masses = np.array([])      # N array of point masses
        self.bonds = []                 # List of (i,j) tuples for bonds
        self.ks = np.array([])          # Spring constants for each bond
        self.fixed = set()              # Set of fixed point indices

    def build_rope(self, n: int, m: float = 1.0, m_end: float = 1000.0, l: float = 1.0):
        """Build a simple rope-like truss"""
        self.bonds = [(i, i+1) for i in range(n-1)]
        self.masses = np.ones(n) * m
        self.masses[0] = self.masses[-1] = m_end
        self.points = np.zeros((n, 3))
        x0 = -l * n / 2
        self.points[:, 0] = np.arange(n) * l + x0 + l/2
        self.ks = np.ones(len(self.bonds))
        self.fixed = {0, n-1}

    def build_grid_2d(self, nx: int, ny: int, m: float = 1.0, m_end: float = 1000.0, l: float = 1.0, k: float = 1.0, k_diag: float = -1.0):
        """Build a 2D grid truss"""
        np_total = (nx + 1) * (ny + 1)
        self.points = np.zeros((np_total, 3))
        self.masses = np.ones(np_total) * m
        self.bonds = []
        self.ks = []
        
        # Create points
        for iy in range(ny + 1):
            for ix in range(nx + 1):
                i = iy * (nx + 1) + ix
                self.points[i] = [ix * l, -iy * l, 0]

        # Create bonds
        for iy in range(ny + 1):
            for ix in range(nx + 1):
                i = iy * (nx + 1) + ix
                # Horizontal bonds
                if ix < nx:
                    self.bonds.append((i, i + 1))
                    self.ks.append(k)
                # Vertical bonds
                if iy < ny:
                    self.bonds.append((i, i + nx + 1))
                    self.ks.append(k)
                # Diagonal bonds
                if k_diag > 0:
                    if ix < nx and iy < ny:
                        self.bonds.append((i, i + nx + 2))
                        self.ks.append(k_diag)
                    if ix > 0 and iy < ny:
                        self.bonds.append((i, i + nx))
                        self.ks.append(k_diag)

        self.ks = np.array(self.ks)
        # Fix corners
        self.fixed = {0, nx}
        self.masses[list(self.fixed)] = m_end

    def ngon_truss(self, p0: np.ndarray, p1: np.ndarray, ax: np.ndarray, 
                   n: int = 8, k: float = 1.0):
        """Create a regular n-gon truss"""
        dir = p1 - p0
        r = np.linalg.norm(dir)
        dir = dir / r
        
        # Make ax orthogonal to dir
        ax = ax - np.dot(ax, dir) * dir
        ax = ax / np.linalg.norm(ax)
        side = np.cross(dir, ax)
        
        # Create points using complex number rotation
        points = []
        bonds = []
        rot = 1 + 0j
        drot = np.exp(2j * np.pi / n)
        
        for i in range(n):
            R = dir * rot.real + side * rot.imag
            points.append(p0 + R)
            bonds.append((i, (i-1) % n))
            rot *= drot
            
        self.points = np.array(points)
        self.bonds = bonds
        self.ks = np.ones(len(bonds)) * k
        self.masses = np.ones(n)

    def wheel(self, p0: np.ndarray, p1: np.ndarray, ax: np.ndarray, width: float, n: int = 8, k_scale: float = 1.0, k_dict: List[float] = None):
        """
        Create a wheel-like truss structure.
        
        Args:
            p0: Center point
            p1: Point defining direction and radius
            ax: Axis direction (will be made orthogonal to p1-p0)
            width: Width of the wheel
            n: Number of segments
            k_scale: Global scaling factor for spring constants
            k_dict: List of spring constants for different bond types [long, perp, zigIn, zigOut]
        """
        if k_dict is None:
            k_dict = [1.0,  # long bonds (perimeter)
                      1.0,   # perpendicular bonds
                      1.0,   # inner zigzag
                      1.0]   # outer zigzag
        
        # Bond types
        KIND_LONG = 0   # Along the perimeter
        KIND_PERP = 1   # Perpendicular elements
        KIND_ZIGIN = 2  # Inner zigzag connections
        KIND_ZIGOUT = 3 # Outer zigzag connections
        
        # Calculate directions
        dir = p1 - p0
        r = np.linalg.norm(dir)
        dir = dir / r
        
        # Make ax orthogonal to dir
        ax = ax - np.dot(ax, dir) * dir
        ax = ax / np.linalg.norm(ax)
        side = np.cross(dir, ax)
        
        # Initialize lists
        points = []
        edges = []
        
        # Complex number for rotation
        rot = 1 + 0j
        drot = np.exp(1j * np.pi / n)  # Half angle rotation
        
        dnp = 4  # Points per segment
        i00 = 0  # Base index (using 0-based indexing)
        i000 = i00
        
        # Generate points and edges
        for i in range(n):
            i01 = i00 + 1
            i10 = i00 + 2
            i11 = i00 + 3
            
            # First two points
            R = dir * rot.real + side * rot.imag
            points.append(p0 + R * (r + width))
            points.append(p0 + R * (r - width))
            
            # Rotate halfway
            rot *= drot
            R = dir * rot.real + side * rot.imag
            
            # Next two points
            points.append(p0 + ax * width + R * r)
            points.append(p0 + ax * -width + R * r)
            
            # Rotate halfway again
            rot *= drot
            
            # Add edges within segment
            edges.extend([
                TrussEdge(i00, i01, KIND_PERP),
                TrussEdge(i10, i11, KIND_PERP),
                TrussEdge(i00, i10, KIND_ZIGIN),
                TrussEdge(i00, i11, KIND_ZIGIN),
                TrussEdge(i01, i10, KIND_ZIGIN),
                TrussEdge(i01, i11, KIND_ZIGIN),
            ])
            
            # Add edges to next segment or back to start
            if i < n - 1:
                edges.extend([
                    TrussEdge(i10, i00 + dnp, KIND_ZIGOUT),
                    TrussEdge(i10, i01 + dnp, KIND_ZIGOUT),
                    TrussEdge(i11, i00 + dnp, KIND_ZIGOUT),
                    TrussEdge(i11, i01 + dnp, KIND_ZIGOUT),
                    TrussEdge(i00, i00 + dnp, KIND_LONG),
                    TrussEdge(i01, i01 + dnp, KIND_LONG),
                    TrussEdge(i10, i10 + dnp, KIND_LONG),
                    TrussEdge(i11, i11 + dnp, KIND_LONG),
                ])
            else:
                # Connect back to the first segment
                edges.extend([
                    TrussEdge(i10, i000 + 0, KIND_ZIGOUT),
                    TrussEdge(i10, i000 + 1, KIND_ZIGOUT),
                    TrussEdge(i11, i000 + 0, KIND_ZIGOUT),
                    TrussEdge(i11, i000 + 1, KIND_ZIGOUT),
                    TrussEdge(i00, i000 + 0, KIND_LONG),
                    TrussEdge(i01, i000 + 1, KIND_LONG),
                    TrussEdge(i10, i000 + 2, KIND_LONG),
                    TrussEdge(i11, i000 + 3, KIND_LONG),
                ])
            
            i00 += dnp
            
        # Convert to numpy arrays
        self.points = np.array(points)
        self.bonds = [(e.i, e.j) for e in edges]
        self.ks = np.array([k_dict[e.kind] * k_scale for e in edges])
        self.masses = np.ones(len(points))

    def get_neighbor_list(self) -> List[Set[int]]:
        """Return list of neighbors for each point"""
        neighbors = [set() for _ in range(len(self.points))]
        for i, j in self.bonds:
            neighbors[i].add(j)
            neighbors[j].add(i)
        return neighbors

    def get_bond_vectors(self) -> np.ndarray:
        """Return vectors for all bonds"""
        vecs = np.zeros((len(self.bonds), 3))
        for idx, (i, j) in enumerate(self.bonds):
            vecs[idx] = self.points[j] - self.points[i]
        return vecs

    def get_bond_lengths(self) -> np.ndarray:
        """Return lengths of all bonds"""
        return np.linalg.norm(self.get_bond_vectors(), axis=1)

    def get_rest_lengths(self) -> np.ndarray:
        """Return rest lengths of all bonds based on initial configuration"""
        return np.linalg.norm(self.get_bond_vectors(), axis=1)
        
    def get_pd_quantities(self):
        """
        Get quantities needed for projective dynamics solver.

        Returns:
            Tuple containing:
            - bonds: Array of (i,j) indices for each bond
            - points: Nx3 array of point positions
            - masses: Array of point masses
            - ks: Array of spring constants
            - fixed: List of fixed point indices
            - l0s: Array of rest lengths
            - neighbs: List of neighboring bond indices for each point
        """
        bonds   = np.array(self.bonds)
        l0s     = self.get_rest_lengths()
        neighs = self.get_neighbor_list()
        return bonds, self.points, self.masses, self.ks, list(self.fixed), l0s, neighs


def solve_vbd_numpy(truss: Truss, dt: float, gravity: np.ndarray, *, niter: int = 10,
                    det_eps: float = 1e-6, velocities: np.ndarray | None = None,
                    fixed_points: List[int] | None = None, verbose: int = 0,
                    track_indices: List[int] | None = None,
                    trajectory: List[np.ndarray] | None = None) -> Tuple[np.ndarray, np.ndarray]:
    """Vertex Block Descent reference solver implemented with NumPy.

    Args:
        truss: Input `Truss` instance providing geometry, masses, bonds, stiffnesses, and rest lengths.
        dt: Time step size.
        gravity: External acceleration (3-vector).
        niter: Number of Gauss–Seidel VBD iterations.
        det_eps: Determinant threshold for skipping ill-conditioned solves.
        velocities: Optional initial velocity array (Nx3). Defaults to zeros.
        fixed_points: Optional iterable overriding `truss.fixed`.
        verbose: When >0 prints per-iteration diagnostics.
        track_indices: Optional list of point indices whose positions should be recorded each iteration.
        trajectory: Optional list to append tracked positions (shape=(ntrack,3)) per iteration; initial state appended before iterations.

    Returns:
        Tuple `(positions, velocities)` with updated positions and implied velocities.
    """
    if dt <= 0:
        raise ValueError("dt must be positive")

    bonds, x0, masses, ks, fixed_default, l0s, _ = truss.get_pd_quantities()
    if len(bonds) == 0:
        raise ValueError("truss has no bonds; VBD update undefined")

    x0 = x0.astype(np.float64, copy=True)
    masses = masses.astype(np.float64, copy=False)
    ks = ks.astype(np.float64, copy=False)
    l0s = l0s.astype(np.float64, copy=False)
    gravity = np.asarray(gravity, dtype=np.float64)

    n_points = x0.shape[0]
    if velocities is None:
        velocities = np.zeros_like(x0)
    else:
        velocities = np.asarray(velocities, dtype=np.float64, order="C")
        if velocities.shape != x0.shape:
            raise ValueError("velocities array must match truss points shape")

    fixed_pts = np.array(sorted(set(fixed_default if fixed_points is None else fixed_points)), dtype=int)
    fixed_mask = np.zeros(n_points, dtype=bool)
    if fixed_pts.size:
        fixed_mask[fixed_pts] = True

    bonds = np.asarray(bonds, dtype=np.int32)
    if bonds.ndim != 2 or bonds.shape[1] != 2:
        raise ValueError("bonds array must be shape (nbond, 2)")

    dt_sq = dt * dt
    inv_h2 = 1.0 / dt_sq
    I3 = np.eye(3, dtype=np.float64)

    x = x0.copy()
    y = x0 + dt * velocities + dt_sq * gravity
    if fixed_pts.size:
        x[fixed_mask] = x0[fixed_mask]
        y[fixed_mask] = x0[fixed_mask]

    eye_mass = np.einsum('i,ab->iab', masses * inv_h2, I3)

    track_idx = None
    if track_indices is not None and len(track_indices) > 0:
        track_idx = np.asarray(track_indices, dtype=int)
        if trajectory is not None:
            trajectory.append(x[track_idx].copy())
    elif trajectory is not None:
        trajectory.clear()

    for itr in range(niter):
        if fixed_pts.size:
            x[fixed_mask] = x0[fixed_mask]

        xi = x[bonds[:, 0]]
        xj = x[bonds[:, 1]]
        d = xi - xj
        L = np.linalg.norm(d, axis=1)
        safe_L = np.maximum(L, 1e-9)
        dir_vec = np.zeros_like(d)
        mask_nonzero = L > 1e-9
        dir_vec[mask_nonzero] = d[mask_nonzero] / safe_L[mask_nonzero, None]

        stretch = L - l0s
        grad_edge = ks[:, None] * stretch[:, None] * dir_vec

        grad = (masses * inv_h2)[:, None] * (x - y)
        np.add.at(grad, bonds[:, 0], grad_edge)
        np.add.at(grad, bonds[:, 1], -grad_edge)

        coeff_iso = ks * (1.0 - l0s / safe_L)
        coeff_dir = ks * (l0s / safe_L)
        coeff_iso[~mask_nonzero] = ks[~mask_nonzero]
        coeff_dir[~mask_nonzero] = ks[~mask_nonzero]
        outer = dir_vec[:, :, None] * dir_vec[:, None, :]
        H_edges = coeff_iso[:, None, None] * I3 + coeff_dir[:, None, None] * outer

        H = eye_mass.copy()
        np.add.at(H, bonds[:, 0], H_edges)
        np.add.at(H, bonds[:, 1], H_edges)

        grad[fixed_mask] = 0.0
        if fixed_pts.size:
            H[fixed_mask] = eye_mass[fixed_mask]

        max_dx = 0.0
        for vid in range(n_points):
            if fixed_mask[vid]:
                continue
            Hi = H[vid]
            gi = grad[vid]
            det = np.linalg.det(Hi)
            if abs(det) < det_eps:
                continue
            dx = -np.linalg.solve(Hi, gi)
            x[vid] += dx
            step = np.linalg.norm(dx)
            if step > max_dx:
                max_dx = step

        if verbose:
            print(f"solve_vbd_numpy itr {itr}: max |dx| = {max_dx:.3e}")

        if track_idx is not None and trajectory is not None:
            trajectory.append(x[track_idx].copy())

    if fixed_pts.size:
        x[fixed_mask] = x0[fixed_mask]

    v = (x - x0) / dt
    return x.astype(np.float64, copy=False), v.astype(np.float64, copy=False)

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import numpy as np
    import plot_utils as pu
    
    # Create a wheel truss with different stiffnesses for different bond types
    truss = Truss()
    p0 = np.array([0., 0., 0.])
    p1 = np.array([2., 0., 0.])
    ax_dir = np.array([0., 0., 1.])
    
    # Different stiffnesses for different bond types
    k_dict = [10000.0,  # long bonds (perimeter)
              5000.0,   # perpendicular bonds
              2000.0,   # inner zigzag
              2000.0]   # outer zigzag
    
    truss.wheel(p0, p1, ax_dir, width=0.4, n=12,  k_scale=1.0, k_dict=k_dict)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))
    
    # Plot with uniform colors
    pu.plot_truss( truss.points, truss.bonds, ax1, edge_color='k', edge_alpha=0.7, point_color='b', point_size=30)
    ax1.set_title('Uniform Edge Colors')
    
    # Plot with colors based on stiffness
    pu.plot_truss(truss.points, truss.bonds, ax2, color_by_stiffness=True, edge_alpha=0.7, point_color='k', point_size=30, cmap='viridis', ks=truss.ks)
    ax2.set_title('Edges Colored by Stiffness')
    
    plt.tight_layout()
    plt.show()
