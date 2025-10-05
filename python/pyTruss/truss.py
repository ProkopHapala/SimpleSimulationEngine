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
        self._rest_lengths = np.array([])

    def build_rope(self, n: int, m: float = 1.0, m_end: float = 1000.0, l: float = 1.0):
        """Build a simple rope-like truss"""
        self._rest_lengths = np.array([])
        self.bonds = [(i, i+1) for i in range(n-1)]
        self.masses = np.ones(n) * m
        self.masses[0] = self.masses[-1] = m_end
        self.points = np.zeros((n, 3))
        x0 = -l * n / 2
        self.points[:, 0] = np.arange(n) * l + x0 + l/2
        self.ks = np.ones(len(self.bonds))
        self.fixed = {0, n-1}
        self._update_rest_lengths()

    def build_grid_2d(self, nx: int, ny: int, m: float = 1.0, m_end: float = 1000.0, l: float = 1.0, k: float = 1.0, k_diag: float = -1.0):
        """Build a 2D grid truss"""
        self._rest_lengths = np.array([])
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
        self._update_rest_lengths()

    def ngon_truss(self, p0: np.ndarray, p1: np.ndarray, ax: np.ndarray, 
                   n: int = 8, k: float = 1.0):
        """Create a regular n-gon truss"""
        self._rest_lengths = np.array([])
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
        self._update_rest_lengths()

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
        self._update_rest_lengths()

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
        if self._rest_lengths.size != len(self.bonds):
            self._update_rest_lengths()
        return self._rest_lengths.copy()

    def _update_rest_lengths(self):
        if len(self.bonds) == 0:
            self._rest_lengths = np.array([])
        else:
            self._rest_lengths = np.linalg.norm(self.get_bond_vectors(), axis=1)
        
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

    def color_graph(self, seed=None) -> Tuple[np.ndarray, List[List[int]]]:
        """
        Partitions the truss vertices into colors for parallel Gauss-Seidel solvers.

        This method implements a parallel randomized coloring algorithm inspired by
        the principles in the "Vivace" paper (e.g., Luby's MIS algorithm). In each
        round, it identifies a "Maximal Independent Set" (MIS) of the remaining
        uncolored vertices and assigns them a new color. An independent set is a
        group of vertices where no two are connected by a bond.
        see: 
            - Marco Fratarcangeli et al. "Vivace: A Practical Gauss-seidel Method for Stable Soft Body Dynamics", ACM Transactions on Graphics (SIGGRAPH Asia), (2016)
            - http://doi.acm.org/10.1145/2980179.2982437
            - https://mfratarcangeli.github.io/publication/sigasia2016/

        Args:
            seed (int, optional): A seed for the random number generator to ensurereproducible colorings. Defaults to None.

        Returns:
            Tuple[np.ndarray, List[List[int]]]:
            - A numpy array where the index is the vertex ID and the value is its color ID.
            - A list of lists, where each inner list contains the vertex IDs for one color.
        """
        if seed is not None:
            np.random.seed(seed)

        n_points = len(self.points)
        if n_points == 0:
            return np.array([]), []
            
        neighbors = self.get_neighbor_list()
        
        vertex_colors = np.full(n_points, -1, dtype=int)
        uncolored_nodes = set(range(n_points))
        current_color_id = 0

        while uncolored_nodes:
            # 1. Assign a random value to each uncolored node.
            # This is the "parallel" step where each node acts independently.
            random_values = {node: np.random.rand() for node in uncolored_nodes}
            
            independent_set = set()
            
            # 2. Identify the Maximal Independent Set (MIS).
            # A node joins the MIS if its random value is higher than all of its
            # uncolored neighbors. This is a local, parallelizable check.
            for node in uncolored_nodes:
                is_local_max = True
                for neighbor in neighbors[node]:
                    if neighbor in uncolored_nodes and random_values.get(neighbor, -1) > random_values[node]:
                        is_local_max = False
                        break
                if is_local_max:
                    independent_set.add(node)
            
            # 3. Assign the current color to all nodes in the MIS.
            for node in independent_set:
                vertex_colors[node] = current_color_id
            
            # 4. Update the set of uncolored nodes and advance the color.
            uncolored_nodes -= independent_set
            current_color_id += 1
            
        # 5. Format the output into partitions.
        num_colors = current_color_id
        partitions = [[] for _ in range(num_colors)]
        for i, color in enumerate(vertex_colors):
            if color != -1:
                partitions[color].append(i)
        
        return vertex_colors, partitions

    def verify_graph_coloring( self, vertex_colors: np.ndarray) -> None:
        if vertex_colors.size != len(self.points):
            raise ValueError(f"color array has length {vertex_colors.size}, expected {len(self.points)}")
        neighbors = self.get_neighbor_list()
        conflicts = []
        for i, neighs in enumerate(neighbors):
            for j in neighs:
                if j <= i:
                    continue
                if vertex_colors[i] == vertex_colors[j]:
                    conflicts.append((i, j, int(vertex_colors[i])))
        if conflicts:
            for i, j, color in conflicts:
                print(f"Coloring conflict: vertices {i} and {j} share color {color}")
            raise ValueError(f"Invalid graph coloring: detected {len(conflicts)} conflicts")


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
