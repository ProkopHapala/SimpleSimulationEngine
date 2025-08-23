import numpy as np
import matplotlib.pyplot as plt
import argparse
from matplotlib.colors import LogNorm

# --- Part 1: Discretization and Data Structure Setup ---

def discretize_surfaces(surfaces, default_width=1.0):
    """
    Takes a list of surfaces and discretizes them into smaller elements.

    Args:
        surfaces (list): A list of tuples, where each tuple defines a surface. Supported forms:
                         1) ((x1, y1), (x2, y2), n, epsilon, heat_in)
                         2) ((x1, y1), (x2, y2), n, epsilon, heat_in, width)
                         3) ((x1, y1), (x2, y2), n, eps_front, eps_back, heat_in, width)
                         4) LEGACY: ((x1, y1), (x2, y2), n, eps_front, eps_back, heat_front, heat_back, width)
                            In legacy (4), heat_in := heat_front + heat_back.
                         Two-sided variants are stored as a single geometric element per segment
                         with per-side emissivities. There is only one heat_in per element.

    Returns:
        dict: A dictionary of NumPy arrays representing all elements.
              Keys: 'p1', 'p2', 'center', 'length', 'normal', 'eps_front', 'eps_back',
                    'rho_front', 'rho_back', 'heat_in', 'width', 'area'.
    """
    all_p1 = []
    all_p2 = []
    all_eps_front = []
    all_eps_back  = []
    all_heat_in   = []
    all_width = []

    for surf in surfaces:
        # Backward-compatible unpacking and feature-detection
        # Recognized lengths: 5,6 (single-sided), 7 (two-sided, single heat), 8 (legacy two heats)
        if len(surf) == 5:
            (p1, p2, n_elements, epsilon_f, heat_in) = surf
            epsilon_b = None
            width = default_width
        elif len(surf) == 6:
            (p1, p2, n_elements, epsilon_f, heat_in, width) = surf
            epsilon_b = None
        elif len(surf) == 7:
            # Two-sided with single heat_in and width
            (p1, p2, n_elements, epsilon_f, epsilon_b, heat_in, width) = surf
        elif len(surf) == 8:
            # Legacy: two heats provided; sum them
            (p1, p2, n_elements, epsilon_f, epsilon_b, heat_f, heat_b, width) = surf
            heat_in = heat_f + heat_b
        else:
            raise ValueError(f"Unsupported surface tuple length {len(surf)} for {surf}")

        p1 = np.array(p1)
        p2 = np.array(p2)
        # Generate points along the line segment
        points = np.linspace(p1, p2, n_elements + 1)
        
        # Create elements from the points
        seg_p1 = points[:-1]
        seg_p2 = points[1:]

        # Do not duplicate elements. Store per-side properties per element.
        all_p1.extend(seg_p1)
        all_p2.extend(seg_p2)
        if epsilon_b is None:
            # Single-sided: mirror emissivity to back; heating is a single scalar per element
            all_eps_front.extend([epsilon_f] * n_elements)
            all_eps_back.extend([epsilon_f] * n_elements)
            all_heat_in.extend([heat_in / n_elements] * n_elements)
        else:
            all_eps_front.extend([epsilon_f] * n_elements)
            all_eps_back.extend([epsilon_b] * n_elements)
            all_heat_in.extend([heat_in / n_elements] * n_elements)
        all_width.extend([width] * n_elements)

    # Convert lists to NumPy arrays for efficient computation
    elements = {
        'p1': np.array(all_p1),
        'p2': np.array(all_p2)
    }
    
    # Calculate derived properties
    elements['center'] = (elements['p1'] + elements['p2']) / 2
    
    vec = elements['p2'] - elements['p1']
    elements['length'] = np.linalg.norm(vec, axis=1)
    
    # Calculate normals (90-degree counter-clockwise rotation)
    # The normal points "out" from the surface.
    # For a vector (dx, dy), the perpendicular is (-dy, dx).
    normals = np.array([-vec[:, 1], vec[:, 0]]).T
    elements['normal'] = normals / np.linalg.norm(normals, axis=1)[:, np.newaxis]
    
    elements['eps_front'] = np.array(all_eps_front)
    elements['eps_back']  = np.array(all_eps_back)
    elements['rho_front'] = 1.0 - elements['eps_front']
    elements['rho_back']  = 1.0 - elements['eps_back']
    elements['heat_in']   = np.array(all_heat_in)
    elements['width'] = np.array(all_width) if len(all_width) > 0 else np.ones_like(elements['length']) * default_width
    elements['area'] = elements['length'] * elements['width']

    return elements

# --- Part 2: View Factor Calculation ---

def compute_view_factors(elements, apply_visibility=True):
    """
    Computes the geometry-only coupling matrix F using a midpoint inverse-square kernel
    with absolute cosine projections (Lambertian geometry term) and emitter area weighting.
    This is robust for arbitrary orientations (parallel, perpendicular, oblique).

    F[i, j] approximates the fraction of energy leaving element i that reaches j.

    Args:
        elements (dict): The dictionary of element properties.
        apply_visibility (bool): Ignored. Kept for backward compatibility.

    Returns:
        np.ndarray: An (N, N) matrix of view factors.
    """
    n = len(elements['length'])

    C = elements['center']        # (n,2)
    N = elements['normal']        # (n,2)
    A = elements['area']          # (n,)

    # Vector from i->j and distances
    R = C[np.newaxis, :, :] - C[:, np.newaxis, :]          # (n,n,2) r_ij = c_j - c_i
    d2 = np.einsum('ijk,ijk->ij', R, R) + 1e-12            # squared distance with epsilon
    d = np.sqrt(d2)
    Rhat = R / d[:, :, np.newaxis]

    # Absolute cosine projections on receiver and emitter sides
    cos_i = np.abs(np.einsum('ik,ijk->ij', N, Rhat))       # |n_i · rhat_ij|
    cos_j = np.abs(np.einsum('jk,ijk->ij', N, -Rhat))      # |n_j · (-rhat_ij)| = |n_j · rhat_ji|

    # Midpoint-based kernel (2D analog of radiosity): base = (cos_i*cos_j)/(pi*r^2) * A_j
    F = (cos_i * cos_j / (np.pi * d2)) * A[np.newaxis, :]

    # Zero self-coupling and clip numerically
    np.fill_diagonal(F, 0.0)
    F = np.clip(F, 0.0, 1.0)

    # DEBUG diagnostics
    row_sums = F.sum(axis=1)
    print(f"   [DEBUG] F_geom:  min={F.min():.3e} max={F.max():.3e} nnz={np.count_nonzero(F)}")
    print(f"   [DEBUG] Row sums (geom): min={row_sums.min():.3e} max={row_sums.max():.3e} mean={row_sums.mean():.3e}")

    return F

# --- Part 3: Solving the Radiosity System ---

def solve_radiosity_system(F, elements):
    """
    Solves the linear system to find the radiosity of each element.
    The system is: (I - F) * B = P_density
    where P_density is the heating power per unit area.

    Args:
        F (np.ndarray): The view factor matrix.
        elements (dict): The dictionary of element properties.

    Returns:
        np.ndarray: A vector B of radiosities for each element.
    """
    n = len(elements['length'])
    I = np.identity(n)

    # Geometry terms at midpoints with signed cosines to split front/back per side
    C = elements['center']
    N = elements['normal']
    A = elements['area']

    R = C[np.newaxis, :, :] - C[:, np.newaxis, :]  # (n,n,2) r_ij = c_j - c_i
    d2 = np.einsum('ijk,ijk->ij', R, R) + 1e-12
    d = np.sqrt(d2)
    Rhat = R / d[:, :, np.newaxis]

    cos_i = np.einsum('ik,ijk->ij', N, Rhat)       # signed cos on receiver i
    cos_j = np.einsum('jk,ijk->ij', N, -Rhat)      # signed cos on emitter j

    cip = np.clip(cos_i, 0.0, None)                # front of i
    cin = np.clip(-cos_i, 0.0, None)               # back of i
    cjp = np.clip(cos_j, 0.0, None)                # front of j
    cjn = np.clip(-cos_j, 0.0, None)               # back of j

    K_ff = (cip * cjp / (np.pi * d2)) * A[np.newaxis, :]
    K_fb = (cip * cjn / (np.pi * d2)) * A[np.newaxis, :]
    K_bf = (cin * cjp / (np.pi * d2)) * A[np.newaxis, :]
    K_bb = (cin * cjn / (np.pi * d2)) * A[np.newaxis, :]

    # Zero self-coupling (numerically stabilize)
    for K in (K_ff, K_fb, K_bf, K_bb):
        np.fill_diagonal(K, 0.0)

    # Receiver-side reflectivity applied per receiving side
    F_front_i = K_ff + K_fb
    F_back_i  = K_bf + K_bb
    W = elements['rho_front'][:, None] * F_front_i + elements['rho_back'][:, None] * F_back_i

    # System matrix
    M = I - W

    # Right-hand side: single heat_in per element -> heating per area
    P_density = elements['heat_in'] / (elements['area'] + 1e-12)

    # Solve the system M * B = P_density
    B = np.linalg.solve(M, P_density)
    
    return B

# --- Part 4: Temperature Calculation ---

def calculate_temperatures(B, F, elements):
    """
    Calculates the final steady-state temperature of each element.
    T^4 is proportional to (HeatingPower/Area + AbsorbedIrradiation) / Epsilon

    Args:
        B (np.ndarray): The solved radiosity vector.
        F (np.ndarray): The view factor matrix.
        elements (dict): The dictionary of element properties.

    Returns:
        np.ndarray: A vector of T^4 values for each element.
    """
    # Recompute directional geometric kernels consistent with the solver (four combinations)
    C = elements['center']
    N = elements['normal']
    A = elements['area']

    R = C[np.newaxis, :, :] - C[:, np.newaxis, :]
    d2 = np.einsum('ijk,ijk->ij', R, R) + 1e-12
    d = np.sqrt(d2)
    Rhat = R / d[:, :, np.newaxis]

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

    # Thin sheet energy balance with one heat_in per element:
    # (eps_f + eps_b)*T^4 = heat_in/area + eps_f*H_front + eps_b*H_back
    num = (elements['heat_in'] / (elements['area'] + 1e-12)) + elements['eps_front'] * H_front + elements['eps_back'] * H_back
    den = (elements['eps_front'] + elements['eps_back'] + 1e-12)
    T4 = num / den
    return T4

# --- Part 5: Visualization ---

def visualize_results(surfaces, elements, T4, labels=None):
    """
    Creates a plot of the geometry and the temperature of each element.
    """
    #plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot the original wireframe geometry (support variable-length tuples)
    for surf in surfaces:
        p1, p2 = surf[0], surf[1]
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color='gray', linewidth=2, zorder=1)

    # Aggregate coincident centers (front/back sides share same center); take max T4 so hot side is visible
    C = elements['center']
    T = T4
    # Unique by rows
    C_unique, inv = np.unique(C, axis=0, return_inverse=True)
    T_agg = np.full(len(C_unique), -np.inf)
    for i, u in enumerate(inv):
        if T[i] > T_agg[u]: T_agg[u] = T[i]
    # Replace -inf with 0 for safety if any unmatched
    T_agg[T_agg == -np.inf] = 0.0

    # Create a scatter plot of unique centers
    scatter = ax.scatter(
        C_unique[:, 0],
        C_unique[:, 1],
        c=T_agg,
        cmap='inferno', # A good colormap for heat
        s=50, # Size of the dots
        zorder=2
    )

    # Add a color bar to show the temperature scale
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label('Proportional to T^4 (Temperature⁴)', fontsize=12)

    # Optional labels
    if labels is not None and labels != 'none':
        if labels == 'T4':
            texts = [f"{v:.2e}" for v in T_agg]
        elif labels == 'temp':
            texts = [f"{max(v,0.0)**0.25:.2f}" for v in T_agg]
        elif labels == 'heat':
            Pden = elements['heat_in'] / (elements['area'] + 1e-12)
            # Aggregate by max for visibility
            P_agg = np.full(len(C_unique), -np.inf)
            for i, u in enumerate(inv):
                if Pden[i] > P_agg[u]: P_agg[u] = Pden[i]
            P_agg[P_agg == -np.inf] = 0.0
            texts = [f"{v:.2e}" for v in P_agg]
        else:
            texts = None
        if texts is not None:
            for (x, y), t in zip(C_unique, texts):
                ax.text(x, y, t, fontsize=8, color='white', ha='center', va='center', zorder=3,
                        bbox=dict(boxstyle='round,pad=0.2', fc='black', ec='none', alpha=0.4))

    ax.set_aspect('equal', adjustable='box')
    ax.set_title('Radiosity Heat Transfer Simulation', fontsize=16)
    ax.set_xlabel('X coordinate', fontsize=12)
    ax.set_ylabel('Y coordinate', fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.3)
    plt.show()

def flip_normals_toward_centroid(elements):
    """Flip per-element normals so they face toward the centroid of all centers. Returns number flipped."""
    c = elements['center'].mean(axis=0)
    v = c - elements['center']
    d = np.einsum('ij,ij->i', elements['normal'], v)
    mask = d < 0
    elements['normal'][mask] *= -1.0
    return int(mask.sum())


def visualize_matrix(F, cutoff=1e-8):
    """Show coupling matrix in linear and logarithmic scales in a separate figure."""
    # Linear scale
    fig1, ax1 = plt.subplots(figsize=(6, 5))
    im1 = ax1.imshow(F, interpolation='nearest', cmap='viridis', origin='lower')
    ax1.set_title('View Factor Matrix (linear)')
    fig1.colorbar(im1, ax=ax1)

    # Log scale with cutoff
    F_clip = np.clip(F, cutoff, None)
    fig2, ax2 = plt.subplots(figsize=(6, 5))
    im2 = ax2.imshow(F_clip, interpolation='nearest', cmap='viridis', origin='lower', norm=LogNorm(vmin=cutoff, vmax=F_clip.max()))
    ax2.set_title(f'View Factor Matrix (log, cutoff={cutoff:g})')
    fig2.colorbar(im2, ax=ax2)
    plt.show()

# --- Main Execution Block ---

if __name__ == "__main__":
    # run like this:
    #   python Radiosity.py --labels temp --show-matrix
    parser = argparse.ArgumentParser(description='2D Radiosity solver')
    parser.add_argument('--show-matrix', action='store_true', help='Show view factor matrix (linear and log)')
    parser.add_argument('--log-cutoff', type=float, default=1e-8, help='Cutoff for log-scale matrix view')
    parser.add_argument('--labels', type=str, default='none', choices=['none','T4','temp','heat'], help='Annotate elements with values')
    parser.add_argument('--width', type=float, default=1.0, help='Default width (z) for surfaces when not specified')
    parser.add_argument('--no-visibility', action='store_true', help='Do not apply visibility mask (debug)')
    parser.add_argument('--flip-inward', action='store_true', help='Flip normals to face scene centroid (enclosure)')
    args = parser.parse_args()

    # Define the scene geometry: A simple box
    # Surface tuple formats:
    #   1) ((x1,y1),(x2,y2), n, eps, heat)
    #   2) ((x1,y1),(x2,y2), n, eps, heat, width)
    #   3) ((x1,y1),(x2,y2), n, eps_front, eps_back, heat_in)
    #   4) ((x1,y1),(x2,y2), n, eps_front, eps_back, heat_inck, width)
    # Notes:
    #  - Two-sided (3/4) creates front/back sides by duplicating segments with opposite normals.
    #  - Heat values are total per surface; they are split uniformly among its n elements.
    #  - width is the z-extent; area = length * width.
    scene_surfaces = [
        # ((x1,y1),(x2,y2), n, eps_front, eps_back, heat_in, width)
        #((-1.0, -1.0), (-1.0,  1.0), 10, 0.8, 0.0 ),                    # Left wall (passive)
        #((-1.0, -1.0), ( 1.0, -1.0), 10, 0.8, 0.0 ),                    # Bottom wall (passive)
        ((-1.0,  1.0), ( 1.0,  1.0), 10, 0.8, 0.8, 0.0,   0.1 ),      # Top wall (two-sided, passive)
        (( 1.0,  1.0), ( 1.0, -1.0), 10, 0.8, 0.8, 100.0, 0.1 )      # Right wall (two-sided, width=1)
    ]
    
    print("1. Discretizing surfaces...")
    elements_data = discretize_surfaces(scene_surfaces, default_width=args.width)
    num_elements = len(elements_data['length'])
    print(f"   ...created {num_elements} elements.")

    if args.flip_inward:
        flipped = flip_normals_toward_centroid(elements_data)
        print(f"   [DEBUG] Flipped {flipped} normals toward centroid.")

    print("2. Computing view factor matrix...")
    view_factor_matrix = compute_view_factors(elements_data, apply_visibility=not args.no_visibility)

    print(f"   ...computed ({num_elements}x{num_elements}) matrix.")
    print(f"   [DEBUG] F summary: min={view_factor_matrix.min():.3e} max={view_factor_matrix.max():.3e} nnz={np.count_nonzero(view_factor_matrix)}")

    print("3. Solving the radiosity system...")
    radiosity_vector = solve_radiosity_system(view_factor_matrix, elements_data)
    print(f"   ...solved for radiosities.")

    print("4. Calculating final temperatures...")
    temperature_vector_T4 = calculate_temperatures(radiosity_vector, view_factor_matrix, elements_data)
    print(f"   ...calculated T^4 values. Min: {temperature_vector_T4.min():.2f}, Max: {temperature_vector_T4.max():.2f}")

    if args.show_matrix:
        print("5. Visualizing matrix...")
        visualize_matrix(view_factor_matrix, cutoff=args.log_cutoff)

    print("6. Visualizing results...")
    visualize_results(scene_surfaces, elements_data, temperature_vector_T4, labels=args.labels)