import numpy as np
from scipy.special import ellipk, ellipe

MU0 = 4 * np.pi * 1e-7


def calc_mutual_inductance(r1, z1, r2, z2):
    """Mutual inductance between two coaxial circular loops.

    Parameters
    ----------
    r1, z1 : float
        Radius and axial position of loop 1.
    r2, z2 : float
        Radius and axial position of loop 2.

    Returns
    -------
    float
        Mutual inductance in Henry.
    """
    if r1 == 0.0 or r2 == 0.0:
        return 0.0

    dz = z2 - z1
    sr = r2 + r1
    k2 = (4.0 * r1 * r2) / (sr * sr + dz * dz)
    k = np.sqrt(k2)

    K = ellipk(k2)
    E = ellipe(k2)

    return 4e-7 * np.pi * np.sqrt(r1 * r2) * ((2.0 / k - k) * K - (2.0 / k) * E)


def calc_self_inductance(R, r_wire=0.01):
    """Self-inductance of a single circular loop of radius R.

    Parameters
    ----------
    R : float
        Loop radius.
    r_wire : float, optional
        Wire radius used in the logarithmic core-regularization.

    Returns
    -------
    float
        Self inductance in Henry.
    """
    return 4e-7 * np.pi * R * (np.log(8.0 * R / r_wire) - 1.75)


def build_inductance_matrix(rs, zs, self_r_wire=0.01):
    """Build full inductance matrix K for a set of coaxial loops.

    Parameters
    ----------
    rs, zs : array_like
        Radii and axial positions of all loops (same length N).
    self_r_wire : float, optional
        Wire radius used for self-inductance.

    Returns
    -------
    K : ndarray, shape (N, N)
        Inductance matrix such that Phi = K @ I.
    """
    rs = np.asarray(rs, dtype=float)
    zs = np.asarray(zs, dtype=float)
    N = rs.size
    K = np.zeros((N, N), dtype=float)

    for i in range(N):
        for j in range(N):
            if i == j:
                K[i, j] = calc_self_inductance(rs[i], r_wire=self_r_wire)
            else:
                K[i, j] = calc_mutual_inductance(rs[i], zs[i], rs[j], zs[j])
    return K


def compute_flux(K, I):
    """Compute flux vector Phi = K @ I."""
    K = np.asarray(K, dtype=float)
    I = np.asarray(I, dtype=float)
    return K @ I


def init_flux(K0, I0, flux_override=None):
    """Compute initial flux Phi0 from K0 and I0 with optional overrides.

    Parameters
    ----------
    K0 : ndarray, shape (N, N)
        Initial inductance matrix.
    I0 : ndarray, shape (N,)
        Initial current vector.
    flux_override : dict or None, optional
        Mapping index -> value for Phi0[i] overrides (e.g. diamagnetic loops).

    Returns
    -------
    Phi0 : ndarray, shape (N,)
        Initial flux vector.
    """
    Phi0 = compute_flux(K0, I0)
    if flux_override is not None:
        for i, val in flux_override.items():
            Phi0[int(i)] = float(val)
    return Phi0


def solve_flux_conserving(K1, Phi0, fixed_mask=None, I_fixed=None):
    """Solve for currents at t1 given flux conservation constraints.

    This solves K1 @ I1 = Phi0 for the subset of coils that conserve flux,
    optionally treating some coils as fixed-current sources.

    Parameters
    ----------
    K1 : ndarray, shape (N, N)
        Inductance matrix at t1.
    Phi0 : ndarray, shape (N,)
        Target flux for each loop (typically Phi(t0)).
    fixed_mask : array_like of bool, optional
        Boolean mask of length N; True where current is fixed.
        If None, all currents are treated as unknown.
    I_fixed : array_like, optional
        Values of fixed currents (same length N, entries for non-fixed are ignored).

    Returns
    -------
    I1 : ndarray, shape (N,)
        Solved current vector at t1.
    """
    K1 = np.asarray(K1, dtype=float)
    Phi0 = np.asarray(Phi0, dtype=float)
    N = Phi0.size

    if fixed_mask is None:
        return np.linalg.solve(K1, Phi0)

    fixed_mask = np.asarray(fixed_mask, dtype=bool)
    if I_fixed is None:
        raise ValueError("I_fixed must be provided when fixed_mask is used")

    I_fixed = np.asarray(I_fixed, dtype=float)

    idx_u = np.where(~fixed_mask)[0]
    idx_f = np.where(fixed_mask)[0]

    if idx_u.size == 0:
        return I_fixed.copy()

    K_uu = K1[np.ix_(idx_u, idx_u)]
    if idx_f.size > 0:
        K_uf = K1[np.ix_(idx_u, idx_f)]
        RHS = Phi0[idx_u] - K_uf @ I_fixed[idx_f]
    else:
        RHS = Phi0[idx_u]

    I_u = np.linalg.solve(K_uu, RHS)

    I1 = I_fixed.copy()
    I1[idx_u] = I_u
    return I1


def magnetic_energy(K, I):
    """Total magnetic energy 0.5 * I^T K I."""
    K = np.asarray(K, dtype=float)
    I = np.asarray(I, dtype=float)
    return 0.5 * float(I.T @ (K @ I))


def generate_parabolic_nozzle(
    n_rings=12,
    r_throat=1.0,
    r_exit=3.0,
    z_start=0.0,
    z_end=2.0,
    sc_current=1.0e6,
    plasma_r_start=0.2,
    plasma_r_end=None,
    plasma_r_end_factor=0.90,
    sc_r=None,
    sc_z=None,
    plasma_z_start=None,
    plasma_z_end=None,
):
    """Generate a simple parabolic nozzle configuration as coil definitions.

    Returns a multi-line string with lines of the form:
        TYPE  R0  Z0  R1  Z1  I0

    Where:
    - SC     : seed superconducting coil (radius/z can be overridden)
    - CAGE   : a set of nozzle cage rings following a parabolic wall
    - PLASMA : plasma armature moving from throat to exit (r,z path configurable)
    """
    if sc_r is None:
        sc_r = r_throat * 1.5
    if sc_z is None:
        sc_z = z_start

    if plasma_z_start is None:
        pl_z_start = z_start
    else:
        pl_z_start = plasma_z_start

    if plasma_z_end is None:
        pl_z_end = z_end
    else:
        pl_z_end = plasma_z_end

    if plasma_r_end is None:
        pl_r_end = r_exit * plasma_r_end_factor
    else:
        pl_r_end = plasma_r_end

    # Parabola parameter A such that (r_exit, z_end) lies on the wall
    A = (z_end - z_start) / ((r_exit - r_throat) ** 2)
    r_step = (r_exit - r_throat) / max(1, n_rings - 1)

    lines = []

    # Seed SC coil
    lines.append(f"SC      {sc_r:.3f}  {sc_z:.3f}   {sc_r:.3f}  {sc_z:.3f}   {sc_current:.1e}")

    # Cage coils along parabolic wall
    for i in range(n_rings):
        r_curr = r_throat + i * r_step
        z_curr = A * r_curr ** 2 + z_start
        lines.append(f"CAGE    {r_curr:.3f}  {z_curr:.3f}   {r_curr:.3f}  {z_curr:.3f}   0.0")

    # Plasma coil trajectory
    lines.append(
        f"PLASMA  {plasma_r_start:.3f}  {pl_z_start:.3f}   {pl_r_end:.3f}  {pl_z_end:.3f}   0.0"
    )

    return "\n".join(lines)


def generate_spherical_plasma_loops(n_plasma=8, r0=0.01, r1=0.4, z0=0.0):
    """Generate PLASMA loop lines approximating a spherical shell in (r,z).

    The output is a multi-line string with lines of the form

        PLASMA  R0  Z0   R1  Z1   I0

    where (R0,Z0) and (R1,Z1) describe initial and final spherical shells
    centered at z0 with radii r0 and r1, respectively, and I0 = 0.
    """
    import numpy as _np

    lines = []
    theta_min = 0.5 * _np.pi / max(4.0, float(n_plasma))
    theta_max = _np.pi - theta_min
    thetas = _np.linspace(theta_min, theta_max, n_plasma)

    for th in thetas:
        sr0 =      r0 * _np.sin(th)
        zz0 = z0 + r0 * _np.cos(th)
        sr1 =      r1 * _np.sin(th)
        zz1 = z0 + r1 * _np.cos(th)
        lines.append(f"PLASMA  {sr0:.6f}  {zz0:.6f}   {sr1:.6f}  {zz1:.6f}   {0.0:.6e}")

    return "\n".join(lines)


def field_loop_rz(a, z0, I, r_grid, z_grid, eps=1e-12):
    """Biot–Savart field of a circular loop at (a, z0) with current I in (r,z) plane.

    Parameters
    ----------
    a : float
        Loop radius.
    z0 : float
        Axial position of loop center.
    I : float
        Current (A).
    r_grid, z_grid : ndarray
        Meshgrid of evaluation points (same shape).
    eps : float
        Small regularization to avoid divide-by-zero at r=0.

    Returns
    -------
    Br, Bz : ndarray
        Radial and axial components on the grid.
    """
    r = np.asarray(r_grid, dtype=float)
    z = np.asarray(z_grid, dtype=float)
    rho = r
    z_rel = z - z0

    # Avoid division by zero on axis
    rho_safe = np.where(rho < eps, eps, rho)

    k2 = (4.0 * a * rho_safe) / ((a + rho_safe) ** 2 + z_rel ** 2)
    k = np.sqrt(k2)
    K = ellipk(k2)
    E = ellipe(k2)

    denom = np.sqrt((a + rho_safe) ** 2 + z_rel ** 2)
    factor = MU0 * I / (2.0 * np.pi * denom)
    common = (a - rho_safe) ** 2 + z_rel ** 2

    Br = factor * z_rel / rho_safe * (
        -K + ((a * a + rho_safe * rho_safe + z_rel * z_rel) / common) * E
    )
    Bz = factor * (K + ((a * a - rho_safe * rho_safe - z_rel * z_rel) / common) * E)

    # Set Br to zero on axis where symmetry enforces it
    Br = np.where(rho < eps, 0.0, Br)
    return Br, Bz
