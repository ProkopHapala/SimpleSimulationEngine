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


def generate_sc_seed_coil(sc_r, sc_z, sc_current=1.0e6, coil_type="SC"):
    """Return a single seed coil definition line (type configurable, default "SC")."""
    t = f"{coil_type:6s}"
    return f"{t}  {sc_r:.3f}  {sc_z:.3f}   {sc_r:.3f}  {sc_z:.3f}   {sc_current:.1e}"


def generate_parabolic_rings(n_rings=12, r_throat=1.0, r_exit=3.0, z_start=0.0, z_end=2.0, coil_type="CAGE"):
    """Generate rings following a parabolic wall (generic coil type).

    By default this is used for CAGE coils of a parabolic nozzle.
    """
    A = (z_end - z_start) / ((r_exit - r_throat) ** 2)  # parabola parameter
    r_step = (r_exit - r_throat) / max(1, n_rings - 1)
    t = f"{coil_type:6s}"
    lines = []
    for i in range(n_rings):
        r_curr = r_throat + i * r_step
        z_curr = A * r_curr ** 2 + z_start
        lines.append(f"{t}  {r_curr:.3f}  {z_curr:.3f}   {r_curr:.3f}  {z_curr:.3f}   0.0")
    return "\n".join(lines)


def generate_parabolic_cage_rings(n_rings=12, r_throat=1.0, r_exit=3.0, z_start=0.0, z_end=2.0):
    """Backward-compatible wrapper: CAGE rings along a parabolic nozzle wall."""
    return generate_parabolic_rings(n_rings=n_rings, r_throat=r_throat, r_exit=r_exit, z_start=z_start, z_end=z_end, coil_type="CAGE")


def generate_parabolic_plasma_armature(plasma_r_start=0.2, plasma_r_end=1.0, plasma_z_start=0.0, plasma_z_end=2.0, coil_type="PLASMA"):
    """Generate only a single armature line between two (r,z) endpoints.

    Historically used for the PLASMA armature of a parabolic nozzle; coil_type
    allows reuse for other types.
    """
    t = f"{coil_type:6s}"
    return (
        f"{t}  {plasma_r_start:.3f}  {plasma_z_start:.3f}   "
        f"{plasma_r_end:.3f}  {plasma_z_end:.3f}   0.0"
    )

def generate_spherical_rings(n_rings=8, r0=0.01, r1=0.4, z0=0.0, coil_type="PLASMA"):
    """Generate loop lines approximating a spherical shell in (r,z).

    coil_type allows using the same geometry for PLASMA, CAGE, etc.
    """
    import numpy as _np

    lines = []
    theta_min = 0.5 * _np.pi / max(4.0, float(n_rings))
    theta_max = _np.pi - theta_min
    thetas = _np.linspace(theta_min, theta_max, n_rings)
    t = f"{coil_type:6s}"

    for th in thetas:
        sr0 =      r0 * _np.sin(th)
        zz0 = z0 + r0 * _np.cos(th)
        sr1 =      r1 * _np.sin(th)
        zz1 = z0 + r1 * _np.cos(th)
        lines.append(f"{t}  {sr0:.6f}  {zz0:.6f}   {sr1:.6f}  {zz1:.6f}   {0.0:.6e}")

    return "\n".join(lines)


def generate_spherical_plasma_loops(n_plasma=8, r0=0.01, r1=0.4, z0=0.0):
    """Backward-compatible wrapper for PLASMA spherical shell loops."""
    return generate_spherical_rings(n_rings=n_plasma, r0=r0, r1=r1, z0=z0, coil_type="PLASMA")


def generate_disk_rings(n_rings=8, r_inner=0.1, r_outer=1.0, z0=0.0, coil_type="PLASMA"):
    """Generate concentric rings in a disk at fixed z0."""
    rs = np.linspace(r_inner, r_outer, n_rings)
    t = f"{coil_type:6s}"
    lines = []
    for r in rs:
        lines.append(f"{t}  {r:.6f}  {z0:.6f}   {r:.6f}  {z0:.6f}   {0.0:.6e}")
    return "\n".join(lines)


def generate_tube_rings(n_rings=8, radius=1.0, z_start=0.0, z_end=2.0, coil_type="PLASMA"):
    """Generate rings forming a straight tube between z_start and z_end."""
    zs = np.linspace(z_start, z_end, n_rings)
    t = f"{coil_type:6s}"
    lines = []
    for z in zs:
        lines.append(f"{t}  {radius:.6f}  {z:.6f}   {radius:.6f}  {z:.6f}   {0.0:.6e}")
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

    Br = factor * z_rel / rho_safe * (  -K + ((a * a + rho_safe * rho_safe + z_rel * z_rel) / common) * E )
    Bz = factor * (K + ((a * a - rho_safe * rho_safe - z_rel * z_rel) / common) * E)

    # Set Br to zero on axis where symmetry enforces it
    Br = np.where(rho < eps, 0.0, Br)
    return Br, Bz


def field_dipole_rz(m, z_m, r_grid, z_grid, eps=1e-12):
    """Magnetic field (Br,Bz) of an axial point dipole in the (r,z) meridional plane.

    Parameters
    ----------
    m : float
        Dipole moment magnitude (A m^2), oriented along +z.
    z_m : float
        Axial position of the dipole (on the symmetry axis r=0).
    r_grid, z_grid : ndarray
        Evaluation grid in cylindrical coordinates.
    eps : float, optional
        Small regularization radius to avoid singularity at the dipole location.

    Returns
    -------
    Br, Bz : ndarray
        Radial and axial components of the dipole field.
    """
    r = np.asarray(r_grid, dtype=float)
    z = np.asarray(z_grid, dtype=float)
    rho = r
    z_rel = z - z_m

    # Regularized distance from dipole
    rho_safe = np.where(rho < eps, eps, rho)
    R2 = rho_safe * rho_safe + z_rel * z_rel
    R = np.sqrt(R2)
    R5 = R2 * R2 * R + eps

    pref = MU0 * m / (4.0 * np.pi)

    Br = pref * 3.0 * rho_safe * z_rel / R5
    Bz = pref * (2.0 * z_rel * z_rel - rho_safe * rho_safe) / R5

    # On-axis limit: Br must vanish by symmetry
    Br = np.where(rho < eps, 0.0, Br)
    return Br, Bz


def Aphi_dipole_rz(m, z_m, r_grid, z_grid, eps=1e-12):
    """Azimuthal vector potential A_phi of an axial point dipole in the (r,z) plane.

    The dipole moment m is oriented along +z and located at (r=0, z=z_m). In
    Coulomb gauge, the vector potential is

        A = mu0 / (4 pi) * (m x R) / |R|^3 ,

    which for an on-axis dipole reduces to a purely azimuthal component

        A_phi = mu0 * m * r / (4 pi * R^3) ,   R^2 = r^2 + (z - z_m)^2.
    """
    r = np.asarray(r_grid, dtype=float)
    z = np.asarray(z_grid, dtype=float)
    rho = r
    z_rel = z - z_m

    rho_safe = np.where(rho < eps, eps, rho)
    R2 = rho_safe * rho_safe + z_rel * z_rel
    R = np.sqrt(R2)
    R3 = R2 * R + eps

    pref = MU0 * m / (4.0 * np.pi)
    Aphi = pref * rho_safe / R3

    # On axis the azimuthal component must vanish by symmetry
    Aphi = np.where(rho < eps, 0.0, Aphi)
    return Aphi


def Aphi_loop_rz(a, z0, I, r_grid, z_grid, eps=1e-12):
    """Azimuthal vector potential A_phi of a circular loop at (a,z0) with current I.

    This uses the standard elliptic-integral expression in cylindrical
    coordinates. With

        k^2 = 4 a r / ((a + r)^2 + (z - z0)^2),

    the azimuthal component for a loop of radius a and current I is

        A_phi(r,z) = mu0 I / (2 pi k) * sqrt(a / r)
                      * [ (1 - 0.5 k^2) K(k^2) - E(k^2) ],

    where K,E are complete elliptic integrals of the first and second kind.
    The expression is regularized near r=0 and on-axis A_phi is forced to 0
    by symmetry.
    """
    r = np.asarray(r_grid, dtype=float)
    z = np.asarray(z_grid, dtype=float)
    rho = r
    z_rel = z - z0

    # Avoid singularity at r=0: A_phi must vanish on axis, but formula has 1/sqrt(r)
    rho_safe = np.where(rho < eps, eps, rho)

    # Elliptic parameter m = k^2
    denom = (a + rho_safe) ** 2 + z_rel * z_rel
    m = (4.0 * a * rho_safe) / denom
    K = ellipk(m)
    E = ellipe(m)

    # A_phi = (mu0 I / pi) * sqrt(a/r) * (1/k) * [ (1 - m/2) K(m) - E(m) ]
    # Algebraically this simplifies to
    # A_phi = (mu0 I / (2 pi r)) * sqrt(denom) * [ (1 - m/2) K(m) - E(m) ]
    prefactor = (MU0 * I) / (2.0 * np.pi * rho_safe) * np.sqrt(denom)
    bracket = (1.0 - 0.5 * m) * K - E
    Aphi = prefactor * bracket

    # Enforce A_phi = 0 on axis by symmetry
    Aphi = np.where(rho < eps, 0.0, Aphi)
    return Aphi


def coil_axis_Bz(a, z0, I, z_eval):
    """Axial Bz of a circular loop at radius a, center z0, on the symmetry axis.

    Parameters
    ----------
    a : float
        Loop radius.
    z0 : float
        Axial coordinate of loop center.
    I : float
        Loop current (A).
    z_eval : float or ndarray
        Axial position(s) where Bz is evaluated (r=0).

    Returns
    -------
    Bz_axis : float or ndarray
        Axial component of B at (r=0, z=z_eval).
    """
    z_eval = np.asarray(z_eval, dtype=float)
    z_rel = z_eval - float(z0)
    denom = (a * a + z_rel * z_rel) ** 1.5
    return MU0 * I * a * a / (2.0 * denom)


def coil_dipole_energy(a, z0, I, m, z_m):
    """Interaction energy between a circular coil and an axial point dipole.

    The dipole is on-axis at z_m with moment m along +z. The energy is

        U = - m * Bz_coil(r=0, z=z_m)

    where Bz_coil is the axial field produced by the loop at the dipole
    location.
    """
    Bz_at_dip = coil_axis_Bz(a, z0, I, z_m)
    return -m * float(Bz_at_dip)
