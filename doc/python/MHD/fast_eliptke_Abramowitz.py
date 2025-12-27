import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ellipk, ellipe


def fast_ellip_ke(m):
    """
    Polynomial approximation for Complete Elliptic Integrals K(m) and E(m).
    Source: Abramowitz and Stegun, formulas 17.3.34 and 17.3.36.
    Error < 2e-8.
    """
    m1 = 1.0 - m
    # Avoid log(0) at singularity
    m1 = np.where(m1 < 1e-12, 1e-12, m1) 

    # Coefficients for K
    a0, a1, a2, a3, a4 = 1.38629436112, 0.09666344259, 0.03590092383, 0.03742563713, 0.01451196212
    b0, b1, b2, b3, b4 = 0.5, 0.12498593597, 0.06880248576, 0.03328355346, 0.00441787012

    # Coefficients for E
    c1, c2, c3, c4 = 0.44325141463, 0.06260601220, 0.04757383546, 0.01736506451
    d1, d2, d3, d4 = 0.24998368310, 0.09200180037, 0.04069697526, 0.00526449639

    log_m1 = np.log(1.0/m1)

    # Horner's method for polynomial eval
    poly_K1 = a0 + m1*(a1 + m1*(a2 + m1*(a3 + m1*a4)))
    poly_K2 = b0 + m1*(b1 + m1*(b2 + m1*(b3 + m1*b4)))
    K = poly_K1 + poly_K2 * log_m1

    poly_E1 = 1.0 + m1*(c1 + m1*(c2 + m1*(c3 + m1*c4)))
    poly_E2 = m1*(d1 + m1*(d2 + m1*(d3 + m1*d4)))
    E = poly_E1 + poly_E2 * log_m1

    return K, E


def field_loop_rz_fast(a, z0, I, r_grid, z_grid, eps=1e-12):
    """Thin-loop Br,Bz using fast elliptic KE approximation."""
    r = np.asarray(r_grid, dtype=float)
    z = np.asarray(z_grid, dtype=float)
    z_rel = z - z0
    r_safe = np.where(r < eps, eps, r)

    denom = (a + r_safe) ** 2 + z_rel ** 2
    m = (4.0 * a * r_safe) / denom
    K, E = fast_ellip_ke(m)

    factor = 4e-7 * np.pi * I / np.sqrt(denom)  # MU0/(2pi)=2e-7*2pi -> 4e-7*pi
    alpha2 = a * a + r_safe * r_safe + z_rel * z_rel
    beta2 = (a - r_safe) ** 2 + z_rel * z_rel

    Br = factor * (z_rel / r_safe) * (-K + (alpha2 / beta2) * E)
    Bz = factor * (K + ((a * a - r_safe * r_safe - z_rel * z_rel) / beta2) * E)
    Br = np.where(r < eps, 0.0, Br)
    return Br, Bz


def field_loop_rz_ref(a, z0, I, r_grid, z_grid, eps=1e-12):
    """Reference thin-loop Br,Bz using scipy ellipk/ellipe (mirrors inductance_core)."""
    r = np.asarray(r_grid, dtype=float)
    z = np.asarray(z_grid, dtype=float)
    z_rel = z - z0
    r_safe = np.where(r < eps, eps, r)

    denom = (a + r_safe) ** 2 + z_rel ** 2
    m = (4.0 * a * r_safe) / denom
    K = ellipk(m)
    E = ellipe(m)

    factor = 4e-7 * np.pi * I / np.sqrt(denom)
    alpha2 = a * a + r_safe * r_safe + z_rel * z_rel
    beta2 = (a - r_safe) ** 2 + z_rel * z_rel

    Br = factor * (z_rel / r_safe) * (-K + (alpha2 / beta2) * E)
    Bz = factor * (K + ((a * a - r_safe * r_safe - z_rel * z_rel) / beta2) * E)
    Br = np.where(r < eps, 0.0, Br)
    return Br, Bz


def demo_compare():
    """Compare fast vs reference loop field: 2D panels + 1D cuts."""
    a = 1.0
    z0 = 0.0
    I = 1.0

    r_max = 2.5
    z_min, z_max = -2.0, 2.0
    Nr, Nz = 200, 220
    r = np.linspace(-r_max, r_max, Nr)
    z = np.linspace(z_min, z_max, Nz)
    Zg, Rg = np.meshgrid(z, r)  # (Nr, Nz)
    Rabs = np.abs(Rg)

    Br_f, Bz_f = field_loop_rz_fast(a, z0, I, Rabs, Zg)
    Br_r, Bz_r = field_loop_rz_ref(a, z0, I, Rabs, Zg)

    # Magnitudes
    Bmag_f = np.sqrt(Br_f * Br_f + Bz_f * Bz_f)
    Bmag_r = np.sqrt(Br_r * Br_r + Bz_r * Bz_r)
    diff = Bmag_f - Bmag_r

    vmin = 0.0
    vmax = np.percentile(Bmag_r, 99.5)
    dlim = np.percentile(np.abs(diff), 99.5)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)
    im0 = axes[0].imshow(Bmag_f, extent=(z_min, z_max, -r_max, r_max), origin="lower", aspect="auto", vmin=vmin, vmax=vmax)
    axes[0].set_title("Fast B magnitude")
    im1 = axes[1].imshow(Bmag_r, extent=(z_min, z_max, -r_max, r_max), origin="lower", aspect="auto", vmin=vmin, vmax=vmax)
    axes[1].set_title("Reference B magnitude")
    im2 = axes[2].imshow(diff, extent=(z_min, z_max, -r_max, r_max), origin="lower", aspect="auto", cmap="coolwarm", vmin=-dlim, vmax=dlim)
    axes[2].set_title("Fast - Ref")
    for ax in axes:
        ax.set_xlabel("z [m]")
        ax.set_ylabel("r [m]")
        ax.axvline(z0, color="k", lw=0.5, ls="--")
        ax.axhline(0.0, color="k", lw=0.5, ls="--")
    fig.colorbar(im0, ax=axes[0], shrink=0.8)
    fig.colorbar(im1, ax=axes[1], shrink=0.8)
    fig.colorbar(im2, ax=axes[2], shrink=0.8, label="Δ|B|")
    fig.suptitle("Loop field: fast elliptic KE vs scipy (a=1, I=1)")
    plt.tight_layout()

    # 1D cuts
    idx_r0 = np.argmin(np.abs(r))
    idx_z0 = np.argmin(np.abs(z - z0))

    cuts = [
        ("On-axis Bz (r=0)", z, Bz_f[idx_r0, :], Bz_r[idx_r0, :]),
        ("Radial Br (z=0)", r, Br_f[:, idx_z0], Br_r[:, idx_z0]),
        ("Radial Bz (z=0)", r, Bz_f[:, idx_z0], Bz_r[:, idx_z0]),
        (f"Axial through r=a ({a})", z, Bz_f[np.argmin(np.abs(r - a)), :], Bz_r[np.argmin(np.abs(r - a)), :]),
    ]

    fig1, axs = plt.subplots(2, 2, figsize=(12, 8))
    axs = axs.ravel()
    for ax, (title, x, f_fast, f_ref) in zip(axs, cuts):
        ax.plot(x, f_ref, "k-", lw=1.5, label="Ref")
        ax.plot(x, f_fast, "r--", lw=1.0, label="Fast")
        ax.plot(x, f_fast - f_ref, "b:", lw=1.0, label="Fast-Ref")
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        ax.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    demo_compare()