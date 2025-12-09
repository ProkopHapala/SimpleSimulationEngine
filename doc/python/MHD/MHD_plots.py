import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb

from inductance_core import field_loop_rz


def plot_coil_geometry( types, R0,  z0, I0, R1, z1, I_final, bg_mode="hsv", title="Coil geometry (initial vs final)" ):
    """Plot initial and final coil positions with current labels in (z,r) plane.

    Parameters
    ----------
    types : array_like of str
        Coil type labels (e.g. "SC", "CAGE", "PLASMA"). Used only for colors.
    R0, z0 : array_like of float
        Initial radii/positions of coils.
    R1, z1 : array_like of float
        Final radii/positions of coils.
    I0 : array_like of float
        Initial currents.
    I_final : array_like of float
        Final currents from solver (same order as types/R0/z0).
    title : str, optional
        Figure title.
    bg_mode : {"hsv", "mag"}, optional
        Background coloring mode: HSV encodes direction and magnitude,
        MAG encodes magnitude only using a colormap.
    """
    types = np.asarray(types)
    R0 = np.asarray(R0, dtype=float)
    z0 = np.asarray(z0, dtype=float)
    R1 = np.asarray(R1, dtype=float)
    z1 = np.asarray(z1, dtype=float)
    I0 = np.asarray(I0, dtype=float)
    I_final = np.asarray(I_final, dtype=float)

    fig, axes = plt.subplots(1, 2, figsize=(10, 10), sharey=True)
    states = [
        (R0, z0, I0, "t=0"),
        (R1, z1, I_final, "t=final"),
    ]

    # Determine global bounds for margins
    all_R = np.concatenate([R0, R1])
    all_Z = np.concatenate([z0, z1])
    r_span = max(1.0, np.max(np.abs(all_R)) - np.min(np.abs(all_R)))
    r_pad = 0.1 * r_span
    r_extent = np.max(np.abs(all_R)) + r_pad
    r_min_plot, r_max_plot = -r_extent, r_extent
    z_min, z_max = all_Z.min(), all_Z.max()
    dz = 0.1 * max(1.0, z_max - z_min)

    Nr, Nz = 120, 240
    r_grid = np.linspace(r_min_plot, r_max_plot, Nr)
    z_grid = np.linspace(z_min - dz, z_max + dz, Nz)
    Zg, Rg = np.meshgrid(z_grid, r_grid)

    def make_field_and_rgb(Rs, Zs, Is):
        Br = np.zeros_like(Rg)
        Bz = np.zeros_like(Rg)
        for a, zc, I in zip(Rs, Zs, Is):
            Br_i_abs, Bz_i = field_loop_rz(a, zc, I, np.abs(Rg), Zg)
            Br += np.sign(Rg) * Br_i_abs
            Bz += Bz_i
        Bmag = np.sqrt(Br * Br + Bz * Bz)
        B_ref = 0.01 * np.max(Bmag)  # damp dominant field to reveal weaker structures
        if B_ref <= 0.0:
            B_ref = 1.0
        Bmag_norm = np.clip(Bmag / B_ref, 0.0, 1.0)
        phi = np.arctan2(Bz, Br)  # direction in (r,z)
        if bg_mode == "hsv":
            hue = (phi + np.pi) / (2 * np.pi)
            sat = np.ones_like(hue)
            val = Bmag_norm
            hsv = np.stack([hue, sat, val], axis=-1)
            rgb = hsv_to_rgb(hsv)
        else:
            cmap = plt.cm.magma
            rgb = cmap(Bmag_norm)[..., :3]
        return Br, Bz, rgb

    Br_init, Bz_init, rgb_initial = make_field_and_rgb(R0, z0, I0)
    Br_final, Bz_final, rgb_final = make_field_and_rgb(R1, z1, I_final)
    field_states = [(Br_init, Bz_init, rgb_initial), (Br_final, Bz_final, rgb_final)]

    colors={'CAGE':'g','PLASMA':'r','SC':'k'}

    for ax, (Rs, Zs, Is, subt), (Br_s, Bz_s, rgb) in zip(axes, states, field_states):
        extent = (z_min - dz, z_max + dz, r_min_plot, r_max_plot)
        ax.imshow(  rgb, extent=extent,origin="lower",aspect="equal",alpha=0.8,)
        ax.streamplot( z_grid, r_grid,Bz_s, Br_s, color="gray",linewidth=0.4,density=1.8,arrowsize=0.6,)
        for typ, r, z, I in zip(types, Rs, Zs, Is):
            I_MA = I / 1e6
            color = colors.get(typ, "k")
            for sign in (+1, -1):
                r_pos = sign * r
                ax.plot(z, r_pos, marker="o", markersize=4, color=color)
            ax.text(z, r, f"{I_MA:.6f} MA", ha="center", va="bottom", fontsize=8)

        ax.set_title(subt)
        ax.set_xlabel("z [m]")
        ax.set_xlim(z_min - dz, z_max + dz)
        ax.set_ylim(r_min_plot, r_max_plot)
        ax.set_aspect("equal", adjustable="box")
        ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.5)

    axes[0].set_ylabel("r [m]")
    fig.suptitle(title)
    plt.tight_layout()
    return fig, axes


def compute_axial_Bz_profile(Rs, Zs, Is, z_min, z_max, n_samples=400):
    """Compute Bz(r=0,z) profile from a set of loops.

    Returns z_axis, Bz_axis.
    """
    Rs = np.asarray(Rs, dtype=float); Zs = np.asarray(Zs, dtype=float); Is = np.asarray(Is, dtype=float)
    z_axis = np.linspace(z_min, z_max, n_samples)
    r_axis = np.zeros_like(z_axis)
    Br_ax = np.zeros_like(z_axis); Bz_ax = np.zeros_like(z_axis)
    for a, zc, I in zip(Rs, Zs, Is):
        Rg, Zg = np.meshgrid(r_axis, z_axis)
        Br_i, Bz_i = field_loop_rz(a, zc, I, np.abs(Rg), Zg)
        Br_ax += Br_i[:, 0]
        Bz_ax += Bz_i[:, 0]
    return z_axis, Bz_ax


def compute_radial_B_profile(Rs, Zs, Is, z0, r_max, n_samples=400):
    """Compute Br, Bz, |B| as functions of radius at fixed z0.

    Returns r, Br(r), Bz(r), |B|(r).
    """
    Rs = np.asarray(Rs, dtype=float); Zs = np.asarray(Zs, dtype=float); Is = np.asarray(Is, dtype=float)
    r = np.linspace(0.0, r_max, n_samples)
    z = np.full_like(r, float(z0))
    Br_rad = np.zeros_like(r); Bz_rad = np.zeros_like(r)
    for a, zc, I in zip(Rs, Zs, Is):
        Rg, Zg = np.meshgrid(r, z)
        Br_i, Bz_i = field_loop_rz(a, zc, I, Rg, Zg)
        Br_rad += Br_i[0, :]
        Bz_rad += Bz_i[0, :]
    Bmag_rad = np.sqrt(Br_rad * Br_rad + Bz_rad * Bz_rad)
    return r, Br_rad, Bz_rad, Bmag_rad


def overlay_B_profiles_on_axes(ax, z0, Rs, Zs, Is, frac_height=0.25, frac_width=0.25):
    """Overlay axial and radial B profiles onto an existing (z,r) axes.

    This is a thin wrapper around :func:`overlay_axial_B_profile_on_axes` and
    :func:`overlay_radial_B_profile_on_axes`.
    """
    (z_axis, Bz_ax)                   = overlay_axial_B_profile_on_axes (ax,     Rs, Zs, Is, frac_height=frac_height)
    (r_rad, Br_rad, Bz_rad, Bmag_rad) = overlay_radial_B_profile_on_axes(ax, z0, Rs, Zs, Is, frac_width=frac_width)
    return (z_axis, Bz_ax), (r_rad, Br_rad, Bz_rad, Bmag_rad)


def overlay_axial_B_profile_on_axes(ax, Rs, Zs, Is, frac_height=0.25):
    """Overlay axial Bz(r=0,z) profile onto an existing (z,r) axes.

    The curve is normalized and scaled to occupy a fraction of the vertical
    span; a zero-reference line at r=0 is drawn for comparison.
    """
    z_min, z_max = ax.get_xlim(); r_min, r_max = ax.get_ylim()
    z_axis, Bz_ax = compute_axial_Bz_profile(Rs, Zs, Is, z_min, z_max)
    Bz_norm = Bz_ax / (np.max(np.abs(Bz_ax)) + 1e-30)
    r_scale = frac_height * (r_max - r_min)
    r_curve = Bz_norm * r_scale
    ax.axhline(0.0, color="k", linestyle=":", linewidth=0.7)
    ax.plot(z_axis, r_curve, color="k", linewidth=1.0)
    return z_axis, Bz_ax


def overlay_radial_B_profile_on_axes(ax, z0, Rs, Zs, Is, frac_width=0.25):
    """Overlay radial Br/Bz/|B| profile at fixed z0 onto an existing (z,r) axes.

    The three components are drawn as horizontal deviations around z0, scaled
    to occupy a fraction of the horizontal span; a zero-reference line at
    z=z0 is drawn for comparison.
    """
    z_min, z_max = ax.get_xlim(); r_min, r_max = ax.get_ylim()
    r_max_rad = max(abs(r_min), abs(r_max))
    r_rad, Br_rad, Bz_rad, Bmag_rad = compute_radial_B_profile(Rs, Zs, Is, z0, r_max_rad)
    Bmax = np.max(Bmag_rad) + 1e-30
    z_scale = frac_width * (z_max - z_min)
    z_Br = z0 + (Br_rad / Bmax) * z_scale
    z_Bz = z0 + (Bz_rad / Bmax) * z_scale
    z_Bm = z0 + (Bmag_rad / Bmax) * z_scale
    ax.axvline(z0, color="k", linestyle=":", linewidth=0.7)
    ax.plot(z_Br, r_rad, color="b", linewidth=1.0)
    ax.plot(z_Bz, r_rad, color="r", linewidth=1.0)
    ax.plot(z_Bm, r_rad, color="k", linewidth=1.0)
    return r_rad, Br_rad, Bz_rad, Bmag_rad


def plot_B_profiles(Rs, Zs, Is, z0, axes=None, show_separate=True):
    """High-level helper to draw B profiles.

    - If ``axes`` is a 2-panel geometry axes array (as returned by plot_coil_geometry),
      overlays profiles on ``axes[1]``.
    - If ``axes`` is None, creates a fresh (z,r) axes for the overlay.
    - Optionally creates a separate figure with two 1D plots (axial and radial).

    Returns
    -------
    (z_axis, Bz_ax), (r_rad, Br_rad, Bz_rad, Bmag_rad)
        Raw profile data used for the plots.
    """
    Rs = np.asarray(Rs, dtype=float); Zs = np.asarray(Zs, dtype=float); Is = np.asarray(Is, dtype=float)

    if axes is None:
        fig, ax_final = plt.subplots(figsize=(6, 6))
        ax_final.set_xlabel("z [m]"); ax_final.set_ylabel("r [m]")
        z_min_plot, z_max_plot = float(np.min(Zs)), float(np.max(Zs))
        r_extent = float(np.max(np.abs(Rs)))
        r_min_plot, r_max_plot = -r_extent, r_extent
        ax_final.set_xlim(z_min_plot, z_max_plot); ax_final.set_ylim(r_min_plot, r_max_plot)
    else:
        # axes may be a 2-panel array (from plot_coil_geometry) or a single Axes
        try:
            # treat as sequence and use final panel
            ax_final = axes[1]
        except (TypeError, IndexError):
            ax_final = axes

    (z_axis, Bz_ax), (r_rad, Br_rad, Bz_rad, Bmag_rad) = overlay_B_profiles_on_axes(ax_final, float(z0), Rs, Zs, Is)

    if show_separate:
        fig_p, ax_p = plt.subplots(1, 2, figsize=(10, 4))
        # Axial profile
        ax_p[0].plot(z_axis, Bz_ax, color="k")
        ax_p[0].axhline(0.0, color="gray", linestyle=":", linewidth=0.7)
        ax_p[0].set_xlabel("z [m]"); ax_p[0].set_ylabel("Bz(r=0) [T]"); ax_p[0].set_title("Axial Bz profile")
        ax_p[0].grid(True)
        # Radial profile
        ax_p[1].plot(r_rad, Br_rad, color="b", label="Br")
        ax_p[1].plot(r_rad, Bz_rad, color="r", label="Bz")
        ax_p[1].plot(r_rad, Bmag_rad, color="k", label="|B|")
        ax_p[1].axhline(0.0, color="gray", linestyle=":", linewidth=0.7)
        ax_p[1].set_xlabel("r [m]"); ax_p[1].set_ylabel("B components [T]"); ax_p[1].set_title("Radial B profile at z0")
        ax_p[1].grid(True); ax_p[1].legend()
        plt.tight_layout()

    return (z_axis, Bz_ax), (r_rad, Br_rad, Bz_rad, Bmag_rad)
