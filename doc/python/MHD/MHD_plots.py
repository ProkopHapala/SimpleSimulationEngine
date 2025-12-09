import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb

from inductance_core import field_loop_rz


def plot_coil_geometry(
    types,
    R0,
    z0,
    I0,
    R1,
    z1,
    I_final,
    title="Coil geometry (initial vs final)",
    bg_mode="hsv",
):
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

    for ax, (Rs, Zs, Is, subt), (Br_s, Bz_s, rgb) in zip(axes, states, field_states):
        ax.imshow(
            rgb,
            extent=(z_min - dz, z_max + dz, r_min_plot, r_max_plot),
            origin="lower",
            aspect="equal",
            alpha=0.8,
        )
        ax.streamplot(
            z_grid,
            r_grid,
            Bz_s,
            Br_s,
            color="k",
            linewidth=0.6,
            density=1.2,
            arrowsize=0.7,
        )
        for typ, r, z, I in zip(types, Rs, Zs, Is):
            I_MA = I / 1e6
            if typ == "SC":
                color = "k"
            elif typ == "CAGE":
                color = "g"
            else:
                color = "r"
            radius_draw = 0.02 * max(1.0, r_extent)
            for sign in (+1, -1):
                r_pos = sign * r
                circle = plt.Circle(
                    (z, r_pos),
                    radius_draw,
                    fill=False,
                    color=color,
                    linewidth=1.0,
                )
                ax.add_patch(circle)
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
