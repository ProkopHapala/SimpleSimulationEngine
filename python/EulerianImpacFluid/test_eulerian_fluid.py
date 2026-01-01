import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from EulerianImpacFluid import EulerianImpactFluid

# EOS params (must match OpenCL)
GAMMA_1 = 1.4
PI_1 = 0.0
GAMMA_2 = 3.0
PI_2 = 40.0e9

# Material properties (SI units)
RHO_H2 = 71.0         # kg/m^3 (Liquid Hydrogen)
RHO_URANIUM = 19050.0 # kg/m^3 (Uranium)

"""
##### RUN LIKE THIS: 
python test_eulerian_fluid.py --fields density,pressure,sound,speed --p0 1e9 -t 1e-8 -n 10 -f 50
"""

EPS = 1e-9

# -------------------------------------------------------------------------
def heaviside(phi, dx):
    epsilon = 1.5 * dx
    x = (phi + epsilon) / (2.0 * epsilon)
    x = np.clip(x, 0.0, 1.0)
    return x * x * (3.0 - 2.0 * x)

# Pretty print initial regions (SI units)
# -------------------------------------------------------------------------
def summarize_region(label, rho1, rho2, vx, vy, p0):
    rho_total = rho1 + rho2
    alpha1 = rho1 / (rho_total + EPS) if rho_total > 0 else 0.0
    alpha2 = 1.0 - alpha1
    g1_m1 = GAMMA_1 - 1.0
    g2_m1 = GAMMA_2 - 1.0
    denom = (alpha1 / g1_m1) + (alpha2 / g2_m1) + EPS
    pinf_term = (alpha1 * GAMMA_1 * PI_1 / g1_m1) + (alpha2 * GAMMA_2 * PI_2 / g2_m1)
    internal_e = p0 * denom + pinf_term
    ke = 0.5 * rho_total * (vx * vx + vy * vy)
    E_val = internal_e + ke
    # Sound speeds per phase
    c1 = np.sqrt(GAMMA_1 * (p0 + PI_1) / (rho1 + EPS)) if rho1 > 0 else 0.0
    c2 = np.sqrt(GAMMA_2 * (p0 + PI_2) / (rho2 + EPS)) if rho2 > 0 else 0.0
    # Mixture frozen sound speed (same structure as kernel)
    gamma_avg_num = (alpha1 * GAMMA_1 / g1_m1) + (alpha2 * GAMMA_2 / g2_m1)
    c_mix = np.sqrt(max((gamma_avg_num / denom) * (p0 + pinf_term / max(gamma_avg_num, EPS)) / (rho_total + EPS), 0.0))
    print(
        f"[init-region] {label}: rho1={rho1:.3e} kg/m^3 rho2={rho2:.3e} kg/m^3 "
        f"alpha1={alpha1:.3f} alpha2={alpha2:.3f} "
        f"p0={p0:.3e} Pa c_mix={c_mix:.3e} m/s c1={c1:.3e} c2={c2:.3e} "
        f"vx={vx:.3f} m/s vy={vy:.3f} m/s E={E_val:.3e} J/m^3"
    )

# CLI -------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Eulerian impact fluid demo (PyOpenCL)")
parser.add_argument("-W", "--width",  type=int, default=256)
parser.add_argument("-H", "--height", type=int, default=128)
parser.add_argument("-x", "--dx",     type=float, default=0.01)
# Stable defaults from tested run: dt=1e-8, 10 substeps per frame
parser.add_argument("-t", "--dt",     type=float, default=1e-7)
parser.add_argument("-n", "--steps",  type=int, default=10, help="steps per frame")
parser.add_argument("-r", "--redist", type=int, default=20, help="frames between redistance")
parser.add_argument("-s", "--snap",   type=int, default=1000, help="frames between saved PNGs")
parser.add_argument("-f", "--frames", type=int, default=200)
parser.add_argument("-p", "--prefix", type=str, default="snapshot_frame")
parser.add_argument("--p0", type=float, default=1e9, help="initial target pressure (Pa)")
parser.add_argument("-F", "--fields", type=str, default="density,phi,pressure,energy,speed,sound", help="Comma-separated fields to plot: density,phi,pressure,energy,speed,sound")
args = parser.parse_args()

width, height = args.width, args.height
dx = args.dx
dt = args.dt
num_steps_per_frame = args.steps
redistance_every = args.redist
snapshot_every = args.snap
snapshot_prefix = args.prefix

# Initialize the solver
solver = EulerianImpactFluid(width, height, dx)

print(f"[config] width={width} height={height} dx={dx} dt={dt} "
      f"steps_per_frame={num_steps_per_frame} redist_every={redistance_every} "
      f"snap_every={snapshot_every} prefix='{snapshot_prefix}' fields='{args.fields}' frames={args.frames} p0={args.p0}")
print("[init-summary] Reporting initial region properties (SI units):")
summarize_region("Background H2", RHO_H2, 0.0, 0.0, 0.0, args.p0)
summarize_region("Projectile U", 0.0, RHO_URANIUM, 5000.0, 0.0, args.p0)

# Initial state setup
x = np.linspace(0, width*dx, width)
y = np.linspace(0, height*dx, height)
X, Y = np.meshgrid(x, y)

# phi > 0 is Uranium, phi < 0 is Hydrogen
cx, cy, r = 0.5, height*dx/2, 0.2
phi = r - np.sqrt((X - cx)**2 + (Y - cy)**2)

# alpha2 = H(phi) => Uranium region (phi > 0)
# alpha1 = 1 - H(phi) => Hydrogen region (phi < 0)
a2 = heaviside(phi, dx)
a1 = 1.0 - a2
is_uranium = phi > 0

rho_a1 = RHO_H2 * a1
rho_a2 = RHO_URANIUM * a2
rho_tot = rho_a1 + rho_a2

# Smooth velocity transition
v_smooth = 5000.0 * a2
ru = rho_tot * v_smooth
rv = np.zeros_like(rho_a1)

# Energy: mixture EOS (isobaric p=p0)
g1_m1 = GAMMA_1 - 1.0
g2_m1 = GAMMA_2 - 1.0
denom = (a1 / g1_m1) + (a2 / g2_m1) + EPS
pinf_term = (a1 * GAMMA_1 * PI_1 / g1_m1) + (a2 * GAMMA_2 * PI_2 / g2_m1)

internal_e = args.p0 * denom + pinf_term
ke = 0.5 * (v_smooth**2) * rho_tot
E = internal_e + ke

solver.initialize(rho_a1, rho_a2, ru, rv, E, phi)
    
# Run one step with dt=0 to populate pressure/sound diagnostics at t=0
solver.step(0.0)
solver.queue.finish() # Ensure GPU is done before reading back for t=0 plot

step_counter = 0

field_list = [f.strip().lower() for f in args.fields.split(",") if f.strip()]
extent = [0, width*dx, 0, height*dx]
data = solver.get_data()
density = np.clip(data['rho_a1'] + data['rho_a2'], 1e-6, None)
ru64 = np.nan_to_num(data['ru'], copy=False).astype(np.float64); rv64 = np.nan_to_num(data['rv'], copy=False).astype(np.float64); density64 = density.astype(np.float64)
vel_sq_init = np.clip((ru64*ru64 + rv64*rv64) / (density64 + 1e-6)**2, 0.0, 1e20)
speed = np.sqrt(vel_sq_init)
pressure_init = data.get('p', None)
if pressure_init is None:
    pressure_init = np.ones_like(density)
sound_init = data.get('c', None)
if sound_init is None:
    sound_init = np.ones_like(density)

# Host-side EOS diagnostic to verify initialization
a2_d = heaviside(data['phi'], dx)
a1_d = 1.0 - a2_d
rho_sum = data['rho_a1'] + data['rho_a2']
    
g1_m1 = GAMMA_1 - 1.0
g2_m1 = GAMMA_2 - 1.0
denom_d = (a1_d / g1_m1) + (a2_d / g2_m1) + EPS
pinf_d = (a1_d * GAMMA_1 * PI_1 / g1_m1) + (a2_d * GAMMA_2 * PI_2 / g2_m1)
    
# Kinetic energy using double precision for diagnostic
ru_d = data['ru'].astype(np.float64)
rv_d = data['rv'].astype(np.float64)
rho_d = rho_sum.astype(np.float64)
ke_d = 0.5 * (ru_d**2 + rv_d**2) / (rho_d + EPS)
    
internal_d = data['E'].astype(np.float64) - ke_d
p_cpu = (internal_d - pinf_d) / denom_d
    
print(f"[diag-init] p_target={args.p0:.3e}")
print(f"[diag-init] CPU-EOS p=[{p_cpu.min():.3e},{p_cpu.max():.3e}] gpu-p=[{pressure_init.min():.3e},{pressure_init.max():.3e}]")
print(f"[diag-init] E=[{data['E'].min():.3e},{data['E'].max():.3e}] ke=[{ke_d.min():.3e},{ke_d.max():.3e}] internal=[{internal_d.min():.3e},{internal_d.max():.3e}]")
    
if is_uranium.any():
    print(f"[diag-init] Uranium mask (phi>0): p_cpu=[{p_cpu[is_uranium].min():.3e},{p_cpu[is_uranium].max():.3e}] "
          f"internal=[{internal_d[is_uranium].min():.3e},{internal_d[is_uranium].max():.3e}] pinf_d=[{pinf_d[is_uranium].min():.3e},{pinf_d[is_uranium].max():.3e}]")

plot_data = []
if "density"  in field_list: plot_data.append(("Density", density, "jet"))
if "phi"      in field_list: plot_data.append(("Level-set phi", data['phi'], "RdBu"))
if "pressure" in field_list: plot_data.append(("Pressure p", pressure_init, "plasma"))
if "energy"   in field_list: plot_data.append(("Energy E", data['E'], "magma"))
if "speed"    in field_list: plot_data.append(("Speed |u|", speed, "viridis"))
if "sound"    in field_list: plot_data.append(("Sound Speed c", sound_init, "plasma"))

nplots = len(plot_data)
fig, axes = plt.subplots(1, nplots, figsize=(4.5 * nplots, 4))
if nplots == 1:
    axes = [axes]

artists = {}
for ax, (title, arr, cmap) in zip(axes, plot_data):
    if title == "Level-set phi":
        im = ax.imshow(arr, extent=extent, origin='lower', cmap=cmap, vmin=-1.0, vmax=1.0)
    else:
        im = ax.imshow(arr, extent=extent, origin='lower', cmap=cmap)
    ax.set_title(title)
    cbar = plt.colorbar(im, ax=ax)
    artists[title] = (im, cbar)

def update(frame):
    global step_counter
    for _ in range(num_steps_per_frame):
        solver.step(dt)
        step_counter += 1
    
    # Redistance every few frames to keep phi stable
    if frame % redistance_every == 0:
        print(f"[redistance] frame={frame} step={step_counter} iterations=5")
        solver.redistance(iterations=5)
        
    data = solver.get_data()
    density = np.clip(data['rho_a1'] + data['rho_a2'], 1e-6, None)  # floor to avoid div/neg issues
    ru64 = np.nan_to_num(data['ru'], copy=False).astype(np.float64); rv64 = np.nan_to_num(data['rv'], copy=False).astype(np.float64); density64 = density.astype(np.float64)
    vel_sq = np.clip((ru64*ru64 + rv64*rv64) / (density64 + 1e-6)**2, 0.0, 1e20)
    speed = np.sqrt(vel_sq)
    pressure = data.get('p', None)
    if pressure is None:
        pressure = np.ones_like(density)
    sound_speed = data.get('c', None)
    if sound_speed is None:
        sound_speed = np.ones_like(density)

    for title, (im, cbar) in artists.items():
        if title == "Density":
            im.set_data(density); im.set_clim(density.min(), density.max()); cbar.update_normal(im)
        elif title == "Level-set phi":
            im.set_data(data['phi']); im.set_clim(-1.0, 1.0); cbar.update_normal(im)
        elif title == "Energy E":
            im.set_data(data['E']); im.set_clim(data['E'].min(), data['E'].max()); cbar.update_normal(im)
        elif title == "Speed |u|":
            im.set_data(speed); im.set_clim(speed.min(), speed.max()); cbar.update_normal(im)
        elif title == "Pressure p":
            im.set_data(pressure); im.set_clim(pressure.min(), pressure.max()); cbar.update_normal(im)
        elif title == "Sound Speed c":
            im.set_data(sound_speed); im.set_clim(sound_speed.min(), sound_speed.max()); cbar.update_normal(im)
    
    # Diagnostics
    # Mass diagnostics per material
    mass_h = float((data['rho_a1'] * (dx * dx)).sum())
    mass_u = float((data['rho_a2'] * (dx * dx)).sum())
    total_mass = mass_h + mass_u
    total_energy = float((data['E'] * (dx * dx)).sum())
    rho_min = float(density.min()); rho_max = float(density.max())
    E_min = float(data['E'].min()); E_max = float(data['E'].max())
    phi_min = float(data['phi'].min()); phi_max = float(data['phi'].max())
    p_min = float(pressure.min()); p_max = float(pressure.max())
    c_min = float(sound_speed.min()); c_max = float(sound_speed.max())
    total_pressure = float((pressure * (dx * dx)).sum())
    finite_ok = np.isfinite(density).all() and np.isfinite(data['E']).all() and np.isfinite(pressure).all()
    print(f"frame={frame} step={step_counter} mass={total_mass:.4e} (H={mass_h:.4e}, U={mass_u:.4e}) energy={total_energy:.4e} p_sum={total_pressure:.4e} "
          f"rho=[{rho_min:.3e},{rho_max:.3e}] E=[{E_min:.3e},{E_max:.3e}] phi=[{phi_min:.3e},{phi_max:.3e}] "
          f"p=[{p_min:.3e},{p_max:.3e}] c=[{c_min:.3e},{c_max:.3e}]")
    if not finite_ok:
        # Identify where NaNs first appear for debugging
        nan_mask = ~np.isfinite(density)
        if nan_mask.any():
            idx = np.argwhere(nan_mask)[0]
            print(f"NaN in density at {idx}, value={density[tuple(idx)]}")
        nan_mask_E = ~np.isfinite(data['E'])
        if nan_mask_E.any():
            idx = np.argwhere(nan_mask_E)[0]
            print(f"NaN in energy at {idx}, value={data['E'][tuple(idx)]}")
        nan_mask_p = ~np.isfinite(pressure)
        if nan_mask_p.any():
            idx = np.argwhere(nan_mask_p)[0]
            print(f"NaN in pressure at {idx}, value={pressure[tuple(idx)]}")
        # Stop animation updates to avoid flooding
        ani.event_source.stop()
    
    if frame % snapshot_every == 0:
        fname = f"{snapshot_prefix}{frame:04d}.png"
        fig.savefig(fname, dpi=150)
        print(f"[snapshot] saved {fname}")
    
    return [im for im, _ in artists.values()]

ani = FuncAnimation(fig, update, frames=args.frames, interval=50, blit=False)
plt.show()
