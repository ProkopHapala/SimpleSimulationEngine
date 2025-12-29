import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Rectangle
from scipy.linalg import solve_banded
from matplotlib import colormaps

# --- CONFIGURATION ---
parser = argparse.ArgumentParser()
parser.add_argument("--cells",        type=int,   default=40,   help="Spatial cells")
parser.add_argument("--groups",       type=int,   default=8,    help="Energy groups")
parser.add_argument("--dt",           type=float, default=1e-3, help="Timestep (us)")
parser.add_argument("--core-radius",  type=float, default=3.0,  help="Core (fuel) outer radius (cm)")
parser.add_argument("--mod-radius",   type=float, default=6.5,  help="Moderator outer radius (cm)")
parser.add_argument("--poison-scale", type=float, default=0.1,  help="Multiplier for poison absorption/scatter")
# Burnup accelerator: Real reactors take months, bombs take nanoseconds. 
#parser.add_argument("--burn-scale",      type=float, default=0.005, help="Burnup speed multiplier")
parser.add_argument("--implosion-speed",     type=float, default=10.0, help="Initial implosion speed (cm/us)")
parser.add_argument("--mod-thickness",       type=float, default=2.0, help="Moderator thickness (cm); overrides --mod-radius via core_radius + mod_thickness")
parser.add_argument("--refl-thickness",      type=float, default=2.0, help="Reflector/tamper thickness (cm); sets outer radius = mod_radius + refl_thickness")
parser.add_argument("--outer-radius",        type=float, default=10.0, help="Outer radius (cm) if refl-thickness is not provided")
parser.add_argument("--fuel-init-density",   type=float, default=0.2, help="Initial relative fuel density (dimensionless, 1.0 = reference)")
parser.add_argument("--nu",                  type=float, default=2.8, help="Neutrons per fission")
parser.add_argument("--target-radius",       type=float, default=2.0, help="Final compressed outer radius target (cm)")
parser.add_argument("--fuel-number-density", type=float, default=4.8e22, help="Atoms per cm^3 for density=1 fuel")
args = parser.parse_args()

# Material property table
# Note: These are Macroscopic Cross Sections (1/cm) for density=1.0 material
MATERIAL_XS = {
    'fuel': {
        'f_fast': 0.25, 
        'f_thrm': 0.80, 
        'a_fast': 0.005, 
        'a_thrm': 0.02, 
        's': 0.25      
    },
    'pois': {'a': 0.10, 's': 0.10},
    'mod':  {'a': 0.0001, 's': 1.2}, 
    'refl': {'a': 0.001,  's': 0.8}, 
}

class ImplosionBurnupSim:
    def __init__(self):
        # 1. Geometry & Units (cm, microseconds)
        self.N = args.cells
        self.G = args.groups
        self.dt = args.dt
        # Allow geometry to be configured by thickness for moderator and reflector
        self.core_radius = args.core_radius
        self.mod_radius = args.mod_radius if args.mod_thickness is None else self.core_radius + args.mod_thickness
        outer_from_refl = None
        if args.refl_thickness is not None:
            outer_from_refl = self.mod_radius + args.refl_thickness
        # Ensure outer radius is at least beyond moderator; fall back to provided outer-radius
        self.R_initial = max(self.mod_radius, outer_from_refl or args.outer_radius)
        
        # Lagrangian Mesh
        self.r_faces = np.linspace(0, self.R_initial, self.N + 1)
        self.update_geometry()
        
        # 2. Physics: Energy Strata (Velocity cm/us)
        # 14 MeV (Fast) -> Thermal
        self.v = np.geomspace(5200.0, 0.22, self.G)
        
        # 3. Material Tracking (Number Densities)
        self.mat_densities = {
            'fuel': np.zeros(self.N),
            'pois': np.zeros(self.N),
            'mod':  np.zeros(self.N),
            'refl': np.zeros(self.N)
        }
        
        # Initial Setup (Layers) from CLI radii
        r_mid = 0.5 * (self.r_faces[:-1] + self.r_faces[1:])
        mask_core = r_mid < self.core_radius
        mask_mod = (r_mid >= self.core_radius) & (r_mid < self.mod_radius)
        mask_refl = r_mid >= self.mod_radius
        # Start intentionally subcritical: allow initial fuel density to be below reference (1.0)
        self.mat_densities['fuel'][mask_core] = args.fuel_init_density
        self.mat_densities['mod'][mask_mod] = 1.0
        self.mat_densities['refl'][mask_refl] = 1.0
        
        # Store initial masks to track layer boundaries for plotting
        self.initial_masks = {'core': mask_core, 'mod': mask_mod, 'refl': mask_refl}
        # Total initial fuel atoms
        self.fuel_atoms0 = np.sum(self.mat_densities['fuel'] * self.V_cells) * args.fuel_number_density
        self.cum_neutrons = 0.0
        self.cum_burn_atoms = 0.0
        print(f"[init] fuel_atoms={self.fuel_atoms0:.3e} | core_r={args.core_radius} cm | mod_r={args.mod_radius} cm")
        # Print key run properties so they are visible in logs
        print("[props] cells={:d} | groups={:d} | dt={:.3e} us | nu={:.2f} | poison_scale={:.3e}".format( self.N, self.G, self.dt, args.nu, args.poison_scale))
        print("[xs] fuel: f_fast={:.3f}, f_thrm={:.3f}, a_fast={:.4f}, a_thrm={:.3f}, s={:.3f}".format(
            MATERIAL_XS['fuel']['f_fast'], MATERIAL_XS['fuel']['f_thrm'],
            MATERIAL_XS['fuel']['a_fast'], MATERIAL_XS['fuel']['a_thrm'],
            MATERIAL_XS['fuel']['s']))
        print("[xs] mod: a={:.5f} (scaled 1/v), s={:.3f} | refl: a={:.3f}, s={:.3f} | poison: a={:.3f}*scale, s={:.3f}*scale".format(
            MATERIAL_XS['mod']['a'], MATERIAL_XS['mod']['s'],
            MATERIAL_XS['refl']['a'], MATERIAL_XS['refl']['s'],
            MATERIAL_XS['pois']['a'], MATERIAL_XS['pois']['s']))

        # 4. Microscopic Cross Sections
        self.micro_xs = self._generate_micro_xs()
        
        # 5. State Variables
        self.phi = np.zeros((self.G, self.N))
        self.phi[0, mask_core] = 1e8 # Initial spark
        self.k_eff = 0.0
        self.time = 0.0
        self.print_every = 20
        self.total_power = 0.0
        self.step_count = 0
        
        # Hydrodynamics
        self.u = np.zeros(self.N + 1) # Velocity

    def _generate_micro_xs(self):
        """ Define physics properties for each material/group """
        xs = {
            'fuel': {'f': np.zeros(self.G), 'a': np.zeros(self.G), 's': np.zeros(self.G)},
            'pois': {'f': np.zeros(self.G), 'a': np.zeros(self.G), 's': np.zeros(self.G)},
            'mod':  {'f': np.zeros(self.G), 'a': np.zeros(self.G), 's': np.zeros(self.G)},
            'refl': {'f': np.zeros(self.G), 'a': np.zeros(self.G), 's': np.zeros(self.G)},
        }
        
        for g in range(self.G):
            # 1/v scaling factor
            vr = np.sqrt(self.v[0] / self.v[g]) 
            
            # FUEL
            if g < self.G // 2: # Fast groups
                xs['fuel']['f'][g] = MATERIAL_XS['fuel']['f_fast']
                xs['fuel']['a'][g] = MATERIAL_XS['fuel']['a_fast']
            else: # Thermal groups
                xs['fuel']['f'][g] = MATERIAL_XS['fuel']['f_thrm']
                xs['fuel']['a'][g] = MATERIAL_XS['fuel']['a_thrm']
            xs['fuel']['s'][g] = MATERIAL_XS['fuel']['s']
            
            # POISON
            xs['pois']['a'][g] = MATERIAL_XS['pois']['a'] * vr * args.poison_scale
            xs['pois']['s'][g] = MATERIAL_XS['pois']['s'] * args.poison_scale
            
            # MODERATOR
            xs['mod']['s'][g] = MATERIAL_XS['mod']['s']
            xs['mod']['a'][g] = MATERIAL_XS['mod']['a'] * vr
            
            # REFLECTOR
            xs['refl']['s'][g] = MATERIAL_XS['refl']['s']
            xs['refl']['a'][g] = MATERIAL_XS['refl']['a']
            
        return xs

    def ensure_finite(self, name, arr):
        if not np.all(np.isfinite(arr)):
            print(f"[NaN/Inf] {name} at t={self.time:.6f}us")
            raise RuntimeError(f"Non-finite values detected in {name}")
        return arr

    def update_geometry(self):
        self.dr = np.diff(self.r_faces)
        self.r_centers = 0.5 * (self.r_faces[:-1] + self.r_faces[1:])
        self.A_faces = 4 * np.pi * self.r_faces**2
        self.V_cells = (4/3) * np.pi * (self.r_faces[1:]**3 - self.r_faces[:-1]**3)

    def step_hydro(self):
        target_R = args.target_radius  # Compress toward configurable radius
        current_R = self.r_faces[-1]
        
        if current_R > target_R:
            compression_speed = args.implosion_speed # cm/us
            self.u = -compression_speed * (self.r_faces / self.R_initial)
        else:
            self.u[:] = 0.0 
            
        self.r_faces += self.u * self.dt
        V_old = self.V_cells.copy()
        self.update_geometry()
        compression_ratio = V_old / self.V_cells
        
        for mat in self.mat_densities:
            self.mat_densities[mat] *= compression_ratio

    def step_physics(self):
        # 1. Build Macroscopic XS
        sig_t = np.zeros((self.N, self.G))
        sig_a_rem = np.zeros((self.N, self.G)) # Removal (Absorb + Fission)
        sig_f = np.zeros((self.N, self.G))
        sig_s = np.zeros((self.N, self.G))
        
        for mat, N_array in self.mat_densities.items():
            dens = N_array[:, None] 
            sig_f += dens * self.micro_xs[mat]['f']
            sig_s += dens * self.micro_xs[mat]['s']
            sig_a_rem += dens * (self.micro_xs[mat]['a'] + self.micro_xs[mat]['f'])
            
        sig_t = sig_a_rem + sig_s + 1e-4

        # 2. Source Calculation (Fission)
        nu = args.nu
        # Fission rate in fissions/us per cell
        fission_rate_cells = np.sum(sig_f * self.phi.T, axis=1) * self.V_cells 
        fission_neutrons_rate_cells = nu * fission_rate_cells 
        
        step_neutrons = np.sum(fission_neutrons_rate_cells * self.dt)
        self.cum_neutrons += step_neutrons
        
        # Burnup Tracking
        fuel_fission_rate_atoms = fission_rate_cells
        fuel_capture_rate_atoms = np.sum((sig_a_rem - sig_f) * self.phi.T, axis=1) * self.V_cells

        # K-Eff Variables
        total_prod = np.sum(fission_neutrons_rate_cells)
        total_loss = 0.0 # Absorbtion + Leakage (NOT Scattering)
        self.total_power = total_prod
        
        # 3. Solve Diffusion
        for g in range(self.G):
            # Source Term (Fission Spectrum + Downscattering)
            S_vol = np.zeros(self.N)
            
            # Fission Spectrum: 90% Grp 0, 10% Grp 1
            if g == 0: S_vol += (fission_neutrons_rate_cells * 0.90) / self.V_cells
            if g == 1: S_vol += (fission_neutrons_rate_cells * 0.10) / self.V_cells
            
            # Scattering IN from previous group
            if g > 0:
                S_vol += sig_s[:, g-1] * self.phi[g-1]
            
            # --- SOLVER LOGIC (Neutrons leaving this group) ---
            # Removal for matrix = Absorption + Fission + Scatter OUT
            rem_g_solver = sig_a_rem[:, g].copy()
            if g < self.G - 1:
                rem_g_solver += sig_s[:, g]
            
            # --- K_EFF LOGIC (Neutrons leaving the SYSTEM) ---
            # Loss for K_eff = Absorption + Fission + Leakage
            # We do NOT include scattering out here, because those neutrons exist in group g+1
            total_loss += np.sum(sig_a_rem[:, g] * self.phi[g] * self.V_cells)
            
            # Diffusion Coeff
            D = 1.0 / (3.0 * np.maximum(sig_t[:, g], 1e-2))
            
            # Matrix Setup
            inertia = self.V_cells / (self.v[g] * self.dt)
            diag = inertia + rem_g_solver * self.V_cells
            
            # Coupling
            D_mid = np.zeros(self.N+1)
            D_mid[1:-1] = 2*D[:-1]*D[1:]/(D[:-1]+D[1:])
            dr_centers = np.diff(self.r_centers) 
            conductance = D_mid[1:-1] * self.A_faces[1:-1] / np.maximum(dr_centers, 1e-6)
            
            # Leakage (Surface Loss)
            leakage_coeff = 0.5 * self.A_faces[-1]
            
            # Add leakage to K_eff loss
            total_loss += leakage_coeff * self.phi[g, -1]
            
            # Banded Solver
            main = diag.copy()
            upper = np.zeros(self.N-1)
            lower = np.zeros(self.N-1)
            
            main[1:] += conductance  
            lower[:] = -conductance
            main[:-1] += conductance 
            upper[:] = -conductance
            main[-1] += leakage_coeff
            
            rhs = self.phi[g] * inertia + S_vol * self.V_cells
            
            ab = np.vstack((np.append(0, upper), main, np.append(lower, 0)))
            self.phi[g] = solve_banded((1, 1), ab, rhs)

        # Update K_eff (avoid div by zero)
        if total_loss > 1e-20: 
            self.k_eff = total_prod / total_loss
        else:
            self.k_eff = 0.0

        # 4. Burnup (Depletion)
        atoms_burned_cells = (fuel_fission_rate_atoms + fuel_capture_rate_atoms) * self.dt
        
        fuel_atoms_cells = self.mat_densities['fuel'] * self.V_cells * args.fuel_number_density
        
        burn_fraction = np.zeros_like(fuel_atoms_cells)
        mask = fuel_atoms_cells > 1e-20
        burn_fraction[mask] = atoms_burned_cells[mask] / fuel_atoms_cells[mask]
        
        actual_burned_atoms_cells = fuel_atoms_cells * burn_fraction
        self.mat_densities['fuel'] *= (1.0 - burn_fraction)
        self.mat_densities['pois'] += actual_burned_atoms_cells / (self.V_cells * args.fuel_number_density + 1e-30)
        
        self.cum_burn_atoms += np.sum(actual_burned_atoms_cells)
        self.mat_densities['fuel'] = np.maximum(self.mat_densities['fuel'], 0)
        
        self.time += self.dt
        self.step_count += 1
        if self.step_count % self.print_every == 0:
            fuel_atoms = np.sum(self.mat_densities['fuel'] * self.V_cells) * args.fuel_number_density
            burnup_frac = 1.0 - (fuel_atoms / (self.fuel_atoms0 + 1e-30))
            print(f"[t={self.time:.4f} us | step={self.step_count}] k_eff={self.k_eff:.3f} | power={self.total_power:.2e} | peak_flux={np.max(self.phi):.2e} | fuel_atoms={fuel_atoms:.3e} | burnup={burnup_frac*100:.3f}%")

# --- VISUALIZATION SETUP ---
sim = ImplosionBurnupSim()

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={'height_ratios': [3, 1]})

ax1.set_yscale('log')
ax1.set_xlim(0, 10.5)
ax1.set_ylim(1, 1e25)
ax1.set_ylabel(r'Flux [n / $cm^2 \mu s$]')
ax1.set_xlabel('Radius [cm]')
ax1.grid(True, which='both', alpha=0.3)

bg_core = Rectangle((0,0), 0, 1, color='red', alpha=0.15, label='Core (Fuel)', transform=ax1.get_xaxis_transform())
bg_mod  = Rectangle((0,0), 0, 1, color='blue', alpha=0.15, label='Moderator', transform=ax1.get_xaxis_transform())
bg_refl = Rectangle((0,0), 0, 1, color='gray', alpha=0.15, label='Tamper', transform=ax1.get_xaxis_transform())
ax1.add_patch(bg_core)
ax1.add_patch(bg_mod)
ax1.add_patch(bg_refl)

cmap = colormaps.get_cmap('plasma').resampled(sim.G)
lines = []
energies_eV = (sim.v / sim.v[0])**2 * 14e6 
for g in range(sim.G):
    ln, = ax1.plot([], [], color=cmap(g), lw=1.5)
    lines.append(ln)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=sim.G-1))
cb = plt.colorbar(sm, ax=ax1, ticks=[0, sim.G-1])
cb.ax.set_yticklabels([f'Fast\n~{energies_eV[0]/1e6:.1f} MeV', f'Thermal\n~{energies_eV[-1]:.2f} eV'])

ax1b = ax1.twinx()
ax1b.set_ylim(0, 5.0) 
ax1b.set_ylabel('Fuel Density (Relative)')
ln_dens, = ax1b.plot([], [], 'k--', lw=1.5, label='Fuel Density')
ax1.legend([ln_dens, bg_core, bg_mod, bg_refl], ['Fuel Density', 'Core', 'Mod', 'Tamper'], loc='upper right', fontsize='small')

info_text = ax1.text(0.02, 0.96, "", transform=ax1.transAxes, va='top', fontfamily='monospace', fontsize=9)

ax2.set_xlim(0, 2.0)
ax2.set_ylim(0, 3.0) # Expanded range for supercritical k
ax2.set_xlabel('Time [$\mu s$]')
ax2.set_ylabel('$k_{eff}$')
ax2.axhline(1.0, color='red', linestyle=':', alpha=0.5)
ax2.grid(True)
ln_k, = ax2.plot([], [], 'g-', lw=2, label='k_eff')
# Secondary axis for neutron population and fuel atoms
ax2b = ax2.twinx()
ax2b.set_ylabel('Total Neutrons / Fuel Atoms')
ax2b.set_yscale('log')
ax2b.set_ylim(1e1, 1e28)
ln_neutrons, = ax2b.plot([], [], 'b:', lw=1.2, label='Total Neutrons')
ln_fuel_atoms, = ax2b.plot([], [], 'r--', lw=1.2, label='Fuel Atoms')
ax2.legend([ln_k, ln_neutrons, ln_fuel_atoms], ['k_eff', 'Total Neutrons', 'Fuel Atoms'], loc='upper left', fontsize='small')

hist_t, hist_k, hist_neutrons, hist_fuel_atoms = [], [], [], []

def update(frame):
    for _ in range(5):
        sim.step_hydro()
        sim.step_physics()
    
    for g in range(sim.G):
        lines[g].set_data(sim.r_centers, sim.phi[g])
    
    ln_dens.set_data(sim.r_centers, sim.mat_densities['fuel'])
    dens_max = np.max(sim.mat_densities['fuel'])
    ax1b.set_ylim(0, dens_max*1.2 if dens_max > 0 else 1.0)
    
    def get_outer_r(mask):
        idx = np.where(mask)[0]
        return sim.r_faces[idx[-1]+1] if len(idx) > 0 else 0
    
    r_c = get_outer_r(sim.initial_masks['core'])
    r_m = get_outer_r(sim.initial_masks['mod'])
    r_t = sim.r_faces[-1]
    
    bg_core.set_width(r_c)
    bg_mod.set_x(r_c)
    bg_mod.set_width(r_m - r_c)
    bg_refl.set_x(r_m)
    bg_refl.set_width(r_t - r_m)
    
    hist_t.append(sim.time)
    hist_k.append(sim.k_eff)
    hist_neutrons.append(np.sum(sim.phi * sim.V_cells))
    hist_fuel_atoms.append(np.sum(sim.mat_densities['fuel'] * sim.V_cells) * args.fuel_number_density)
    ln_k.set_data(hist_t, hist_k)
    ln_neutrons.set_data(hist_t, hist_neutrons)
    ln_fuel_atoms.set_data(hist_t, hist_fuel_atoms)
    
    if sim.time > ax2.get_xlim()[1]:
        ax2.set_xlim(0, sim.time * 1.5)
        ax2b.set_xlim(0, sim.time * 1.5)
        
    info_text.set_text(f"t={sim.time:.4f} us\nk_eff={sim.k_eff:.3f}\npower={sim.total_power:.2e}")
    ax1.set_title(f"Time: {sim.time:.3f} $\mu s$ | $k_{{eff}}$: {sim.k_eff:.3f}")

    return lines + [ln_dens, bg_core, bg_mod, bg_refl, ln_k, ln_neutrons, ln_fuel_atoms, info_text]

ani = animation.FuncAnimation(fig, update, frames=500, interval=50, blit=False)
plt.show()