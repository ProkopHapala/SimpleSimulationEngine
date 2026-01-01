import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. Physics Engine Constants & State
# ==========================================

import argparse

# Parameters (overridden by CLI)
dt      = 0.01
n_steps = 2000
mass    = 1.0
inertia = 0.5  # Moment of inertia
k_bond  = 50.0
l_bond  = 1.0
k_align = 2.55 # Strength of orthogonality constraint

# Initial State (overridden by CLI)
pos   = np.array([[0.0, 0.0], [1.0, 0.0]])
vel   = np.zeros_like(pos)
theta = np.array([np.pi/2, np.pi/2 + 0.3])
omega = np.array([0.0, 0.0])

# ==========================================
# 2. Force Evaluation (The Physics)
# ==========================================

def get_forces(p, v, th, w):
    """
    Computes Forces and Torques.
    Includes coupling between linear and angular parts.
    """
    F = np.zeros_like(p)
    T = np.zeros_like(th)
    PE = 0.0
    
    # 1. Bond Force (Spring)
    d     = p[1] - p[0]
    dist  = np.linalg.norm(d)
    dir_u = d / dist
    
    # F = k * (r - r0)
    f_mag  = k_bond * (dist - l_bond)
    F_bond = f_mag * dir_u
    
    F[0] += F_bond
    F[1] -= F_bond
    PE += 0.5 * k_bond * (dist - l_bond)**2
    
    # 2. Alignment Force (Sigma-Pi)
    # The stick h should be orthogonal to bond direction u.
    # Energy = k * (h . u)^2   (0 when perpendicular)
    
    for i in range(2):
        # Orientation vector h
        h = np.array([np.cos(th[i]), np.sin(th[i])])
        
        # Dot product (should be 0)
        # We use 'dir_u' for bond direction. 
        # Note: dir_u points 0->1. For particle 1, bond direction is reversed.
        u_local = dir_u if i == 0 else -dir_u
        
        c = np.dot(h, u_local)
        
        # --- Torque on Angle ---
        # dE/dTheta. h' = (-sin, cos)
        h_prime = np.array([-np.sin(th[i]), np.cos(th[i])])
        dE_dth = 2 * k_align * c * np.dot(h_prime, u_local)
        
        T[i] -= dE_dth # Torque = -Gradient
        
        # --- Recoil Force on Position ---
        # This is the tricky part often missed!
        # Rotating the bond vector changes the energy.
        # V = k * (h . u)^2
        # F = - grad_r(V)
        # grad_r(u) is complex. 
        # Force transverse to bond: F_trans = - (2 * k * c) * (h_perp_component) / r
        
        # Analytic gradient of (h . u) w.r.t position of *this* atom
        # u = (p_target - p_me) / |r|
        # ... standard derivation ...
        # Force = -2 * k * c * [ (h - c*u) / r ]
        
        force_recoil = 2 * k_align * c * (h - c * u_local) / dist
        
        # Apply to self
        F[i] += force_recoil
        # Apply opposite to neighbor
        neighbor = 1 - i
        F[neighbor] -= force_recoil
        
        PE += k_align * c**2

    return F, T, PE

# ==========================================
# 3. Propagators
# ==========================================

def semi_implicit_init(p, v, th, w):
    return {}

def semi_implicit_step(p, v, th, w, aux):
    F, T, pe = get_forces(p, v, th, w)
    a = F/mass
    alpha = T/inertia
    v += a*dt
    w += alpha*dt
    p += v*dt
    th += w*dt
    return p, v, th, w, pe, aux

def velocity_verlet_init(p, v, th, w):
    return {}

def velocity_verlet_step(p, v, th, w, aux):
    F, T, pe = get_forces(p, v, th, w)
    a = F/mass
    alpha = T/inertia
    v += 0.5*a*dt
    w += 0.5*alpha*dt
    p += v*dt
    th += w*dt
    F_new, T_new, pe_new = get_forces(p, v, th, w)
    a_new = F_new/mass
    alpha_new = T_new/inertia
    v += 0.5*a_new*dt
    w += 0.5*alpha_new*dt
    return p, v, th, w, pe_new, aux

def split_verlet_init(p, v, th, w):
    F_init, T_init, _ = get_forces(p, v, th, w)
    v_half = v - 0.5*(F_init/mass)*dt
    w_half = w - 0.5*(T_init/inertia)*dt
    return (v_half, w_half, v.copy(), w.copy())

def split_verlet_step(p, v, th, w, aux):
    v_half, w_half, v_pred, w_pred = aux
    F, T, pe = get_forces(p, v_pred, th, w_pred)
    a = F/mass
    alpha = T/inertia
    v_now = v_half + 0.5*a*dt
    w_now = w_half + 0.5*alpha*dt
    v_half += a*dt
    w_half += alpha*dt
    p += v_half*dt
    th += w_half*dt
    v_pred = v_half + 0.5*a*dt
    w_pred = w_half + 0.5*alpha*dt
    return p, v_now, th, w_now, pe, (v_half, w_half, v_pred, w_pred)

def stored_force_verlet_init(p, v, th, w):
    F, T, _ = get_forces(p, v, th, w)
    return (F, T)

def stored_force_verlet_step(p, v, th, w, aux):
    F_old, T_old = aux
    v_half = v + 0.5*(F_old/mass)*dt
    w_half = w + 0.5*(T_old/inertia)*dt
    p_new = p + v_half*dt
    th_new = th + w_half*dt
    F_new, T_new, pe = get_forces(p_new, v_half, th_new, w_half)
    v_new = v_half + 0.5*(F_new/mass)*dt
    w_new = w_half + 0.5*(T_new/inertia)*dt
    return p_new, v_new, th_new, w_new, pe, (F_new, T_new)

# ==========================================
# 4. Solvers
# ==========================================

def run_simulation(init_fn, step_fn, verbose_every=0, label=""):
    p  = pos.copy()
    v  = vel.copy()
    th = theta.copy()
    w  = omega.copy()
    aux = init_fn(p, v, th, w)
    history = {'E':[], 'L':[], 'P':[], 'pos':[], 'com':[]}

    for step in range(n_steps):
        p, v, th, w, pe, aux = step_fn(p, v, th, w, aux)
        ke_trans = 0.5*mass*np.sum(v**2)
        ke_rot = 0.5*inertia*np.sum(w**2)
        total_E = ke_trans + ke_rot + pe
        L_total = 0
        for i in range(2):
            L_total += np.cross(p[i], mass*v[i]) + inertia*w[i]
        P_total = mass*np.sum(v, axis=0)
        com = np.mean(p, axis=0)
        history['E'].append(total_E)
        history['L'].append(L_total)
        history['P'].append(P_total)
        history['pos'].append(p.copy())
        history['com'].append(com)
        if verbose_every>0 and step%verbose_every==0:
            print(f"[{label or 'sim'}] step={step:6d}  E={total_E: .9e}  L={L_total: .9e}  Px={P_total[0]: .9e}  Py={P_total[1]: .9e}")

    return history

def _parse_vec2(text):
    parts = text.split(',')
    if len(parts)!=2:
        raise argparse.ArgumentTypeError("Expected 'x,y'")
    return np.array([float(parts[0]), float(parts[1])])

def _parse_angle(text):
    return float(text)

def plot_panel(ax, series, labels, title, ylabel, styles=None):
    styles = styles or {}
    for idx, data in enumerate(series):
        style = styles.get(idx, {})
        ax.plot(data, **style, label=labels[idx])
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.legend()

def plot_results(res_semi, res_vv, res_split, res_stored):
    fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, figsize=(8, 18))
    pos_split = np.array(res_split['pos'])
    ax1.plot(pos_split[:,0,0], pos_split[:,0,1], ls='-', lw=0.5, label='Atom 1 (Split)')
    ax1.plot(pos_split[:,1,0], pos_split[:,1,1], ls='-', lw=0.5, label='Atom 2 (Split)')
    ax1.set_title("Trajectory (Split Verlet)")
    ax1.set_aspect('equal')
    ax1.legend()
    labels = ["Semi-Implicit (Leapfrog)", "Std Velocity Verlet (2x Force)", "Split Verlet (1x Force)", "Stored-Force Verlet (1x Force)"]
    styles = {0:{'ls':'-','lw':0.5},1:{'ls':':','lw':1.5},2:{'ls':'-','lw':0.5},3:{'ls':'-','lw':1.0}}
    plot_panel(ax2, [res_semi['E'], res_vv['E'], res_split['E'], res_stored['E']], labels, "Total Energy (Hamiltonian)", "Energy", styles)
    plot_panel(ax3, [res_semi['L'], res_vv['L'], res_split['L'], res_stored['L']], labels, "Total Angular Momentum", "L", styles)
    P_split = np.array(res_split['P'])
    plot_panel(ax4, [P_split[:,0], P_split[:,1]], ["Px (Split)", "Py (Split)"], "Total Linear Momentum (Split Verlet)", "P", {0:{'ls':'-','lw':0.7},1:{'ls':'--','lw':0.7}})
    com_split = np.array(res_split['com'])
    plot_panel(ax5, [com_split[:,0], com_split[:,1]], ["x_COM", "y_COM"], "Center of Mass Position (Split Verlet)", "COM", {0:{'ls':'-','lw':0.7},1:{'ls':'--','lw':0.7}})
    plt.tight_layout()
    plt.show()

# USER LIKE THIS:
# python rigid_body_proof.py --dt 0.002 --n-steps 5000 --pos0 0,0 --pos1 1,0.1 --theta0 1.57 --theta1 1.57 --omega0 5 --omega1 0 --verbose-every 100

def main():
    global dt, n_steps, mass, inertia, k_bond, l_bond, k_align, pos, vel, theta, omega
    parser = argparse.ArgumentParser(description="Rigid body toy: compare propagators")
    parser.add_argument("--dt", type=float, default=dt)
    parser.add_argument("--n-steps", type=int, default=n_steps)
    parser.add_argument("--mass", type=float, default=mass)
    parser.add_argument("--inertia", type=float, default=inertia)
    parser.add_argument("--k-bond", type=float, default=k_bond)
    parser.add_argument("--l-bond", type=float, default=l_bond)
    parser.add_argument("--k-align", type=float, default=k_align)
    parser.add_argument("--pos0", type=_parse_vec2, default="{},{}".format(*pos[0]))
    parser.add_argument("--pos1", type=_parse_vec2, default="{},{}".format(*pos[1]))
    parser.add_argument("--vel0", type=_parse_vec2, default="{},{}".format(*vel[0]))
    parser.add_argument("--vel1", type=_parse_vec2, default="{},{}".format(*vel[1]))
    parser.add_argument("--theta0", type=_parse_angle, default=theta[0])
    parser.add_argument("--theta1", type=_parse_angle, default=theta[1])
    parser.add_argument("--omega0", type=float, default=omega[0])
    parser.add_argument("--omega1", type=float, default=omega[1])
    parser.add_argument("--verbose-every", type=int, default=100, help="Print E/L/P every N steps (0 disables)")
    parser.add_argument("--no-plot", action="store_true", help="Skip matplotlib plots")
    args = parser.parse_args()

    dt = args.dt
    n_steps = args.n_steps
    mass = args.mass
    inertia = args.inertia
    k_bond = args.k_bond
    l_bond = args.l_bond
    k_align = args.k_align
    pos = np.array([args.pos0, args.pos1], dtype=float)
    vel = np.array([args.vel0, args.vel1], dtype=float)
    theta = np.array([args.theta0, args.theta1], dtype=float)
    omega = np.array([args.omega0, args.omega1], dtype=float)

    res_semi   = run_simulation( semi_implicit_init,   semi_implicit_step, verbose_every=args.verbose_every, label="semi")
    res_vv     = run_simulation( velocity_verlet_init, velocity_verlet_step, verbose_every=args.verbose_every, label="vv")
    res_split  = run_simulation( split_verlet_init,    split_verlet_step, verbose_every=args.verbose_every, label="split")
    res_stored = run_simulation( stored_force_verlet_init, stored_force_verlet_step, verbose_every=args.verbose_every, label="stored")

    print(f"Initial total linear momentum (input state): {mass*np.sum(vel,axis=0)}")

    # Plotting
    if args.no_plot:
        return
    drift_semi = (res_semi['E'][-1] - res_semi['E'][0])
    drift_split = (res_split['E'][-1] - res_split['E'][0])
    print(f"Energy Drift Semi-Implicit: {drift_semi:.6f}")
    print(f"Energy Drift Split-Verlet:  {drift_split:.6f}")
    plot_results(res_semi, res_vv, res_split, res_stored)

if __name__ == "__main__":
    main()