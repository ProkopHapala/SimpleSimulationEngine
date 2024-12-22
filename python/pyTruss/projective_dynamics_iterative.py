import numpy as np
from projective_dynamics import build_neighbor_list, make_pd_matrix, make_pd_rhs, update_velocity

def jacobi_iteration(A, b, x):
    Aii = np.diag(A)
    Ax = A @ x
    return (b - (Ax - Aii*x) ) / Aii

def estimate_spectral_radius(A, x0, n_iter=10):
    """Estimate spectral radius using power iteration"""
    x      = x0.copy()
    for _ in range(n_iter):
        x      = A @ x
        x_next = A @ x
        eK1    = np.dot(x_next, x)
        eK2    = np.dot(x, x)
        rho    = np.sqrt( abs( eK1 / eK2 ) )
    return min(rho, 0.9992)  # Cap at 0.9992 as suggested in paper

def calculate_omega(k, rho, prev_omega, delay_start=10):
    """Calculate Chebyshev weight"""
    if k < delay_start:
        return 1.0
    elif k == delay_start:
        return 2.0 / (2.0 - rho**2)
    else:
        return 4.0 / (4.0 - rho**2 * prev_omega)

def solve_pd_jacobi(points, velocity, bonds, masses, ks, dt=0.1, n_iter=100, 
                   inner_iter=50, gravity=np.array([0, -9.81, 0]), 
                   fixed_points=None, call_back=None, damping=0.01):
    """Solve the system using projective dynamics with Jacobi iterations"""
    n_points = len(points)
    neighbs = build_neighbor_list(bonds, n_points)
    A, Mt = make_pd_matrix(neighbs, bonds, masses, dt, ks)
    
    # Initialize
    pos = points.copy()
    pos_cor = points.copy()
    pos_pred = points.copy()
    l0s = np.array([np.linalg.norm(points[j] - points[i]) for i, j in bonds])
    
    # Main simulation loop
    for itr in range(n_iter):
        velocity *= (1-damping)
        velocity[:,:] += gravity[None,:] * dt
        pos_pred[:] = pos + velocity * dt
        
        # Apply fixed point constraints
        if fixed_points is not None:
            pos_pred[fixed_points] = points[fixed_points]
        
        # Build right-hand side
        b = make_pd_rhs(neighbs, bonds, masses, dt, ks, pos, l0s, pos_pred)
        
        # Solve system for each coordinate using Jacobi iteration
        for i in range(3):
            x = pos_pred[:, i].copy()  # Use predicted position as initial guess
            for j in range(inner_iter):
                x = jacobi_iteration(A, b[:, i], x, fixed_points)
            pos_cor[:, i] = x
        
        # Apply fixed point constraints
        if fixed_points is not None:
            pos_cor[fixed_points] = points[fixed_points]
        
        # Update velocity and position
        velocity = update_velocity(pos_cor, pos, velocity, dt)
        pos[:] = pos_cor[:]
        
        if call_back is not None:
            call_back(pos)
        
        # Print debug info
        v_norm = np.linalg.norm(velocity)
        pos_diff = np.linalg.norm(pos - pos_pred)
        print(f"iter:{itr} dt={dt:.3e} |v|={v_norm:.3e} dp={pos_diff:.3e}")
    
    return pos, velocity

def solve_pd_chebyshev(points, velocity, bonds, masses, ks, dt=0.1, n_iter=100, 
                      inner_iter=50, gravity=np.array([0, -9.81, 0]), 
                      fixed_points=None, call_back=None, damping=0.01):
    """Solve the system using projective dynamics with Chebyshev-accelerated Jacobi iterations"""
    n_points = len(points)
    neighbs = build_neighbor_list(bonds, n_points)
    A, Mt = make_pd_matrix(neighbs, bonds, masses, dt, ks)
    
    # Initialize
    pos = points.copy()
    pos_cor = points.copy()
    pos_pred = points.copy()
    l0s = np.array([np.linalg.norm(points[j] - points[i]) for i, j in bonds])
    
    # Estimate spectral radius for Chebyshev acceleration
    rho = estimate_spectral_radius(A, np.random.rand(n_points))
    
    # Main simulation loop
    for itr in range(n_iter):
        velocity *= (1-damping)
        velocity[:,:] += gravity[None,:] * dt
        pos_pred[:] = pos + velocity * dt
        
        # Apply fixed point constraints
        if fixed_points is not None:
            pos_pred[fixed_points] = points[fixed_points]
        
        # Build right-hand side
        b = make_pd_rhs(neighbs, bonds, masses, dt, ks, pos, l0s, pos_pred)
        
        # Solve system for each coordinate using Chebyshev-accelerated Jacobi
        for i in range(3):
            x = pos_pred[:, i].copy()  # Current solution
            x_prev = x.copy()          # Previous solution
            omega = 1.0                # Initial omega
            
            for k in range(inner_iter):
                # Perform Jacobi iteration
                x_new = jacobi_iteration(A, b[:, i], x, fixed_points)
                
                # Apply Chebyshev acceleration
                omega_next = calculate_omega(k, rho, omega)
                x_accel = omega_next * (x_new - x) + x
                
                # Update for next iteration
                x_prev = x.copy()
                x = x_accel
                omega = omega_next
            
            pos_cor[:, i] = x
        
        # Apply fixed point constraints
        if fixed_points is not None:
            pos_cor[fixed_points] = points[fixed_points]
        
        # Update velocity and position
        velocity = update_velocity(pos_cor, pos, velocity, dt)
        pos[:] = pos_cor[:]
        
        if call_back is not None:
            call_back(pos)
        
        # Print debug info
        v_norm = np.linalg.norm(velocity)
        pos_diff = np.linalg.norm(pos - pos_pred)
        print(f"iter:{itr} dt={dt:.3e} |v|={v_norm:.3e} dp={pos_diff:.3e}")
    
    return pos, velocity

# Example usage
if __name__ == "__main__":
    from projective_dynamics import build_grid_2d
    
    # Create a simple 5x5 grid
    nx, ny = 5, 5
    bonds, points, masses, ks, fixed = build_grid_2d(nx, ny)
    velocity = np.zeros_like(points)
    
    # Run simulation with Jacobi solver
    print("Running with plain Jacobi solver...")
    new_points_jacobi, new_velocity_jacobi = solve_pd_jacobi(
        points, velocity.copy(), bonds, masses, ks, fixed_points=fixed
    )
    
    # Run simulation with Chebyshev-accelerated solver
    print("\nRunning with Chebyshev-accelerated Jacobi solver...")
    new_points_cheby, new_velocity_cheby = solve_pd_chebyshev(
        points, velocity.copy(), bonds, masses, ks, fixed_points=fixed
    )
