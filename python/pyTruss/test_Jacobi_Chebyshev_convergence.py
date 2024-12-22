import numpy as np
import matplotlib.pyplot as plt
from projective_dynamics import solve_pd, build_grid_2d, make_pd_matrix, make_pd_rhs
from projective_dynamics_iterative import jacobi_iteration, calculate_omega, estimate_spectral_radius

def solve_iterative(A, b, x0, n_iter, bChebyshev=False, bInertial=False, beta=0.5  ):
    """Solve system using Chebyshev-accelerated Jacobi iteration"""
    
    errors = []
    x_direct = np.linalg.solve(A, b)
    
    if bChebyshev:
        rho = estimate_spectral_radius(A, np.random.rand(len(x0)))
    omega_old = 1.0
    
    x_old  = x0.copy()
    x      = x0.copy()
    v_old  = x*0.0 
    for itr in range(n_iter):
        x_pred  = jacobi_iteration(A, b, x)
        
        if bChebyshev:
            omega   = calculate_omega(itr, rho, omega_old )
            x_new = omega * (x_pred - x) + x_old
            x_old = x
            x     = x_new
            omega_old = omega
        elif bInertial:
            v_new = x_pred - x
            v = beta*v_new + (1-beta)*v_old
            x += v
            v_old = v  
        else:
            x = x_pred
        
        error = np.linalg.norm(x - x_direct) / np.linalg.norm(x_direct)
        errors.append(error)
        print(  f" {itr} omega: {omega_old} error: {error} " )
    
    return x, errors

def test_convergence():
    # Create test problem
    nx, ny = 10, 10
    bonds, points, masses, ks, fixed = build_grid_2d(nx, ny)
    dt = 0.1

    # Get RHS for one coordinate (x-coordinate)
    l0s = np.array([np.linalg.norm(points[j] - points[i]) for i, j in bonds])
    l0s += l0s*np.random.randn(len(bonds))*0.3
    ks  += ks*np.random.randn(len(bonds))*5.0

    
    # Build system matrix and RHS
    neighbs = [[] for _ in range(len(points))]
    for i, (i_, j_) in enumerate(bonds):
        neighbs[i_].append(i)
        neighbs[j_].append(i)
    
    A, Mt = make_pd_matrix(neighbs, bonds, masses, dt, ks)
    
    # Create a random initial state and velocity for testing
    velocity = np.random.randn(len(points), 3) * 0.1
    points_pred = points + velocity * dt
    

    b   = make_pd_rhs(neighbs, bonds, masses, dt, ks, points, l0s, points_pred)
    
    # Test both methods for x-coordinate
    n_iter = 20
    x0 = points_pred[:, 0].copy()  # Initial guess

    # add some random 
    x0*=1.1
    x0 += np.random.randn(len(x0))*5
    
    # Solve using both methods
    _, errors_jacobi = solve_iterative(A, b[:, 0], x0, n_iter, bChebyshev=False)
    #_, errors_cheby  = solve_iterative(A, b[:, 0], x0, n_iter, bChebyshev=True)
    _, errors_cheby  = solve_iterative(A, b[:, 0], x0, n_iter, bInertial=True, beta=0.985 )

    
    # Plot results
    plt.figure(figsize=(10, 6))
    plt.semilogy(range(1, n_iter + 1), errors_jacobi, 'b-', label='Jacobi')
    plt.semilogy(range(1, n_iter + 1), errors_cheby, 'r-', label='Chebyshev')
    plt.grid(True)
    plt.xlabel('Iteration')
    plt.ylabel('Relative Error')
    plt.title('Convergence Comparison: Jacobi vs Chebyshev-accelerated Jacobi')
    plt.legend()
    plt.show()
    #plt.savefig('convergence_comparison.png')
    #plt.close()

if __name__ == "__main__":
    test_convergence()
    print("Test completed. Results saved in 'convergence_comparison.png'")
