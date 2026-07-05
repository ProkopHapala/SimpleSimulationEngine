# === AUTO-DOC BEGIN ===
"""
@brief Convergence test for Jacobi vs Chebyshev/inertial on a Projective Dynamics system.

Builds a PD system matrix from a 10x10 grid truss with randomized rest lengths and stiffness,
then solves one coordinate using plain Jacobi, Chebyshev-accelerated Jacobi, and inertial
(heavy-ball momentum) Jacobi. Compares relative error vs direct solve. Uses shared
`jacobi_iteration` from `sparse.py` and `estimate_spectral_radius` / `chebyshev_omega_update`
from `IterativeLinearSolvers.py`. The local `solve_iterative` is a thin wrapper that handles
the acceleration policy selection.
"""
# === AUTO-DOC END ===

import numpy as np
import matplotlib.pyplot as plt
import argparse
from projective_dynamics import make_pd_matrix, make_pd_rhs
from sparse import build_neighbor_list, jacobi_iteration
from IterativeLinearSolvers import estimate_spectral_radius, chebyshev_omega_update, momentum_accelerator
from truss import Truss

def solve_iterative(A, b, x0, n_iter, bChebyshev=False, bInertial=False, beta=0.5):
    """Solve system using Jacobi with optional Chebyshev or momentum acceleration."""
    errors = []
    x_direct = np.linalg.solve(A, b)
    
    if bChebyshev:
        rho = estimate_spectral_radius(A)
    omega_old = 1.0
    
    x_old  = x0.copy()
    x      = x0.copy()
    v_old  = x*0.0 
    for itr in range(n_iter):
        x_pred, _ = jacobi_iteration(A, b, x)
        
        if bChebyshev:
            omega   = chebyshev_omega_update(itr, rho, S=10, prev_omega=omega_old)
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

def test_convergence(args):
    # Create test problem
    nx, ny = 10, 10
    truss = Truss()
    truss.build_grid_2d(nx=nx, ny=ny, m=1.0, m_end=1000.0, l=1.0, k=10000.0, k_diag=1000.0)
    bonds, points, masses, ks, fixed, l0s, neighbs = truss.get_pd_quantities()
    dt = 0.1

    l0s = l0s + l0s*np.random.randn(len(bonds))*0.3
    ks   = ks + ks*np.random.randn(len(bonds))*5.0
    
    # Build system matrix and RHS
    neighbs = build_neighbor_list(bonds, len(points))
    A = make_pd_matrix(neighbs, bonds, masses, dt, ks)
    
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
    if not args.noshow: plt.show()
    if args.savefig: plt.savefig(args.savefig, dpi=150)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Test Jacobi vs Chebyshev/inertial convergence.")
    parser.add_argument("--noshow", action="store_true", help="Skip plt.show() for headless execution.")
    parser.add_argument("--savefig", type=str, default="", help="Path to save the plot.")
    args = parser.parse_args()
    test_convergence(args)
    print("Test completed.")
