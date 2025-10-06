import numpy as np

def get_iteration_matrix(A, method='jacobi'):
    """
    Calculates the iteration matrix T for a given method (Jacobi or Gauss-Seidel).
    The update rule is x_k+1 = T * x_k + c.
    Convergence depends on the spectral radius of T.
    """
    D = np.diag(np.diag(A))
    L = np.tril(A, k=-1)
    U = np.triu(A, k=1)

    if method.lower() == 'jacobi':
        # T_J = -D^-1 * (L + U)
        D_inv = np.linalg.inv(D)
        T     = -np.dot(D_inv, L + U)
    elif method.lower() == 'gauss-seidel':
        # T_GS = -(D + L)^-1 * U
        D_plus_L_inv = np.linalg.inv(D + L)
        T            = -np.dot(D_plus_L_inv, U)
    else:
        raise ValueError(f"Unknown method: {method}. Choose 'jacobi' or 'gauss-seidel'.")

    return T

# ====================================
#  Modular Solver Components
# ====================================

def jacobi_step(A, b, x):
    """Performs a single Jacobi iteration step."""
    D_inv = 1.0 / np.diag(A)
    R = A - np.diag(np.diag(A))
    x_new = (b - np.dot(R, x)) * D_inv
    return x_new

def gauss_seidel_step(A, b, x):
    """Performs a single Gauss-Seidel iteration step."""
    x_new = x.copy()
    for i in range(A.shape[0]):
        s1 = np.dot(A[i, :i], x_new[:i])
        s2 = np.dot(A[i, i + 1:], x[i + 1:])
        x_new[i] = (b[i] - s1 - s2) / A[i, i]
    return x_new

def chebyshev_accelerator(x_tilde_k1, x_k, x_k_minus_1, omega):
    """Applies the Chebyshev acceleration formula."""
    return x_k_minus_1 + omega * (x_tilde_k1 - x_k_minus_1)

def momentum_accelerator(x_tilde_k1, x_k, x_k_minus_1, bmix):
    """Applies the Heavy Ball momentum acceleration formula."""
    inertia = x_k - x_k_minus_1
    return x_tilde_k1 + bmix * inertia

def chebyshev_omega_update(k, rho, S, prev_omega, **kwargs):
    """
    Calculates the dynamic omega parameter for Chebyshev acceleration.

    Args:
        k (int): Current iteration number.
        rho (float): Spectral radius of the base iteration matrix.
        S (int): Delayed start iteration.
        prev_omega (float): The omega value from the previous iteration (k-1).

    Returns:
        float: The new omega value for iteration k.
    """
    if k < S: return 1.0
    if k == S: return 2.0 / (2.0 - rho**2) if rho**2 < 2.0 else 1.0
    return 4.0 / (4.0 - rho**2 * prev_omega) if (4.0 - rho**2 * prev_omega) != 0 else 1.0

def momentum_bmix_update(k, **kwargs):
    """
    Calculates the b_mix parameter for momentum acceleration.
    This follows the schedule from the C++ implementation.
    """
    niter = kwargs.get('max_iters', 20)
    istart = kwargs.get('b_istart', 3)
    iend = kwargs.get('b_iend', niter - 1)
    b_start = kwargs.get('b_start', 0.2)
    b_end = kwargs.get('b_end', 0.2)
    b_last = kwargs.get('b_last', 0.0)

    if k == 0 or k >= niter - 1:
        return b_last
    if k < istart:
        return 0.0
    if k >= iend:
        return b_end
    
    span = iend - istart
    if span <= 0:
        return b_end
    t = (k - istart) / span
    return b_start + t * (b_end - b_start)

def solve_iterative(A, b, step_func, accelerator_func=None, omega_update_func=None, max_iters=20, tol=1e-9, **kwargs):
    """
    A general iterative solver framework.

    Args:
        step_func (callable): A function for a single solver step, e.g., jacobi_step.
        accelerator_func (callable, optional): An acceleration function, e.g., chebyshev_accelerator.
        omega_update_func (callable, optional): A function to update the omega parameter dynamically.
    """
    size = A.shape[0]
    x_k = np.zeros(size)
    x_k_minus_1 = np.zeros(size)
    x_true = np.linalg.solve(A, b)

    residuals = []
    errors = []

    # --- Accelerator specific setup ---
    rho = kwargs.get('rho', 0.0)
    S = kwargs.get('S', 0)
    kwargs['max_iters'] = max_iters # Pass max_iters to omega_update_func
    omega = 1.0

    # --- Print Header ---
    header = f"{'Iter':<5} {'Residual':<15} {'Error':<15}"
    if accelerator_func:
        header += f" {'omega_k':<10}"
    print(header)
    print("-" * len(header))

    for k in range(max_iters):
        # 1. Perform the base solver step
        x_tilde_k1 = step_func(A, b, x_k)

        # 2. Apply acceleration if provided
        if accelerator_func:
            # Create a clean kwargs dict for the update function to avoid passing duplicate arguments
            update_kwargs = kwargs.copy()
            update_kwargs.pop('rho', None)
            update_kwargs.pop('S', None)
            omega = omega_update_func(k, rho=rho, S=S, prev_omega=omega, **update_kwargs) if omega_update_func else 1.0
            x_k1 = accelerator_func(x_tilde_k1, x_k, x_k_minus_1, omega)
        else:
            x_k1 = x_tilde_k1

        # 3. Update state for the next iteration
        x_k_minus_1 = x_k
        x_k = x_k1

        # 4. Logging and convergence check
        residual = np.linalg.norm(b - A @ x_k)
        error = np.linalg.norm(x_k - x_true)
        residuals.append(residual)
        errors.append(error)

        log_line = f"{k:<5} {residual:<15.6e} {error:<15.6e}"
        if accelerator_func:
            log_line += f" {omega:<10.4f}"
        print(log_line)

        if residual < tol:
            print("Converged!")
            break

    return residuals, errors

# ====================================
#  User-Facing Wrapper Functions
# ====================================

def solve_jacobi(A, b, max_iters=20, tol=1e-9):
    """Solves Ax=b using the standard Jacobi method."""
    print("\n--- Running Standard Jacobi ---")
    return solve_iterative(A, b, jacobi_step, max_iters=max_iters, tol=tol)

def solve_gauss_seidel(A, b, max_iters=20, tol=1e-9):
    """Solves Ax=b using the standard Gauss-Seidel method."""
    print("\n--- Running Standard Gauss-Seidel ---")
    return solve_iterative(A, b, gauss_seidel_step, max_iters=max_iters, tol=tol)

def solve_jacobi_chebyshev(A, b, rho, max_iters=20, S=0, tol=1e-9):
    """Solves Ax=b using Jacobi accelerated with the Chebyshev method."""
    print(f"\n--- Running Jacobi with Chebyshev Acceleration (rho={rho:.4f}, S={S}) ---")
    return solve_iterative(A, b, jacobi_step, chebyshev_accelerator, chebyshev_omega_update, max_iters=max_iters, tol=tol, rho=rho, S=S)

def solve_gauss_seidel_chebyshev(A, b, rho, max_iters=20, S=0, tol=1e-9):
    """Solves Ax=b using Gauss-Seidel accelerated with the Chebyshev method."""
    print(f"\n--- Running Gauss-Seidel with Chebyshev Acceleration (rho={rho:.4f}, S={S}) ---")
    return solve_iterative(A, b, gauss_seidel_step, chebyshev_accelerator, chebyshev_omega_update, max_iters=max_iters, tol=tol, rho=rho, S=S)

def solve_jacobi_momentum(A, b, max_iters=20, tol=1e-9, **kwargs):
    """Solves Ax=b using Jacobi accelerated with the Heavy Ball momentum method."""
    print(f"\n--- Running Jacobi with Momentum Acceleration ---")
    return solve_iterative(A, b, jacobi_step, momentum_accelerator, momentum_bmix_update, max_iters=max_iters, tol=tol, **kwargs)

def solve_gauss_seidel_momentum(A, b, max_iters=20, tol=1e-9, **kwargs):
    """Solves Ax=b using Gauss-Seidel accelerated with the Heavy Ball momentum method."""
    print(f"\n--- Running Gauss-Seidel with Momentum Acceleration ---")
    return solve_iterative(A, b, gauss_seidel_step, momentum_accelerator, momentum_bmix_update, max_iters=max_iters, tol=tol, **kwargs)