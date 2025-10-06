import numpy as np
import argparse

from IterativeLinearSolvers import get_iteration_matrix, solve_jacobi, solve_jacobi_chebyshev, solve_gauss_seidel, solve_gauss_seidel_chebyshev

# unlimited line lengh for numpy array print
np.set_printoptions(linewidth=np.inf)

def generate_sdd_matrix(size, diag_dominance_factor=1.1):
    """
    Generates a Symmetric, Diagonally Dominant (SDD) matrix.
    Such matrices are guaranteed to be positive definite, and the Jacobi method
    is guaranteed to converge for them.

    Args:
        size (int): The dimension of the square matrix.
        diag_dominance_factor (float): How much larger the diagonal elements are
            than the sum of the absolute values of the other elements in their row.
            Must be > 1.0 for strict diagonal dominance.

    Returns:
        numpy.ndarray: The generated SDD matrix.
    """
    # Start with a random symmetric matrix
    A = np.random.rand(size, size)
    A = (A + A.T) / 2.0
    np.fill_diagonal(A, 0) # Zero out the diagonal for now

    # Calculate the sum of absolute values for each row
    row_sum = np.sum(np.abs(A), axis=1)

    # Set the diagonal elements to be dominant
    diagonal = row_sum * diag_dominance_factor
    np.fill_diagonal(A, diagonal)

    return A


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare Jacobi and Chebyshev-accelerated Jacobi solvers.")
    parser.add_argument('--size',           type=int,   default=10,   help='Size of the matrix to generate.')
    parser.add_argument('--diag-dominance', type=float, default=1.1,  help='Diagonal dominance factor (>1). Smaller values make convergence harder.')
    parser.add_argument('--seed',           type=int,   default=None, help='Random seed for matrix generation.')
    parser.add_argument('--method',         type=str,   default='jacobi', choices=['jacobi', 'gauss-seidel'], help='Base iterative method to use.')
    parser.add_argument('--load-matrix',    type=str,   default=None, help='Path to a .npy file to load the matrix A from.')
    parser.add_argument('--print-matrix',   type=int,   default=1,    help='Print the matrix A before solving.')
    parser.add_argument('--max-iters',      type=int,   default=15,   help='Maximum number of iterations for the solvers.')
    parser.add_argument('--cheby-delay',    type=int,   default=2,    help='Delayed start (S) for Chebyshev acceleration.')
    args = parser.parse_args()

    # --- 1. System Setup: Generate or load the matrix A ---
    if args.load_matrix:
        print(f"Loading matrix from {args.load_matrix}...")
        A = np.load(args.load_matrix)
        if A.ndim != 2 or A.shape[0] != A.shape[1]:
            raise ValueError("Loaded matrix must be a square 2D array.")
        matrix_size = A.shape[0]
        print(f"Loaded a {matrix_size}x{matrix_size} matrix.")
    else:
        if args.seed is not None:
            np.random.seed(args.seed)
        matrix_size = args.size
        A = generate_sdd_matrix(matrix_size, diag_dominance_factor=args.diag_dominance)
        print(f"Generated a {matrix_size}x{matrix_size} matrix with diagonal dominance factor = {args.diag_dominance}")

    if args.print_matrix:
        print("Matrix A:\n", A)

    b = np.ones(matrix_size)

    # --- 2. Analyze the Iteration Matrix for the chosen method ---
    T = get_iteration_matrix(A, method=args.method)
    eigenvalues = np.linalg.eigvals(T)
    spectral_radius_actual = np.max(np.abs(eigenvalues))

    print("="*50)
    print(f"ACTUAL Spectral Radius of {args.method.title()} Iteration Matrix: {spectral_radius_actual:.6f}")
    print("This is the value we should use for 'rho' in Chebyshev.")
    print("="*50)

    # --- 3. Run and Compare the Solvers ---
    # We will use the true spectral radius for optimal performance
    rho_param = spectral_radius_actual
    
    if args.method == 'jacobi':
        # Jacobi without acceleration
        res_base, err_base = solve_jacobi(A, b, max_iters=args.max_iters)
        # Jacobi with Chebyshev acceleration
        res_cheby, err_cheby = solve_jacobi_chebyshev(A, b, rho=rho_param, max_iters=args.max_iters, S=args.cheby_delay)
    else: # gauss-seidel
        # Gauss-Seidel without acceleration
        res_base, err_base = solve_gauss_seidel(A, b, max_iters=args.max_iters)
        # Gauss-Seidel with Chebyshev acceleration
        res_cheby, err_cheby = solve_gauss_seidel_chebyshev(A, b, rho=rho_param, max_iters=args.max_iters, S=args.cheby_delay)

    # --- 4. Plotting the results ---
    try:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(10, 5))
        plt.semilogy(err_base, 'o-', label=f'Standard {args.method.title()} Error')
        plt.semilogy(err_cheby, 's-', label=f'{args.method.title()} + Chebyshev Error')
        plt.title('Convergence Comparison')
        plt.xlabel('Iteration')
        plt.ylabel('Error (log scale)')
        plt.grid(True, which='both', linestyle='--')
        plt.legend()
        plt.show()
    except ImportError:
        print("\nMatplotlib not found. Please install it (`pip install matplotlib`) to see the convergence plot.")