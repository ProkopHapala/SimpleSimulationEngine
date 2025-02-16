# Lingebra.h

This file provides a collection of functions and a class (`LinSolver`) for performing linear algebra operations, including matrix manipulation, linear system solving, and eigenvalue decomposition. It's designed to provide efficient and reusable tools for solving linear algebra problems in various applications.

## Includes

- `<math.h>`: Provides standard mathematical functions, such as `sqrt` and `fabs`.
- `<cstdlib>`: Provides general utility functions, including memory allocation and random number generation.
- `<stdio.h>`: Provides standard input/output functions, such as `printf`.
- `fastmath.h`: Provides optimized mathematical functions.
- `VecN.h`: Defines the `VecN` class for N-dimensional vector operations.
- `Lingebra.h`: (typo in original template) Self-inclusion, likely unintentional.

## Free functions

### namespace `Lingebra`

The `Lingebra` namespace encapsulates a set of free functions for performing various linear algebra operations. It provides tools for matrix creation, manipulation, linear system solving, and eigenvalue decomposition.

- `from_continuous` - Creates a double pointer matrix from a continuous memory block.
- `new_matrix` - Allocates memory for a new double pointer matrix.
- `delete_matrix` - Deallocates memory for a double pointer matrix.
- `transpose` - Transposes a matrix.
- `symCopy` - Copies the upper or lower triangle of a symmetric matrix to the other triangle.
- `dot` - Performs a matrix-vector multiplication.
- `dotT` - Performs a matrix-vector multiplication with the transpose of the matrix.
- `mmul_ik_kj` - Performs matrix-matrix multiplication using the i-k-j loop order.
- `mmul_ik_jk` - Performs matrix-matrix multiplication using the i-k-jk loop order.
- `mmul_ki_kj` - Performs matrix-matrix multiplication using the k-i-kj loop order.
- `mmul_ki_jk` - Performs matrix-matrix multiplication using the k-i-jk loop order.
- `random_matrix` - Generates a matrix with random values within a specified range.
- `print_matrix` - Prints a matrix to the console.
- `makeQuadricFormMatrix` - Creates a matrix representing a quadric form from a set of directions and coefficients.
- `evalQudraticForm` - Evaluates a quadratic form given a matrix and a vector.
- `evalQudraticFormDirs` - Evaluates a quadratic form given a set of directions, coefficients, and a vector.
- `GaussElimination` - Performs Gaussian elimination with partial pivoting on a matrix.
- `linSolve_gauss` - Solves a linear system of equations using Gaussian elimination.
- `linSolve_CG` - Solves a linear system of equations using the Conjugate Gradient method.
- `linSolve_BCG` - Solves a linear system of equations using the Bi-Conjugate Gradient method.
- `leastSquareFit_Gauss` - Solves a least squares fitting problem using Gaussian elimination.
- `get_diag_vector` - Extracts the diagonal elements of a matrix into a vector.
- `fill_identity` - Fills a matrix with the identity matrix.
- `jacobi_rot` - Performs a Jacobi rotation on two elements of a matrix.
- `jacobi_rot_cs` - Performs a Jacobi rotation on two elements of a matrix, given cosine and sine.
- `jacobi_rotation` - Performs a Jacobi rotation on a symmetric matrix to reduce off-diagonal elements.
- `jacobi_rotation_small` - Performs a Jacobi rotation on a symmetric matrix, optimized for small matrices.
- `getMaxIndex` - Returns the index of the maximum element in a vector, using a provided function to evaluate the elements.
- `updateRowMax` - Updates the maximum off-diagonal element in each row of a matrix.
- `printmatrix` - Prints a matrix to the console using a specified format.
- `eig_Jacobi_init` - Initializes the Jacobi eigenvalue algorithm.
- `eig_Jacobi_step` - Performs a single step of the Jacobi eigenvalue algorithm.
- `eig_Jacobi` - Computes the eigenvalues and eigenvectors of a symmetric matrix using the Jacobi method.
- `aproxOrthoNormStep` - Applies an approximate orthonormalization step to a set of vectors.
- `orthtoForce` - Applies a force to a set of vectors to make them more orthogonal.

---
## Types (classes and structs)
---

### class `LinSolver`

The `LinSolver` class provides an abstract base class for solving linear systems of equations. It defines a virtual `dotFunc` method that must be implemented by derived classes to perform the matrix-vector multiplication, allowing for different matrix storage formats and multiplication algorithms.

#### properties

- `n`:`int`:`public:` - The size of the linear system (number of equations and unknowns).
- `istep`:`int` - The current iteration step in the iterative solver.
- `r`:`double *` - A pointer to an array storing the residual vector.
- `r2`:`double *` - A pointer to an array storing a temporary residual vector.
- `p`:`double *` - A pointer to an array storing the search direction vector.
- `Ap`:`double *` - A pointer to an array storing the matrix-vector product of the matrix and the search direction.
- `rho`:`double` - A scalar value used in the Conjugate Gradient method.
- `alpha`:`double` - A scalar value used in the Conjugate Gradient method.
- `x`:`double*` - A pointer to an array storing the solution vector.
- `b`:`double*` - A pointer to an array storing the right-hand side vector.
- `M`:`double*` - A pointer to the matrix data (optional, may be managed externally).

#### methods

- `realloc` - Reallocates memory for the internal arrays if the size of the linear system has changed.
- `dealloc` - Deallocates memory for the internal arrays.
- `setLinearProblem` - Sets the linear system to be solved, providing pointers to the solution vector, right-hand side vector, and matrix data.
- `step_GD` - Performs a single step of the Gradient Descent method.
- `step_CG` - Performs a single step of the Conjugate Gradient method.
- `step0_CG` - Performs the initial setup for the Conjugate Gradient method.
- `step_CG_simple` - Performs a simplified step of the Conjugate Gradient method.
- `solve_CG` - Solves the linear system using the Conjugate Gradient method for a specified number of iterations or until a specified error tolerance is reached.
- `solve_CG_` - Solves the linear system using the Conjugate Gradient method, with a slightly different implementation.