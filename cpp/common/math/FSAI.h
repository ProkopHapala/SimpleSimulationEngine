
#ifndef  FSAI_h
#define  FSAI_h
#include <vector>
#include <algorithm>
#include <cmath>

/*

Factored Sparse Approximate Inverse (FSAI) algorithm. (done by Claude 3.5 Sonet)

fsai function:
Iterates through each row of the matrix using your neighbor list structure.
Determines the sparsity pattern for each row using the indexed array.
Sets up and solves a local dense problem using Gaussian elimination.
Updates the G matrix (stored as a dense matrix for simplicity) with the solutio

The Least Squares Problem in FSAI:

In the FSAI method, we're trying to find a sparse lower triangular matrix G such that G*G^T approximates A^(-1). This is equivalent to minimizing ||AG - I||_F under the constraint that G has a prescribed sparsity pattern.

For each row i of G, we solve a separate least squares problem:

minimize ||A_i * g_i - e_i||_2

where:
- A_i is the i-th row of A, restricted to the columns in the sparsity pattern of G
- g_i is the i-th row of G (which we're solving for)
- e_i is the i-th unit vector

Constraints:
1. Sparsity: g_i is only allowed to have non-zero entries where the sparsity pattern of G allows.
2. Lower triangular: We only consider columns j ≤ i.

This formulation indeed aims to minimize the error in (G*G^T * A) = I, but we do it row by row, and only for the prescribed sparsity pattern.

Conjugate Gradient Normal Equation (CGNE):

CGNE is an iterative method for solving least squares problems, especially useful for large, sparse systems. Here's how it fits into our problem:

1. Our least squares problem for row i can be written as:
   minimize ||A_i * g_i - e_i||_2

2. The normal equations for this problem are:
   (A_i^T * A_i) * g_i = A_i^T * e_i

3. CGNE solves these normal equations iteratively, without explicitly forming A_i^T * A_i.

The algorithm works as follows:

1. Initialize g_i = 0, r = e_i, p = A_i^T * r
2. Repeat until convergence:
   a. α = (r^T * r) / (p^T * A_i^T * A_i * p)
   b. g_i = g_i + α * p
   c. r_new = r - α * A_i * A_i^T * r
   d. β = (r_new^T * r_new) / (r^T * r)
   e. p = A_i^T * r_new + β * p
   f. r = r_new

This method allows us to solve the least squares problem efficiently while maintaining sparsity, as we only update the entries of g_i that correspond to the chosen sparsity pattern.

Benefits of this approach:

1. Sparsity Preservation: We only compute entries of G where we allow non-zeros, maintaining the desired sparsity pattern.

2. Flexibility: We can control the trade-off between accuracy and sparsity by adjusting the sparsity pattern and the number of CGNE iterations.

3. Scalability: This method is well-suited for large, sparse matrices as it doesn't require explicit formation of A^T * A.

4. Approximation Quality: By solving these least squares problems, we're finding the best possible approximation of A^(-1) given the sparsity constraints.

In summary, this approach allows us to construct a sparse approximate inverse that balances accuracy and sparsity. We're essentially trying to make G*G^T as close to A^(-1) as possible, under the constraint that G maintains a prescribed sparsity pattern. The CGNE method provides an efficient way to solve these constrained least squares problems, allowing us to handle large, sparse matrices effectively.
*/

class SparseMatrix { public:
    int n;
    int* lens;
    int* i0s;
    double* values;
    int* indexed;

    void dot_dens_vector(const double* v, double* out) const {
        for (int i = 0; i < n; i++) {
            int ni = lens[i];
            int i0 = i0s[i];
            const int* indi = indexed + i0;
            const double* Ai = values + i0;
            double sum = 0;
            for (int k = 0; k < ni; k++) {
                int j = indi[k];
                sum += Ai[k] * v[j];
            }
            out[i] = sum;
        }
    }

};

std::vector<double> fsai(const SparseMatrix& A, int p) {
    int n = A.n;
    std::vector<double> G(n * n, 0.0);  // Initialize G as a dense matrix for simplicity

    for (int i = 0; i < n; ++i) {
        // Step 1: Determine the sparsity pattern for row i
        std::vector<int> pattern;
        int ni = A.lens[i];
        int i0 = A.i0s[i];
        const int* indi = A.indexed + i0;
        for (int k = 0; k < ni; ++k) {
            if (indi[k] <= i) {
                pattern.push_back(indi[k]);
            }
        }

        // Ensure we don't exceed p non-zero elements
        if (pattern.size() > p) {
            pattern.resize(p);
        }

        int pattern_size = pattern.size();

        // Step 2: Set up the local problem
        std::vector<std::vector<double>> A_local(pattern_size, std::vector<double>(pattern_size, 0.0));   // basically just square matrix A_loc[pattern_size,pattern_size]
        std::vector<double> e_i(pattern_size, 0.0);

        // Fill A_local
        for (int k = 0; k < pattern_size; ++k) {
            int row = pattern[k];
            int row_ni = A.lens[row];
            int row_i0 = A.i0s[row];
            const int* row_indi = A.indexed + row_i0;
            const double* row_values = A.values + row_i0;
            
            for (int l = 0; l < pattern_size; ++l) {
                int col = pattern[l];
                for (int m = 0; m < row_ni; ++m) {
                    if (row_indi[m] == col) {
                        A_local[k][l] = row_values[m];
                        break;
                    }
                }
            }
            if (pattern[k] == i) {
                e_i[k] = 1.0;
            }
        }

        // Step 3: Solve the local problem using Gaussian elimination
        std::vector<double> g_i = e_i;
        for (int k = 0; k < pattern_size; ++k) {
            for (int l = k + 1; l < pattern_size; ++l) {
                double factor = A_local[l][k] / A_local[k][k];
                for (int m = k; m < pattern_size; ++m) {
                    A_local[l][m] -= factor * A_local[k][m];
                }
                g_i[l] -= factor * g_i[k];
            }
        }
        for (int k = pattern_size - 1; k >= 0; --k) {
            for (int l = k + 1; l < pattern_size; ++l) {
                g_i[k] -= A_local[k][l] * g_i[l];
            }
            g_i[k] /= A_local[k][k];
        }

        // Step 4: Update G with the solution
        for (int k = 0; k < pattern_size; ++k) {
            G[i * n + pattern[k]] = g_i[k];
        }
    }

    return G;
}

// Helper function to apply the preconditioner
std::vector<double> apply_preconditioner(const std::vector<double>& G, const std::vector<double>& r, int n) {
    std::vector<double> z(n, 0.0);
    std::vector<double> temp(n, 0.0);

    // Apply G
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            temp[i] += G[i * n + j] * r[j];
        }
    }

    // Apply G^T
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            z[i] += G[j * n + i] * temp[j];
        }
    }

    return z;
}


    void CGNE( int n, double* A, double* g_i, double* b, int niter=100 ){

        // // Step 3: Solve the least squares problem using Conjugate Gradient Normal Equation (CGNE)
        // std::vector<double> g_i(pattern_size, 0.0);
        
        std::vector<double> r   (n);
        std::vector<double> p   (n);
        std::vector<double> Ap  (n);
        std::vector<double> AtAp(n);

        for (int iter=0; iter<niter; ++iter) {
            // Compute A^T * A * p
            for (int k = 0; k < n; ++k) {
                double Api = 0.0;
                for (int l = 0; l < n; ++l) {
                    Api += A[l*n+k] * p[l];
                }
                Ap[k] = Api;
            }
            for (int k = 0; k < n; ++k) {
                double AtApi=0.0;
                for (int l = 0; l < n; ++l) {
                    AtApi += A[k*n+l] * Ap[l];
                }
                AtAp[k]=AtApi;
            }

            double alpha = 0.0;
            double pAtAp = 0.0;
            for (int k = 0; k < n; ++k) {
                alpha += r[k] * r[k];
                pAtAp += p[k] * AtAp[k];
            }
            alpha /= pAtAp;

            for (int k = 0; k < n; ++k) {
                g_i[k] += alpha * p[k];
                r  [k] -= alpha * AtAp[k];
            }

            double beta = 0.0;
            for (int k = 0; k < n; ++k) {
                beta += r[k] * r[k];
            }
            if (sqrt(beta) < 1e-10) break;  // Convergence check

            beta /= alpha * pAtAp;

            for (int k = 0; k < n; ++k) {
                p[k] = r[k] + beta * p[k];
            }
        }
    }


void sparse_fsai( const SparseMatrix& A, int niter=100 ) {
    int n = A.n;
    std::vector<double> G(n*n, 0.0);  // Initialize G as a dense matrix for simplicity

    for (int i=0; i<n; ++i){ // loop over rows
        // Step 1: Determine the sparsity pattern for row i
        std::vector<int> pattern;
        int ni = A.lens[i];
        int i0 = A.i0s[i];
        const int* indi = A.indexed + i0;
        for (int k = 0; k < ni; ++k) {
            if (indi[k] <= i) {    // only upper triangular (?)
                pattern.push_back(indi[k]);
            }
        }

        //if (pattern.size() > ni_max) {  pattern.resize(ni_max); }   // Ensure we don't exceed p non-zero elements
        //int pattern_size = pattern.size();

        // Step 2: Set up the least squares problem
        std::vector<double> A_i(ni*ni,0.0);
        std::vector<double> b_i(ni, 0.0);
        std::vector<double> g_i(ni, 0.0);

        // Build local linear system Ai*bi=gi  to fit one row of approximate matrix G   
        for (int k=0; k<ni; ++k){
            int row = pattern[k];
            int row_ni = A.lens[row];
            int row_i0 = A.i0s[row];
            const int*    row_indi   = A.indexed + row_i0;
            const double* row_values = A.values  + row_i0;
            for (int l=0; l<ni; ++l){
                int col = pattern[l];
                for (int m = 0; m < row_ni; ++m){
                    if (row_indi[m] == col){
                        A_i[ k*ni + l ] = row_values[m];
                        break;
                    }
                }
            }
            if(row==i){  b_i[k]=1.0;  }  // basis vectors like b4 = {0,0,0,1,0,0}
        }

        
        CGNE( ni, A_i.data(), g_i.data(), b_i.data(), niter );

        for (int k=0; k<ni; ++k) {
            int j = pattern[k];
            G[i*n+j] = g_i[k];
        }

    }  // for i     // loop over rows
}

std::vector<double> sparse_preserving_fsai(const SparseMatrix& A, int p) {
    int n = A.n;
    std::vector<double> G(n * n, 0.0);  // Initialize G as a dense matrix for simplicity

    for (int i = 0; i < n; ++i) {
        // Step 1: Determine the sparsity pattern for row i
        std::vector<int> pattern;
        int ni = A.lens[i];
        int i0 = A.i0s[i];
        const int* indi = A.indexed + i0;
        for (int k = 0; k < ni; ++k) {
            if (indi[k] <= i) {
                pattern.push_back(indi[k]);
            }
        }

        // Ensure we don't exceed p non-zero elements
        if (pattern.size() > p) {
            pattern.resize(p);
        }

        int pattern_size = pattern.size();

        // Step 2: Set up the least squares problem
        std::vector<std::vector<double>> A_local(pattern_size, std::vector<double>(pattern_size, 0.0));
        std::vector<double> b(pattern_size, 0.0);

        // Fill A_local and b
        for (int k = 0; k < pattern_size; ++k) {
            int row = pattern[k];
            int row_ni = A.lens[row];
            int row_i0 = A.i0s[row];
            const int* row_indi = A.indexed + row_i0;
            const double* row_values = A.values + row_i0;
            
            for (int l = 0; l < pattern_size; ++l) {
                int col = pattern[l];
                for (int m = 0; m < row_ni; ++m) {
                    if (row_indi[m] == col) {
                        A_local[k][l] = row_values[m];
                        break;
                    }
                }
            }
            if (row == i) {
                b[k] = 1.0;
            }
        }


        // Step 3: Solve the least squares problem using Conjugate Gradient Normal Equation (CGNE)
        std::vector<double> g_i(pattern_size, 0.0);
        std::vector<double> r = b;
        std::vector<double> p = r;
        std::vector<double> Ap(pattern_size);

        for (int iter = 0; iter < 100; ++iter) {  // Max 100 iterations for CGNE
            // Compute A^T * A * p
            for (int k = 0; k < pattern_size; ++k) {
                Ap[k] = 0.0;
                for (int l = 0; l < pattern_size; ++l) {
                    Ap[k] += A_local[l][k] * p[l];
                }
            }
            std::vector<double> AtAp(pattern_size, 0.0);
            for (int k = 0; k < pattern_size; ++k) {
                for (int l = 0; l < pattern_size; ++l) {
                    AtAp[k] += A_local[k][l] * Ap[l];
                }
            }

            double alpha = 0.0, pAtAp = 0.0;
            for (int k = 0; k < pattern_size; ++k) {
                alpha += r[k] * r[k];
                pAtAp += p[k] * AtAp[k];
            }
            alpha /= pAtAp;

            for (int k = 0; k < pattern_size; ++k) {
                g_i[k] += alpha * p[k];
                r[k] -= alpha * AtAp[k];
            }

            double beta = 0.0;
            for (int k = 0; k < pattern_size; ++k) {
                beta += r[k] * r[k];
            }
            if (sqrt(beta) < 1e-10) break;  // Convergence check

            beta /= alpha * pAtAp;

            for (int k = 0; k < pattern_size; ++k) {
                p[k] = r[k] + beta * p[k];
            }
        }

        // Step 4: Update G with the solution
        for (int k = 0; k < pattern_size; ++k) {
            G[i * n + pattern[k]] = g_i[k];
        }
    }

    return G;
}


// Example usage in a simplified conjugate gradient method
std::vector<double> conjugate_gradient(const SparseMatrix& A, const std::vector<double>& b, 
                                       const std::vector<double>& G, int max_iter, double tol) {
    int n = A.n;
    std::vector<double> x(n, 0.0);
    std::vector<double> r = b;  // Assuming initial x is zero
    std::vector<double> z = apply_preconditioner(G, r, n);
    std::vector<double> p = z;
    std::vector<double> Ap(n);

    for (int iter = 0; iter < max_iter; ++iter) {
        // Compute Ap
        A.dot_dens_vector(p.data(), Ap.data());

        // Compute alpha
        double r_dot_z = 0.0, p_dot_Ap = 0.0;
        for (int i = 0; i < n; ++i) {
            r_dot_z += r[i] * z[i];
            p_dot_Ap += p[i] * Ap[i];
        }
        double alpha = r_dot_z / p_dot_Ap;

        // Update x and r
        for (int i = 0; i < n; ++i) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        // Check for convergence
        double r_norm = 0.0;
        for (int i = 0; i < n; ++i) {
            r_norm += r[i] * r[i];
        }
        if (std::sqrt(r_norm) < tol) {
            break;
        }

        // Compute new z
        std::vector<double> z_new = apply_preconditioner(G, r, n);

        // Compute beta
        double r_dot_z_new = 0.0;
        for (int i = 0; i < n; ++i) {
            r_dot_z_new += r[i] * z_new[i];
        }
        double beta = r_dot_z_new / r_dot_z;

        // Update p
        for (int i = 0; i < n; ++i) {
            p[i] = z_new[i] + beta * p[i];
        }

        z = z_new;
    }

    return x;
}

#endif

