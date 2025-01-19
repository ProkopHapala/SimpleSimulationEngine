#pragma once

//#ifndef  CGNE_h
//#define  CGNE_h

#include <vector>
#include <algorithm>
#include <cmath>

/*

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

*/

inline void CGNE(int n, double* A, double* g_i, double* b, int niter=100) {
    std::vector<double> r(n);
    std::vector<double> p(n);
    std::vector<double> Ap(n);

    // Initialize r and p
    for (int k = 0; k < n; ++k) {
        r[k] = b[k];
        for (int l = 0; l < n; ++l) {
            r[k] -= A[k*n + l] * g_i[l];
        }
        p[k] = r[k];
    }

    double lr2 = 0.0;
    for (int k = 0; k < n; ++k) {
        lr2 += r[k] * r[k];
    }

    for (int iter = 0; iter < niter; ++iter) {
        // Compute A * p
        for (int k = 0; k < n; ++k) {
            Ap[k] = 0.0;
            for (int l = 0; l < n; ++l) {
                Ap[k] += A[k*n + l] * p[l];
            }
        }

        double pAp = 0.0;
        for (int k = 0; k < n; ++k) {
            pAp += p[k] * Ap[k];
        }

        double alpha = lr2 / pAp;

        double lr2_new = 0.0;
        for (int k = 0; k < n; ++k) {
            g_i[k] += alpha * p[k];
            r[k] -= alpha * Ap[k];
            lr2_new += r[k] * r[k];
        }

        if (sqrt(lr2_new) < 1e-10) break;  // Convergence check

        double beta = lr2_new / lr2;
        lr2 = lr2_new;

        for (int k = 0; k < n; ++k) {
            p[k] = r[k] + beta * p[k];
        }
    }
}


inline void CGNE_bak( int n, double* A, double* g_i, double* b, int niter=100 ){

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


//#endif

