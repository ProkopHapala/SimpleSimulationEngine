# SparseMatrix2.h

This file defines the `SparseMatrix2` class, a template class that implements a sparse matrix data structure with a varying number of neighbors per row. This structure is optimized for storing triangular matrices, such as those resulting from sparse Cholesky factorization. The storage is densely packed for memory efficiency, but inserting new elements requires recreating the matrix.

## Includes

- `arrayAlgs.h`: Provides utility functions for array manipulation, including `binSearchBetween`.

---
## Types (classes and structs)
---

### class `SparseMatrix2`

The `SparseMatrix2` class implements a sparse matrix data structure where each row can have a different number of non-zero elements (neighbors). This design is particularly well-suited for storing triangular matrices, where the number of non-zero elements in each row varies.

#### properties

- `n`:`int`:`public:` - The dimension of the matrix (number of rows).
- `ntot`:`int` - The total number of non-zero elements in the matrix.
- `nngs`:`int*` - An array of size `n` storing the number of non-zero elements in each row.
- `i0s`:`int*` - An array of size `n` storing the starting index of each row in the `vals` and `inds` arrays.
- `vals`:`T*` - An array of size `ntot` storing the non-zero values of the matrix in a folded format.
- `inds`:`int*` - An array of size `ntot` storing the column indices (j) for each non-zero value in the `vals` array.

#### methods

- `dot_dens_vector` - Performs a sparse matrix-dense vector multiplication.
- `get` - Retrieves the value at a given row and column index (i, j). It uses binary search to efficiently check if a non-zero value exists at that location. Returns the value if found, otherwise returns 0.
- `fwd_subs_mi` - Performs forward substitution for a single row `i` in a linear system, assuming the matrix represents a lower triangular matrix. It operates on multi-component vectors (size `ns`).
- `fwd_subs_m_` - Performs forward substitution for the entire linear system, calling `fwd_subs_mi` for each row. It operates on multi-component vectors (size `ns`).
- `fwd_subs_T_mi` - Performs forward substitution for a single row `i` in a linear system, assuming the matrix represents the transpose of a lower triangular matrix. It operates on multi-component vectors (size `ns`).
- `fwd_subs_T_m_` - Performs forward substitution for the entire linear system, calling `fwd_subs_T_mi` for each row. It operates on multi-component vectors (size `ns`).
- `fwd_subs_m` - Performs forward substitution for the entire linear system, assuming the matrix represents a lower triangular matrix. It operates on multi-component vectors (size `ns`) and uses `#pragma omp simd` for potential SIMD parallelization.
- `fwd_subs_T_m` - Performs forward substitution for the entire linear system, assuming the matrix represents the transpose of a lower triangular matrix. It operates on multi-component vectors (size `ns`) and uses `#pragma omp simd` for potential SIMD parallelization.
- `fromDense` - Creates a sparse matrix from a dense matrix, storing only the elements with an absolute value greater than a specified tolerance (`tol`).
- `fromFwdSubT_` - Helper function for `fromFwdSubT`, populating the `inds` and `vals` arrays from a forward substitution matrix.
- `fromFwdSubT` - Creates a sparse matrix from a forward substitution matrix (lower triangular matrix).
- `fprint_inds` - Prints the column indices of the non-zero elements to a file.
- `fprint_vals` - Prints the values of the non-zero elements to a file.