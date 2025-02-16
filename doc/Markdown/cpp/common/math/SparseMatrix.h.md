# SparseMatrix.h

This file defines the `SparseMatrix` class, a template class that implements a sparse matrix data structure with a constant maximum number of neighbors per row. This structure is optimized for storing matrices arising from nearest neighbor interactions, trusses, and molecular structures, where each node typically interacts with a limited number of other nodes. It is less efficient for storing triangular matrices.

## Includes

- `arrayAlgs.h`: Provides utility functions for array manipulation, including `binarySearch_ignor` and `insertElement`.
- `CGNE.h`: Likely provides an implementation of the Conjugate Gradient Normal Error (CGNE) algorithm, used in the `sparse_fsai` function.

---

## Free functions

- `sparse_fsai` - Computes an approximate inverse of a sparse matrix using the FSAI (Factorized Sparse Approximate Inverse) technique, employing the CGNE algorithm to solve the least squares problem.

---
## Types (classes and structs)
---

### class `SparseMatrix`

The `SparseMatrix` class implements a sparse matrix data structure where each row has a limited number of non-zero elements (neighbors). This design is optimized for scenarios where the matrix represents interactions between a limited number of entities, such as in nearest neighbor problems or structural mechanics.

#### properties

- `n`:`int`:`public:` - The dimension of the matrix (number of rows).
- `m`:`int` - The maximum number of neighbors (non-zero elements) allowed in each row.
- `nng`:`int*` - An array of size `n` storing the number of neighbors (non-zero values) in each row.
- `vals`:`T*` - An array of size `m*n` storing the non-zero values of the matrix in a folded format.
- `inds`:`int*` - An array of size `m*n` storing the column indices (j) for each non-zero value in the `vals` array. A value of -1 indicates an empty slot.

#### methods

- `realloc` - Reallocates memory for the internal arrays (`nng`, `inds`, and `vals`) to accommodate a new matrix size (`n_`) and maximum number of neighbors (`m_`). It initializes `nng` to 0 and `inds` to -1.
- `get` - Retrieves the value at a given row and column index (i, j). It uses binary search to efficiently check if a non-zero value exists at that location. Returns the value if found, otherwise returns 0.
- `set` - Sets the value at a given row and column index (i, j). It uses binary search to find the element. If the element exists, its value is updated. If it doesn't exist, a new element is inserted into the sorted neighbor list, shifting existing elements if necessary. Returns `true` if a new element was inserted, `false` if an existing element was updated.
- `dot_dens_vector` - Performs a sparse matrix-dense vector multiplication.
- `dot_dens_vector_m` - Performs a sparse matrix-dense vector multiplication where the dense vector has multiple components per element.
- `checkRow` - Checks if a given row of the sparse matrix matches a reference dense vector within a specified tolerance.
- `checkDens` - Checks if the entire sparse matrix matches a reference dense matrix within a specified tolerance.