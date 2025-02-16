# Buckets.h

This file defines the `Buckets` class, which is used for managing a collection of objects distributed among a set of cells or buckets. It provides functionalities for assigning objects to cells, tracking the number of objects in each cell, and efficiently retrieving the objects within a given cell. This is useful in spatial partitioning and collision detection algorithms, as well as other applications where grouping objects into buckets based on some criteria is needed.

## Includes

- `macroUtils.h`: Provides utility macros, such as `_realloc`, used for memory management within the `Buckets` class.

---
## Types (classes and structs)
---

### class `Buckets`

The `Buckets` class manages the assignment of objects to cells (buckets). It provides methods for updating the cell occupancy based on an object-to-cell mapping, resizing the cell and object arrays, and retrieving the objects contained within a specific cell.  It's useful in scenarios where spatial partitioning or grouping objects into buckets is required.

#### properties

- `ncell`:`int` - The number of cells (buckets) managed by the `Buckets` object.
- `nobj`:`int`:`public:` - The number of objects managed by the `Buckets` object.
- `cellNs`:`int*` - An array of size `ncell` storing the number of objects contained within each cell. `cellNs[i]` represents the number of objects in the i-th cell.
- `cellI0s`:`int*` - An array of size `ncell` storing the starting index of the objects contained in each cell within the `cell2obj` array. `cellI0s[i]` is the index in `cell2obj` where the objects belonging to the i-th cell begin.
- `cell2obj`:`int*` - An array of size `nobj` storing the indices of the objects, organized by cell. The objects belonging to cell `i` are stored in the range `cell2obj[cellI0s[i] ... cellI0s[i+1]-1]`.
- `maxInBucket`:`int` - The maximum number of objects found in any single cell. This value is updated during the `updateOffsets` method.
- `nobjSize`:`int` - The allocated size of the `cell2obj` array. This may be larger than `nobj` to avoid frequent reallocations.
- `obj2cell`:`int*` - An array of size `nobj` mapping each object to its corresponding cell. `obj2cell[i]` represents the index of the cell that contains object `i`.

#### methods

- `clean` - Resets the number of objects in each cell (`cellNs`) to zero.
- `cleanO2C` - Resets the `obj2cell` array, setting each object's cell assignment to a default value (usually -1, indicating no cell).
- `count` - Counts the number of objects in each cell based on the provided `obj2cell` mapping and stores the counts in the `cellNs` array.
- `updateOffsets` - Updates the `cellI0s` array based on the counts in `cellNs`, calculating the starting index of each cell in the `cell2obj` array.  It also determines and stores the maximum number of objects in a single cell in `maxInBucket`.
- `objectsToCells` - Assigns objects to their corresponding cells in the `cell2obj` array, using the `obj2cell` mapping and the offsets stored in `cellI0s`.
- `updateCells` - Updates the cell occupancy information based on the object-to-cell mapping. It performs cleaning, counting, offset updates, and object assignment in sequence.
- `resizeCells` - Resizes the `cellNs` and `cellI0s` arrays if the given `ncell_` is larger than the current number of cells (`ncell`).  Uses the `_realloc` macro from `macroUtils.h`.
- `resizeObjs` - Resizes the `cell2obj` and optionally the `obj2cell` arrays if the given `nobj_` is larger than the current allocated object size (`nobjSize`). Uses the `_realloc` macro from `macroUtils.h`.
- `bindObjs` - Binds a pre-existing `obj2cell_` array to the `Buckets` object.  This avoids allocating a new `obj2cell` array within the class.  It also resizes the object arrays.
- `realloc` - Resizes both the cell and object arrays using `resizeCells` and `resizeObjs`.  Provides a convenient way to reallocate all the internal storage of the `Buckets` object.
- `getInCell` - Copies the indices of objects in a given cell (`icell`) to the output array `out`. Returns the number of objects in the cell.
- `inCell` - Returns a pointer to the beginning of the object indices for a given cell (`icell`) within the `cell2obj` array. This allows direct access to the objects in the cell without copying.
- `printObjCellMaping` - Prints the mapping between objects and cells, useful for debugging.  It iterates through the `obj2cell` and `cell2obj` arrays within a specified range.
- `printCells` - Prints information about the cells, including the number of objects in each cell and the starting index of each cell in the `cell2obj` array.  Can optionally print the contents of each cell (object indices) if `verb` is greater than 0.