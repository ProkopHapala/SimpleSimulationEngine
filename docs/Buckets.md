# Buckets Data Structure Documentation

## General Philosophy

The `Buckets` class is a specialized data structure designed for efficiently mapping objects to cells (buckets) and performing fast lookups. It's particularly useful for:

- Storing adjacency relationships (like edges connected to vertices in a mesh)
- Spatial partitioning (like grid cells in spatial hashing)
- Any scenario requiring fast "find all objects in cell X" operations

The key advantages are:
- O(1) access to all objects in a given cell
- Memory efficiency through compact storage
- Fast insertion and updates when using the `updateCells()` method

## Core Data Members

```cpp
int ncell, nobj;       // Number of cells and objects
int* cellNs = 0;       // [ncell] Number of objects per cell
int* cellI0s = 0;      // [ncell] Starting index of each cell in cell2obj
int* cell2obj = 0;     // [nobj] Contiguous array of object indices
int* obj2cell = 0;     // [nobj] Optional reverse mapping (object to cell)
```

## Key Methods

### Initialization

```cpp
void realloc(int ncell_, int nobj_, bool bO2C=false)
```
Allocates memory for `ncell_` cells and `nobj_` objects. Set `bO2C=true` if you need reverse mapping.

### Building the Structure

```cpp
void bindObjs(int nobj_, int* obj2cell_)
```
Binds an external object-to-cell mapping array.

```cpp
void updateCells(int nobj_=-1, int* obj2cell_=0)
```
Rebuilds the entire structure in 4 steps:
1. Clears cell counts
2. Counts objects per cell
3. Computes cell offsets
4. Distributes objects to cells

### Querying

```cpp
int* inCell(int icell)
```
Returns pointer to objects in cell `icell`.

```cpp
int getInCell(int icell, int* out)
```
Copies objects from cell `icell` to `out` array.

## Example Usage

### Mesh Edge Adjacency Example

```cpp
// In MeshBuilder2.cpp
void Builder2::build_edgesOfVerts(){
    int nv = verts.size();  // Number of cells (vertices)
    int ne = edges.size();  // Number of objects (edges)
    
    // Create object-to-cell map (each edge appears twice - once per vertex)
    int* obj2cell = new int[ne*2];
    for(int ie=0; ie<ne; ie++){
        Vec2i e = edges[ie].lo;
        obj2cell[2*ie]   = e.x;  // Edge belongs to vertex e.x
        obj2cell[2*ie+1] = e.y;  // Edge belongs to vertex e.y
    }
    
    // Initialize and build the Buckets structure
    edgesOfVerts.realloc(nv, ne*2, true);
    edgesOfVerts.bindObjs(ne*2, obj2cell);
    edgesOfVerts.updateCells();
    
    delete[] obj2cell;
}
```

## Operations Guide

### Adding Objects
1. Ensure sufficient capacity with `realloc()`
2. Update your `obj2cell` mapping
3. Call `updateCells()`

### Removing Objects
1. Update your `obj2cell` mapping (set cell to -1)
2. Call `updateCells()`

### Querying Objects
```cpp
// Get edges connected to vertex iv
int n = edgesOfVerts.cellNs[iv];
int* edgesInCell = edgesOfVerts.inCell(iv);
```

## Performance Notes
- `updateCells()` is O(N) but reconstructs the whole structure
- Queries are O(1) after construction
- Ideal for static or batch-updated data
- For dynamic updates, consider alternative strategies

## Debugging

```cpp
void printCells(int verb=0);
void printObjCellMaping(int i0=-1, int i1=-1, bool bBoundary=true);
```
Use these to inspect the internal state of the data structure.
