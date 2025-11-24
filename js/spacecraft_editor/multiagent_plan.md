# Multi-Agent Mesh Generation Workflow

This document outlines the workflow for parallelizing the porting of C++ mesh generation functions to JavaScript.

## Goal
Port various mesh generation functions from `constructionBlockApp.cpp` (and other C++ files) to JavaScript using `MeshBuilder.js`. Each function should be implemented in a separate JS script (or a shared library if appropriate) and verified using the Node.js test environment.

## Workflow for Each Agent

1.  **Select a Function**: Pick a function from the "Functions to Port" list below.
2.  **Create/Update Script**:
    *   If it's a new category, create a new test script (e.g., `test_shapes.js`).
    *   If it fits an existing category, add it to the relevant script.
3.  **Implement Logic**:
    *   Translate the C++ logic to JavaScript.
    *   Use `MeshBuilder` methods (`vert`, `edge`, `chunk`, etc.).
    *   Use `Vec3` for vector math.
4.  **Verify**:
    *   Run the script using Node.js: `node js/spacecraft_editor/test_your_script.js`
    *   Check the output `.obj` file in a 3D viewer (e.g., online OBJ viewer or Blender) or visually inspect the file content.
    *   Check the logs for correctness.

## Functions to Port (from `constructionBlockApp.cpp`)

| Function Name | Status | Assigned To | Notes |
| :--- | :--- | :--- | :--- |
| `parabola` | [x] | Agent | ✅ Ported to `MeshesUV.js` as `Parabola_Wire_new` |
| `panel` | [x] | Agent | ✅ Ported to `MeshesUV.js` as `QuadPanel` |
| `QuadSlab` | [x] | Agent | ✅ Ported to `MeshesUV.js` |
| `SlabTube` | [x] | Agent | ✅ Ported to `MeshesUV.js` |
| `QuadSheet` | [x] | Agent | ✅ Ported to `MeshesUV.js` |
| `TubeSheet` | [x] | Agent | ✅ Ported to `MeshesUV.js` |
| `TorusSheet` | [x] | Agent | ✅ Ported to `MeshesUV.js` |
| **Validation** | [x] | Agent | ✅ Added `validateMesh()` with NaN/Inf/range checking |
| `bevel` | [ ] | | Beveling edges - requires adjacency data |
| `skelet` | [ ] | | Skeleton to truss conversion |
| `blocks` | [ ] | | Block network generation |
| `extrude_octahedron` | [ ] | | Octahedron extrusion |
| `oct_nodes` | [ ] | | Octahedron nodes |
| `cube_nodes` | [ ] | | Cube nodes |

## Test Environment Setup

We use a simple Node.js setup to run the scripts without a browser.

### `test_utils.js`
This script sets up the global environment (`window`, `logger`, `Vec3`, `MeshBuilder`).

```javascript
// test_utils.js
const fs = require('fs');
const path = require('path');

// Mock window and logger
global.window = {
    logger: {
        info: (msg) => console.log(`[INFO] ${msg}`),
        debug: (msg) => console.log(`[DEBUG] ${msg}`),
        warn: (msg) => console.warn(`[WARN] ${msg}`),
        error: (msg) => console.error(`[ERROR] ${msg}`),
    }
};

// Load dependencies
// We assume these files are in the correct relative paths
const { Vec3 } = require('../common_js/Vec3.js');
global.Vec3 = Vec3;

const { MeshBuilder } = require('./js/MeshBuilder.js');
global.MeshBuilder = MeshBuilder;

module.exports = {
    saveObj: (filename, meshBuilder) => {
        const content = meshBuilder.toObjString();
        fs.writeFileSync(filename, content);
        console.log(`Saved ${filename} (${content.length} bytes)`);
    }
};
```

## Example Task Script

```javascript
// test_parabola.js
const { saveObj } = require('./test_utils.js');

// Initialize
const mesh = new new MeshBuilder();

// Implement Parabola Logic (Ported from C++)
function makeParabola(mb) {
    mb.logOperation("makeParabola", {});
    // ... implementation ...
    // mb.vert(...)
    // mb.edge(...)
}

// Run
makeParabola(mesh);

//# Export
saveObj('output_parabola.obj', mesh);

## Implementation Status & Notes

### Refactoring (2024-11-24)

All UV-based mesh generation functions have been refactored into a separate module:
- **New Module**: `js/MeshesUV.js` contains all UV parametric surface functions
- **Extension Pattern**: Uses `extendMeshBuilder()` to dynamically add methods to `MeshBuilder.prototype`
- **Node.js & Browser Compatible**: Both environments supported via conditional exports

### Key Bug Fixes

1. **QuadUVfunc NaN Issue** (Critical):
   - **Problem**: `Vec3.setMul(vec3, scalar)` was being called incorrectly - `setMul` expects two Vec3 objects for component-wise multiplication
   - **Fix**: Changed to use `setLincomb3()` for proper bilinear interpolation with scalar coefficients
   - **Impact**: Affected QuadPanel, QuadSlab, and QuadSheet - all now working correctly

2. **Validation System**:
   - Added `validateMesh()` method to MeshBuilder
   - Checks for NaN, Inf, and calculates coordinate bounds
   - Reports degenerate geometries (collapsed dimensions)
   - Automatically called before OBJ export via `saveObj()`

### Test Files

- `tests/test_parabola.js` - Tests `Parabola_Wire_new`
- `tests/test_panel.js` - Tests `QuadPanel` and `QuadSlab`  
- `tests/test_sheets.js` - Tests `QuadSheet`, `TubeSheet`, `TorusSheet`, `SlabTube`

All tests generate valid OBJ files with comprehensive validation logging.

