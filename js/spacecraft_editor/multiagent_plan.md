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
| **Edge Adjacency** | [x] | Agent | ✅ Added `edgesOfVerts` system for bevel support |
| `bevel` | [~] | Agent | ⚠️ **Implemented but needs debugging** - geometry not symmetric |
| `skelet` | [~] | Agent | ⚠️ **Implemented but needs verification** - simplified `skeleton()` generates girder meshes |
| `blocks` | [~] | Agent | ⚠️ **Implemented but needs verification** - simplified `block()` generates subdivided boxes |
| `extrude_octahedron` | [ ] | | **PRIORITY 2**: Standalone platonic solid extrusion |
| `oct_nodes` | [ ] | | **DEFERRED**: Requires full SpaceCraft system (high-level assembly) |
| `cube_nodes` | [ ] | | **DEFERRED**: Requires full SpaceCraft system (high-level assembly) |

## Implementation Priorities

### Priority 1: Simplified Block/Skeleton Functions (CURRENT)
These are more approachable as we can implement simplified versions focused on mesh generation:
- **`blocks`** - Basic block mesh without full BlockBuilder topology
- **`skelet`** - Simplified skeleton-to-truss conversion

### Priority 2: Standalone Extrusion
- **`extrude_octahedron`** - Self-contained platonic solid extrusion

### Deferred: High-Level SpaceCraft Assembly  
These require the complete SpaceCraft class system and are beyond basic mesh generation:
- **`oct_nodes`**, **`cube_nodes`** - Defer until SpaceCraft system is ported

## Known Issues

### Mesh Geometry Verification Needed (⚠️ Pending Gemini 3)

The following implementations generate meshes but produce unexpected geometry that needs verification against C++ reference:

1. **Bevel (`bevel_vert`)**
   - **Issue**: Beveled meshes not symmetric as expected
   - **Possible Cause**: `sortVertEdgesByNormal()` or corner position calculation
   
2. **Skeleton (`skeleton`)**  
   - **Issue**: Generated girder network geometry looks incorrect
   - **Possible Cause**: Up vector calculation or girder positioning

3. **Blocks (`block`)**
   - **Issue**: Subdivided box mesh geometry appears incorrect  
   - **Possible Cause**: Face grid generation or vertex positioning

**Debug Strategy**: 
- Wait for Gemini 3 availability
- Generate detailed C++ execution logs
- Generate detailed JS execution logs  
- Compare logs line-by-line to find discrepancies
- Focus on vertex positions, edge connections, and geometric calculations

**Status**: Implementations functional but geometry needs validation before production use.

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

