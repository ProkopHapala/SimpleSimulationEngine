# Spacecraft Editor WebGL Implementation Plan

## Overview
This document maps the existing C++ Spacecraft Editor functionality to a new WebGL-based implementation. The goal is to create a performant, scriptable editor that runs in the browser, leveraging the architecture of the existing Molecular Editor (`MolGUI_web`).

## Progress Checklist

- [x] **Phase 1: Infrastructure & Scripting**
    - [x] Basic WebGL Boilerplate (Three.js, OrbitControls)
    - [x] [SpaceCraftEngine](file:///home/prokophapala/git/SimpleSimulationEngine/js/spacecraft_editor/js/SpaceCraftEngine.js#1-98) (Data Management)
    - [x] `SpaceCraftWorker` (Scripting Sandbox)
    - [x] [SpaceCraftRenderer](file:///home/prokophapala/git/SimpleSimulationEngine/js/spacecraft_editor/js/SpaceCraftRenderer.js#2-212) (InstancedMesh Setup)
    - [ ] **GUI Refactor** (Single Side Panel)
    - [ ] **Mouse Fix** (Disable Damping/Inertia)
    - [ ] **Debug Rendering** (Fix invisible ship, add verbose logs)

- [ ] **Phase 2: Truss Rendering & Physics**
    - [ ] Implement Truss rendering (instanced spheres/cylinders using texture buffers)
    - [ ] Implement basic physics/simulation loop

- [ ] **Phase 3: Interaction**
    - [ ] Gizmo Controls (Translate/Rotate)
    - [ ] Selection Logic (Raycasting)

## Efficiency & Coding Rules

### 1. Performance First
*   **Zero Allocation in Loops**: Never create new objects (`new Vector3`, `new Array`) inside the render loop or tight simulation loops. Reuse pre-allocated temporary objects.
*   **TypedArrays**: Use `Float32Array` and `Int32Array` for all bulk data (positions, indices, states). Avoid standard JS Arrays for geometry data.
*   **Texture-Based State**: The "Source of Truth" for positions on the GPU is a [DataTexture](file:///home/prokophapala/git/SimpleSimulationEngine/js/spacecraft_editor/js/SpaceCraftRenderer.js#55-74). Do not update attributes per-vertex if possible; update the texture instead.

### 2. Scripting Architecture
*   **Shadow IDs**: The worker predicts IDs to allow synchronous-looking syntax (`n1 = Node()`).
*   **Batching**: Commands are sent from Worker to Main thread in batches, not individually, to minimize message passing overhead.

## Reference Mapping

| Feature | C++ Implementation | WebGL Target Implementation |
| :--- | :--- | :--- |
| **Core Logic** | [SpaceCraft.h](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/Orbital/SpaceCraft.h), [SpaceCraftComponents.h](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/Orbital/SpaceCraftComponents.h) | [SpaceCraftEngine.js](file:///home/prokophapala/git/SimpleSimulationEngine/js/spacecraft_editor/js/SpaceCraftEngine.js) (Main Thread) |
| **Scripting** | [EditSpaceCraft.h](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/Orbital/EditSpaceCraft.h) (Lua) | [SpaceCraftWorker.js](file:///home/prokophapala/git/SimpleSimulationEngine/js/spacecraft_editor/js/SpaceCraftWorker.js) (JS Web Worker) |
| **Rendering** | `DrawOGL3.h` (Immediate/Display Lists) | [SpaceCraftRenderer.js](file:///home/prokophapala/git/SimpleSimulationEngine/js/spacecraft_editor/js/SpaceCraftRenderer.js) (Three.js + InstancedMesh) |
| **Mesh Gen** | [MeshBuilder2.h](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/geometry/MeshBuilder2.h) | `TrussRenderer.js` (Procedural Geometry / Impostors) |

## Architecture

The application will follow a **Pipeline Architecture**:
`Script -> Abstract SpaceCraft -> Concrete Mesh -> Renderer`

### 1. Abstract SpaceCraft (Logical Layer)
*   **Role:** Stores the high-level component definitions (Nodes, Girders, Rings) exactly as defined in [SpaceCraft.h](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/Orbital/SpaceCraft.h).
*   **Data:** JS Objects/Arrays (not optimized for GPU yet).
    *   [Node](file:///home/prokophapala/git/SimpleSimulationEngine/js/spacecraft_editor/js/SpaceCraft.js#4-11): `{ pos: Vec3, type: int, ... }`
    *   [Girder](file:///home/prokophapala/git/SimpleSimulationEngine/js/spacecraft_editor/js/SpaceCraftWorker.js#23-28): `{ nodeA: Node, nodeB: Node, segments: int, ... }`
*   **Responsibility:**
    *   Maintains the logical structure of the ship.
    *   Handled by [SpaceCraft](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/Orbital/SpaceCraft.h#23-571) class in JS.

### 2. Mesh Builder (Generation Layer)
*   **Role:** Converts Abstract SpaceCraft into a concrete Truss Mesh.
*   **Reference:** [SpaceCraft2Mesh2.h](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/Orbital/SpaceCraft2Mesh2.h) and [MeshBuilder2.h](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/geometry/MeshBuilder2.h).
*   **Process:**
    1.  Iterate over all abstract components.
    2.  Generate vertices and edges (sticks) for each component.
    3.  Store mapping: `Component -> [VertexRange, EdgeRange]`.
*   **Output:** A [Truss](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/Orbital/SpaceCraft2Mesh2.h#592-614) object containing:
    *   [verts](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/geometry/MeshBuilder2.cpp#1892-1899): Array of positions.
    *   [edges](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/geometry/MeshBuilder2.cpp#1634-1638): Array of indices `[a, b]`.

### 3. The Renderer (Visual Layer)
*   **Role:** Visualizing the concrete Truss Mesh.
*   **Technique:**
    *   **Nodes (Vertices):** `THREE.InstancedMesh` (Spheres) at `truss.verts`.
    *   **Sticks (Edges):** `THREE.InstancedMesh` (Cylinders) connecting `truss.edges`.
*   **Optimization:**
    *   The [Truss](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/Orbital/SpaceCraft2Mesh2.h#592-614) data is uploaded to GPU buffers (`Float32Array`).
    *   We can still use the "Shared Buffer" technique, but now the texture contains the *generated* truss vertices, not just the abstract nodes.

## Data Structures

### Abstract Components (JS)
```javascript
class Node { pos, type, ... }
class Girder { nodeA, nodeB, nseg, ... }
class SpaceCraft { nodes=[], girders=[], ... }
```

### Concrete Mesh (JS/TypedArray)
```javascript
class MeshBuilder {
    verts = []; // Float32Array [x,y,z, x,y,z...]
    edges = []; // Int32Array [a,b, a,b...]
    
    addVert(pos) { ... }
    addEdge(a, b) { ... }
    
    // Generation methods
    buildGirder(p0, p1, nseg) { ... }
    buildRing(center, axis, radius, nseg) { ... }
}
```

## Implementation Roadmap

### Phase 1: Infrastructure & Scripting (Done)
*   (We will refactor the existing Engine to support the new pipeline)

### Phase 2: Mesh Generation (New Focus)
1.  **Implement [MeshBuilder](file:///home/prokophapala/git/SimpleSimulationEngine/js/common_js/MeshBuilder.js#2-114)**: A class to manage vertices and edges, with methods like [line](file:///home/prokophapala/git/SimpleSimulationEngine/js/common_js/MeshBuilder.js#49-55), [ngon](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/geometry/MeshBuilder2.cpp#2376-2395), [girder](file:///home/prokophapala/git/SimpleSimulationEngine/js/common_js/MeshBuilder.js#79-104).
2.  **Implement [SpaceCraft](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/Orbital/SpaceCraft.h#23-571)**: A class to hold the abstract components.
3.  **Implement `SpaceCraft2Mesh`**: The logic to convert [SpaceCraft](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/Orbital/SpaceCraft.h#23-571) -> [MeshBuilder](file:///home/prokophapala/git/SimpleSimulationEngine/js/common_js/MeshBuilder.js#2-114).
4.  **Update Renderer**: Render the *generated mesh* instead of the abstract nodes.

### Phase 3: Interaction
*   Gizmos will manipulate the *Abstract SpaceCraft* components.
*   On change, we re-run `SpaceCraft2Mesh` (or partially update) to regenerate the visual mesh.

## Open Questions / Future Work
*   **Mesh Generation:** [MeshBuilder2.h](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/geometry/MeshBuilder2.h) has complex logic for "skinning" the truss (plates, hulls). For Phase 1, we will just render the skeleton (Truss). Later, we can implement a Compute Shader or WASM module to generate the hull mesh if needed.
*   **Physics:** The C++ engine has physics. We might eventually need a WASM port of the C++ physics engine if JS performance isn't enough, but for a "Game Editor", JS physics (Cannon.js or custom simple Verlet) might suffice.

---

## Volumetric Truss Generation Architecture

### Core Concept
The system generates **volumetric trusses** (not just surface meshes) for dynamic physics simulation. Each structural element consists of:
- **Mass points** (vertices)
- **Sticks** (edges) connecting the mass points
- The truss has both structural integrity and can be simulated dynamically

### Node-to-Cube-to-Bridge Pipeline

#### 1. Node → Cube Conversion
Each node in the skeleton is converted to a cube (or other platonic solid):
- The cube has 6 faces, each face can have different types:
  - **Type 1**: Single edge attachment → [snapBoxFace](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/geometry/MeshBuilder2.cpp#1358-1373)
  - **Type 2**: Fork (2 edges) → [snapPrismFace](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/geometry/MeshBuilder2.cpp#1416-1438)
  - **Type 3**: Fork (3+ edges) → [snapFrustrumFace](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/geometry/MeshBuilder2.cpp#1374-1395)
- Face geometry creates quad strips (4 vertices per face)

#### 2. Face Identification
When connecting two nodes:
- System finds which faces should be bridged using [findFace(Vec3d d)](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/geometry/ConstructionBlock.h#220-235)
- Determines best face based on direction vector
- Faces stored as "chunks" (quad strips) with vertex indices

#### 3. Bridge Generation ([bridge_quads](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/geometry/MeshBuilder2.cpp#1988-2065))
Takes two quad faces and bridges them with a volumetric truss:
- **Parameters**:
  - `q1, q2`: Two quads (4 vertex indices each)
  - `nseg`: Number of segments along bridge
  - `stickTypes`: Edge types `{longitudinal, ring, spiral, internal}`
  - `mask`: Controls which diagonal/internal edges to generate
  - `bAlling`: Whether to align quads rotationally

- **Structure**:
  - **Ring edges**: Connect vertices within each cross-section
  - **Longitudinal edges**: Connect corresponding vertices between segments
  - **Spiral edges**: Diagonal bracing (mask.x, mask.y)
  - **Internal edges**: Cross-bracing within segments (mask.z, mask.w)

#### Example: [bridge_quads](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/geometry/MeshBuilder2.cpp#1988-2065) with nseg=3
```
q1 = {A1,B1,C1,D1}  →  Ring 0: A1─B1
                                │  │
                                D1─C1
                        
                        Ring 1: A'─B'  (interpolated)
                                │  │
                                D'─C'
                        
                        Ring 2: A''─B'' (interpolated)
                                │  │
                                D''─C''
                        
q2 = {A2,B2,C2,D2}  →  Ring 3: A2─B2
                                │  │
                                D2─C2

Plus: Longitudinal (A1→A'→A''→A2), Spiral, Internal edges
```

### Implementation Strategy for JavaScript

1. **Use `Vec3.js`** for all vector math (performance-critical)
2. **Port key functions**:
   - [alling_polygons](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/geometry/MeshBuilder2.cpp#1944-1987): Rotationally align two quads
   - [bridge_quads](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/geometry/MeshBuilder2.cpp#1988-2065): Generate volumetric bridge
   - [snapBoxFace](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/geometry/MeshBuilder2.cpp#1358-1373), [snapPrismFace](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/geometry/MeshBuilder2.cpp#1416-1438), [snapFrustrumFace](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/geometry/MeshBuilder2.cpp#1374-1395): Face geometry
3. **Data structures**:
   - [verts](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/geometry/MeshBuilder2.cpp#1892-1899): `[{pos: Vec3, nor: Vec3, uv: Vec2}, ...]`
   - [edges](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/geometry/MeshBuilder2.cpp#1634-1638): `[{x: int, y: int, type: int}, ...]`
   - `chunks`: `[{x: stripStart, y: edgeStart, z: count, w: type}, ...]`
   - `strips`: Flat array of vertex indices for faces

### Performance Considerations
- Pre-allocate arrays where possible
- Use Vec3.js methods to minimize temporary object creation
- Batch geometry updates before GPU upload
- Use TypedArrays for large vertex/edge buffers

## C++ Reference Notes

### Key Functions to Port
The following functions from `MeshBuilder2.cpp` and `ConstructionBlock.h` are critical for the "bridge two quads" functionality:

1.  **`alling_polygons(n, ivs1, ivs2, ipiv)`**
    *   Aligns `ivs2` to `ivs1` rotationally to minimize twist.
    *   Uses UV projection onto a plane orthogonal to the axis between COGs.
    *   **Fix Needed:** Ensure `dmin` initialization handles negative dot products correctly (use `-Infinity`).

2.  **`bridge_quads(q1, q2, nseg, stickTypes, mask, bAlling)`**
    *   Generates the volumetric truss between two quads.
    *   Interpolates vertices for intermediate rings.
    *   Creates longitudinal, ring, spiral, and internal edges.

3.  **`polygonNormal(ich)`**
    *   Calculates the normal vector of a polygon chunk.
    *   Essential for determining which face is pointing towards the other block.

4.  **`getChunkCOG(ich)`**
    *   Calculates the Center of Gravity of a chunk.
    *   Used in `findMostFacingNormal` for distance weighting.

5.  **`findMostFacingNormal(hray, chrange, cosMin, bTwoSide, distWeight, ray0)`**
    *   Iterates over a range of chunks (faces).
    *   Computes dot product of face normal with `hray` (direction vector).
    *   Selects the face that best aligns with the direction.

6.  **`bridgeFacingPolygons(p1, p2, chr1, chr2, nseg, stickTypes, maks)`**
    *   High-level function to bridge two blocks.
    *   Calculates direction `hray = p2 - p1`.
    *   Calls `findMostFacingNormal` for both blocks to find the best faces.
    *   Calls `bridge_quads` to connect them.

### Implementation Details
*   **`ConstructionBlock`**: Represents a block with faces. In JS, we might simplify this or map it to `MeshBuilder` chunks.
*   **`BlockBuilder`**: Manages a collection of blocks.
*   **`ConstructionBlockToMeshBuilder`**: Bridges the logical blocks to the mesh generation.

We will implement these directly in `MeshBuilder.js` to enhance the procedural generation capabilities.
