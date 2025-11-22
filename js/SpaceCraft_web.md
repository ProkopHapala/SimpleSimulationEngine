# Spacecraft Editor WebGL Implementation Plan

## Overview
This document maps the existing C++ Spacecraft Editor functionality to a new WebGL-based implementation. The goal is to create a performant, scriptable editor that runs in the browser, leveraging the architecture of the existing Molecular Editor (`MolGUI_web`).

## Reference Mapping

| Feature | C++ Implementation | WebGL Target Implementation |
| :--- | :--- | :--- |
| **Core Logic** | [SpaceCraft.h](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/Orbital/SpaceCraft.h), [SpaceCraftComponents.h](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/Orbital/SpaceCraftComponents.h) | `SpaceCraftEngine.js` (Main Thread) |
| **Scripting** | [EditSpaceCraft.h](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/Orbital/EditSpaceCraft.h) (Lua) | `SpaceCraftWorker.js` (JS Web Worker) |
| **Rendering** | `DrawOGL3.h` (Immediate/Display Lists) | `SpaceCraftRenderer.js` (Three.js + InstancedMesh) |
| **Mesh Gen** | [MeshBuilder2.h](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/geometry/MeshBuilder2.h) | `TrussRenderer.js` (Procedural Geometry / Impostors) |

## Architecture

The application will follow a **Command-Query Separation (CQS)** pattern with a **Shared GPU State**.

### 1. The Engine (Main Thread)
*   **Role:** Source of Truth. Manages the physics/structural state.
*   **Data:** Stores components in TypedArrays for performance and easy GPU upload.
    *   `Nodes`: `Float32Array` (Positions: x, y, z, w=type)
    *   `Girders`: `Int32Array` (Topology: nodeA, nodeB, type, state)
    *   `Ropes`: `Int32Array` (Topology: nodeA, nodeB, type, state)
*   **Responsibility:**
    *   Receives commands from the Script Worker (`CREATE_NODE`, `CREATE_GIRDER`).
    *   Updates the `DataTexture` (GPU State).
    *   Handles physics simulation (future).

### 2. The Scripting Sandbox (Web Worker)
*   **Role:** User Logic Execution.
*   **Security:** Runs in a Web Worker with no DOM access.
*   **Pattern:** "Shadow ID" / Deterministic Execution.
    *   The worker predicts IDs (0, 1, 2...) locally to allow synchronous-like scripting syntax (`g1 = Girder(n1, n2)`).
    *   Sends batched commands to the Engine.
*   **API:** Mirrors [EditSpaceCraft.h](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/Orbital/EditSpaceCraft.h) but in JS.
    *   `api.Node(pos)`
    *   `api.Girder(n1, n2, type)`
    *   `api.Rope(n1, n2, thick)`

### 3. The Renderer (WebGL / Three.js)
*   **Role:** Visualizing the state.
*   **Technique:** **Texture-Fetch Instancing** (Derived from `MolGUI_web`).
    *   **Nodes:** `THREE.InstancedMesh` (Sphere Impostors).
        *   Vertex Shader fetches position from `uPosTex` using `gl_InstanceID`.
    *   **Girders/Ropes:** `THREE.InstancedMesh` (Cylinder Impostors or Low-Poly Mesh).
        *   **Crucial Difference from Molecules:** Bonds/Girders connect two dynamic points.
        *   Vertex Shader fetches **two** positions (`posA`, `posB`) from `uPosTex` based on `attr_NodeIDs` (vec2 attribute per instance).
        *   Shader computes the center, orientation (quaternion/matrix), and length on the fly.
        *   **Benefit:** Zero CPU cost for updating girder positions when nodes move.

## Data Structures & GPU Layout

### Texture Layout (`uPosTex`)
A `DataTexture` (Float32, RGBA) holding node states.
*   Size: `N x 1` (where N is max nodes, e.g., 4096).
*   Pixel `i`: `[x, y, z, type]` of Node `i`.

### Component Buffers (CPU -> GPU Attributes)

#### Nodes (InstancedMesh)
*   **Geometry:** Quad (Billboard) or Low-poly Sphere.
*   **Instance Attributes:**
    *   `aNodeID` (int): Index into `uPosTex`.

#### Girders (InstancedMesh)
*   **Geometry:** Cylinder (oriented Y-up).
*   **Instance Attributes:**
    *   `aNodeIDs` (vec2 int): Indices of start and end nodes `[idA, idB]`.
    *   `aType` (int): Material/Type index.
    *   `aParams` (vec4): Extra params (thickness, deformation, etc.).

## Implementation Roadmap

### Phase 1: Infrastructure & Scripting Bridge
1.  Setup `SpaceCraftEngine` class.
2.  Setup `SpaceCraftWorker` with the "Shadow ID" logic.
3.  Implement basic commands: [Node](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/Orbital/SpaceCraftComponents.h#386-413), [Girder](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/Orbital/SpaceCraftComponents.h#449-526).
4.  Verify: Run a script that generates a simple cube truss and log the Engine state.

### Phase 2: Truss Rendering (The "Hard" Part)
1.  Implement `uPosTex` management in `SpaceCraftRenderer`.
2.  Implement **Node Renderer**:
    *   Shader that reads `uPosTex`.
3.  Implement **Girder Renderer**:
    *   Shader that reads `uPosTex[idA]` and `uPosTex[idB]`.
    *   Computes `midpoint = (A+B)*0.5`.
    *   Computes rotation matrix to align Y-axis with `B-A`.
    *   Scales Y-axis by `length(B-A)`.

### Phase 3: Interaction (Gizmos)
1.  Port `EditorController` from MolGUI.
2.  Adapt it to write to `uPosTex` (or the Engine's CPU buffer -> Texture update).
3.  Since Girders are procedurally rendered from Nodes, moving a Node with a Gizmo will **instantly** update all connected Girders without extra logic.

## Open Questions / Future Work
*   **Mesh Generation:** [MeshBuilder2.h](file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/geometry/MeshBuilder2.h) has complex logic for "skinning" the truss (plates, hulls). For Phase 1, we will just render the skeleton (Truss). Later, we can implement a Compute Shader or WASM module to generate the hull mesh if needed.
*   **Physics:** The C++ engine has physics. We might eventually need a WASM port of the C++ physics engine if JS performance isn't enough, but for a "Game Editor", JS physics (Cannon.js or custom simple Verlet) might suffice.
