# JS Module Refactor Plan

This document tracks the planned refactor from global/script-based JS to proper ES modules, and the consolidation of shared geometry code.

## 0. Scope

- Shared JS utilities in:
  - `js/common_js`
- Two main browser apps:
  - `js/molgui_web`
  - `js/spacecraft_editor`

Goals:

- Use **ES modules** (`import` / `export`) instead of implicit globals from `<script src="...">`.
- Centralize shared geometry / selection code in `js/common_js`.

---

## 1. Move shared geometry code into `js/common_js`

### 1.1 Files to move

From `js/spacecraft_editor/js/` into `js/common_js/`:

- `MeshesUV.js`  →  `js/common_js/MeshesUV.js`
- `MeshBuilder.js`  →  `js/common_js/MeshBuilder.js`
- `Selection.js`  →  `js/common_js/Selection.js`

### 1.2 After moving, update all references

**1.2.1 SpaceCraft editor side**

- In `js/spacecraft_editor/index.html` (or equivalent entry point):
  - Replace references like:
    - `js/spacecraft_editor/js/MeshBuilder.js`
    - `js/spacecraft_editor/js/MeshesUV.js`
    - `js/spacecraft_editor/js/Selection.js`
  - With new module imports (see §2 once ES modules are enabled) or, temporarily, new script paths:
    - `../common_js/MeshBuilder.js`
    - `../common_js/MeshesUV.js`
    - `../common_js/Selection.js`

- In `js/spacecraft_editor/js/SpaceCraft.js`, `MeshGenTestGUI.js`, `MeshGenerators.js`, and any other files that currently assume global `MeshBuilder`, `MeshesUV`, `Selection`:
  - In **global-script phase** (before ES modules):
    - Just ensure the new common_js scripts are loaded before these files in HTML.
  - In **ES-module phase** (after §2):
    - Replace implicit globals with explicit imports, e.g.:
      ```js
      import { MeshBuilder } from '../../common_js/MeshBuilder.js';
      import { MeshesUV }   from '../../common_js/MeshesUV.js';
      import { Selection, SelectionBanks } from '../../common_js/Selection.js';
      ```

**1.2.2 MolGUI side**

- In `js/molgui_web/index.html` and any JS (e.g. `MoleculeRenderer.js`, `Editor.js`) that will use these shared geometry utilities:
  - After the move, reference shared modules from `../common_js/…` instead of duplicating logic.
  - Once ES modules are enabled (see §2), use imports such as:
    ```js
    import { MeshBuilder } from '../common_js/MeshBuilder.js';
    import { Selection, SelectionBanks } from '../common_js/Selection.js';
    ```

**1.2.3 Keep C++ ↔ JS architectural mapping

- Mirror C++ layout:
  - `cpp/common/geometry/MeshBuilder2.*` ↔ `js/common_js/MeshBuilder.js`
  - `cpp/common/geometry/Selection.*`   ↔ `js/common_js/Selection.js`
  - `cpp/common/geometry/SDfuncs.h`     ↔ `js/common_js/SDfuncs.js`

Document key API parallels directly in JS files as short comments (no long docs change required).

---

## 2. Convert to ES modules (browser)

### 2.1 Entry points: `index.html`

For **MolGUI Web** (`js/molgui_web/index.html`):

- Change current script bootstrapping:

  ```html
  <script src="../common_js/Logger.js?v=5"></script>
  <script src="js/MMParams.js?v=5"></script>
  <script src="js/MoleculeSystem.js?v=5"></script>
  <!-- ... many more ... -->
  <script src="js/main.js?v=5"></script>
  ```

- To a single ES-module entry:

  ```html
  <script type="module" src="js/main.js"></script>
  ```

- Remove `?v=…` cache busting or handle it in the module URL if still desired.

For **SpaceCraft Editor** (`js/spacecraft_editor/index.html` or analogous):

- Apply the same pattern: replace the chain of `<script src>` tags with a single:

  ```html
  <script type="module" src="js/main.js"></script>
  ```

- Ensure `main.js` takes responsibility for importing all other modules.

### 2.2 Refactor `js/common_js` to ES modules

For each file under `js/common_js` (e.g. `Vec3.js`, `Logger.js`, `Draw3D.js`, `MeshRenderer.js`, `SDfuncs.js`, `MeshBuilder.js`, `MeshesUV.js`, `Selection.js`):

1. **Exports**
   - Replace global assignments like:
     ```js
     if (typeof window !== 'undefined') {
         window.Vec3 = Vec3;
     }
     ```
   - With ES exports:
     ```js
     export class Vec3 { /* ... */ }
     // or
     export { Vec3 };
     ```

2. **Imports**
   - Replace implicit globals with explicit imports:
     - Example (`SDfuncs.js`):
       ```js
       import { Vec3 } from './Vec3.js';
       ```
   - Use **relative module paths** consistent with where the file lives.

3. **Optional Node support**
   - If Node/testing support is still required, consider a dual build or keep a small CommonJS wrapper; for the first phase, prioritize browser ES modules.

### 2.3 Refactor `js/molgui_web` to ES modules

For each file in `js/molgui_web/js` (`main.js`, `Editor.js`, `GUI.js`, `ShortcutManager.js`, `MoleculeSystem.js`, `MoleculeRenderer.js`, `LabelRenderer.js`, `IO.js`, etc.):

1. **Stop using `window.*` as the primary export**.
   - Example: instead of `window.app = { ... }`, prefer explicit exports or a created singleton exported from `main.js`.

2. **Add imports at the top** for dependencies:
   - `Vec3`, `SDfuncs`, `MeshBuilder`, `Selection`, etc. from `../common_js/`.
   - Three.js and controls:
     - Either continue to use CDN/global `THREE`, or switch to bundling/importing three.js separately in a later phase.

3. **Wire startup in `main.js`**:
   - `main.js` becomes the central place that constructs `Editor`, `GUI`, `MoleculeSystem`, `MoleculeRenderer`, etc., using imports.

### 2.4 Refactor `js/spacecraft_editor` to ES modules

Similar steps as §2.3:

- For `SpaceCraft.js`, `MeshGenTestGUI.js`, `MeshGenerators.js`, `MeshesUV.js` (after move), `MeshBuilder.js` (after move), `Selection.js` (after move):
  - Import from `../common_js/...` instead of assuming globals.
  - Export classes/functions as ES modules instead of attaching to `window`.
- Centralize app bootstrapping in `js/spacecraft_editor/js/main.js` loaded via `type="module"`.

---

## 3. Migration Strategy

1. **Phase A: Centralize shared geometry (no module system yet)**
   - Move `MeshesUV.js`, `MeshBuilder.js`, `Selection.js` into `js/common_js/`.
   - Update HTML script order in both apps so these are loaded before any code that uses them.
   - Verify no regressions.

2. **Phase B: Introduce ES modules in `common_js`**
   - Convert `Vec3.js`, `SDfuncs.js`, `Selection.js`, `MeshBuilder.js`, `MeshesUV.js`, etc., to `export`/`import` style.
   - Temporarily adapt one app (e.g. `spacecraft_editor`) to ES modules first while leaving MolGUI on globals, or switch both at once if convenient.

3. **Phase C: Convert MolGUI Web and SpaceCraft editor**
   - Replace multiple `<script>` tags with a single `type="module"` entry in each `index.html`.
   - Update all JS files under `js/molgui_web/js` and `js/spacecraft_editor/js` to use ES imports/exports.
   - Remove legacy `window.*` exports where no longer needed.

4. **Phase D: Cleanup**
   - Remove leftover CommonJS / global guards in `common_js` once browsers and tooling rely fully on ES modules.
   - Optionally introduce a bundler (Vite/Rollup/Webpack) if desired, but this is not required for basic ES module support in modern browsers.

---

## 4. Notes / Open Questions

### 4.1 Three.js, CDN, and ES modules

- Three.js is currently loaded via **CDN as a global `THREE`**. This is **compatible** with converting our own code to ES modules:
  - We can keep the existing CDN `<script src="..."></script>` for Three.js and controls.
  - In our ES-module code (`main.js`, etc.) we simply **use the global**:
    - `const scene = new THREE.Scene();` (no `import` yet).
  - This is the **initial plan for Phase B/C**: modularize our code, keep Three.js as-is.

- Later, if desired, Three.js itself can be loaded as an ES module:
  - From a CDN that serves ESM bundles, e.g.
    - `import * as THREE from 'https://unpkg.com/three@0.128.0/build/three.module.js';`
  - Or from a **local copy** checked into the repo, e.g.
    - `import * as THREE from '../../vendor/three/build/three.module.js';`
  - This choice (CDN vs local) is **orthogonal** to using ES modules:
    - CDN gives shared caching and no repo bloat but requires network.
    - Local copy enables fully offline use and version control of the exact Three.js file.

- We can make Three.js source (CDN URL vs local file) configurable later by centralizing the import in `main.js` or a tiny `three_entry.js`.

### 4.2 Misc

- Versioned query strings (`?v=5`) on scripts can be kept initially or replaced with a more systematic cache-busting strategy later.
- Node/test harnesses for geometry code may eventually want a dedicated entry (e.g. `common_js/index.js`) that re-exports `Vec3`, `SDfuncs`, `MeshBuilder`, `Selection`, etc.
