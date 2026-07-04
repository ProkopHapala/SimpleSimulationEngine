# Review Untracked Files and Folders

**Created**: 2025-07-04
**Context**: After gitignore cleanup (including artifacts like .png, .svg, .obj, .truss, .xyz, .xls, OUT*, .log), only code and .md files remain untracked. Need to review each and decide: commit, archive, or delete.

---

## Untracked Folders (new experiment projects)

### Python
- [ ] **`python/GLCL/`** — GLCL.py, NBody_glcl.py, NBody_glsl.py. OpenCL/GL integration experiment?
- [ ] **`python/GUI/`** — GUITemplate.py. GUI template for Python apps?
- [ ] **`python/terrain_ocl_2/`** — terrain.cl, terrain.py, test_terrain.py, verify_logic.cl, tmp.md. Iteration on terrain_ocl?
- [ ] **`python/pySymGLSL/shaders/shader_toy_ref/`** — BuffA/B/C/D.glslf, common.glslf, Image.glslf. ShaderToy reference shaders.
- [ ] **`python/terrain_ocl/`** — test_errosion.py, test_errosion_2.py, tmp.md. (OUT is gitignored)

### JS
- [ ] **`js/ppstm_web/`** — index.html, ppstm.js, stm_shader.glslf. Web-based STM experiment?
- [ ] **`js/LandCraft_web/data/`** — elements.csv. Data for LandCraft web.
- [ ] **`js/FlowField/`** — FlowField.md, FlowField_orbs*.HTML. Flow field rendering.

### Docs
- [ ] **`doc/js/`** — WebGL/WebGPU notes, JS Vec3, user scripting notes from various AI tools (Gemini, Grok, Kimi).
- [ ] **`doc/Markdown/cpp/common/OCL/`** — ProjectiveDynamicsOCL.md. OCL documentation.

### C++
- [ ] **`cpp/common_resources/Visualizer/`** — Various .glslf shader files (BelousovZhabotinsky, FluidDrift, Kaleidoscope, etc.)
- [ ] **`cpp/common_resources/shaders/ShaderToy/`** — SDF shader files.

---

## Untracked Files (individual)

### Python experiments
- [ ] **`python/OrbitalWar/`** — trajectoryOptmization.{py,cl,md}, test_trajectoryOptmization.py. Orbital trajectory optimization.
- [ ] **`python/pyTruss/`** — ARCHITECTURE.md, IMPLEMENTATION_SUMMARY.md, README_new_solver.md, example_ocl_new.py, run_vbd_cloth_old.py. New OCL truss solver docs.
- [ ] **`python/EulerianImpacFluid/`** — tmp.md.
- [ ] **`python/GLCL2/`** — cl/soldiers.cl, shaders/*.glslf. GLCL2 shaders.
- [ ] **`python/pyMolecular/MMFF-ref.py`** — MMFF reference?
- [ ] **`python/pyShaderToy/fragment_code.glslf`**
- [ ] **`python/pySymGLSL/shaders/fluid/solveFluid_mini.glslf`**
- [ ] **`python/p`** — unclear what this is (single file named "p").
- [ ] **`python/test.json`**

### JS experiments
- [ ] **`js/RigidBodyAero/`** — RigidBody.js, RigidBodyAero_plan.md. Rigid body aerodynamics.
- [ ] **`js/spacecraft_editor/`** — MeshGenerators---.js, tests/diag_wheel.js, test_slider_paths.js, test_slider_perf.js, test_wheel_generators.js.
- [ ] **`js/index-.html`**

### Doc/python experiments (Burn1D, MHD, PotentialFlow)
- [ ] **`doc/python/Burn1D/`** — ~20 .py files (EulerianTube, PipeValvelesEngine, SphericalImplosion, wave_in_tube, etc.) + .md. 1D combustion/flow experiments.
- [ ] **`doc/python/MHD/`** — check_B_consistency.py, check_B_kernels.py, demo_dipole_gemini.py, tmp.md. MHD validation scripts.
- [ ] **`doc/python/PotentialFlow/`** — PotentialFlow.py, SlenderBodyTheory*.py, airfoil*.py, simple_vlm.py. Aerodynamics experiments.

### Other
- [ ] **`package.json` / `package-lock.json`** — Node.js package files for three.js. Should these be committed?
- [ ] **`classify_md.py`** — already staged.
- [ ] **`tests_bash/sketches_SDL/test_2D.sh`** — test script.

---

## Action Items

1. **Artifacts gitignored** — `.png`, `.svg`, `.obj`, `.truss`, `.xyz`, `.xls`, `OUT*`, `OUT`, `.log`, `LDLT_*.txt` added to .gitignore
2. **Experimental Python/JS projects** → review each, either commit or move to archive branch
3. **Doc/python experiments** (Burn1D, MHD, PotentialFlow) → likely worth committing as reference
4. **Shader files** (Visualizer, ShaderToy) → likely worth committing
5. **`python/p`** → investigate what this is, likely delete
