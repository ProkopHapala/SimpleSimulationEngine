# Soldiers 2D Army Simulation (GLCL2)

Design document for a 2D soldier simulation using OpenCL (compute) and OpenGL (render) within the GLCL2 framework. Focuses on small-scale MVP (10 vs 10) with flocking/line formation behavior and clear path to extensions.

- Codebase: `python/GLCL2/`
- Reference examples: `scripts/nbody.py`, `cl/nbody.cl`
- Core engine: `GLCLBrowser.py`, `GLCLGUI.py`, `OGLsystem.py`

---

## 1. Objective (MVP)

Simulate two opposing soldier groups as particles in a 2D top-down plane. Each soldier has:
- Position `p=(x,y)`
- Orientation unit vector `u=(ux,uy)` (weapon/face direction)
- Velocity `v=(vx,vy)` (independent of `u`)
- Mass `m`, radius `r`
- Team and type identifiers

Behavior (within cutoff radius):
- Friends: align orientations; strong lateral preference (±90° relative to `u`) to encourage line formation; mild cohesion; separation to avoid overlap.
- Enemies: turn to face nearest enemy; regulate forward spacing.

Rendering MVP: points positioned by `p`, oriented glyphs later.

---

## 2. Data Layout and Buffers

To match GLCL’s single-attribute vertex path (see `GLCLGUI.py` geometry baking), we pack `p,u` in one `vec4` buffer used for rendering. Remaining state in CL-only buffers.

Buffers (names as used in GLCL config and kernel):

- `state_pos_dir`: float4 per soldier
  - `.xy = p` (position in NDC space for MVP)
  - `.zw = u` (orientation unit vector)
  - Used as OpenGL vertex buffer (layout(location=0) vec4)

- `state_vel_mr`: float4 per soldier
  - `.xy = v`
  - `.z  = m`
  - `.w  = r`

- `state_team_type`: float4 per soldier
  - `.x = team_id` (as float)
  - `.y = type_id` (as float)
  - `.zw = reserved` (future per-soldier params)

Optional (future):
- `type_params_A`: float4 per type [d_front, d_side, max_speed, max_turn]
- `type_params_B`: float4 per type [w_sep, w_align, w_coh, w_enemy]

Rationale:
- Keeps GL rendering simple (one attribute buffer), avoids GLSL `float8`.
- CL has coalesced and compact access patterns.

---

## 3. Parameters (GUI-exposed)

Defined in `config["parameters"]` so GLCL can build controls and update uniforms/args automatically:

- Simulation:
  - `particle_count: int` (default 20)
  - `dt: float`

- Cutoffs and spacings:
  - `r_cut: float` (neighbor radius)
  - `d_front: float` (preferred enemy distance ahead)
  - `d_side: float` (preferred lateral spacing for line)

- Weights:
  - `w_sep: float` (separation/repulsion)
  - `w_align: float` (orientation alignment among friends)
  - `w_coh: float` (cohesion within line, mild)
  - `w_enemy: float` (turn toward nearest enemy)

- Kinematics:
  - `max_speed: float`
  - `max_turn: float` (rad/s cap for orientation change)
  - `friction: float` (velocity damping)

- Anisotropy shaping:
  - `side_bias: float` (boost for ±90° neighbors to form lines)
  - `front_bias: float` (boost/suppress front/back neighbors)

---

## 4. Kernel Design (OpenCL)

File: `python/GLCL2/cl/soldiers.cl`
Kernel: `soldiers_step(__global float4* pos_dir, __global float4* vel_mr, __global float4* team_type, float dt, int n, float r_cut, float d_front, float d_side, float w_sep, float w_align, float w_coh, float w_enemy, float max_speed, float max_turn, float friction, float side_bias, float front_bias)`

Algorithm per soldier i (O(N^2) MVP):
- Read state: `p_i, u_i, v_i, m_i, r_i, team_i`
- Loop neighbors j:
  - `d = p_j - p_i`, `r = length(d)`, `h = d / (r+eps)`
  - If `r <= r_cut`:
    - Separation (all): accumulate `sep += -h * softcore(r, r_i, r_j)`
    - If friend (`team_j == team_i`):
      - Alignment: accumulate `u_align += u_j * w_lat`, where `w_lat = (1 - abs(dot(h, u_i))) * side_bias`
      - Lateral cohesion: accumulate `coh += dot(d, n_i) * n_i` with `n_i = perp(u_i)=(−uy,ux)`, targeting `|coh| -> d_side`
    - If enemy: track nearest direction `h_enemy` and distance `r_enemy`

- Orientation target: `u_des = normalize( w_enemy*h_enemy + w_align*normalize(u_align) + w_coh*normalize(project(coh,n_i)) )`
  - Rotate `u_i` toward `u_des` with angular step ≤ `max_turn*dt`; normalize.

- Velocity update: `dv = w_sep*sep + v_lane_adjust + enemy_spacing_adjust`
  - Friction: `v_i = (1 - friction*dt)*v_i + dv*dt`
  - Optionally bias motion along facing `u_i`; clamp `|v_i| ≤ max_speed`.

- Integrate: `p_i += v_i * dt`
- Bounds: wrap or bounce inside [-1,1]^2 (MVP: simple clamp/bounce)
- Write back to buffers.

Notes:
- Anisotropy is driven by `w_lat = 1 - |dot(h, u_i)|` favoring lateral neighbors to form lines.
- Keep epsilons to avoid NaNs; minimize branches.
- O(N^2) is fine for small N; grid neighbor search later.

---

## 5. Rendering

MVP: GL_POINTS with custom VS.

- Vertex shader: `python/GLCL2/shaders/soldier_points.glslv`
  - `layout(location=0) in vec4 pos_dir;`
  - `gl_Position = vec4(pos_dir.xy, 0.0, 1.0);`
  - Optionally `gl_PointSize` from uniform.

- Fragment shader: reuse `python/GLCL2/shaders/monocolor.glslf` initially.

### Formation rendering (instanced rectangles)

- Goal: render each formation as an oriented rectangle sized by formation widths.
- Data source: `formation_params` (`float8` per formation): `(o.x,o.y,d.x,d.y,w_par,w_perp,k_par,k_perp)`.
- Shader: `python/GLCL2/shaders/pose2d_rect.glslv`
  - Attributes
    - `location=0` per-vertex: `vec2 vpos` corners of a unit quad in model space (use `quad_verts` with 4 verts for TRIANGLE_STRIP)
    - `location=1` per-instance: `vec4 inst_pose = (px,py,ux,uy)` sourced from `formation_params[0:4]`
    - `location=2` per-instance: `vec4 inst_misc = (w_par,w_perp,k_par,k_perp)` sourced from `formation_params[4:8]` (only `.xy` used)
  - Transform: `world = p + u*(vpos.x*w_par) + v*(vpos.y*w_perp)`, where `v=perp(u)`
- GL config:
  - Vertex buffer: `quad_verts` with layout `[2]`
  - Instancing: `instance_buffer="formation_params"`, `instance_attribs=[4,4]`, `instance_divisor=1`
  - Draw mode: `TRIANGLE_STRIP`, instances = `formation_count`

Notes:
- This keeps a clear separation from soldier rendering, which stays as points for now.
- Buffer layout matches OpenCL: no repacking; GLSL reads two `vec4`s from one `float8` buffer.

Phase 1.5: Oriented point sprites (no engine change)
- FS uses `gl_PointCoord` and per-vertex `u` (passed from VS) to draw a thin oriented line or a T-shape glyph in screen space.

Team colors (later): separate passes per team or extend attribute handling.

### Instanced rendering (general primer)

- __Concept__: The vertex shader runs `vertexCount * instanceCount` times. Per-vertex attributes advance with the vertex index, per-instance attributes advance with the instance index.
- __API__: `glDrawArraysInstanced(mode, first, vertexCount, instanceCount)` (or `glDrawElementsInstanced`) triggers the nested loop over vertices × instances.
- __Divisors__: `glVertexAttribDivisor(loc, d)` controls how often an attribute advances. `d=0` (default) = per-vertex; `d=1` = per-instance; `d>1` = advance once per `d` instances.
- __Typical layout__:
  - location 0: per-vertex mesh data (e.g., `vec2 vpos` corners of a unit quad), divisor 0
  - location 1..k: per-instance transforms/params (e.g., `vec4 pos_dir`, widths), divisor 1
- __Setup sketch__:
  - Bind VAO; bind mesh VBO; `glVertexAttribPointer(0, ...)`; `glVertexAttribDivisor(0, 0)`
  - Bind instance VBO; `glVertexAttribPointer(1, ...)`; `glVertexAttribDivisor(1, 1)` (and similarly for other instance attrs)
  - Call `glDrawArraysInstanced(..., vertexCount, instanceCount)`
- __Shader side__: Use separate inputs for vertex vs instance attributes; `gl_InstanceID` is also available for indexing per-instance resources (SSBO/texture) if needed.
- __VBOs__: Mesh and instance data can live in separate VBOs (common) or a single interleaved VBO; only the divisor determines stepping behavior.

---

## 6. Script Wiring (GLCL config)

File: `python/GLCL2/scripts/soldiers.py`

- `config` keys:
  - `parameters`: as in section 3.
  - `buffers`:
    - `"state_pos_dir": ("particle_count", 4, "f4")`
    - `"state_vel_mr":  ("particle_count", 4, "f4")`
    - `"state_team_type": ("particle_count", 4, "i4")`
    - `"formation_params": ("formation_count", 8, "f4")`  // float8 as used by CL
    - `"quad_verts": (4, 2, "f4")`  // local quad geometry used by instanced rectangle
    - Optional for visualization: `"formation_lines"`, `"formation_points"`
  - `opencl_source`: `["../cl/soldiers.cl"]`
  - `kernels`: define `soldiers_step` with appropriate buffer/parameter binding
  - `kernel_pipeline`: `["soldiers_step"]`
  - `opengl_shaders`: `{ "soldier_render": ("../shaders/soldier_points.glslv", "../shaders/monocolor.glslf", ["color", "point_size"]), "pose2d_rect": ("../shaders/pose2d_rect.glslv", "../shaders/monocolor.glslf", ["color"]) }`
  - `render_pipeline`: `[ ("soldier_render", "particle_count", "state_pos_dir", None), ("pose2d_rect", "formation_count", "quad_verts", None, {mode: TRIANGLE_STRIP, attribs:[2], instance_buffer:"formation_params", instance_attribs:[4,4], instance_divisor:1}) ]`

### Update (instancing for formations)

- Buffers added:
  - `"formation_params": ("formation_count", 8, "f4")`  // float8 as used by CL
  - `"quad_verts": (4, 2, "f4")`  // local quad geometry used by instanced rectangle
  - Optional for visualization: `"formation_lines"`, `"formation_points"`
- Shaders added:
  - `"pose2d_rect": ("../shaders/pose2d_rect.glslv", "../shaders/monocolor.glslf", ["color"])`
- Render pipeline additions:
  - `( "line_color",    "formation_point_count", "formation_lines", None, {mode: LINES, attribs:[2,4]} )`  // goal lines
  - `( "pose2d_rect",  "formation_count",       "quad_verts",       None, {mode: TRIANGLE_STRIP, attribs:[2], instance_buffer:"formation_params", instance_attribs:[4,4], instance_divisor:1} )`
- Soldier oriented quads (`instanced_quad`) are intentionally disabled in the pipeline for clarity while focusing on formations; soldiers render as points.

 `init()`:
  - Build 10+10 soldiers: two lines with small jitter, opposing `u` (e.g., team 0 at y=-0.4 facing +x; team 1 at y=+0.4 facing -x)
  - Zero velocity, set reasonable radii (0.01–0.02)

---

{{ ... }}

- Team 0: positions along x in [-0.6, -0.1], y=-0.4 with jitter; `u=(+1,0)`
- Team 1: positions along x in [0.1, 0.6],  y=+0.4 with jitter; `u=(-1,0)`
- `v=(0,0)`, `m=1`, `r=0.015`

---

## 8. Debugging and Stability

- Start with very small timestep and weights; gradually increase.
- Clamp `max_speed`, `max_turn`; always normalize `u`.
- Add temporary `printf` for a few particles in the kernel (commented out in commit).
- Prefer fail-loudly inside kernel during development.

---

## 9. Performance Roadmap

- MVP O(N^2) neighbor loop is fine for N<=1k.
- Phase 3: Uniform grid / hashed cells in CL for O(N) neighbor queries.
- Data-oriented memory layout (SoA-ish across buffers) already in place.

---

## 10. Phases and Files

Phase 1 (MVP):
- Create:
  - `python/GLCL2/cl/soldiers.cl` (kernel)
  - `python/GLCL2/shaders/soldier_points.glslv` (VS)
  - `python/GLCL2/scripts/soldiers.py` (script wiring + init)
- Use existing `python/GLCL2/shaders/monocolor.glslf` (FS)

Phase 1.5:
- Add oriented sprite FS `soldier_points_oriented.glslf` using `gl_PointCoord` and `u`.

Phase 2:
- Team color rendering (split passes or attribute extension).

Phase 3:
- Spatial hashing / uniform grid neighbor search.

Phase 4:
- Leader behavior, morale/stamina/health buffers, terrain influence texture; optional type-specific params via `type_params_*` buffers.

---

## 11. Assumptions and Compatibility

- Rendering path relies on single-attribute VBO per buffer (see `GLCLGUI.py::bake_render_objects()`), hence packing `p,u` into `vec4`.
- GL state: `GL_PROGRAM_POINT_SIZE` enabled by `GLCLGUI.initializeGL()`.
- Browser handles syncing buffer data from CL to GL for buffers used in render pipeline.
- Parameters map to kernel arguments and shader uniforms by name via `GLCLBrowser`.

---

## 12. Quick Tuning Defaults (initial)

- `dt = 0.01`
- `r_cut = 0.25`
- `d_front = 0.2`
- `d_side = 0.06`
- `w_sep = 1.0`, `w_align = 0.8`, `w_coh = 0.3`, `w_enemy = 1.2`
- `max_speed = 0.5`
- `max_turn = 2.0` (rad/s)
- `friction = 0.5`
- `side_bias = 1.5`
- `front_bias = 0.5`

These are starting points; adjust interactively.

---

## 13. Testing

- Start with 10 vs 10.
- Verify: line stabilization, facing toward nearest enemy, spacing maintenance.
- Instrument: frame-by-frame observation; temporary prints for aggregates for a few frames.

---

## 14. Next Actions

1) Implement Phase 1 files.
2) Run in GLCLBrowser and tune weights until visible line formation.
3) Add oriented glyph FS (Phase 1.5).

---

## 15. Implementation Progress (2025-08-17)

- __Files created__
  - `python/GLCL2/cl/soldiers.cl` – Phase 1 kernel implementing separation, alignment (lateral anisotropy), cohesion, enemy facing, kinematics, and bounds.
  - `python/GLCL2/scripts/soldiers.py` – `config` with parameters, buffers, kernel pipeline, shaders, render pipeline; `init()` builds 10+10 soldiers in two lines with jitter; sets `m=1`, `r≈0.015`, initial `u=±x`, `v=0`.
  - `python/GLCL2/shaders/soldier_points.glslv` – simple VS reading `pos_dir` and setting `gl_Position`; default `gl_PointSize=5` for visibility.
  - `python/GLCL2/shaders/pose2d_rect.glslv` – instanced rectangle VS using `formation_params` split into two `vec4` instance attributes.

- __Run status (GLCLBrowser)__
  - Successfully loads config, allocates CL buffers (`state_pos_dir`, `state_vel_mr`, `state_team_type`).
  - Renders via GL_POINTS using `state_pos_dir` as VBO; FS `monocolor.glslf` used.
  - Initial kernel launch failed with `INVALID_WORK_GROUP_SIZE` when `local_size=(32,)` and `global=(20,)`.
  - __Fix__: set `local_size=None` in `config["kernels"]["soldiers_step"]` in `python/GLCL2/scripts/soldiers.py`, letting driver choose. Kernel now executes each frame.
  - Formation rendering added: instanced rectangles draw one per formation from `formation_params`; soldier oriented quads pass disabled.
  - __Update (2025-08-18) – Instanced quads enabled and fixed__:
    - Problem: soldier instanced quads didn’t move although soldier points did. Two issues:
      1) VAO collision: both formations and soldiers used `quad_verts`, causing VAO reuse and wrong instance buffer bound.
      2) Buffer sync: `GLCLGUI.update_buffer_data()` updated only the vertex VBO owner, not an instance VBO owner when a buffer was shared (e.g., `state_pos_dir`).
    - Fixes:
      - Added separate mesh VBO `quad_verts_soldier` and used it in soldier quad pass to get a dedicated VAO.
      - Modified `python/GLCL2/GLCLGUI.py::update_buffer_data()` to update both the vertex VBO (`gl_objects[buf]`) and the instance VBO owner (`instance_owners[buf]`).
    - Result: soldier instanced quads now follow `state_pos_dir` positions/orientations (match GL points). Formation rectangles remain correct.

- __Data/Sync__
  - CL→GL sync: `state_pos_dir` is flagged for rendering and synchronized each frame by GLCL.
  - GLobject baked with 20 elements; device chosen: NVIDIA GeForce GTX 1650 (via platform selection).
  - Instance buffer sync: `formation_params` is referenced by `render_pipeline` instancing options and is auto-synchronized CL→GL by `GLCLBrowser`.

- __Parameters exposed__
  - Matches Section 3: `particle_count, dt, r_cut, d_front, d_side, w_sep, w_align, w_coh, w_enemy, max_speed, max_turn, friction, side_bias, front_bias`.

- __Defaults and tuning__
  - Start values as in Section 12. For visible line formation: increase `w_align` and `side_bias`, modestly raise `w_coh`, adjust `d_side`, and tune `w_enemy` with `d_front`.

- __Known limitations (Phase 1)__
  - O(N^2) neighbor loop; simple bounds; orientation not yet visualized as glyph; single color per pass.

- __Next steps__
  - Tune parameters interactively (t5). Add additional GUI params if needed (t6). Implement oriented glyph rendering (t7). Plan uniform grid neighbor search (t8).

### Task Checklist

- [x] Draft design doc for GLCL Soldiers (2D army) simulation
- [x] Create OpenCL kernel `soldiers.cl` implementing 2D dynamics with lateral anisotropy
- [x] Create user script `soldiers.py` wiring buffers, params, kernel, and render pipeline
- [x] Generate small init (10+10 soldiers) with two teams; place lines and random jitter
- [ ] Run in GLCLBrowser and tune weights for visible line formation (in progress)
- [ ] Add GUI parameters for weights and cutoffs; expose dt, enemy distance and line spacings
- [ ] Plan rendering upgrade to oriented glyphs (T-shape)
- [ ] Performance follow-up: implement uniform grid neighbor search for O(N)
- [ ] Document future extensions: leaders, morale/stamina buffers, terrain influence

__Instancing debug/fix (2025-08-18)__

- [x] Separate soldier and formation VAOs by introducing `quad_verts_soldier`
- [x] Enable soldier instanced quads using `state_pos_dir` as instance buffer
- [x] Fix shared buffer sync in `GLCLGUI.update_buffer_data()` so shared buffers update both vertex and instance VBOs

---

## 16. Formation Goals (Desire Potentials)

__Motivation__
- Soldiers currently spread from initial positions due to separation and lack of strong global objectives.
- Introduce per-formation “goal lines” acting as a harmonic desire potential. They are not physical forces; they bias desired orientation and acceleration subject to physical limits.

__Data additions__
- Buffer `formation_params`: float8 per formation
  - `.xy = o` (origin)
  - `.zw = d` (normalized direction)
  - `.s0 = w_par` (longitudinal half-width)
  - `.s1 = w_perp` (lateral half-width)
  - `.s2 = k_par` (longitudinal stiffness)
  - `.s3 = k_perp` (lateral stiffness)
  - Note: two formations (one per team) for MVP.
- Mapping soldier→formation:
  - Use `state_team_type.z` as `formation_id` (float), or derive from `team_id` (x) for 1:1 team→formation.

__Kernel integration (`python/GLCL2/cl/soldiers.cl`)__
- For soldier i with position `p` and formation (o,d,w_par,w_perp,k_par,k_perp):
  - Nearest point on line: `q = o + dot(p - o, d) * d`
  - Error vector: `e = p - q`, decompose: `e_par = dot(e, d) * d`, `e_perp = e - e_par`
  - Harmonic desire acceleration (clamped by widths):
    - `a_goal = -k_par * clamp_len(e_par, w_par) - k_perp * clamp_len(e_perp, w_perp)`
  - Blend into steering and velocity update:
    - Orientation target term: add `w_goal * normalize(-e_perp + beta * d * sign(-dot(e, d)))` to `u_des` mix to keep soldiers on the line and roughly along it.
    - Acceleration term: add `desire_gain * a_goal` to `dv` (before friction and clamps), then cap by `max_accel` and `max_speed`.
  - Keep collision radius `r` much smaller than neighbor/align radius `r_cut` (e.g., `r ≈ r_cut/5`).

__Parameters (GUI)__
- `desire_gain: float` (scales a_goal into dv)
- `w_goal: float` (weight into orientation blend)
- `beta: float` (balance longitudinal vs lateral bias in orientation target)
- `max_accel: float` (cap on |dv/dt|)
- Per-formation constants are in `formation_params` (widths, stiffness);
  optional global multipliers: `k_goal_mul_par`, `k_goal_mul_perp`.

__Rendering goal lines__
- Add simple GL line pass to visualize two formations:
  - Build buffer `formation_lines` with 2 vertices per line: `o - L*d`, `o + L*d` (L≈1.2 for full-screen span).
  - Shader `line2d.glslv` + `monocolor.glslf`; set `gl_Position` from positions; use `glLineWidth` for visibility.
  - Optionally draw lateral width as translucent band later.

__Initial defaults (MVP)__
- `w_par=0.20`, `w_perp=0.06`, `k_par=1.0`, `k_perp=2.0`
- `desire_gain=1.0`, `w_goal=0.8`, `beta=0.5`, `max_accel=3.0`
- `r_cut=0.25`, `r≈0.05` (set per soldier radius), `friction=0.5`

__Notes__
- Desire affects orientation and acceleration targets but final motion remains constrained by `max_turn`, `max_speed`, `max_accel`, and collisions/separation.
- Keep numerical stability: clamp/soften with widths to avoid large instantaneous accelerations when far from the line.

__Implementation plan__
1) Add buffers and params in `python/GLCL2/scripts/soldiers.py` (`formation_params`, optional `formation_lines` for rendering).
2) Initialize two formations (one per team) in `init()`; fill buffers and render pipeline for lines.
3) Extend kernel signature to take `__global float8* formation_params` and per-soldier `formation_id` (or compute from `team_id`).
4) Implement `a_goal` and steering mix; add `desire_gain`, `w_goal`, `beta`, `max_accel` params.
5) Tune interaction: ensure `r` (collision) << `r_cut` (alignment); adjust `k_par,k_perp` and gains.
