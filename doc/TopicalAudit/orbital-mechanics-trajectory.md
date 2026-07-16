---
type: TopicalAudit
title: Orbital Mechanics & Trajectory Optimization
tags: [topic, cross-language, orbital, trajectory, optimization, spline, n-body, ode]
---

## Summary

Orbital mechanics with Keplerian elements, N-body gravitational simulation, and trajectory optimization via spline-based variational methods. C++ `SpaceBodies.h` defines `OrbitalElements`, `Orbit`, `SpaceBody` with Keplerian propagation and astorb database loading. `SpaceWorld.h` implements N-body integration via RKF45. Trajectory optimization uses cubic B-spline control points optimized by FIRE (`DynamicOpt.h`) with `TrajectoryVariation.h` objective. `SplineManager.h` for non-uniform hermite spline evaluation. Application: `OrbitalWar`/`SpaceTactics` game. Python Jupyter notebooks for combat calculations.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/common/Orbital/SpaceBodies.h` | active | `OrbitalElements` (semi-major, eccentricity, inclination, arg_periapsis, node_longitude, mean_anomaly). `Orbit` with vis-viva equation, true anomaly solver, `pointAtEpoch()`, `toPoints()`. `SpaceBody` with trajectory arrays (`trjPos`, `trjThrust`), hermite spline interpolation. `SpaceBodyIntegrator` for single-body RKF45. `fromString_astorb()` for asteroid database loading. |
| C++ | `cpp/common/Orbital/SpaceWorld.h` | active | N-body simulation: `getDerivODE()` computes gravitational forces + ship thrust. `predictTrjs()` integrates all bodies. `ODEintegrator_RKF45` with adaptive timestep. Dual representation: Keplerian orbits for planets + N-body for ships. `load_astorb()` for real asteroid data. |
| C++ | `cpp/common/dynamics/ODEintegrator.h` | active | `ODEintegrator` base (Euler) + `ODEintegrator_RKF45` (Runge-Kutta-Fehlberg 4(5) with adaptive step). Supports both function pointer and virtual object derivatives. |
| C++ | `cpp/common/math/spline_hermite.h` | active | Hermite spline basis functions: `basis()`, `dbasis()`, `ddbasis()`, `curve_point()` |
| C++ | `cpp/common/math/SplineManager.h` | active | Non-uniform hermite spline manager: control points + optional explicit derivatives. `eval()`, `evalUniform()`, `insertPoint()`, `removePoint()`. Binary search for time index. |
| C++ | `cpp/common/math/CubicBSpline.h` | active | Cubic B-spline basis: `val()`, `dval()`, `ddval()`, `basis()`, `dbasis()`, `ddbasis()` |
| C++ | `cpp/common/Orbital/TrajectoryVariation.h` | active | `getTrajectoryVariation()` — objective function for trajectory optimization: minimizes thrust integral along spline control points |
| C++ | `cpp/common/dynamics/DynamicOpt.h` | active | FIRE (Fast Inertial Relaxation Engine) optimizer: `bindArrays()`, `move_FIRE()`. Used for trajectory optimization. |
| C++ | `cpp/apps/OrbitalWar/test_OptContinuousThrust.cpp` | active | Trajectory optimization demo: 20 B-spline control points, FIRE optimizer, `getTrajectoryVariation` objective. Visualizes control points + forces. |
| C++ | `cpp/apps/OrbitalWar/spaceTactics.cpp` | active | Space combat game: N-body physics, trajectory splines, time scrubbing, weapon tests |
| C++ | `cpp/common/Orbital/spaceCombat.h` | active | Combat assembly: guns, projectiles, hit detection, damage evaluation |
| C++ | `cpp/common/Orbital/asteroidEngineering.h` | active | Asteroid composition (mineral/rock types), orbital maneuver cost (`manuever_planeChange`) |
| C++ | `cpp/common/dynamics/appliedPhysics.h` | active | `centralGravityForce()`, physical constants (G, AU) |
| Python | `projects/SpaceCombat/ch1_basics.ipynb` | active | Orbital mechanics basics notebook |
| Python | `projects/SpaceCombat/Manueverability.ipynb` | active | Spacecraft maneuverability analysis |
| Python | `projects/SpaceCombat/Weapons.ipynb` | active | Weapon performance calculations |
| Python | `projects/SpaceCombat/MISO_paper_draft.ipynb` | active | MISO (Multi-Input Single-Output) paper draft for space combat |
| Python | `projects/SpaceCombat/Industry.ipynb` | active | Space industry/economy model |
| Python | `projects/SpaceCombat/py/` | active | Python support modules |
| Doc | `docs/SpaceTactics/SpaceTactics.md` | active | Comprehensive code map: N-body, weapons, propulsion, shields |
| Doc | `docs/SpaceTactics/SpaceBodies.md` | active | Celestial body definitions |
| Doc | `docs/SpaceTactics/SpaceWorld.md` | active | World model documentation |
| Doc | `docs/SpaceTactics/SplineManager.md` | active | Spline trajectory management |
| Doc | `docs/SpaceTactics/appliedPhysics.md` | active | Applied physics models |
| Doc | `docs/SpaceTactics/asteroidEngineering.md` | active | Asteroid mining/engineering |
| Doc | `docs/SpaceTactics/spaceCombat.md` | active | Combat model documentation |

## Sub-topics

### Keplerian Orbital Mechanics

- `OrbitalElements` → `Orbit` conversion via Euler angles (`fromEuler_orb`)
- True anomaly from mean anomaly: iterative Kepler equation solver (`true_anomaly()`)
- Vis-viva equation: `v = sqrt(GM*(2/r - 1/a))`
- Orbital point at epoch: `pointAtEpoch()` → 3D position from true anomaly + orbit rotation
- Epoch shifting: `shift_epoch()` updates mean anomaly
- Real asteroid loading from astorb database format

### N-Body Gravitational Simulation

- `SpaceWorld::getDerivODE()` — gravitational forces between all bodies + ship thrust
- RKF45 integrator with adaptive timestep (SAFETY=0.2, PGROW=1.2, PSHRINK=0.7)
- Trajectory storage: `trjPos[]` arrays per body, hermite spline interpolation
- Dual representation: planets on Keplerian orbits, ships on N-body trajectories
- `predictTrjs()` — integrate forward, store positions at each timestep

### Trajectory Optimization

- **B-spline control points** represent trajectory
- **FIRE optimizer** (`DynamicOpt.h`) minimizes objective
- **`getTrajectoryVariation()`** — thrust integral objective (minimize total thrust)
- `test_OptContinuousThrust.cpp` — demo with 20 control points, real-time visualization
- Commented-out spline editor code suggests interactive trajectory editing was planned

## Parity Status

- **C++ N-body ↔ Python notebooks**: Notebooks perform analytical calculations (delta-v, maneuverability) but don't simulate N-body dynamics. No direct parity.
- **Keplerian orbits ↔ N-body integration**: Dual representation in `SpaceWorld`. Planets use analytical Keplerian, ships use numerical integration. No consistency check found.
- **`test_OptContinuousThrust` ↔ `SpaceTactics` game**: Optimization demo is standalone; game uses spline trajectories but optimization not integrated into gameplay.

## Open Issues

- `SplineManager.insertPoint()` is a placeholder (`return -1`)
- `SpaceBodyIntegrator` thrust interpolation commented out (`//f.add()` in `getDerivODE`)
- `test_OptContinuousThrust.cpp` has commented-out interactive editing (mouse handling, event handling)
- No automated test for Keplerian ↔ N-body consistency
- `appliedPhysics.h` location inconsistent — included from `cpp/common/dynamics/` but referenced as `cpp/common/Orbital/`
- Trajectory optimization uses FIRE but no comparison with other optimizers (CG, L-BFGS)
- No multi-body trajectory optimization (only single spacecraft)
- `projects/SpaceCombat/` notebooks are standalone calculations, not integrated with C++ simulation
- See `doc/Markdown/cpp/common/Orbital/` for auto-generated API docs

---

## Algorithm Review: Orbital Mechanics & Trajectory Optimization

This section reviews algorithms for orbital propagation, trajectory representation, and trajectory optimization used in the codebase. For each algorithm we describe: merit, use cases, strengths, weaknesses, GPU parallelization potential, and implementation status across platforms (C++, Python).

### A. Orbital Propagation Methods

#### A1. Keplerian Analytical Propagation (Closed-Form from Orbital Elements)

- **Merit**: Given the 6 classical orbital elements (semi-major axis $a$, eccentricity $e$, inclination $i$, argument of periapsis $\omega$, longitude of ascending node $\Omega$, mean anomaly $M$), compute the 3D position at any epoch analytically. No numerical integration needed — just solve the Kepler equation and apply the orbit rotation matrix.
- **Use cases**: Propagating planets, asteroids, and any body on a known Keplerian orbit. Loading real asteroid databases (astorb format). Rendering orbit ellipses. Predicting planet positions for combat games (no N-body needed for planets).
- **Strengths**: Exact (to Kepler equation solver tolerance). O(1) per evaluation — no time-stepping. Energy and angular momentum conserved by construction. Can jump to any epoch without integrating from start. Simple: `pointAtEpoch()` → `true_anomaly()` → polar-to-Cartesian with rotation matrix.
- **Weaknesses**: Only valid for two-body problem (one central mass). No perturbations (J2 oblateness, drag, third-body). No thrust. Kepler equation ($M = E - e\sin E$) requires iterative solver (fixed-point, ~5–10 iterations for $e < 0.9$). Fails for parabolic/hyperbolic orbits ($e \geq 1$) — current code uses `atan(sqrt((1+e)/(1-e))*tan(E/2))` which is singular at $e=1$.
- **GPU parallelization**: Excellent — each body is independent. Kepler solver is a fixed-iteration loop. Already trivially parallel for thousands of asteroids. No GPU implementation exists but would be straightforward.
- **Implemented**:
  - C++: `OrbitalElements` struct in `cpp/common/Orbital/SpaceBodies.h:145-169`. `Orbit` struct at `:172-279` with `pointAtEpoch()`, `toPoints()`, `toPoints_epochs()`, `speed()` (vis-viva), `evalAtPhi()`. `true_anomaly()` fixed-point solver at `:119-131`. `mean_anomaly_at_epoch()` at `:138-143`. `semi_major2period()` (Kepler's third law) at `:133-135`. `fromString_astorb()` for real asteroid loading at `:343-401`.
  - Python: Notebooks in `projects/SpaceCombat/` perform analytical orbit calculations but do not implement a full Keplerian propagator.
  - Status: **Active, C++ only.** The primary propagation method for planets and asteroids in `SpaceWorld`.

#### A2. N-Body Gravitational Integration (RKF45)

- **Merit**: Integrate Newton's equations of motion for N gravitationally interacting bodies: $\ddot{\mathbf{r}}_i = \sum_{j \neq i} G m_j (\mathbf{r}_j - \mathbf{r}_i) / |\mathbf{r}_j - \mathbf{r}_i|^3$. Uses RKF45 (Runge-Kutta-Fehlberg 4(5)) with adaptive timestep. State vector: 6 DOF per body (position + velocity), packed as flat `double[]` array.
- **Use cases**: Spacecraft trajectory simulation with thrust. Combat dynamics where ships maneuver under gravity from multiple planets. Any scenario where thrust or perturbations make Keplerian propagation invalid.
- **Strengths**: Handles arbitrary number of gravitational sources. Adaptive timestep (SAFETY=0.2, PGROW=1.2, PSHRINK=0.7) — shrinks dt near close encounters, grows during coasting. Dual callback: function pointer or virtual `ODEderivObject`. Stores trajectory arrays (`trjPos[]`) for later spline interpolation and visualization.
- **Weaknesses**: O(N²) force evaluation per step (all pairs). Not symplectic — energy drifts over long simulations. RKF45 is 6-stage (expensive per step). No regularization for close encounters (KS-transform, Encke). Fixed `dt_max=0.1` may be too large for close approaches. Thrust interpolation is linear (`trjThrust[i]*(1-du) + trjThrust[i+1]*du`) — not smooth.
- **GPU parallelization**: Good — force computation is O(N²) parallel (each (i,j) pair independent). Adaptive step control is sequential. For N < 100, CPU is adequate. For N > 1000 (asteroid belt), GPU would help. No GPU implementation exists.
- **Implemented**:
  - C++: `SpaceWorld::getDerivODE()` in `cpp/common/Orbital/SpaceWorld.h:122-172`. `predictTrjs()` at `:290-304` — main integration loop. `ODEintegrator_RKF45` from `cpp/common/dynamics/ODEintegrator.h`. `centralGravityForce()` from `cpp/common/dynamics/appliedPhysics.h`.
  - Python: No N-body simulator. Notebooks do analytical calculations only.
  - Status: **Active, C++ only.** Used for ship trajectories in `OrbitalWar`/`SpaceTactics`.

#### A3. Dual Representation (Keplerian + N-Body Hybrid)

- **Merit**: Planets propagated analytically via Keplerian elements (A1); ships propagated numerically via N-body integration (A2) with planets as fixed gravitational centers. Combines the exactness and speed of analytical propagation for unperturbed bodies with the flexibility of numerical integration for maneuvering spacecraft.
- **Use cases**: Space combat games where planets follow fixed orbits but ships maneuver freely. Any scenario with a clear separation between "heavy bodies on known orbits" and "light bodies that maneuver."
- **Strengths**: Planets never drift (analytical). Ships get full thrust modeling. O(N_planets × N_ships) force evaluation — not O((N_p + N_s)²). Planet positions available at any epoch without integration (for targeting, rendering).
- **Weaknesses**: No mutual perturbation between planets (assumes Keplerian). No planet-ship gravitational back-reaction (ships don't affect planets — fine for mass ratios). Inconsistency risk: if a ship's numerical trajectory diverges from what Keplerian would predict for the same initial conditions, there's no check. Planet positions in `getDerivODE()` use `pvs` (ODE state array) not `orbit->pointAtEpoch()` — so planets ARE numerically integrated, not analytically propagated. The dual representation exists in the data structures but `predictTrjs()` integrates everything numerically.
- **GPU parallelization**: Same as A2 — force computation is parallel.
- **Implemented**:
  - C++: `SpaceWorld` in `cpp/common/Orbital/SpaceWorld.h`. `planets[]` and `ships[]` vectors. Both are integrated in `predictTrjs()` via `ode.step()`. `SpaceBody::pointAtEpoch()` exists for analytical propagation but is used only for `pickPlanet()` (ray-picking), not for trajectory integration.
  - Status: **Partially implemented.** Data structures support dual representation, but `predictTrjs()` integrates all bodies numerically. True hybrid (analytical planets + numerical ships) is not used in the integration loop.

#### A4. Single-Body Multi-Center Integration (SpaceBodyIntegrator)

- **Merit**: Integrate a single body under gravity from multiple fixed centers (planets on pre-computed trajectories). Simpler than full N-body: one body, multiple gravitational sources. The body's thrust can be interpolated from `trjThrust[]`.
- **Use cases**: Re-integrating a single ship's trajectory with different parameters (e.g., thrust profile) while keeping planet trajectories fixed. Testing trajectory variations.
- **Strengths**: O(N_centers) per step — much cheaper than full N-body when only one body needs re-integration. Can use pre-computed planet positions (from spline interpolation). Inherits RKF45 adaptive stepping.
- **Weaknesses**: Thrust interpolation commented out (`//f.add()` in `getDerivODE`) — not functional. Planet positions use `centers[i]` directly but interpolation from trajectory splines is not implemented (marked as TODO: `// d = ( interpolate spline )`). Only single-body — no ship-ship interactions.
- **GPU parallelization**: Poor for single body. Could parallelize over multiple trajectory variations (each variation is independent).
- **Implemented**:
  - C++: `SpaceBodyIntegrator` in `cpp/common/Orbital/SpaceBodies.h:411-446`. Inherits both `ODEderivObject` and `ODEintegrator_RKF45`. `evalTrj()` at `:435-444`.
  - Status: **Incomplete.** Thrust and planet position interpolation not implemented. Framework exists but not used in production.

---

### B. Time Integration Schemes

#### B1. RKF45 (Runge-Kutta-Fehlberg 4(5) with Adaptive Step)

- **Merit**: 6-stage embedded Runge-Kutta method producing both 4th and 5th order solutions. The difference gives an error estimate for adaptive timestep control. Widely used for non-stiff ODEs where accuracy is important.
- **Use cases**: Orbital mechanics (primary integrator). Any non-stiff dynamical system requiring adaptive accuracy.
- **Strengths**: Adaptive — automatically shrinks dt when error is large (close encounters, high thrust) and grows dt when error is small (coasting). Error estimate is "free" (embedded method — no extra function evaluations). 4th-order accuracy with 5th-order error estimate. Well-established (Fehlberg, 1969).
- **Weaknesses**: Not symplectic — energy drifts over long orbital simulations (thousands of orbits). 6 function evaluations per step (vs. 1 for Euler, 2 for Leapfrog). No angular momentum conservation. Adaptive dt can cause issues with trajectory storage (non-uniform time grid). `MAX_STEPS=10000` hardcoded — may be too low for long simulations. `invMaxYerr[]` must be set by user — not auto-initialized.
- **GPU parallelization**: Poor for single trajectory (sequential). Good for ensemble (many independent trajectories — e.g., Monte Carlo trajectory optimization). No GPU implementation.
- **Implemented**:
  - C++: `ODEintegrator_RKF45` in `cpp/common/dynamics/ODEintegrator.h`. Used by `SpaceWorld` and `SpaceBodyIntegrator`. Parameters: SAFETY=0.2, PGROW=1.2, PSHRINK=0.7, MAX_ADAPT=10, dt_min=0.001, dt_max=0.1.
  - Python: No RKF45 implementation. Notebooks use analytical formulas.
  - Status: **Active.** The only ODE integrator used for orbital mechanics. See also `ode-integration.md`.

#### B2. Fixed-Step Euler (Base Class)

- **Merit**: Simplest integrator: $\mathbf{Y}_{n+1} = \mathbf{Y}_n + \dot{\mathbf{Y}} \cdot dt$. One function evaluation per step. Base class for `ODEintegrator_RKF45`.
- **Use cases**: Quick prototyping. Non-critical simulations. Testing derivative functions.
- **Strengths**: Trivial. Fastest per step. No error estimation overhead.
- **Weaknesses**: First-order — error O(dt). Energy drifts rapidly. Unstable for oscillatory systems unless dt << orbital period. Not used for any production orbital simulation.
- **GPU parallelization**: Excellent — trivially parallel for independent trajectories.
- **Implemented**:
  - C++: `ODEintegrator` base class in `cpp/common/dynamics/ODEintegrator.h`. `step()` method.
  - Status: **Exists but not used** for orbital mechanics. RKF45 is always preferred.

---

### C. Trajectory Representation & Interpolation

#### C1. Hermite Spline Interpolation (4-Point, Non-Uniform)

- **Merit**: Interpolate a smooth curve through discrete trajectory points using cubic Hermite basis. 4-point stencil: `curve_point(u, p[i-1], p[i], p[i+1], p[i+2])`. Produces C¹ continuous trajectory from discrete `trjPos[]` samples.
- **Use cases**: Rendering smooth ship trajectories from RKF45 output. Interpolating positions between stored timesteps for combat calculations (hit detection, targeting). Time-scrubbing in the game UI.
- **Strengths**: Smooth (C¹). Local support — changing one point only affects nearby intervals. No global solve needed. Works with non-uniform time grid (via `findNextInd()` binary search). Fast evaluation — 4 multiplications + 3 additions per basis.
- **Weaknesses**: Not C² (acceleration discontinuous at knots — may matter for thrust interpolation). 4-point stencil needs 2 extra points at boundaries (clamped in `validTrjIndex()`). No derivative control — tangents implicit from neighboring points. `findNextInd()` is linear scan, not binary search (O(N) per lookup).
- **GPU parallelization**: Good — each interpolation is independent. Already used in real-time rendering.
- **Implemented**:
  - C++: `Spline_Hermite::curve_point()` in `cpp/common/math/spline_hermite.h`. Used by `SpaceBody::getTrjPos()` in `SpaceBodies.h:314-318`. `nonUni2spline()` at `:85-94` for non-uniform grid.
  - Status: **Active.** Primary trajectory interpolation method.

#### C2. Cubic B-Spline (Uniform, Local Support)

- **Merit**: Represent a trajectory as a sum of cubic B-spline basis functions weighted by control points: $\mathbf{x}(u) = \sum_i \mathbf{c}_i B_i(u)$. Each basis has local support (4 spans). Second derivative is continuous (C²). Control points don't lie on the curve — they "pull" it.
- **Use cases**: Trajectory optimization — control points are the optimization variables. The smooth C² property ensures continuous thrust (thrust ∝ second derivative of position). `getTrajectoryVariation()` uses B-spline second derivatives to compute thrust.
- **Strengths**: C² continuous — smooth thrust profile. Local support — changing one control point affects only 4 spans. Well-conditioned basis. Second derivative (acceleration/thrust) is also a B-spline (lower order). Standard in CAD and trajectory optimization.
- **Weaknesses**: Control points don't pass through data points (no interpolation — only approximation). Needs boundary condition handling (code has ad-hoc correction at `i<5` and `i>ncoef-3`). Uniform parameterization only — no non-uniform B-spline. Not used for trajectory rendering (Hermite is used instead).
- **GPU parallelization**: Good — basis evaluation is local and independent per span.
- **Implemented**:
  - C++: `CubicBSpline` in `cpp/common/math/CubicBSpline.h`. `basis()`, `dbasis()`, `ddbasis()` for position, velocity, acceleration. Used by `getTrajectoryVariation()` in `TrajectoryVariation.h:58-143`. Boundary correction at `:92-108`.
  - Status: **Active.** Used only for trajectory optimization, not for rendering.

#### C3. Linear Interpolation (Thrust)

- **Merit**: Simplest interpolation: $\mathbf{T}(u) = \mathbf{T}_i (1-u) + \mathbf{T}_{i+1} \cdot u$. O(1), no extra points needed.
- **Use cases**: Interpolating thrust vectors between stored timesteps in `SpaceWorld::getDerivODE()`.
- **Strengths**: Trivial. Preserves thrust magnitude bounds. No overshoot.
- **Weaknesses**: C⁰ only — thrust discontinuities at timesteps. Not suitable for optimization (non-smooth objective). Not physically realistic (real thrust changes smoothly).
- **GPU parallelization**: Trivial.
- **Implemented**:
  - C++: `SpaceBody::getThrust()` in `SpaceBodies.h:332-339`. Used in `SpaceWorld::getDerivODE()` at `:152`.
  - Status: **Active but crude.** Thrust interpolation is linear — a known limitation.

---

### D. Trajectory Optimization

#### D1. Variational Thrust Minimization on B-Spline Control Points

- **Merit**: Formulate trajectory optimization as a variational problem: minimize total thrust $\int_0^T |\mathbf{T}(t)|^2 dt$ where $\mathbf{T} = \ddot{\mathbf{x}} - \mathbf{g}(\mathbf{x})$ (thrust = acceleration minus gravity). The trajectory $\mathbf{x}(u)$ is a cubic B-spline with control points $\mathbf{c}_i$ as optimization variables. Gravity $\mathbf{g}(\mathbf{x})$ is evaluated at B-spline-interpolated positions. Integration via 6-point Gauss-Legendre quadrature (exact for polynomials of degree ≤ 11). The gradient $\partial J / \partial \mathbf{c}_i$ is computed analytically by differentiating through the B-spline basis and gravity Jacobian.
- **Use cases**: Finding low-thrust trajectories between orbital states. Continuous-thrust spacecraft trajectory design. The demo (`test_OptContinuousThrust.cpp`) optimizes a 20-control-point trajectory.
- **Strengths**: Smooth objective (C² B-spline → continuous thrust). Analytic gradient — no finite differences. Gauss-Legendre quadrature is exact for the polynomial part (gravity is the only non-polynomial term). Local support means gradient computation is O(N × N_gauss × 4) per segment. Boundary conditions handled (ad-hoc correction at endpoints).
- **Weaknesses**: Gravity Jacobian (`addGravity()`) uses `GM = -1.0` hardcoded — TODO says "read mass of body." Gravity is evaluated relative to origin (no multi-body — `p` not `p - body_pos`). No constraints (no collision avoidance, no thrust limits, no final-state constraints). No multi-body gravity (single central body only). Boundary condition correction is ad-hoc (not derived from B-spline theory). `DEBUG_PLOT = true` hardcoded — OpenGL calls in compute function.
- **GPU parallelization**: Good — Gauss-Legendre quadrature points are independent. Gradient accumulation is a reduction. For Monte Carlo over multiple trajectories, each is independent.
- **Implemented**:
  - C++: `getTrajectoryVariation()` in `cpp/common/Orbital/TrajectoryVariation.h:58-143`. `addGravity()` at `:27-55` with gravity Jacobian. 6-point Gauss-Legendre quadrature at `:9-12`. `CubicBSpline::basis()`, `ddbasis()` for position and acceleration evaluation.
  - Status: **Active (prototype).** Demo works but gravity is hardcoded and single-body. No constraints.

#### D2. FIRE Optimizer for Trajectory Optimization

- **Merit**: FIRE (Fast Inertial Relaxation Engine) is a dynamics-based optimizer that adapts timestep and damping based on the alignment of force and velocity. Applied to trajectory optimization: control points are "particles," gradient of thrust integral is "force," FIRE finds the minimum-thrust trajectory.
- **Use cases**: Optimizing B-spline control points to minimize thrust integral. General nonlinear optimization where gradients are available.
- **Strengths**: Faster convergence than plain gradient descent. Adaptive dt (grows when force-velocity aligned, shrinks on misalignment). Velocity damping schedule. No line search needed. Well-tested in molecular mechanics. `bindArrays()` for zero-copy. Force/velocity clamping for stability.
- **Weaknesses**: No constraint handling (box constraints, equality constraints). No second-order information (Hessian) — may be slow near minimum. No comparison with L-BFGS or conjugate gradient. `optimize()` convergence loop uses simple force-norm threshold.
- **GPU parallelization**: Good — force/velocity are arrays, operations are element-wise. No GPU implementation for trajectory optimization.
- **Implemented**:
  - C++: `DynamicOpt` in `cpp/common/dynamics/DynamicOpt.h`. `move_FIRE()`, `move_GD()`, `move_MD()`, `move_LeapFrog()`. Used in `test_OptContinuousThrust.cpp` with 20 B-spline control points.
  - Status: **Active.** See also `optimization.md` for full FIRE algorithm details.

#### D3. Monte Carlo Random Search (Python)

- **Merit**: Randomly perturb trajectory parameters, accept if the objective improves. No gradient needed. Biased variant adds momentum direction to perturbation.
- **Use cases**: Global optimization of orbital transfer trajectories where the objective landscape is non-convex (multiple local minima). Orbital transfer design.
- **Strengths**: No gradient needed. Escapes local minima. Adaptive step size (grow on hit, shrink on miss ratio). Simple to implement. Can run many independent trajectories in parallel.
- **Weaknesses**: Slow convergence (especially near minimum). No exploitation of local structure. Python 2 print statements (needs porting). No covariance adaptation (CMA-ES would be better).
- **GPU parallelization**: Excellent — each trajectory is independent. Embarrassingly parallel for Monte Carlo.
- **Implemented**:
  - Python: `MC_Run()` and `MCBias_Run()` in `python/pySimE/space/exp/OrbitalTransferOpt/Random_optimization.py`.
  - C++: Not implemented for orbital mechanics. `OptRandomWalk.h` exists in `cpp/common/dynamics/` but is not used for trajectory optimization.
  - Status: **Active (Python only).** No C++ counterpart for orbital trajectory MC optimization.

---

### E. Orbital Maneuver Calculations

#### E1. Kepler Equation Solver (Fixed-Point Iteration)

- **Merit**: Solve Kepler's equation $M = E - e \sin E$ for eccentric anomaly $E$ given mean anomaly $M$ and eccentricity $e$. Fixed-point iteration: $E_{n+1} = M + e \sin E_n$. Converges for $e < 1$.
- **Use cases**: Converting mean anomaly to true anomaly for Keplerian propagation. Required for `pointAtEpoch()`.
- **Strengths**: Simple. Converges for moderate eccentricity ($e < 0.9$). No derivatives needed. Error easily monitored ($|E - e\sin E - M|$).
- **Weaknesses**: Slow for high eccentricity ($e \to 1$ — many iterations). No convergence guarantee for $e \geq 1$ (hyperbolic orbits need different formulation). `errConv = 1e-6` is hardcoded — may not be sufficient for high-precision. No Newton-Raphson fallback for difficult cases.
- **GPU parallelization**: Excellent — each body is independent. Fixed-iteration loop.
- **Implemented**:
  - C++: `true_anomaly()` in `cpp/common/Orbital/SpaceBodies.h:119-131`. Fixed-point iteration with `errConv` parameter.
  - Status: **Active.** Only Kepler equation solver in the codebase.

#### E2. Vis-Viva Equation

- **Merit**: Compute orbital speed at any radius: $v = \sqrt{GM(2/r - 1/a)}$. Fundamental relation from energy conservation in two-body problem.
- **Use cases**: Computing velocity at periapsis/apoapsis. Initial velocity setup for N-body integration. Maneuver delta-v calculations.
- **Strengths**: Exact. O(1). No iteration. Works for elliptic and hyperbolic orbits.
- **Weaknesses**: Only valid for two-body problem. No perturbation corrections.
- **GPU parallelization**: Trivial.
- **Implemented**:
  - C++: `Orbit::speed()` in `SpaceBodies.h:197`. `orbitPoint()` at `:51-60` (radial + tangential velocity components). `appliedPhysics.h` has `centralGravityForce()`.
  - Python: Used in notebooks for delta-v calculations.
  - Status: **Active.**

#### E3. Plane Change Energy

- **Merit**: Compute energy cost of changing orbital plane: $\Delta E = \frac{1}{2} |\mathbf{L}_2 - \mathbf{L}_1|^2 / r^2$ where $\mathbf{L}$ is specific angular momentum vector and $r$ is radius at maneuver point. Related to $\Delta v = 2v \sin(\theta/2)$ for pure plane change.
- **Use cases**: Evaluating maneuver cost for asteroid interception. Spacecraft mission design. `asteroidEngineering.h` uses this for mining mission planning.
- **Strengths**: Simple. Uses angular momentum invariants (conserved quantities). O(1).
- **Weaknesses**: Only pure plane change (no combined plane+altitude change). Assumes impulsive maneuver. No optimization of maneuver location.
- **GPU parallelization**: Trivial.
- **Implemented**:
  - C++: `E_orbitPlaneChange()` in `SpaceBodies.h:45-49`. `asteroidEngineering.h` for maneuver cost estimation.
  - Status: **Active.** Used for asteroid mining mission planning.

---

### F. Implementation Gap Summary

| Method | C++ | Python | GPU | Doc | Gap |
|--------|-----|--------|-----|------|-----|
| Keplerian Propagation | ✅ | ❌ | ❌ | ✅ | Python + GPU gap |
| N-Body Integration (RKF45) | ✅ | ❌ | ❌ | ✅ | Python + GPU gap |
| Dual Representation (hybrid) | ✅ (partial) | ❌ | ❌ | ❌ | Not fully used — planets integrated numerically |
| Single-Body Multi-Center | ✅ (incomplete) | ❌ | ❌ | ❌ | Thrust + planet interpolation not implemented |
| RKF45 Adaptive Integration | ✅ | ❌ | ❌ | ✅ | See `ode-integration.md` |
| Symplectic Integrator | ❌ | ❌ | ❌ | ❌ | Not implemented — needed for long-term stability |
| Hermite Spline Interpolation | ✅ | ❌ | ❌ | ✅ | Python gap |
| Cubic B-Spline Trajectory | ✅ | ❌ | ❌ | ✅ | Python gap |
| Variational Thrust Minimization | ✅ (prototype) | ❌ | ❌ | ❌ | Single-body gravity, no constraints |
| FIRE on Control Points | ✅ | ❌ | ❌ | ✅ | No comparison with other optimizers |
| Monte Carlo Trajectory Search | ❌ | ✅ | ❌ | ❌ | No C++ counterpart |
| Kepler Equation Solver | ✅ | ❌ | ❌ | ✅ | Fails for e ≥ 1 |
| Vis-Viva Equation | ✅ | ✅ (analytical) | ❌ | ✅ | — |
| Plane Change Energy | ✅ | ❌ | ❌ | ✅ | — |
| Encke Formulation | ❌ | ❌ | ❌ | ❌ | Not implemented — would improve long-term accuracy |
| Lambert's Problem | ❌ | ❌ | ❌ | ❌ | Not implemented — needed for rendezvous/targeting |
| Porkchop Plots | ❌ | ❌ | ❌ | ❌ | Not implemented — needed for transfer window analysis |
| Multi-Body Trajectory Optimization | ❌ | ❌ | ❌ | ❌ | Only single spacecraft |
| Constrained Optimization | ❌ | ❌ | ❌ | ❌ | No thrust limits, collision avoidance, or final-state constraints |

### G. Key Discrepancies & Open Opportunities

1. **Dual representation not fully utilized**: `SpaceWorld` stores `Orbit*` for each `SpaceBody` but `predictTrjs()` integrates ALL bodies (planets + ships) numerically via RKF45. The analytical `pointAtEpoch()` is only used for `pickPlanet()` (ray-picking). True hybrid propagation (analytical planets + numerical ships) would eliminate planet energy drift and reduce the ODE system size. The infrastructure exists — just need to use `orbit->pointAtEpoch(t)` for planet positions in `getDerivODE()` instead of integrating them.

2. **No symplectic integrator**: RKF45 is not symplectic — energy drifts over long orbital simulations. For game scenarios (short duration), this is acceptable. For long-term orbital mechanics (asteroid belt evolution, station-keeping), a symplectic integrator (Verlet, leapfrog, or higher-order Yoshida) would conserve energy. `DynamicOpt::move_LeapFrog()` exists but is not used for orbital propagation.

3. **Trajectory optimization is single-body and unconstrained**: `getTrajectoryVariation()` computes gravity from a single body at the origin (`GM = -1.0` hardcoded, `p` not `p - body_pos`). No multi-body gravity. No constraints (thrust limits, collision avoidance, final-state targeting). No comparison with other optimizers (CG, L-BFGS). The variational formulation is sound but the implementation is a prototype.

4. **No Lambert's problem solver**: Lambert's problem (given two positions and a time-of-flight, find the orbit) is fundamental for orbital rendezvous, targeting, and transfer trajectory design. Not implemented. Would enable: ship-to-ship intercept calculation, orbital transfer design, porkchop plots for launch windows. Standard algorithms exist (Universal Variable, p-iteration, Izzo's method).

5. **Thrust interpolation is linear**: `SpaceBody::getThrust()` uses linear interpolation between `trjThrust[i]` and `trjThrust[i+1]`. This produces C⁰ thrust (discontinuous acceleration at timesteps). For optimization, thrust should be C¹ (continuous acceleration) or C² (continuous jerk). The B-spline trajectory representation (C2) is used for optimization but not for simulation — the simulation uses linear interpolation of discrete thrust values. This mismatch means the optimized trajectory and the simulated trajectory may differ.

6. **No multi-spacecraft trajectory optimization**: `getTrajectoryVariation()` optimizes a single spacecraft's trajectory. For fleet operations (formation flying, cooperative maneuvering, multi-ship intercept), multi-body trajectory optimization is needed. The variational formulation extends naturally but the implementation doesn't support it.

7. **`findNextInd()` is O(N) linear scan**: Trajectory interpolation lookup uses linear scan (`findNextInd()` in `SpaceBodies.h:73-83`). For long trajectories (thousands of timesteps), this should be binary search (O(log N)). The `j` parameter is passed by reference to amortize across sequential calls, but random access is still O(N).

8. **No GPU implementation**: All orbital mechanics is CPU-only. For Monte Carlo trajectory optimization (thousands of independent trajectories), GPU would provide massive speedup. For N-body with many asteroids (N > 1000), GPU force evaluation would help. The RKF45 adaptive step is sequential per trajectory but parallel across trajectories.

9. **Python notebooks are standalone**: `projects/SpaceCombat/` notebooks perform analytical calculations (delta-v, maneuverability, weapon performance) but don't interface with the C++ simulation. No ctypes binding exists for orbital mechanics (unlike molecular dynamics which has `pyMolecular/`). A Python wrapper would enable: Jupyter-based trajectory design, automated testing, and cross-validation.

10. **`SpaceBodyIntegrator` is incomplete**: The single-body multi-center integrator exists but has commented-out thrust interpolation and no planet position interpolation from splines. Completing this would enable fast trajectory re-evaluation (e.g., for Monte Carlo optimization) without re-integrating the full N-body system.
