# PotentialFlow.h — Analytical Building Blocks for Potential and Vortex Flows

This document is a didactic atlas of the analytical functions implemented in `cpp/common/math/PotentialFlow.h`. It summarizes what each function computes, the key formulas they implement, and relationships between them (derivatives/integrals). It also references two example applications that use these functions:

- Electromagnetics: `cpp/sketches_SDL/3D/test_Electromagnetic.cpp`
  - Simulates motion of charged particles in electrostatic fields and approximates magnetic fields from coils using polygonal segments (Biot–Savart law).
- Aerodynamics (Vortex Lattice Method): `cpp/sketches_SDL/3D/test_VortexLattice.cpp`
  - Computes induced velocity from a horseshoe vortex system to study lift and induced drag.

All formulas below use dimensionless constants unless noted. The constant `VortexQuantum = 1/(4π)` scales Biot–Savart-like fields.


## Function Index (All functions)

- `pointSource(R)`: Coulomb-like point source; returns $R/|R|^3$ (up to constants), i.e., the gradient of `1/r`.
- `sourceDipol(R, coefs)`: Generalized monopole+dipole field using `coefs.e` (monopole) and `coefs.f` (dipole p); scales like `1/r^3` with correct angular dependence.

- `dBiotSawart(R, dI)`: Differential contribution of a filament/current element; `VortexQuantum · (dI × R)/|R|^3`.
- `ILineFinite(R, hL, l)`: Induced field of a finite straight segment of length `l` along unit `hL`; uses closed-form `(cosθ2 − cosθ1)/|hL × R|^2` factor.
- `ILineSemiInf(R, hL)`: Semi-infinite straight filament extending along `+hL`; closed form proportional to `(1 − cosθ1)/|hL × R|^2`.
- `ILineSemiInfDecay(R, hL, w2)`: Semi-infinite filament with Lorentzian core smoothing; denominator `|hL × R|^2 + (hL·R)^2 w2` regularizes the singular core.
- `horseshoe(B, R, p0, p1, hDir, strength)`: Canonical VLM horseshoe vortex (bound segment `p0→p1` with trailing legs ±`hDir`); adds induced field to `B` with circulation `strength`.
- `horseshoeDecay(B, R, p0, p1, hDir, strength, w2)`: Horseshoe vortex with smoothed trailing legs controlled by `w2` to stabilize near-field induction; adds into `B`.
- `ISemiInfSheet(R, a, b, l)`: Field of a semi-infinite sheet strip formed by integrating a semi-infinite line along `b` over length `l`; direction `n = normalize(a × b)` with logarithmic magnitude.

- `sourceXY(x, y, Fx, Fy)`: Endpoint primitive for constant line density `ρ=1`; writes `(Fx, Fy)` at a single endpoint used to build finite-segment differences.
- `sourceLineF_pow0/1/2/3(x, y, Fx, Fy)`: Endpoint primitives for `ρ(x) = 1, x, x^2, x^3`; fill `(Fx, Fy)` at one endpoint for assembling finite segments.
- `sourceLineF_linear(x, y, Fx, Fy, C0, C1)`: Linear combination `C0 + C1 x` built from `pow0` and `pow1` primitives; writes `(Fx, Fy)` at one endpoint.
- `sourceLineF_poly3(x, y, Fx, Fy, C0, C1, C2, C3)`: Polynomial density up to cubic by combining `pow0..pow3`; writes `(Fx, Fy)` at one endpoint.
- `sourceLine_const(R, fw, L, funcXY)`: 3D assembler for a straight finite segment along `fw` of length `L`; evaluates endpoint primitives `funcXY` at `x` and `x+L` and returns their difference projected along `fw` and the transverse direction.

- `CurrentLoop(dp, R)`: Placeholder for off-axis circular current loop field (requires elliptic integrals); currently returns zero and awaits a tailored implementation.


## 1) Point and Dipole Sources

Let `R ∈ ℝ^3`, `r = |R|`, `r^2 = R·R`.

- Field of a point source (Coulomb-like):
  - Implementation: `pointSource(R) = R / r^3`.
  - Inverse-square radial field, gradient of `1/r` (up to constants).

- Generalized monopole + dipole field:
  - Parameters: `coefs.e` (scalar monopole), `coefs.f` (Vec3 dipole moment `p`).
  - Implementation:
    - `ir2 = 1/r^2`, `ir3 = 1/r^3`
    - `F(R) = R * ((p·R)*ir2 + e) * ir3 + p * ir3`
  - Relation to ideal electric dipole (no monopole, `e=0`): proportional to `(3 ⟨p|r̂⟩ r̂ − p)/r^3`.


## 2) Biot–Savart Primitives

Conventions:
- `R`: vector from the element to the evaluation point.
- `hL`: unit vector along the line/filament direction.
- `l`: length of the finite element.
- `VortexQuantum = 1/(4π)` absorbed into outputs.

- Differential element (Biot–Savart):
  - `dB = VortexQuantum * (dI × R) / r^3`.

- Finite straight segment (from R to R + l hL):
  - Let `a = |hL × R|` (minimum distance), `c1 = cosθ1 = (hL·R)/|R|`, `c2 = cosθ2 = (hL·(R + l hL))/|R + l hL|`.
  - Field:
    - `B = VortexQuantum * (hL × R)/a^2 * (c2 − c1)`
  - Implementation directly computes `(hL × R) * ((c2 − c1) / |hL × R|^2)`.

- Semi-infinite straight segment (starting at the observation’s projection and going along +hL):
  - With `c1 = cosθ1 = (hL·R)/|R|`:
  - `B = VortexQuantum * (hL × R)/|hL × R|^2 * (1 − c1)`.

- Semi-infinite with Lorentzian core smoothing:
  - With `cL = hL·R`, replace denominator by `|hL × R|^2 + (cL)^2 w2`:
  - `B = VortexQuantum * (hL × R) * (1 − c1) / ( |hL × R|^2 + (cL)^2 w2 )`.
  - Parameter `w2` controls core radius/smoothing; `w2 → 0` recovers the singular filament.

- Horseshoe vortex (Vortex Lattice Method canonical element):
  - Bound segment from `p0` to `p1` with strength Γ, trailing semi-infinite legs from `p0`, `p1` in direction `hDir` with strengths `−Γ`, `+Γ`.
  - Implementation sums: `ILineSemiInf(R−p0, hDir) * (−Γ) + ILineSemiInf(R−p1, hDir) * (+Γ) + ILineFinite(R−p0, (p0−p1)/|p0−p1|, |p0−p1|) * Γ`.
  - `horseshoeDecay` uses the smoothed semi-infinite legs.

- Semi-infinite sheet (analytical integral across one dimension):
  - Geometry: orthonormal `a`, `b` span the sheet, evaluation at `R`.
  - Let `x = a·R`, `y = b·R`, and integrate along `b` over `[0, l]`.
  - Field direction is normal `n = normalize(a × b)`.
  - Closed form (logarithmic):
    - `B = VortexQuantum * n * ln( (√(x^2 + (y + l)^2) + x) / (√(x^2 + y^2) + x) )`


## 3) Line-Source Integrals With Polynomial Density

Interpretation:
- We integrate a 3D Coulomb-like kernel along a straight line parallel to x at fixed offset y: integrand ∝ ρ(x) (x, y)/r^3 with `r = √(x^2 + y^2)`.
- Results are expressed as endpoint primitives `F(x, y)`, which are used in differences to form finite segments.
- Useful identity: `asinh(x/y) = ln(x/y + √(1 + (x/y)^2)) = ln((x + r)/y)`.

Closed forms implemented (Fx, Fy = primitive values at a single endpoint):
- ρ(x) = 1:
  - `Fx = − 1/r`
  - `Fy = (x/y) * (1/r)`
- ρ(x) = x:
  - `Fx = − x/r + asinh(x/y)`
  - `Fy = − y/r`
- ρ(x) = x^2:
  - `Fx = (x^2 + 2 y^2)/r`
  - `Fy = y * (− x/r + asinh(x/y))`
- ρ(x) = x^3:
  - `Fx = ( x (x^2 + 3 y^2) − 3 y^2 r asinh(x/y) ) / (2 r)`
  - `Fy = y (x^2 + 2 y^2) / r`

Linear and cubic combinations simply combine these primitives with coefficients.

Relationships (derivatives/integrals):
- With respect to x, higher-order ρ(x) primitives are generated by integrating products like `x^n / r^3`.
- Notably, `d/dx asinh(x/y) = 1/√(x^2 + y^2) = 1/r`, which explains appearances of `asinh(x/y)` in the antiderivatives.
- The sequence {ρ = 1, x, x^2, x^3} shows a ladder where each primitive contributes to the next via algebraic multiplication by x plus a corrective term involving `asinh(x/y)`.


## 4) 3D Assembly Helper for Finite Line Segments

Purpose:
- Builds the 3D field of a straight finite segment of length `L` along unit direction `fw` by evaluating the 2D endpoint primitives `funcXY(x, y)` at the two segment endpoints in local coordinates and differencing them.

Procedure:
- Decompose `R` into `x = fw·R` and `up = R − fw x` with `y = |up|`.
- Evaluate `(Fx0, Fy0) = funcXY(x, y)` and `(Fx1, Fy1) = funcXY(x+L, y)`.
- Return `fw*(Fx1 − Fx0) + up*(Fy1 − Fy0)`.

This is the core pattern used to turn the endpoint primitives from Section 3 into finite-segment fields in 3D.


## 5) Usage in the Example Applications

- `test_Electromagnetic.cpp` (magnetic bottle / plasma nozzle)
  - `coilField(pos, h, R, n)` builds a polygonal approximation of a circular current loop and sums `ILineFinite` contributions along its edges to compute `B(pos)`.
  - `ILineSemiInf` and the 2D grid are also used for visualization and acceleration.
  - Electrostatic wires use a simple inverse-square-like field superposition `E += q (r)/(r^2 + Rwire^2)` (smoothing), separate from `PotentialFlow.h` primitives.

- `test_VortexLattice.cpp` (Vortex Lattice Method)
  - The free-stream velocity `vInf` is combined with the induced velocity from `horseshoeDecay` for each panel/segment between control points `CPs[i] → CPs[i+1]`.
  - `horseshoeDecay` uses smoothed trailing legs (Lorentzian core) to avoid singularities in the induced velocity.
  - Additional numerical integration helpers verify endpoint-primitive formulas against discretized sums.


## 6) Current Loop (Coil) — Elliptic Integrals Note

Background:
- The magnetic field of a circular current loop at an off-axis point can be expressed in closed form using complete elliptic integrals (e.g., with parameter `k^2 = 4α/Q` where `α, β, Q` are functions of geometry as in the code stub).
- The repository currently includes `cpp/common/math/elliptic_integral.c` (with a matching header) providing Carlson forms and elliptic integrals; this is a general and somewhat heavy implementation.

Project preference:
- For coil fields, a specialized numerical approximation of the required elliptic integrals (or direct quadrature with adaptive sampling) may be preferable for simplicity/performance.
- Action item: Keep `CurrentLoop(...)` as a stub and, when time permits, implement a tailored numerical approximation of the two needed complete elliptic integrals to obtain `B_r` and `B_z` efficiently and robustly.

References:
- Elliptic integrals overview and numerics: see references already listed in `PotentialFlow.h` and the general notes in `elliptic_integral.c`.


## Summary of Key Identities and Tips

- `asinh(x/y) = ln((x + √(x^2 + y^2))/y)` is central to the 1D primitives.
- For Biot–Savart segments, the vector direction is always `(hL × R)`, scaled by an angular factor `(cosθ2 − cosθ1)` and inverse distance-square-like denominator `|hL × R|^2` (or its smoothed variant).
- Finite segments are built as endpoint differences of analytical primitives; this pattern underlies both Coulomb-line and vortex-filament assemblies.


 


This documentation is intended to be pasted back as inline comments near each function for future maintainability, once reviewed and adjusted for any unit conventions and sign conventions in your application.
