---
type: TopicalAudit
title: Fast Math Approximations
tags: [topic, cpp, math, approximation, chebyshev, sin-cos, atan2, erf]
---

## Summary

Speed-optimized approximations of transcendental functions using Chebyshev polynomial fits and table lookups. `gonioApprox.h` covers sin/cos (Chebyshev polynomial in x²), atan2 (multiple variants: polynomial, NVIDIA, quadrant-switched), acos (table-based Hermite interpolation, NVIDIA polynomial). `fastmath.h` provides erf (4-term and 6-term polynomial), `erfx_e6` (erf(kx)/x with derivative), utility functions (clip, clamp, fastFloor, dangle). `functions.h` provides test functions for optimization (Rosenbrock, Mandelbrot, Lorenz, harmonic, spiral).

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/common/math/gonioApprox.h` | active | Sin/cos: Chebyshev polynomial in x² with order 6/8/10. `cos_sin()` — range reduction to [-π/4, π/4], polynomial eval, quadrant rotation. Atan2: 3 variants (`atan2_a1`, `atan2_a2`, `atan2_a3`) + NVIDIA version. Tan: polynomial in x² order 6-12. Acos: table-based (66-entry Hermite) + NVIDIA polynomial. |
| C++ | `cpp/common/math/fastmath.h` | active | `erf_4()`, `erf_6()` — polynomial approximations. `erfx_e6()` — erf(k·x)/x with derivative, ~1e-6 accuracy, used for Gaussian-blurred Coulomb potentials. Utility: `sq()`, `clip()`, `clamp()`, `clamp_abs()`, `fastFloor()`, `fastFract()`, `fastModf()`, `dangle()`, `x2grid()`. Constants: GOLDEN_RATIO, DEG2RAD, RAD2DEG, M_TWO_PI. |
| C++ | `cpp/common/math/functions.h` | active | Test functions for optimization benchmarks: `sigmoideAbs()`, `sigmoideSqrt()`, `lorenz()`, `x1period()`, `x2period()`, `harmonic()`, `rosenbrok()`, `sinValey()`, `spiral()`, `mandelbort()`, `particles()`, `cross_valey()`, `gridValy()`. All with analytic derivatives. |
| C++ | `C/math/fastmath.h` | active | C-language version (shared with C++ via include) |

## Sub-topics

### Sin/Cos Chebyshev Approximation

Range reduction: x → [-π/4, π/4] via integer quotient + quadrant bits. Then evaluate polynomial in x²:
- `cos_xx_8(xx) = 1 + xx·(-0.308425... + xx·(0.015854... + xx·(-0.000326... + xx·3.54e-6)))`
- `sin_xx_6(xx) = 0.785398... + xx·(-0.080745... + xx·(0.002490... + xx·(-3.596e-5)))`
- Result: `sin = x·sin_xx(xx)`, `cos = cos_xx(xx)`
- Quadrant correction via bit manipulation (`ix&2`, `ix&4`)

Configurable order version: `cos_sin(x, ca, sa, order, nsplit)` — supports arbitrary polynomial order and range splitting.

### Atan2 Variants

1. `atan2_a1()` — Volkan Salma method: reduce to octant, polynomial in a²
2. `atan2_a2()` — quadrant-switched: 8 cases based on sign(x), sign(y), |x|vs|y|
3. `atan2_a3()` — fully branchless 8-case switch with `atan_poly()`
4. `atan2_nvidia()` — NVIDIA Cg library polynomial (float precision)
5. Templated `atan2_t<poly>()` — generic with polynomial strategy

### Erf Approximation

- `erf_4()`: 4-term polynomial, `(1 + x·(0.278 + x·(0.230 + x·(0.001 + x·0.078))))⁴`
- `erf_6()`: 6-term polynomial, higher accuracy
- `erfx_e6(x, k, dy)`: Special form `erf(k·x)/x` with derivative — used for Gaussian-blurred Coulomb potentials in molecular mechanics. Even/odd polynomial decomposition. Asymptotic branch for x>4.5.

### Acos Table Lookup

66-entry table with Hermite interpolation: `TABLE_ACOS[66]` stores (value, derivative) pairs at uniform spacing dx=0.059375. Cubic Hermite interpolation for |x|<20·dx, sqrt fallback for larger |x|.

## Parity Status

- **`gonioApprox.h` sin/cos ↔ std::sin/cos**: No formal parity test. Chebyshev coefficients derived from known fits (referenced: lolengine.net blog).
- **`atan2_a1` ↔ `atan2_a2` ↔ `atan2_a3` ↔ `atan2_nvidia`**: Multiple variants exist but no comparative accuracy benchmark in repo.
- **`erf_4` ↔ `erf_6` ↔ `erfx_e6`**: Different accuracy levels, no automated test.

## Open Issues

- No automated accuracy benchmark for any approximation
- `cos_sin()` uses mixed orders (cos_xx_8 for cos, sin_xx_6 for sin) — inconsistent
- `acos_table()` has a bug: `i<<1;` on line 341 has no effect (missing assignment)
- `functions.h` test functions use global variables (`VALEY_TIGHTNESS`, `FREQ`, `N_MANDELBROT`) — not thread-safe
- `erfx_e6` comment notes need for further fitting of `dy/x` for Gauss:Coulomb()
- No float32-specific fast path for sin/cos (all double)
- `C/math/fastmath.h` — C version not audited in detail, may diverge from C++ version
