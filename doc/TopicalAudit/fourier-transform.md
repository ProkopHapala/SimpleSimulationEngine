---
type: TopicalAudit
title: Fourier Transform
tags: [topic, cpp, fft, fourier, spectral, poisson]
---

## Summary

FFT (Fast Fourier Transform) implementation following Numerical Recipes in C (Cooley-Tukey radix-2). Supports 1D, 2D, and 3D transforms with forward/inverse modes. Also includes Fourier series coefficient evaluation/evaluation utilities for periodic functions. Used in electromagnetic field simulation (Poisson solver via FFT) and music visualization (power spectrum).

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/common/math/Fourier.h` | active | `FFT()` — 1D radix-2 Cooley-Tukey (Numerical Recipes). `FFT_2D()` — separable 2D via transpose + 1D. `FFT_3D()` — separable 3D via two transpose + 1D. `fourier_coef()` / `fourier_eval()` — Fourier series coefficient accumulation/evaluation from complex phase. `genSinCosArray()`, `fourier_coef_array()`, `fourier_eval_array()`. |
| C++ | `cpp/sketches_SDL/3D/test_Electromagnetic.cpp` | active | `solvePoissonFFT()` — 2D Poisson solver using `FFT_2D()` on density field. Inverse filtering in frequency domain (commented out). `solvePoisson()` — iterative Jacobi-type Poisson solver as alternative. |
| C++ | `cpp/apps/MusicVizualizer/MusicUtils.h` | active | `update()` — FFT on audio waveform, power spectrum computation, spectrum histogram dynamics. Uses `FFT()` from `Fourier.h`. |
| C++ | `cpp/sketches_SDL/music/test_visualize.cpp` | active | Audio visualization: FFT on SDL_mixer stream buffer, spectrum histogram with dynamics. |
| Java | `java/Common/` | not found | No Java FFT implementation found |

## Sub-topics

### FFT (Cooley-Tukey Radix-2)

From Numerical Recipes in C (pp. 507-508):
- Input: `data[2*i] = Re`, `data[2*i+1] = Im`, `nn = 2^m`
- Bit-reversal permutation (Danielson-Lanczos)
- Butterfly operations with trigonometric recurrence (`wpr`, `wpi`)
- Forward: `isign = +1`, Inverse: `isign = -1` (with `1/nn` normalization)
- Complexity: O(N log N)

### 2D/3D FFT

Separable approach: apply 1D FFT along each axis sequentially, with explicit transpose between axes.
- `FFT_2D(nnx, nny, data, isign)`: FFT along x → transpose → FFT along y → transpose back
- `FFT_3D(nx, ny, nz, data, isign)`: FFT along x → transpose(x↔y) → FFT along y → transpose(y↔z) → FFT along z → transpose back

### Fourier Series (Non-FFT)

For periodic functions `f(φ)` given complex phase `z = exp(iφ)`:
- `fourier_coef(ca, sa, val, n, coefs)`: accumulate coefficients using recurrence `cos(nφ) + i·sin(nφ)` from `cos(φ)`, `sin(φ)`
- `fourier_eval(ca, sa, n, coefs)`: evaluate series sum
- Used for angular-dependent quantities (e.g. aerodynamic coefficients, molecular potentials)

### Poisson Solver via FFT

In `test_Electromagnetic.cpp`:
- Forward FFT of density → multiply by Green's function `1/k²` in frequency domain → inverse FFT
- Currently incomplete: Green's function multiplication is commented out
- Alternative iterative solver (`solvePoisson()`) uses Jacobi-type relaxation

## Parity Status

- **`FFT()` ↔ Numerical Recipes reference**: Direct port. No formal parity test but algorithm is well-established.
- **`FFT_2D()` ↔ `FFT_3D()`**: Both use same separable approach. No formal parity test.
- **Poisson FFT solver ↔ iterative Poisson solver**: Both exist in `test_Electromagnetic.cpp` but no comparison run.

## Open Issues

- FFT requires power-of-2 sizes (`nn = 2^m`) — no mixed-radix or Bluestein support
- 2D/3D FFT uses explicit transpose with `new`/`delete` — not optimized for cache, no in-place variant
- Poisson FFT solver incomplete — Green's function multiplication commented out
- No GPU/OpenCL FFT implementation
- No Python FFT wrapper — Python uses `numpy.fft` directly
- `fourier_coef()` / `fourier_eval()` not used outside of potential test sketches
- No windowing functions (Hann, Hamming) for spectral analysis
- Music visualizer code may be stale (SDL_mixer dependency)
