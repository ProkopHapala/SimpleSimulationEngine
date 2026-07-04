---
type: TopicalAudit
title: Noise & Procedural Generation
tags: [topic, javascript, wgsl, python, cpp, noise, perlin, simplex, fbm, terrain, value-noise, gradient-noise]
---

## Summary

Procedural noise generation for terrain synthesis. WGSL compute shaders (LandCraft WebGPU) with 5 noise flavors (sin, value, simplex, gradient, orbit) and FBM/multifractal terrain generators. Python NumPy implementation with value noise FBM and domain warping. C++ simplex noise terrain. All implementations are GPU-friendly (point-wise evaluation, no neighbor dependency).

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| WGSL | `js/LandCraft_web/noise-lib.js` | active | 5 noise functions as WGSL strings: `sinNoise` (sin-hash + bilinear), `valueNoise` (integer hash + smoothstep), `simplexNoise` (2D simplex with hash gradient), `gradNoise` (gradient noise with 8-direction hash), `orbitNoise` (Nuttall-windowed orbit hash, 4×4 stencil). `NOISE_FLAVORS` preset dict. |
| WGSL | `js/LandCraft_web/generators.js` | active | Terrain generators: `sinNoise`/`valueNoise`/`simplexNoise`/`gradNoise`/`orbitNoise` (raw debug), `fbm` (octave loop with matrix transform), `elevated` (IQ-style elevated terrain with derivative-based attenuation), `sirenian` (ridged multifractal with smoothstep). `buildAnalyticSinNoiseFn()` — analytic sin sum (5 directional components). `buildOctaveMatrixSnippet()` — rotation+shear+stretch+anisotropy per octave. |
| Python | `python/terrain_ocl/terrain.py` | active | `fbm_value_noise()` — NumPy FBM with rotation per octave, domain warping support. `value_noise_2d()` — integer hash + smoothstep interpolation. `fbm_value_noise_2()` — simpler FBM variant. OpenCL kernel for GPU terrain generation. ~680 lines. |
| C++ | `cpp/common/maps/TerrainSimplex.h` / `.cpp` | active | Simplex noise terrain generation for C++ applications |
| JS | `js/LandCraft_web/main.js` | active | LandCraft main: WebGPU setup, generator dispatch, texture output |
| Doc | `js/LandCraft_web/doc/FractalTerrain.md` | doc | Discussion: midpoint displacement vs FBM, GPU-friendly algorithms, multifractal terrain, domain warping |

## Sub-topics

### Noise Functions (WGSL)

**sinNoise**: `sinHash(p) = fract(sin(dot(p, vec2(12.9898, 78.233))) * 43758.5453)` — classic sin-hash. Bilinear interpolation with smoothstep.

**valueNoise**: Integer hash `i.x * 374761393 + i.y * 668265263`, bit-mix `(n<<13)^n`, polynomial `n*(n*n*15731+789221)+1376312589`. Returns `[-1, 1]`. Bilinear interpolation.

**simplexNoise**: 2D simplex with skew factor `K1=0.366`, unskew `K2=0.211`. Triangle cell with 3 corner contributions. Gradient hash via `fract(sin(h)*43758.5453)`. Returns `[-1, 1]` scaled by 70.

**gradNoise**: 8-direction gradient hash. Integer hash → 3 bits → select from 8 gradient directions (axis-aligned and diagonal). Dot product with fractional offset. Bilinear interpolation.

**orbitNoise**: Nuttall-windowed orbit hash. 4×4 stencil (wider support). `orbitHash` → 4 random numbers (position offset + gradient). Nuttall window for smooth falloff. Returns `[0, 1]`.

### FBM (Fractal Brownian Motion)

Standard FBM loop:
```
v = 0; a = 0.5; x = p;
for i in octaves:
    v += noise(x) * a;
    x = octMat * x;  // rotation + shear + stretch + anisotropy
    a *= 0.5;
```

`octMat` = rotation × shear × stretch (per octave):
- `rotation`: angle in degrees
- `grow`: frequency multiplier (default 1.8)
- `anisotropy`: Y-axis stretch factor
- `skew`: shear factor

### Elevated Terrain (IQ-style)

`terrainH()` from `elevated` generator:
- Uses `noised()` — noise + analytic derivative
- Accumulates derivative `d` for terrain attenuation: `a += b * n.x / (1 + dot(d,d))`
- Prevents excessive detail at high slopes
- 250×120 height scale

### Sirenian (Ridged Multifractal)

- `abs()` on derivatives for ridge effect
- `smoothstep(-0.5, 1.5, n.y)` for valley smoothing
- Alternating sign `z *= zscl; zscl *= 0.8` for varied octaves
- Radial falloff: `h -= smoothstep(0.3, 0.5, dist)`

### Domain Warping (Python)

`fbm_value_noise()` with `warp != 0`:
- `wx = value_noise_2d(rx + 13.1, ry + 7.7, seed + 101 + i)`
- `wy = value_noise_2d(rx - 9.2, ry + 21.3, seed + 211 + i)`
- `rx = rx + warp * (wx - 0.5)`
- `ry = ry + warp * (wy - 0.5)`

## Parity Status

- **WGSL `valueNoise` ↔ Python `value_noise_2d`**: Same algorithm (integer hash + smoothstep bilinear). Different hash functions. No formal parity test.
- **WGSL `simplexNoise` ↔ C++ `TerrainSimplex`**: Both 2D simplex. Implementation details may differ. No parity test.
- **WGSL `fbm` ↔ Python `fbm_value_noise`**: Same FBM structure (octave loop, amplitude halving). WGSL adds matrix transform per octave. Python adds domain warping. No formal parity test.

## Open Issues

- No formal parity tests between WGSL, Python, and C++ noise implementations
- `sinNoise` hash `sin(dot(p, vec2(12.9898, 78.233))) * 43758.5453` — known to have poor distribution at certain frequencies
- `orbitNoise` 4×4 stencil — 16× cost per evaluation vs single-cell noise
- `TerrainSimplex.h/.cpp` not fully audited
- Python `terrain.py` re-implements `value_noise_2d` locally ("external deps were breaking") — DRY violation
- No 3D noise implementation (only 2D terrain)
- No Perlin noise (only simplex, value, gradient variants)
- `generators.js` `elevated` and `sirenian` have hardcoded magic numbers (250×120, 0.00002, etc.)
