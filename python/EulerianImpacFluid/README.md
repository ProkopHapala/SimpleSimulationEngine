2D Eulerian multi-material compressible flow solver with level-set interface tracking, implemented in PyOpenCL. Simulates high-velocity impact between two stiffened-gas materials (e.g. uranium projectile into liquid hydrogen) using a 5-equation diffuse-interface model.

- **EulerianImpacFluid.cl** — OpenCL kernels: `update_fluid` (Lax-Friedrichs flux + non-conservative phi advection, 8x8 tiled with local-memory halo) and `redistance_phi` (Godunov-type HJ signed-distance reinitialization). Mixture stiffened-gas EOS with volume fractions from level-set phi.
- **EulerianImpacFluid.py** — PyOpenCL wrapper class: buffer allocation, ping-pong state swapping, kernel dispatch, host-device transfer. NVIDIA GPU preferred via `_pick_nvidia_gpu()`.
- **test_eulerian_fluid.py** — CLI demo driver with argparse: sets up uranium-in-hydrogen impact scenario, animates diagnostic fields (density, pressure, sound speed, velocity, phi) with matplotlib FuncAnimation. Run: `python test_eulerian_fluid.py --fields density,pressure,sound,speed --p0 1e9 -t 1e-8 -n 10 -f 50`
- **EulerianImpacFluid.md** — Detailed derivation and design notes for the 5-equation model, EOS, flux scheme, and level-set tracking (124K, comprehensive reference).
- **EulerianImpacFluid_debug_presure.md** — Retrospective on pressure debugging: 6 root causes identified and fixed (volume-fraction vs mass-fraction mix-up, conservative phi advection, init inconsistency, missing KE, GPU diagnostics, float32 precision with large PI_2).
- **tmp.md** — Full code dump snapshot of all 3 source files (AI-generated reference copy).
