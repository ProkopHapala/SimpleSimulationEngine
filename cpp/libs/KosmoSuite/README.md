# KosmoSuite

Physics simulation suite for space and nuclear phenomena. Provides shared library with C API for Python/ctypes access. Covers N-body gravity, orbital mechanics, spacecraft launch trajectories, electromagnetic fields, nuclear fission pulses, neutron transport, and spherical shock waves.

## Physics Modules

- **N-body gravity** — RKF45 integration of gravitational N-body system
- **Orbital mechanics** — Keplerian orbital elements from position/velocity, eccentric/true anomaly conversion
- **Spacecraft launch** — multi-stage launch trajectory optimization with control points
- **Ship acceleration** — constant-thrust trajectory with gravity
- **Electromagnetic fields** — Biot-Savart law for coil magnetic fields, Lorentz force on charged particles
- **Nuclear fission pulse** — fission yield model (cross-sections, generation rate, energy per fission)
- **Neutron transport** — Monte Carlo neutron tracing through nuclear materials (fission/capture cross-sections)
- **Spherical shock waves** — compressible material EOS, adiabatic shock wave propagation

## Files

- **KosmoSuite.cpp** — main library: gravity force, Lorentz force, Biot-Savart, exported C API functions
- **cpp/Nbody.h** — N-body gravitational simulation with RKF45 integrator
- **cpp/OrbitalUtils.h** — orbital element conversion (Keplerian elements from state vectors)
- **cpp/SpaceLaunchODE.h** — multi-stage launch trajectory ODE with thrust control points
- **cpp/ShipAccel.h** — constant-acceleration trajectory along velocity direction
- **cpp/trajectory.h** — trajectory interpolation and binary search utilities
- **cpp/elmag.h** — electromagnetic field sampling, Biot-Savart for coils, Lorentz force integration
- **cpp/fissionPulse.h** — nuclear fission pulse model: cross-sections, yield, energy release
- **cpp/neutronTrace.h** — neutron transport: nuclear material properties, mean free path, fission/capture ratios
- **cpp/neutronTraceHomogeneous.h** — simplified neutron transport in homogeneous medium
- **cpp/shockWavesSpherical.h** — spherical shock wave propagation in compressible materials
- **cpp/sphericalShockWaves.h** — spherical shock wave utilities (stub)
- **cpp/CMakeLists.txt** — subdirectory build config
- **CMakeLists.txt** — build target: `KosmoSuite` (SHARED)

## Related Documentation

- [docs/SpaceTactics/appliedPhysics.md](../../docs/SpaceTactics/appliedPhysics.md) — applied physics models for space combat
- [docs/SpaceTactics/SpaceBodies.md](../../docs/SpaceTactics/SpaceBodies.md) — celestial body definitions
- [docs/BiotSavart/](../../docs/BiotSavart/) — potential flow, vortex methods, MHD plasma
