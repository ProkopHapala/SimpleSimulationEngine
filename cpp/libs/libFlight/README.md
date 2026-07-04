# libFlight

Shared library wrapping the flight simulation world for Python/ctypes access. Provides C API for aircraft aerodynamics: `FlightWorld` with `AeroCraft`, control surfaces (pitch/yaw/roll/throttle), gravity, wind, and projectile combat.

## Physics

- Aerodynamic forces from `AeroSurf.h` wing panels and `AeroCraft.h` aircraft model
- Rigid body dynamics via `Body.h`
- Control state with rate-limited actuator model (pitch, yaw, roll, throttle, trigger)
- Projectile combat via `Projectile3D.h` and `ShooterCommon.h`

## Files

- **libFlight.cpp** — `Flight` library: `FlightWorld` singleton, file loading, array pruning utility, exported C API
- **libFlight.h** — `FlightWorld` class: aircraft, control states, gravity, wind, projectiles; `ControlState` actuator model
- **CMakeLists.txt** — build target: `Flight` (SHARED, links `Body`, `AeroSurf`, `AeroCraft` objects)

## Related Documentation

- [docs/AeroCombat/AeroCombat.md](../../docs/AeroCombat/AeroCombat.md) — aerodynamic combat model overview
- [docs/AeroCombat/AeroSurf.md](../../docs/AeroCombat/AeroSurf.md) — aerodynamic surface force calculation
