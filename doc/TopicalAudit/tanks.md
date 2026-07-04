---
type: TopicalAudit
title: Tanks
tags: [topic, cpp, military, simulation, tank, armor, penetration, wheels, turret, projectile, raytrace, 3d, rigid-body]
---

## Summary

3D tank simulation with rigid-body physics, per-polygon armor model, ray-based hit detection, tracked-wheel terrain interaction, turret rotation, and projectile ballistics. `Tank(RigidBody, Warrior3D)` has `VehicleBlock` hull and turret with `ArmorPlate` per polygon. `Wheel` spring-damper suspension. Ray casting against hull+turret mesh returns effective armor thickness (accounting for impact angle). `Shooter` world manages warriors, projectiles, terrain. Armor penetration notes document the physics model.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/apps/Tanks/Tank.h` | active | `Wheel` (spring-damper), `ArmorPlate` (thickness, mass, material), `VehicleBlock(KinematicBody, OMesh)` (armor per polygon, ray casting, armor mass calc), `Tank(RigidBody, Warrior3D)` (wheels, hull, turret, power_gear, turret rotation, terrain interaction, ray casting). |
| C++ | `cpp/apps/Tanks/Tank.cpp` | active | Stub — only includes header. Implementation in `Tank.h` (header-only). |
| C++ | `cpp/apps/Tanks/Tanks_main.cpp` | active | `Tanks_single(AppSDL2OGL_3D)`: `Shooter world`, tank creation, terrain, objects, rendering (armor color, captions, wheels, ray hit), camera, input (WASD drive, mouse turret, click fire). |
| C++ | `cpp/apps/Tanks/TankHelpers.h` | active | `armorColorScale()`, `renderArmor()`, `renderArmorCaptions()`, `drawTankWheels()`, `prepareTerrain()`. |
| C++ | `cpp/apps/Tanks/notes/ArmorPenetration.md` | doc | Armor penetration physics notes. |
| C++ | `cpp/common/CombatModels/Shooter.h` | active | `Shooter` world: warriors, projectiles, terrain, fire, update. |
| C++ | `cpp/common/CombatModels/Projectile3D.h` | active | Projectile ballistics. |
| C++ | `cpp/common/CombatModels/Warrior3D.h` | active | Base warrior class: gun_rot, trigger, move_warrior. |
| C++ | `cpp/common/CombatModels/Turret.h` | active | Turret rotation. |
| C++ | `cpp/common/physics/appliedPhysics.h` | active | `KinematicBody`, `RigidBody` — rigid body dynamics. See `projective-dynamics.md`. |
| C++ | `cpp/common/geometry/Object3D.h` | active | `Object3D`, `OMesh` — 3D object with mesh, ray casting. |
| C++ | `cpp/common/geometry/raytrace.h` | active | Ray-tracing utilities. |
| C++ | `cpp/common/terrain/Terrain25D.h` | active | 2.5D terrain heightmap with `eval(pos, gradient)`. |
| C++ | `cpp/common/geometry/Mesh.h` | active | Mesh data structure. |
| C++ | `cpp/common/physics/balistics.h` | active | Ballistic calculations. |

## Sub-topics

### VehicleBlock & Armor

`VehicleBlock(KinematicBody, OMesh)`:
- `armor` vector of `ArmorPlate` (one per polygon)
- `ray(ray0, hRay, normal, imin)` — ray-mesh intersection in local coordinates
- `getMaxArmor()` — max plate thickness
- `getArmorMass(density)` — `D = thickness*1e-3`, `S = polygonArea(i)`, `m = D*S*density`
- `loadArmor(fname)` — per-polygon thickness from file
- `fromString(str)` — load OBJ mesh + armor file

### Tank Ray Casting (Effective Armor)

`Tank::ray(ray0, hRay, ipl, block, effthick, normal)`:
1. Transform ray to tank-local coordinates via `rotMat`
2. Ray vs hull mesh and turret mesh — pick closer hit
3. Get polygon index `ipl` and armor plate
4. Compute `effthick = thick / -dot(normal, hRay)` — LOS thickness accounting for impact angle
5. Return hit distance `t`

### Wheel Suspension

`Wheel`:
- Spring-damper: `k` (spring), `mass`, `damp` (damping)
- `getSpringForce() = k*(lpos.y - h)` — compression force
- `move(dt)`: `vh += (vh*damp + fh + springForce)*dt`, `lpos.z += vh*dt`
- `lpos` — local position, `h` — terrain height, `vh` — vertical velocity, `fh` — external force

### Terrain Interaction

`Tank::interactTerrain(terrain)`:
- For each wheel: compute world position, sample terrain height
- If wheel penetrates terrain (`dh > 0`):
  - Normal force: `fnormal = dh * k` (damped if moving up)
  - Friction: `-dv * fnormal` (lateral), `-5*gv.dot(c)` (forward grip)
  - Drive: `5*(power_gear[i&1]*maxPower - gv.dot(a)*wheelLock[i&1])` — traction + rolling resistance
  - `wheelLock = 1 - clamp(|power_gear|, 0, 0.7)` — less lock when accelerating
- Apply force at wheel contact point

### Turret Control

`Tank::rotateTurret(angle)` — rotate turret around its B axis
`Tank::rotateTurretToward(dir)`:
- Compute global turret rotation
- Clamp rotation rate to 0.01 rad/step
- Yaw: align `gun_rot` toward `dir` in turret's A-B plane
- Elevation: rotate turret C/A around B by clamped angle

### Drive Control

`power_gear[2]` — left/right track power:
- W: both +1 (forward), S: both -1 (reverse)
- A: left -1, right +1 (turn left), D: left +1, right -1 (turn right)
- `wheelLock[i&1]` — even wheels = left track, odd = right track

### Rendering

- `renderArmor(block, maxThickness)` — color polygons by thickness ratio
- `renderArmorCaptions(block, sz)` — text labels: `"i:thicknessmm mass ton"`
- `drawTankWheels(tank)` — point crosses at wheel positions
- Ray hit visualization: normal vector, polygon border, effective thickness text

## Parity Status

- **`Tank` armor ↔ `LTUnitType` armor**: `Tank` uses per-polygon 3D armor with ray-traced effective thickness. `LTUnitType` uses 5 directional values (front/side/back/top/bottom). No formal parity.
- **`Tank` combat ↔ `DamagePhysicsModel.py`**: `Tank` computes effective thickness but doesn't resolve penetration. Python model does full energy-based penetration. No integration.
- **`Tank` ↔ `LandTactics`**: Different domains (3D tank vs 2D top-down). No parity.

## Open Issues

- `Tank.cpp` is a stub — all implementation in `Tank.h` (header-only)
- `Tank::ray()` has `effthick = thick/-cdot` — if `cdot ≈ 0` (grazing), effective thickness → infinity (correct physics but may cause numerical issues)
- `VehicleBlock::loadArmor()` has `printf` with wrong format: `"%i %s %lf"` but passes `mat_name` (string) for `%s` and `&thickness` for `%lf` — `mat_name` not used
- `Tank::interactTerrain()` uses `mrot.c` for forward grip but `mrot.a` for drive direction — axis convention unclear
- No projectile firing implementation visible in `Tanks_main.cpp` (trigger set but `world.fireProjectile` commented out)
- `Tank::rotateTurretToward()` clamps to 0.01 rad/step — frame-rate dependent
- `Wheel::move()` updates `lpos.z` but spring force uses `lpos.y` — possible axis inconsistency
- No armor penetration resolution — `effthick` computed but no damage model applied
- `ArmorPenetration.md` not fully audited
- `Shooter` world not fully audited — projectile/warrior interaction unknown
