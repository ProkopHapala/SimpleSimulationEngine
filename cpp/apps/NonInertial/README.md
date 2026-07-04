# NonInertial

Real-time combat simulation (like Liero or Worms) set in a rotating non-inertial frame of reference (e.g. a spinning space station). Warriors fight with projectiles while experiencing Coriolis and centrifugal forces.

## Physics

- Non-inertial reference frame with angular velocity `omega`
- Coriolis force and centrifugal force on projectiles and warriors
- Collision detection between warriors and environment
- Projectile ballistics with air drag and restitution

## Files

- **NonInert_main.cpp** — main application: 2D combat in rotating frame, warrior control, firing
- **NonInertWorld.cpp / .h** — world model: rotating frame physics, warriors, projectiles, collisions
- **Warrior2D.cpp / .h** — warrior: position, health, movement, aiming
- **Projectile2D.cpp / .h** — projectile: ballistics in non-inertial frame
- **CMakeLists.txt** — build target: `NonInertial_main`
