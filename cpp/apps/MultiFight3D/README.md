# MultiFight3D

3D multiplayer combat simulation. Multiple fighters (warriors) engage in melee and projectile combat in a 3D arena. Uses rigid body dynamics for movement and collision detection for hit detection.

## Files

- **MultiFight3D_main.cpp** — main application: 3D arena, fighter control, combat
- **MultiFight3DWorld.cpp / .h** — world model: fighters, projectiles, collision, combat resolution
- **CMakeLists.txt** — build target: `MultiFight3D_main`

## Related Documentation

- [docs/RigidBodyDynamics/RigidBodyProperly.md](../../docs/RigidBodyDynamics/RigidBodyProperly.md) — rigid body dynamics formulation and integration
- [docs/RigidBodyDynamics/rigid_body_proof.py](../../docs/RigidBodyDynamics/rigid_body_proof.py) — Python proof/verification of rigid body equations
