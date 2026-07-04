# SwordPlay

Soft-body sword fighting simulation — martial combat system inspired by Mortal Kombat and Blade Symphony. Uses soft-body physics for realistic sword-on-sword collision and character movement. The sword arm is modeled as a soft polyline chain.

## Physics

- Soft polyline dynamics for sword arm (`SoftPolyLine2D.h`)
- Muscle-skeleton model for character animation
- Combat move definitions (attacks, blocks, combos)

## Files

- **SwordPlay_main.cpp** — main application: 2D combat, sword physics, character control
- **SwordArm.h** — sword arm: soft polyline chain with sword attached
- **CombatMoves.h** — combat move definitions: attacks, blocks, combos
- **MusculeSkelet.h** — muscle-skeleton model for character body
- **CMakeLists.txt** — build target: `SwordPlay_main`

## Related Documentation

- [docs/FormationTactics/PolearmsPhysicsModel.md](../../docs/FormationTactics/PolearmsPhysicsModel.md) — polearms physics model: reach, swing dynamics, hit detection
- [docs/RigidBodyDynamics/RigidBodyProperly.md](../../docs/RigidBodyDynamics/RigidBodyProperly.md) — rigid body dynamics formulation
