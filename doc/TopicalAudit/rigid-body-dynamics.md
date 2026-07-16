---
type: TopicalAudit
title: Rigid Body Dynamics
tags: [topic, cross-language, physics, rigid-body]
---

## Summary

Rigid body dynamics in 2D (`Body2D.h`) and 3D (`Body.h`), used in molecular dynamics (`MolecularWorld.h`), flight simulation (`AeroCraft.h`), tank simulation, and sailing. Pairwise force kernels in `Forces.h` (LJ, Morse, Coulomb, spring, angular). Constraint solving for multi-body assemblies in `KinematicSolver.h` (Newton-Raphson with slider joints). Configuration space exploration in `RBodyConfDyn.h`. Python proof/derivation in `docs/RigidBodyDynamics/RigidBodyProperly.md` with companion `rigid_body_proof.py`. Applications span multiple game/sim apps.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/common/dynamics/Body2D.h` | active | 2D rigid body: position, velocity, rotation, angular velocity, mass, inertia |
| C++ | `cpp/common/dynamics/Body.h` | active | 3D rigid body: quaternions, full 6-DOF |
| C++ | `cpp/common/dynamics/Forces.h` | active | Pairwise force kernels: LJ, Morse, Coulomb, spring, repulsion, angular terms. Smoothstep cutoffs (r^2-based). Used by molecular dynamics and soft body sims |
| C++ | `cpp/common/dynamics/KinematicSolver.h` | active | Newton-Raphson (Levenberg-Marquardt) constraint solver for rigid body assemblies with anchor and slider joints. B-spline slider paths |
| C++ | `cpp/common/dynamics/RBodyConfDyn.h` | active | 6-DOF configuration space exploration via MD with repulsion from visited states. For docking/contact pose sampling |
| C++ | `cpp/common/dynamics/AeroCraft.h` | active | Flight dynamics extension: aerodynamic surfaces, control inputs |
| C++ | `cpp/common/dynamics/MolecularWorld.h` | active | Rigid-body molecular dynamics (coarse-grained) |
| C++ | `cpp/apps/Tanks/` | active | Vehicle on terrain |
| C++ | `cpp/apps/SailWar/` | active | 2D sail-ship simulator |
| C++ | `cpp/apps/AeroCombat/` | active | 3D subsonic aircraft combat |
| Python | `docs/RigidBodyDynamics/rigid_body_proof.py` | active | Numerical proof/verification of rigid body equations |
| Doc | `docs/RigidBodyDynamics/RigidBodyProperly.md` | active | LLM chat transcript: derivation of correct rigid body dynamics (58KB) |

## Parity Status

- **C++ `Body.h` ↔ Python `rigid_body_proof.py`**: Python serves as numerical verification of C++ 3D rigid body equations. No automated parity test found.
- 2D and 3D implementations are independent (not a port of each other).

## Open Issues

- No automated regression test linking Python proof to C++ implementation
- `Body.h` and `Body2D.h` share concepts but no shared code — possible duplication of physics logic
- `AeroCraft.h` aerodynamic model documented in `aerodynamics-hydrodynamics.md`
- Motion integration notes: `cpp/common/engine/Notes/MotionIntegration.md`
