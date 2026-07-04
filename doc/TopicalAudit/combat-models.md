---
type: TopicalAudit
title: Combat Models (Melee, Damage, Aero, Land)
tags: [topic, cross-language, combat, damage, games, physics]
---

## Summary

Multiple combat physics models across different game applications: melee weapon physics (polearms, damage penetration), aerial combat aerodynamics, land tactical simulation, formation battle simulation. Each model is grounded in physics (kinetic energy, momentum, aerodynamics) but tuned for gameplay. C++ implementations with Python design/derivation docs (LLM chats).

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/common/CombatModels/AirCombatModel.h` | active | Aerial combat: drag, lift, turn rate, climb rate, energy management, slip-away tactics |
| C++ | `cpp/apps/AeroCombat/` | active | 3D subsonic aircraft combat simulator |
| C++ | `cpp/apps/FormationTactics/` | active | 2D historical battle simulator (TotalWar-like, ~30000 soldiers) |
| C++ | `cpp/apps/LandTactics/` | active | WWII and modern era tactical simulator |
| C++ | `cpp/apps/Tanks/` | active | 3D tank simulator |
| C++ | `cpp/apps/SailWar/` | active | 2D sail-ship simulator |
| Python | `docs/FormationTactics/DamagePhysicsModel.py` | active | Damage model implementation: piercing, cutting, crushing |
| Python | `docs/FormationTactics/PolearmsPhysicsModel.py` | active | Polearm physics: thrust, swing, moment of inertia, biomechanics |
| Doc | `docs/FormationTactics/DamagePhysicsModel.md` | active | LLM chat: armor penetration model. KE vs momentum, inelastic collision `m/(m+M)`, three damage types (piercing/cutting/crushing) |
| Doc | `docs/FormationTactics/PolearmsPhysicsModel.md` | active | LLM chat: polearm combat model. Thrust (linear), swing (angular), moment of inertia, swing time, biomechanics (40KB) |
| Doc | `docs/LandTactics/AirCombatModel.md` | active | Air combat model documentation: aerodynamics, turn performance, tactical maneuvering |
| Doc | `docs/LandTactics/LandTactics.md` | active | Land tactical sim implementation report |
| Doc | `docs/LandTactics/LTFaction.md` | active | Faction model |
| Doc | `docs/LandTactics/LTUnit.md` | active | Unit model |
| Doc | `docs/LandTactics/LTUnitType.md` | active | Unit type definitions |
| Doc | `docs/LandTactics/LTShelter.md` | active | Shelter/building model |
| Doc | `docs/LandTactics/LTWorld.md` | active | World model |
| Doc | `docs/LandTactics/LTcommon.md` | active | Common definitions |
| Doc | `docs/AeroCombat/AeroCombat.md` | active | AeroCombat app documentation |
| Doc | `docs/AeroCombat/AeroSurf.md` | active | Aerodynamic surface model |
| Doc | `cpp/apps/Tanks/notes/ArmorPenetration.md` | active | Tank armor penetration model |
| Doc | `cpp/apps/AeroCombat/Notes/Controls.md` | active | AeroCombat controls |
| Doc | `cpp/apps/AeroCombat/Notes/PanelInteraction.md` | active | Panel interaction |
| Doc | `cpp/common/engine/Notes/ShotHitDetection.md` | active | Shot hit detection algorithm |

## Sub-topics

### Damage Physics Model

Three damage mechanisms derived from physics:
- **Piercing**: KE-driven armor penetration. `p*S*L = E` (pressure × area × depth = energy). Pointy weapons: area increases with depth.
- **Cutting**: Edge length matters. Wound depth + width. Material behind armor is much weaker.
- **Crushing (Blunt)**: Momentum-driven. Ignores armor. Energy transfer fraction: `m_head / (m_head + M_armor)` from inelastic collision.

### Polearm Physics Model

- **Thrust**: `v = sqrt(2*F*d/m_total)`, `t = sqrt(2*d*m_total/F)`. Acceleration distance: 1m one-handed, 1.5m two-handed.
- **Swing**: Moment of inertia `I = (1/3)*m_shaft*L² + m_head*L²`. Angular velocity from torque × time.
- **Swing time**: Determined by soldier strength, weapon length, and mass distribution.

### Aerial Combat Model

- Drag: `Fd = Cd * (ρ/2) * v²`
- Lift: `Fl = Cl * (ρ/2) * v²`
- Turn rate, climb rate, energy management (KE vs PE trade-off)
- Slip-away tactics: fast vs slow aircraft analysis
- Stall and energy limitations

## Parity Status

- **Python `DamagePhysicsModel.py` ↔ C++ `FormationTactics`**: Python is derivation/prototyping, C++ is game integration. No automated parity test.
- **Python `PolearmsPhysicsModel.py` ↔ C++ `FormationTactics`**: Same — Python derivation, C++ integration.
- **C++ `AirCombatModel.h` ↔ `AeroCombat` app**: `AirCombatModel.h` is the shared model, `AeroCombat` app uses it.
- **AeroCombat ↔ LandTactics**: Independent apps, shared engine infrastructure (`cpp/common/`).

## Open Issues

- Damage model Python implementations are in `docs/` not `python/` — unusual placement
- `LandTactics` and `FormationTactics` overlap conceptually (both tactical sims) — relationship unclear
- `cpp/apps/Tanks/notes/ArmorPenetration.md` — separate armor penetration model for tanks, not unified with FormationTactics damage model
- `cpp/common/engine/Notes/ShotHitDetection.md` — hit detection algorithm, not integrated into damage model audit
- No automated tests for any combat model physics
- Multiple TODO files in app directories: `cpp/apps/AeroCombat/Notes/TODO.md`, `cpp/apps/FormationTactics/Notes/ToDo.md`, `cpp/apps/LandTactics/notes/TODO.md`, `cpp/apps/OrbitalWar/Notes/TODO.md`, `cpp/apps/Tanks/TODO.md`
