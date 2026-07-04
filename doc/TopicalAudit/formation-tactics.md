---
type: TopicalAudit
title: Formation Tactics
tags: [topic, cpp, python, military, simulation, formation, soldier, melee, damage, polearm, battle-line, faction, collision-buffer]
---

## Summary

Formation Tactics simulates pre-modern infantry combat with soldiers arranged in formations (rectangular grids). Each soldier is a `Body2D` with melee/ranged combat, stamina, moral, shields, and armor. Formations have bounding boxes, target movement, will-force application, and leave-men-behind logic. Interactions accelerated by `TileBuffer2D` collision buffer. `BattleLine` coordinates multiple formations along a line. Python damage physics model (`DamagePhysicsModel.py`) provides detailed energy-based armor penetration (crushing/cutting/piercing). Polearm physics model (`PolearmsPhysicsModel.py`) computes swing/thrust kinematics.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/apps/FormationTactics/FormationTactics_main.cpp` | active | `FormationTacticsApp(AppSDL2OGL)`: GUI, world update, rendering, faction/formation selection, tiled view. |
| C++ | `cpp/apps/FormationTactics/FormationWorld.h` | active | `FormationWorld`: formations, factions, soldierTypes, collision buffer (`TileBuffer2D<Soldier*,16,16,64>`), terrain (`TerrainCubic`), sim params (dt, damping). |
| C++ | `cpp/apps/FormationTactics/FormationWorld.cpp` | active | `simulationStep()`: eliminate invalids, clean temp, update bbox, formation interactions (brute or buffered), update formations. `formationInteractions_buff()` — collision buffer accelerated. |
| C++ | `cpp/apps/FormationTactics/Formation.h` | active | `Formation`: nrows/ncols grid of soldiers, bbox, center, target movement, will-force, deploy, interact. `Salvo` class for ranged fire. |
| C++ | `cpp/apps/FormationTactics/Formation.cpp` | active | `interact()`: BBox overlap check, enemy/friend interactions. `deploySoldiers()`: grid layout with hex offset. `setupSoldiers()`: type assignment. |
| C++ | `cpp/apps/FormationTactics/Soldier.h` | active | `Soldier(Body2D)`: moral, stamina, time_buf, shield, impair_mask, opponent. `enemy_interaction()` (melee push, opponent selection), `friend_interaction()` (collision push), `attack_melee()`, `getHit()`, `getShot()`, `getDamage()`. |
| C++ | `cpp/apps/FormationTactics/Soldier.cpp` | active | Implementation of soldier update, combat, interactions. |
| C++ | `cpp/apps/FormationTactics/SoldierType.h` | active | `SoldierType`: mass, radius, speed, stamina_regain, melee params (skill/period/range/damage/penetration/fStamina/push), defence params, fire params (period/range/damage/penetration/spread/ammo), armor (Fw/Bk), shield params. `fromString()` semicolon-delimited parser. |
| C++ | `cpp/apps/FormationTactics/Faction.h` | active | `Faction`: name, color, battleLines, formations. `initFaction()` — create formations along battle line, deploy soldiers. `getFormationAt()` — nearest formation to point. |
| C++ | `cpp/apps/FormationTactics/BattleLine.h` | active | `BattleLine`: formations along a line, `setTargetLine()`, spacing. |
| C++ | `cpp/apps/FormationTactics/FormationTacticsCommon.h` | active | Common definitions, constants. |
| Python | `docs/FormationTactics/PolearmsPhysicsModel.py` | active | Polearm performance: swing (torque→ω→v_head→KE), 1H thrust (force→v→KE), 2H thrust. Weapons: Spear, Pike, Glaive, Halberd, Poleaxe. |
| Python | `docs/FormationTactics/DamagePhysicsModel.md` | active | Energy-based damage model: crushing (blunt), cutting (edge length vs toughness), piercing (cone step-through). Targets: head/torso/limb × unarmored/cloth/steel. |
| Python | `docs/FormationTactics/DamagePhysicsModel.py` | active | Python implementation of damage model from the markdown. |
| Doc | `docs/FormationTactics/PolearmsPhysicsModel.md` | doc | Polearm physics model documentation. |

## Sub-topics

### Soldier Combat

`Soldier::enemy_interaction(sj, melee)`:
- If melee: check `melee_range + Rj`, compute direction, select opponent by `c_dir² * invr` score, push force based on facing
- If not melee: standard collision push
- `attack_melee(enemy)`: costs `melee_period` time and `melee_fStamina` stamina, calls `enemy->getHit(damage, penetration, skill, rot)`

`Soldier::getHit(bare_damage, penetration, skill, dir)`:
- Currently simplified: `getDamage(bare_damage)` directly (defense/armor code commented out)
- `getDamage(damage)`: `xdmg = damage/damage_tolerance`, probabilistic injury via `impair_mask` bits

`Soldier::getShot(bare_damage, penetration, dir)`:
- Shield check: if `rot.dot(dir) > shield_cos`, shield absorbs
- Armor: `armor_ramp(penetration, armorFw)` — ratio-based

`Soldier::update(dt, tAttack)`:
- Stamina regeneration, time buffer accumulation
- If `time_buf > tAttack` and has opponent: `attack_melee()`
- Attention direction normalization, velocity damping, force integration

### Formation Movement

`Formation`:
- `setTarget()` / `setEnds()` / `setCenterRot()` — position formation
- `moveToTarget()` — interpolate toward target
- `applyWillForce()` — soldiers steer toward formation center/target
- `leaveMenBehind()` — drop soldiers too far from CoG
- `checkMenBehind()` — BBox check: `r2box > maxBbox2 * (width² + length²)`
- `deploySoldiers()` — grid layout: `dlf = length/ncols`, `dfw = width/nrows`, hex offset on odd rows

### Collision Buffer Acceleration

`FormationWorld::formationInteractions_buff()`:
- For each formation `fi`: setup `Ruler2DFast` around BBox, insert soldiers into `TileBuffer2D`
- For each other formation `fj`: query buffer for nearby soldiers, call `enemy_interaction()` or `friend_interaction()`
- Self-interaction: query buffer for own soldiers, `friend_interaction()`
- Comment: "PROBLEM: may interact several-times with the same neighbor"

### Polearm Physics (Python)

`PolearmsPhysicsModel.py`:
- Swing: `I = (1/3*m_shaft + m_head)*L²`, `ω = sqrt(2*Work/I)`, `v_head = ω*L`, `KE = 0.5*m_head*v_head²`
- 1H thrust: `Work = F*d`, `v = sqrt(2*Work/(m_total + m_dead))`, `KE = 0.5*m_total*v²`
- 2H thrust: same formula with higher force/distance/dead mass
- Weapons: Spear (2.5m, 0.3kg head), Pike (5.0m), Glaive (2.0m, 1.5kg head), Halberd (2.5m, 1.2kg), Poleaxe (1.8m, 2.2kg)

### Damage Physics Model (Python)

`DamagePhysicsModel.py`:
- Crushing: `KE * (m_impact / (m_impact + blunt_mass))` — always applied
- Cutting: `energy_density = KE / edge_length` vs `toughness = yield_strength * thickness²`
- Piercing: step-through cone model, `area = π*(depth*tan(angle/2))²`, energy per step
- Final: `(crushing + penetrating) * damage_multiplier`
- Targets: head/torso/limb × unarmored/cloth/steel

## Parity Status

- **C++ `Soldier::getHit` ↔ Python `DamagePhysicsModel`**: C++ is simplified (direct damage), Python is detailed (crushing/cutting/piercing with armor toughness). No formal parity.
- **C++ `SoldierType` melee ↔ Python `PolearmsPhysicsModel`**: C++ uses abstract damage/penetration values. Python computes from physics (torque, force, mass, length). No formal parity.
- **C++ `FormationTactics` ↔ C++ `LandTactics`**: Both military sims. `LandTactics` uses `LTUnit` with guns/armor. `FormationTactics` uses `Soldier` with melee/ranged. Different combat models. No parity.

## Open Issues

- `Soldier::getHit()` has defense/armor code commented out — simplified to direct damage
- `Soldier::getShot()` computes `damage` but doesn't call `getDamage()` — returns `true` without applying damage
- `formationInteractions_buff()` comment: "may interact several-times with the same neighbor" — no dedup
- `Formation::tAttack` is `bool` but should be `double` (used as `tAttack` in `Soldier::update`)
- `SoldierType::toStrCaptioned()` line 104: prints `min_speed` for `max_speed` — copy-paste bug
- `Formation::~Formation()` uses `delete soldiers` instead of `delete[] soldiers` — array delete mismatch
- `Soldier::Soldier(type_, formation_)` doesn't set `formation` member (parameter name shadowed)
- No ranged combat implementation in `Soldier` (`shoot()` incomplete)
- `Salvo` class declared but not fully implemented
- Python damage model not integrated with C++ simulation
- No unit tests for combat resolution
