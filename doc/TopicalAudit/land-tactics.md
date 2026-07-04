---
type: TopicalAudit
title: Land Tactics
tags: [topic, cpp, military, simulation, units, weapons, armor, damage, combat, tactics, terrain, formation]
---

## Summary

Land Tactics is a 2D military simulation with units (infantry, guns, tanks, assault guns, helicopters), weapon types (kinetic, HEAT, fragmentation, fire), armor modeling (directional: front/side/back/top/bottom), and combat resolution via hit probability + penetration + damage ramp. Units are `RigidBody2D` with goal-seeking movement, suppression mechanics, and turret rotation. Includes `AirCombatModel` for aircraft performance (turn rate, climb rate, max speed from aerodynamic polar). `LTWorld` holds factions, squads, shelters, and terrain.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/apps/LandTactics/LandTactics_main.cpp` | active | `FormationTacticsApp(AppSDL2OGL)`: GUI, world, rendering, factions, squads. Includes `AirCombatModel.h`. |
| C++ | `cpp/apps/LandTactics/LTUnitType.h` | active | `LTGunType` (rps, nburst, vMuzzle, caliber, pMass, AP, ExDMG, dmgType, spread, ballistic coef). `LTUnitType` (kind, size, mass, maxSpeed, enginePower, armor[5], HP, guns[nGunMax]). `UnitKind` enum: inf/gun/tank/stug/heli. `DmgTypes` enum: KE/HEAT/Frag/Fire. |
| C++ | `cpp/apps/LandTactics/LTUnitType.cpp` | active | `fromString()` parsing for gun and unit types from semicolon-delimited format. `updateAux()` computes crossArea, ballisticCoef, szAreas. |
| C++ | `cpp/apps/LandTactics/LTUnit.h` | active | `LTGun` (type, n, ammo, timer, attachedTo). `LTUnit(RigidBody2D)`: alive, wound, suppressed, turret_dir, willForce, goal_pos, target. Methods: `move_to_goal()`, `fireGun()`, `getShot()`, `getProjectedArea()`, `update()`. |
| C++ | `cpp/apps/LandTactics/LTUnit.cpp` | active | Implementation: movement (maxSpeed toward goal), firing (hit probability from spread vs projected area, velocity decay, kinetic damage), damage resolution (penetration vs armor, damage ramp, suppression). |
| C++ | `cpp/apps/LandTactics/LTShelter.h` | active | Shelter/cover objects. |
| C++ | `cpp/apps/LandTactics/LTFaction.h` | active | Faction management: squads, units. |
| C++ | `cpp/apps/LandTactics/LTWorld.h` | active | World state: factions, terrain, squads. |
| C++ | `cpp/apps/LandTactics/LTcommon.h` | active | Common definitions, enums, constants. |
| C++ | `cpp/apps/LandTactics/LTrender.h` | active | Rendering functions for units, terrain, UI. |
| C++ | `cpp/apps/LandTactics/LTdraw.h` | active | Drawing helpers. |
| C++ | `cpp/common/CombatModels/AirCombatModel.h` | active | `CombatAirCraft`: mass, power, wing/hull area, aspect, CD0, CLmax. `maxSpeed_simple()`, `trunRate_simple()` (optimal CL for max lift), `climbRate_simple()`, `climbRate_CLmax()`. `getPoperTrust()` (propeller thrust from kinetic energy). `makeClimbCurves()`. |
| C++ | `cpp/apps/OrbitalWar/spaceTactics.cpp` | active | `SpaceTactics(AppSDL2OGL_3D)`: space warfare simulation with planets, ships, weapons, trajectories. See `formation-tactics.md`. |
| C++ | `cpp/apps/AeroCombat/` | active | Aero combat app. See `tanks.md`. |
| C++ | `cpp/apps/BlockHouseTactics/` | active | Block house tactics app. |
| Doc | `docs/AeroCombat/AeroCombat.md` | doc | Aero combat design. |
| Doc | `docs/AeroCombat/AeroSurf.md` | doc | Aerodynamic surface model. |
| Doc | `encyclopedia/space_warfare/01_Tactics.md` | doc | Space warfare tactics: long-range vs close-range, maneuvering, bluff, internal/external maneuvers, ship reflexes, broadside fire. |

## Sub-topics

### Unit Types & Weapons

`LTGunType`:
- `rps` (rounds/sec), `nburst`, `vMuzzle`, `caliber`, `pMass`
- `crossArea = π*caliber²/4`, `ballisticCoef = crossArea/pMass`
- `AP` (armor penetration [mm]), `ExDMG` (explosive damage [J])
- `dmgType`: KE / HEAT / Frag / Fire
- `spread` (tg angular spread)
- `getVelocityDecay(dist) = exp(-ballisticCoef*dist)` — aerodynamic deceleration
- `getKineticDamage(dist, dHeight) = pMass * (vMuzzle²*decay² - dHeight*g)`
- `getSpread(dist) = (spread*dist)²`

`LTUnitType`:
- `kind`: inf / gun / tank / stug / heli
- `sz` (3D bounding box), `szAreas` (projected areas per face)
- `mass`, `maxSpeed`, `enginePower`
- `armorFront/Side/Back/Top/Bottom` — directional armor
- `HP`, `recovery_rate`
- `guns[nGunMax]` (up to 3)

### Combat Resolution

`LTUnit::fireGun(i, target)`:
1. Compute range `r` and direction `hdir`
2. `hit_prob = 1/(1 + (spread(r)/projectedArea)²)` — hit probability from spread vs target area
3. `velocity = getVelocityDecay(r)`, `Ek = getKineticDamage(r, dh)`
4. Call `target.getShot(hdir, nburst, area, ExDMG, Ek, AP)`

`LTUnit::getShot(from, nburst, area, damage_const, damage_kinetic, penetration)`:
1. `pass = (penetration - armorFront) / penetration` — penetration ratio
2. If `pass < 0`: no damage (armor too thick)
3. `damage = pass * damage_kinetic + damage_const`
4. `hit_prob = hit_area / (hit_area + area)` — hit probability from area ratio
5. `kill_prob = hit_prob * damage_ramp(damage, HP)`
6. `suppressed -= kill_prob` — suppression effect
7. For each burst round: `randf() > kill_prob` → skip; else `wound = 1`

### Movement & Suppression

`LTUnit::move_to_goal(dt)`:
- Direction toward `goal_pos`, speed = `maxSpeed * dt`
- If close enough: snap to goal
- `update(dt)`: suppressed decays at `recovery_rate`, effective dt reduced by suppression

`LTUnit::update(dt)`:
- `suppressed -= recovery_rate * dt; if < 0: reset`
- `dt *= (1 - suppressed)` — suppressed units act slower
- Gun timers increment by `rps * dt`
- `move_to_goal(dt)`

### Air Combat Model

`CombatAirCraft`:
- `maxSpeed_simple()`: `v = (P / (0.5*ρ*area*CD0))^(1/3)` — max speed ignoring gravity
- `trunRate_simple()`: optimize CL for max lift → `CL = sqrt(3*CD0*area_tot / (area_wing/(π*aspect)))`, clamped to CLmax. Then `speed`, `Flift`, `accel`, `Rturn = speed²/accel`, `ω = speed/Rturn`
- `climbRate_simple()`: `vy = c1 - c2*v² - c3/v²` where `c1=P/(mg)`, `c2=CD0*A2/(mg)`, `c3=mg/(A1*π*aspect)`
- `climbRate_CLmax()`: climb rate at stall speed (CL=CLmax)
- `getPoperTrust()`: propeller thrust from `P = 0.5*area*ρ*v0*(2*v0*Δv + Δv²)`

### Space Warfare

`SpaceTactics` (`spaceTactics.cpp`):
- `SpaceWorld`: planets, ships, combatants, guns, projectiles
- Orbital mechanics: planet trajectories, ship trajectories
- Weapons: rail guns, lasers, whipple shields, ablation rockets
- `CombatAssembly`: target + gun + collision simulation
- `SpaceShipMobility`: engine Isp, propellant, delta-V
- See `encyclopedia/space_warfare/` for tactical doctrine

## Parity Status

- **`LTUnit` combat ↔ `DamagePhysicsModel.py`**: Different domains. LTUnit is top-down 2D with simplified penetration. DamagePhysicsModel is detailed energy-based with armor toughness, cutting/piercing/crushing. No formal parity.
- **`AirCombatModel` ↔ `AeroSurf`**: `AirCombatModel` uses simplified analytical formulas. `AeroSurf` provides aerodynamic polar data. `AirCombatModel` can use `AeroPolar` in `makeClimbCurves()`.
- **`SpaceTactics` ↔ `LandTactics`**: Different domains (space vs ground). No parity.

## Open Issues

- `LTUnit::getShot()` uses `armorFront` only — ignores directional armor (side/back/top/bottom)
- `LTUnit::getProjectedArea()` uses `szAreas.dot(from)` — but `from` is not normalized to unit sphere
- `LTUnit::fireGun()` has `if(gun.timer > 0.0)` — should be `>= 1.0` for ready-to-fire (timer counts up by `rps*dt`)
- `LTGunType::fromString()` has `sscanf` with 5 `%lf` but only 4 variables in format string (`rps, nburst, spread, vMuzzle, pMass` — `nburst` is `%li`)
- `LTUnitType::fromString()` uses `exit(-1)` on gun not found — should handle gracefully
- `AirCombatModel::trunRate_simple()` — typo in name (should be `turnRate`)
- `AirCombatModel::climbRate_simple()` comment says "ToDo: Consider limited lift coeff CL (stall)" — not fully resolved
- `makeClimbCurves()` computes forces but doesn't output anything — incomplete
- `SpaceTactics` constructor has `DEBUG` macros and `exit(0)` calls — not production-ready
- No unit test framework for combat resolution — damage/hit probability untested
- `LandTactics_main.cpp` class is named `FormationTacticsApp` — naming inconsistency with app name
