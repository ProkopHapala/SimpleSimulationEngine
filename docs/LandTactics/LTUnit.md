# LTUnit.h - Individual Military Unit System

## Purpose
Represents individual military units with realistic physics, combat capabilities, and tactical behavior for squad-based tactical simulation.

## Core Components

### 1. Unit State Management

#### LTUnit Class
- **Inherits**: RigidBody2D - Provides realistic physics simulation
- **Key State Variables** (from `cpp/apps/LandTactics/LTUnit.h`):
  - `LTUnitType* type`
  - `LTGun guns[nGunMax]` with fields per gun: `type`, `n`, `ammo`, `timer`, `attachedTo`
  - `bool alive`
  - `int wound`
  - `double suppressed`
  - `Vec2d turret_dir`
  - `Vec2d willForce`
  - `Vec2d goal_pos`
  - `LTUnit* target`

#### Weapon Systems
- **LTGun Array**: Up to 3 weapons per unit
- **Ammo tracking**: Individual weapon ammunition counts
- **Fire timers**: Rate-of-fire limitations
- **Target assignment**: Current engagement target

### 2. Combat Mechanics

#### Damage System
- **getShot(const Vec3d& from, int nburst, double area, double damage_const, double damage_kinetic, double penetration)**
  - Computes damage application from a shot burst with simple armor/kinetic model.
  - Uses projected area via `getProjectedArea()` to model hit probability; armor and penetration handling are internal to this method [need-check: exact internal formula not shown in header].

#### Fire Control
- **fireGun(int i, LTUnit& target)**: Fires gun `i` (0..nGunMax-1) at target; respects gun timers/ammo.
- **Target selection**: Manual via `target` pointer; no autonomous targeting in header [need-check].
- **Accuracy**: Likely influenced by `suppressed`, weapon `spread`, and range (see `LTGunType`) [need-check].
- **Damage application**: Calls into `getShot()`; armor model specifics in `LTUnitType` and gun parameters.

### 3. Tactical Movement

#### Movement API (local)
- **goal_pos**: Waypoint position stored on unit.
- **willForce**: Accumulated steering/force vector toward goal.
- **move_to_goal(double dt)**: Advances position/velocity toward `goal_pos` using `RigidBody2D` integration.
- Obstacle avoidance / formation keeping are not present in this header [need-check].

#### Physics Integration
- **update(double dt)**: Integrates per-frame state, dispatches fire/move.
- Backed by `RigidBody2D` for kinematics; collision/terrain coupling not defined here [need-check].

## Technical (from code)
Grounded in `cpp/apps/LandTactics/LTUnit.h`:

### Class Structure
- `class LTUnit : public RigidBody2D` - physics-based unit
- `class LTGun` - individual weapon system

### Key Methods
- `getProjectedArea(Vec3d from)` - visible/targetable projected area
- `fireGun(int i, LTUnit& target)` - weapon firing mechanics
- `move_to_goal(double dt)` - goal-directed movement
- `setType(LTUnitType* type_)` - initializes guns from type
- `update(double dt)` - frame-by-frame state updates

### State Tracking
- `LTUnitType* type` - references unit type definition
- `LTGun guns[nGunMax]` - weapon array (nGunMax from LTUnitType)
- `LTUnit* target` - current engagement target

**Notes:**
- Units inherit from RigidBody2D for realistic physics
- Weapon systems are initialized from unit type definitions via `setType()`/ctor
- Suppression system is present as a scalar; its effect on accuracy/movement is implied but not defined here [need-check].
