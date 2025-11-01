# `spaceCombat.h`

This header file is the core of the combat mechanics in `SpaceTactics`. It defines the data structures for various weapon systems, projectiles, defenses, and the logic for simulating combat encounters.

## Weapon and Defense Structures

### `whippleShieldType` Struct

Models a Whipple shield, a type of hypervelocity impact shield used to protect spacecraft.

*   `int n`: Number of layers in the shield.
*   `double layerDens`: Areal density of each layer (kg/m²).
*   `double spacing`: Spacing between layers (m).
*   `double critEdens`: Critical energy density (J/m²) required to defeat a layer.
*   **`void impact(double& m, double& R, double& v, double& vT)`**: Simulates the impact of a projectile on the shield. It takes the projectile's mass, radius, normal velocity, and tangential velocity as references and modifies them as the projectile passes through the shield layers, fragmenting and losing energy.

### `ProjectileType` Struct

Defines the basic properties of a kinetic projectile.

*   `double mass`: Mass of the projectile (kg).
*   `double caliber`: Caliber (diameter) of the projectile (m).

### `SpaceGunType` Struct

Defines the performance characteristics of a kinetic weapon (e.g., a railgun).

*   `double length`: Length of the accelerator barrel (m).
*   `double maxForce`: Maximum force the gun can exert on the projectile (N).
*   `double maxPower`: Maximum power the gun can deliver (W).
*   `double scatter`: Angular uncertainty of the gun's aim (in radians, as tan(α)).
*   `double fireRate`: The rate of fire in Hz (shots per second).
*   **`double getMuzzleVelocity(ProjectileType* shotType, double& t)`**: Calculates the final muzzle velocity of a projectile. It models a two-stage acceleration: a force-limited phase followed by a power-limited phase. It also returns the time `t` the projectile spends in the barrel.

### `SpaceGun` Struct

Represents one or more instances of a specific gun type on a combatant.

*   `int n`: The number of guns of this type.
*   `ProjectileType* shotType`: Pointer to the projectile type used by this gun.
*   `SpaceGunType* gunType`: Pointer to the gun's performance characteristics.

### `ProjectedTarget` Struct

Represents a target as seen from a specific direction, modeling its cross-sectional area and damage tolerance.

*   `double area`: The target's cross-sectional area (m²).
*   `double HPs`: The total hit points (in Joules) of the target.
*   `double HPexponent`: An exponent for the damage calculation, allowing for non-linear damage responses.
*   `double damageTolerance`: The amount of energy a target can absorb from a single hit before taking damage.
*   `whippleShieldType* wshield`: A pointer to the target's Whipple shield.
*   `double health`: The remaining health of the target (0.0 to 1.0).
*   `double damage`: The accumulated damage in Joules.
*   `void getDamage(double E_Damage, double nhit)`: Calculates and applies damage to the target based on the energy of incoming hits.

## Combat Simulation Classes

### `SpaceSalvo` Struct

Represents a volley of projectiles fired from a gun.

*   `int n`: The number of projectiles in the salvo.
*   `double speed`: The muzzle velocity of the projectiles.
*   `double scatter0`: The initial angular scatter of the salvo.
*   `double delay`: The time the projectile spends in the barrel before launch.
*   `double getSpread(double dist, double aDelta, double& t)`: Calculates the radius of the projectile spread at a given `dist`, considering the initial scatter and the maneuvering acceleration `aDelta` of the firing platform. It also returns the total time `t` to target.

### `CombatAssembly` Class

Orchestrates a single combat event between multiple combatants.

*   `std::vector<SpaceSalvo> salvos`: A list of all salvos fired in the encounter.
*   `std::vector<ProjectedTarget> targets`: A list of all targets involved.
*   `void fireGun(SpaceGun& guns, double burstTime)`: Adds a salvo to the simulation, calculated from a specific `SpaceGun` firing for a `burstTime`.
*   `void addTarget(ProjectedTarget& t)`: Adds a target to the encounter.
*   `void collide(double dist, double aDelta)`: The main simulation method. It calculates the spread of each salvo, determines hit probabilities for each target, simulates the effect of Whipple shields, and applies damage.

## Propulsion Structures

### `PulseEngine` Struct

Models a pulsed propulsion system, such as a nuclear pulse engine.

*   `double fuel_energy_density`: Energy per kg of the fuel (J/kg).
*   `double burnUp`: The fraction of fuel that is consumed in each pulse.
*   `double nozzle_efficiency`: The efficiency of converting thermal energy to kinetic energy.
*   `double pulse_fuel`: Mass of fuel per pulse (kg).
*   `double pulse_inert`: Mass of inert propellant per pulse (kg).
*   `double rate`: The pulse rate in Hz.
*   **`void evalPulse(double& ve, double& thrust)`**: Calculates the effective exhaust velocity (`ve`) and average thrust for the engine.

### `SpaceShipMobility` Struct

Models the overall mobility of a spaceship, potentially blending two different engine types (e.g., high-Isp and low-Isp).

*   `PulseEngine engine_hiIsp`, `engine_loIsp`: The two engine types.
*   `double mass_propletant`, `mass_fuel`, `mass_empty`: The mass breakdown of the ship.
*   **`double evalDeltaV(double time, double fIsp)`**: Calculates the total delta-v the ship can achieve in a given `time` by blending the two engines with a factor `fIsp`.
