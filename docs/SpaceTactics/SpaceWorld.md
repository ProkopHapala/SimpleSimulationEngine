# `SpaceWorld.h`

This header file defines the central `SpaceWorld` class, which acts as the main container and manager for the entire space combat simulation. It holds all simulation objects, manages the physics integration, and orchestrates combat encounters.

## `Time` Namespace

A utility namespace for handling time-related constants and conversions.

### Constants

*   `minute`, `hour`, `day`, `week`, `month`, `year`: Double-precision floating-point values representing standard time durations in seconds.
*   `niceUnits`: An array of "nice" time units for formatting time values in a human-readable way.

### Functions

*   **`toStr(double nsec, char* out)`**: Converts a time value in seconds into a human-readable string (e.g., "2.5 days", "1.5 years").
*   **`niceUnitAbove(double val)`**: Returns the smallest "nice" time unit that is larger than the given value.

## `SpaceCombatant` Struct

Represents an active participant in a combat scenario.

*   `int faction`: The faction identifier for the combatant.
*   `SpaceBody* body`: A pointer to the `SpaceBody` representing the combatant's physical presence.
*   `std::vector<SpaceGun> guns`: A list of weapon systems available to the combatant.
*   `double accel`: The maximum maneuvering acceleration of the combatant in m/s².
*   `std::vector<ProjectedTarget> targets`: A temporary list of projected targets. This is likely to be replaced by a more sophisticated 3D targeting system.

## `SpaceWorld` Class

The core class that orchestrates the simulation. It inherits from `ODEderivObject` to provide the derivative function for the ODE solver.

### Public Members

*   `std::vector<SpaceBody> planets`: A vector holding all planetary bodies in the simulation.
*   `std::vector<SpaceBody> ships`: A vector holding all ships in the simulation.
*   `std::vector<SpaceCombatant> combatants`: A vector of all active combatants.
*   `std::unordered_map<std::string, SpaceGunType*> gunTypes`: A map of available gun types, accessible by name.
*   `std::unordered_map<std::string, ProjectileType*> projectileTypes`: A map of available projectile types, accessible by name.
*   `ODEintegrator_RKF45 ode`: The Runge-Kutta-Fehlberg 4(5) ODE solver used for physics integration.
*   `int trj_n`: The number of points in the pre-calculated trajectories.
*   `double trj_dt`: The time step between trajectory points.
*   `double* masses`: An array of masses for all objects, used by the ODE solver.

### Public Methods

*   `void intertialTansform(Vec3d p0, Vec3d v0)`: Applies an inertial transform to all objects in the world, shifting their positions by `p0` and velocities by `v0`.
*   `void addPlanet(...)`: Adds a new planet to the simulation.
*   `void addShip(...)`: Adds a new ship to the simulation.
*   `virtual void getDerivODE(double t, int n, double* Ys, double* dYs)`**: The core derivative function for the ODE solver. It calculates the gravitational forces and thrust for all objects at a given time `t`.
*   `void makeGun(SpaceGun& g, int n, const char* gunName, const char* prjName)`**: A helper function to create a `SpaceGun` instance from predefined types.
*   `double evalDamage(SpaceCombatant& target, double duration, int iTrj, double du)`**: Simulates a combat encounter and calculates the damage received by a `target`.
*   `void allocateODE()`: Allocates memory for the ODE solver based on the number of objects.
*   `void allocateTrjs(int n)`: Allocates memory for the trajectory arrays of all objects.
*   `void objects2ode()`: Copies the state (position, velocity) of all simulation objects into the ODE solver's state array.
*   `void ode2objects()`: Copies the state from the ODE solver back to the simulation objects.
*   `void ode2trjs(int j)`: Copies the state from the ODE solver to a specific point `j` in the trajectory arrays.
*   `void predictTrjs(int n, double dt)`: Runs the simulation for `n` steps with time step `dt` to pre-calculate the trajectories of all objects.
*   `int load_astorb(char* fname, int n_reserve)`: Loads asteroid data from a file in the `astorb` format.
*   `int pickPlanet(const Vec3d& ro, const Vec3d& rd, double epoch)`: A picking function to identify a planet based on a ray cast at a specific epoch.
