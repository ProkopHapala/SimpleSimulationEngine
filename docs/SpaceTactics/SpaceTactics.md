# SpaceTactics: Code Map and Documentation

This document provides a comprehensive overview and code map of the `SpaceTactics` project, a space combat simulation game with a focus on realistic orbital mechanics and physics.

## 1. Overall Concept & Gameplay

`SpaceTactics` is a strategic simulation where players command fleets of spaceships. The gameplay abstracts away from detailed ship micromanagement and instead emphasizes the planning and execution of complex orbital maneuvers.

### Core Simulation Concepts:

*   **N-Body Physics:** The simulation is built on an N-body gravitational model, managed by the `SpaceWorld` class. All celestial bodies (planets, moons) and ships exert gravitational forces on each other, requiring players to account for complex gravitational interactions.
*   **Trajectory-Based Movement:** Ship movement is not based on direct, continuous thrust. Instead, players define trajectories using splines (`SplineManager`). This allows for planning long-term maneuvers like:
    *   Low-thrust orbital transfers.
    *   Gravitational slingshot (fly-by) maneuvers.
    *   Complex orbital insertions.
    The `trjThrust` array in the `SpaceBody` class stores the thrust vectors applied over time along these splines.
*   **Physically Realistic Weapons:** The game models various weapon systems with a basis in real-world physics concepts:
    *   **Railguns/Mass Drivers:** Modeled via `SpaceGunType` and `ProjectileType`. The simulation calculates muzzle velocity based on the gun's length, maximum force, and power (`getMuzzleVelocity`). Projectile damage is determined by kinetic energy upon impact.
    *   **Lasers:** The code includes functions (`print_Laser`, `difractionLimit_spot`) to calculate laser spot size and intensity over vast distances, considering diffraction limits based on aperture and wavelength.
    *   **Whipple Shields:** The `whippleShieldType` struct models a realistic defense against kinetic projectiles. It simulates the impact on multiple spaced layers, calculating how the projectile fragments and loses energy (`impact` method).
*   **Asteroid Engineering:** A unique feature, hinted at in `asteroidEngineering.h`, involves the potential for utilizing asteroids. This could include mining for resources, converting them into mobile habitats, or using them as kinetic weapons.

## 1.1. Key Features and Capabilities

Beyond the core concepts, the simulation incorporates several advanced features that provide a rich and realistic strategic environment.

### Advanced Combat Physics

The combat model goes beyond simple hitpoint attrition, incorporating detailed physics for both offense and defense.

*   **Dynamic Weapon Accuracy:** The accuracy of kinetic weapons is not static. The spread of a projectile salvo is dynamically calculated based on the distance to the target and, crucially, the maneuvering of the firing ship. Firing while accelerating or turning will increase the spread, creating a realistic trade-off between evasion and accuracy.
*   **Layered Shield Simulation:** Defense is modeled with `whippleShieldType`, which simulates multi-layer Whipple shields. When a projectile hits, the model calculates its fragmentation and energy loss as it passes through each layer. This allows for nuanced armor design and weapon effects, where some projectiles might be stopped completely while others punch through with reduced energy.
*   **Realistic Railgun Performance:** The performance of kinetic accelerators (railguns) is not a simple constant. The `getMuzzleVelocity` function models a two-stage acceleration process: an initial phase limited by the maximum force the projectile can withstand, followed by a phase limited by the weapon's maximum power output. This captures the real-world engineering constraints of such weapons.

### Sophisticated Propulsion Models

Ship movement and strategic mobility are governed by detailed and realistic propulsion models.

*   **Advanced Pulse Engines:** The simulation includes models for `PulseEngine` systems, such as nuclear pulse propulsion. These are defined by physical parameters like fuel energy density, burn-up fraction, and nozzle efficiency.
*   **Hybrid Propulsion Systems:** The `SpaceShipMobility` model allows for blending two different engine types (e.g., a high-thrust, low-efficiency engine and a low-thrust, high-efficiency one). This enables ships to have flexible mission profiles, using different engine settings for rapid maneuvers versus long-duration, efficient cruises.

### Realistic Celestial Mechanics

The simulation environment is built on a robust and flexible model of celestial mechanics.

*   **Dual Orbit Representation:** The engine can represent celestial bodies in two ways: through classical `OrbitalElements` for predictable, analytical Keplerian orbits, or through fully integrated N-body physics (`trjPos`). This allows for a mix of stable, long-term orbits and dynamically simulated trajectories for ships and other objects.
*   **Real-World Data Integration:** The simulation can load real-world asteroid data directly from the `astorb` database format. This feature allows for the creation of highly realistic scenarios populated with thousands of known asteroids, each with its correct orbital parameters.

### Strategic Resource Management

The framework includes features for in-depth strategic resource management, particularly concerning asteroids.

*   **Detailed Asteroid Composition:** Asteroids are not just generic resource nodes. The `asteroidEngineering.h` module defines a system of `MineralType` and `RockType`, allowing asteroids to have a detailed geological composition. This creates a basis for a complex mining and resource-refining economy.
*   **Physics-Based Asteroid Moving:** The cost and difficulty of moving an asteroid are modeled realistically. The `manuever_planeChange` function calculates the immense delta-v and propellant mass required to alter an asteroid's orbit, making "asteroid wrangling" a significant, high-investment strategic objective rather than a trivial task.

## 2. Code Map & Module Breakdown

The project's functionality is modularized into several key files. Here is a breakdown of the most relevant components.

### `cpp/apps/OrbitalWar/spaceTactics.cpp`

This is the main application file.

*   **Role:** Entry point of the game, handles the main game loop, rendering (using SDL2 and OpenGL), and user input.
*   **Key Class:** `SpaceTactics` inherits from `AppSDL2OGL_3D` and orchestrates the entire simulation.
*   **Functionality:**
    *   **Initialization:** Sets up the game world by creating celestial bodies (`world.addPlanet`) and ships (`world.addShip`). The initial positions and velocities create a scenario (e.g., ships around Jupiter's moons).
    *   **Game Loop (`draw`, `drawHUD`):** Renders the state of the `SpaceWorld`. It draws the trajectories of all bodies and displays ships and planets at their positions corresponding to the current simulation time (`timeCur`).
    *   **Time Control:** The HUD includes a timeline that the player can scrub through to view the predicted trajectories at different points in time.
    *   **Weapon Tests:** The `weapon_tests` function and various `print_` functions serve as a sandbox for testing and balancing the physics of different weapon systems.

### `cpp/common/Orbital/SpaceWorld.h`

This file defines the core of the simulation environment.

*   **Role:** Manages the state of the entire simulation, including all celestial bodies, ships, and their interactions.
*   **Key Class:** `SpaceWorld`
*   **Functionality:**
    *   **Object Management:** Holds `std::vector`s of `planets` and `ships`.
    *   **Physics Integration:** Implements the `ODEderivObject` interface. The `getDerivODE` method calculates the gravitational forces between all bodies and applies thrust from ship trajectories. It uses a Runge-Kutta-Fehlberg (RKF45) integrator (`ode`) to step the simulation forward in time.
    *   **Trajectory Prediction (`predictTrjs`):** This is a crucial function that runs the physics simulation for a specified number of steps (`n`) with a given time delta (`dt`) to pre-calculate the future paths of all objects. The results are stored in the `trjPos` arrays of each `SpaceBody`.
    *   **Combat Simulation (`evalDamage`):** This function simulates a combat encounter at a specific point in time. It gathers all combatants, calculates weapon firing solutions, projectile flight time, and resulting damage based on the `spaceCombat.h` module.

### `cpp/common/Orbital/SpaceBodies.h`

This file defines the fundamental objects that exist within the `SpaceWorld`.

*   **Role:** Defines the data structures for all physical objects in the simulation.
*   **Key Class:** `SpaceBody` (inherits from `RigidBody`)
*   **Functionality:**
    *   **Attributes:** Stores physical properties like `mass`, `radius`, `pos` (position), and `vel` (velocity).
    *   **Trajectory Storage:** Contains pointers to arrays that store the pre-calculated trajectory data:
        *   `trjPos`: A series of `Vec3d` points representing the body's position over time.
        *   `trjThrust`: A series of `Vec3d` vectors representing the thrust applied at each step of the trajectory.
    *   **Orbital Representation:** Can also represent objects via classical `OrbitalElements` and an `Orbit` class, allowing for setup from astronomical data (e.g., `fromString_astorb` for asteroids).

### `cpp/common/math/SplineManager.h`

This file provides the tools for creating and managing the splines used for ship trajectories.

*   **Role:** Manages Hermite splines for defining smooth paths.
*   **Key Class:** `SplineManager`
*   **Functionality:**
    *   **Control Points:** Stores control points (`CPs`) and their corresponding times (`ts`).
    *   **Evaluation:** The `eval` and `evalIt` methods allow for sampling the spline at any given time `t` to get the position, velocity (`dval`), and acceleration (`ddval`). This is how a ship's planned trajectory is translated into concrete positions and thrusts for the physics engine.
    *   The function `nonUni2spline` is used in `spaceTactics.cpp` to convert a series of discrete thrust commands over time into a smooth `trjThrust` profile for a ship.

### `cpp/common/Orbital/spaceCombat.h`

This file contains the detailed implementation of weapon physics and combat mechanics.

*   **Role:** Defines the building blocks of combat, from weapon types to damage calculation.
*   **Key Structs & Classes:**
    *   `ProjectileType`: Defines a projectile's basic properties (`mass`, `caliber`).
    *   `SpaceGunType`: Defines the performance characteristics of a kinetic weapon (`length`, `maxForce`, `maxPower`, `scatter`).
    *   `whippleShieldType`: Models a multi-layer shield and its interaction with projectiles.
    *   `ProjectedTarget`: Represents a target's cross-section and hit points from a specific direction.
    *   `SpaceSalvo`: Represents a volley of projectiles fired from a `SpaceGun`. It calculates the spread of the salvo over distance.
    *   `CombatAssembly`: A temporary class that orchestrates a single combat event. It takes salvos and targets, calculates hit probabilities, and applies damage.

### `cpp/common/dynamics/appliedPhysics.h`

This file is a collection of utility functions for various physics calculations.

*   **Role:** Provides a library of fundamental physics formulas used throughout the simulation.
*   **Functionality:**
    *   **Gravity:** `centralGravityForce` calculates the force vector between two masses.
    *   **Rocketry:** Includes the Tsiolkovsky rocket equation (`tsielkovsky_speed`, `tsielkovsky_payload`) and formulas for exhaust velocity (`exhaustVelocity`).
    *   **Constants:** Defines a wide range of physical constants (`const_Graviational`, `const_AU`, etc.).

### `cpp/common/Orbital/asteroidEngineering.h`

This file lays the groundwork for mechanics related to asteroids.

*   **Role:** Defines data structures for asteroid composition and potential economic/strategic actions.
*   **Key Structs:**
    *   `RockType`, `MineralType`: Define the geological makeup of an asteroid.
    *   `Asteroide`: A class that could represent a game object with specific mineral deposits and an orbit.
*   **Functionality:**
    *   The function `manuever_planeChange` calculates the energy (`deltaV` and propellant mass) required to change an asteroid's orbital plane, hinting at the strategic possibility of moving asteroids.

## 3. Path to Python Integration

The current structure is well-suited for being refactored into a C++ library with Python bindings.

*   **Core Library:** The classes in `SpaceWorld`, `SpaceBodies`, `spaceCombat`, and `appliedPhysics` form a strong foundation for a core simulation library.
*   **Scripting Interface:** Python could be used to:
    *   **Define Scenarios:** Set up initial conditions (planets, ships, orbits) without recompiling C++ code.
    *   **Control AI:** Script the strategic decisions for AI opponents.
    *   **Data Analysis & Visualization:** Use libraries like `numpy` and `matplotlib` to analyze simulation output, plot trajectories, and visualize combat results, which is currently done with basic `printf` statements and simple OpenGL rendering.
    *   **Rapid Prototyping:** Quickly test new game mechanics, weapon parameters, or ship designs in Python before implementing them in C++.

## 4. Python Ctypes Interface

A C-style interface has been added to allow controlling the `SpaceTactics` simulation from Python using the `ctypes` library. This is useful for scripting scenarios, AI, and for data analysis.

The C functions are exposed in `cpp/libs/OrbitalWar/SpaceTacticsLib.cpp`. A Python wrapper can be created similar to `python/LandCraft/landcraft.py` to interact with the simulation library.

### 4.1. Loading the Library

The compiled shared library (e.g., `libSpaceTacticsLib.so` or `SpaceTacticsLib.dll`) needs to be loaded in Python. The `ctypes` library provides the necessary tools. An example pattern can be found in `python/pyMeta/cpp_utils.py`.

```python
import ctypes
import numpy as np

# Assuming the library is named 'libSpaceTactics.so' and is in the path
lib = ctypes.CDLL('libSpaceTactics.so')
```

### 4.2. API Functions

Here is a description of the available functions.

#### World Management

*   **`sw_init()`**
    *   Initializes the simulation world. Must be called before any other simulation functions.
    *   **Python prototype:** `lib.sw_init.restype = None`

*   **`sw_clear()`**
    *   Clears all objects from the simulation world and resets it to a clean state.
    *   **Python prototype:** `lib.sw_clear.restype = None`

#### Object Management

*   **`sw_add_planet(name, mass, radius, pos, vel)`**
    *   Adds a planet to the simulation.
    *   `name` (c_char_p): Name of the planet.
    *   `mass` (c_double): Mass of the planet.
    *   `radius` (c_double): Radius of the planet.
    *   `pos` (POINTER(c_double)): 3D position vector `[x, y, z]`.
    *   `vel` (POINTER(c_double)): 3D velocity vector `[vx, vy, vz]`.
    *   Returns (c_int): The index of the newly created planet.
    *   **Python prototype:**
        ```python
        from ctypes import c_char_p, c_double, c_int, POINTER
        lib.sw_add_planet.argtypes = [c_char_p, c_double, c_double, POINTER(c_double), POINTER(c_double)]
        lib.sw_add_planet.restype = c_int
        ```

*   **`sw_add_ship(name, mass, radius, pos, vel)`**
    *   Adds a ship to the simulation.
    *   Parameters are the same as for `sw_add_planet`.
    *   Returns (c_int): The index of the newly created ship.
    *   **Python prototype:**
        ```python
        lib.sw_add_ship.argtypes = [c_char_p, c_double, c_double, POINTER(c_double), POINTER(c_double)]
        lib.sw_add_ship.restype = c_int
        ```

*   **`sw_get_n_planets()`**
    *   Returns the current number of planets in the simulation.
    *   Returns (c_int): Number of planets.
    *   **Python prototype:** `lib.sw_get_n_planets.restype = c_int`

*   **`sw_get_n_ships()`**
    *   Returns the current number of ships in the simulation.
    *   Returns (c_int): Number of ships.
    *   **Python prototype:** `lib.sw_get_n_ships.restype = c_int`

#### Simulation and Trajectories

*   **`sw_allocate_trjs(n)`**
    *   Allocates memory for storing `n` steps of trajectory data for all bodies. This must be called before predicting trajectories.
    *   `n` (c_int): The number of trajectory points to allocate.
    *   **Python prototype:** `lib.sw_allocate_trjs.argtypes = [c_int]`

*   **`sw_predict_trjs(n, dt)`**
    *   Runs the N-body simulation for `n` steps with a time step of `dt` to pre-calculate trajectories.
    *   `n` (c_int): Number of steps to simulate.
    *   `dt` (c_double): Time step for the simulation.
    *   **Python prototype:** `lib.sw_predict_trjs.argtypes = [c_int, c_double]`

*   **`sw_set_ship_thrust(ship_idx, n_points, ts, thrusts)`**
    *   Defines a thrust profile for a ship over time using a spline.
    *   `ship_idx` (c_int): The index of the ship to modify.
    *   `n_points` (c_int): The number of control points in the thrust profile.
    *   `ts` (POINTER(c_double)): An array of time values for each control point.
    *   `thrusts` (POINTER(c_double)): A flattened array of 3D thrust vectors `[tx1, ty1, tz1, tx2, ty2, tz2, ...]`.
    *   **Python prototype:**
        ```python
        lib.sw_set_ship_thrust.argtypes = [c_int, c_int, POINTER(c_double), POINTER(c_double)]
        ```

#### Data Access

*   **`sw_get_trj_pos(body_type, body_idx, out_pos, max_points)`**
    *   Retrieves the pre-calculated trajectory positions for a specific body.
    *   `body_type` (c_int): `0` for a planet, `1` for a ship.
    *   `body_idx` (c_int): The index of the body.
    *   `out_pos` (POINTER(c_double)): A pre-allocated numpy array to store the output positions. The shape should be `(max_points, 3)`.
    *   `max_points` (c_int): The maximum number of points to copy into `out_pos`.
    *   Returns (c_int): The number of points actually copied.
    *   **Python prototype:**
        ```python
        lib.sw_get_trj_pos.argtypes = [c_int, c_int, POINTER(c_double), c_int]
        lib.sw_get_trj_pos.restype = c_int
        ```

### 4.3. Example Python Usage

```python
import numpy as np
from ctypes import c_double, c_int, POINTER

# --- Setup ---
# (Load library and set argtypes as shown above)

# --- Initialize simulation ---
lib.sw_init()

# --- Add bodies ---
earth_pos = np.array([1.0, 0.0, 0.0])
earth_vel = np.array([0.0, 0.0, 0.0])
lib.sw_add_planet(b"Earth", 1.0, 0.1, earth_pos.ctypes.data_as(POINTER(c_double)), earth_vel.ctypes.data_as(POINTER(c_double)))

ship_pos = np.array([1.1, 0.0, 0.0])
ship_vel = np.array([0.0, 0.1, 0.0])
ship_idx = lib.sw_add_ship(b"Voyager", 0.001, 0.01, ship_pos.ctypes.data_as(POINTER(c_double)), ship_vel.ctypes.data_as(POINTER(c_double)))

# --- Predict trajectory ---
n_steps = 1000
dt = 0.01
lib.sw_allocate_trjs(n_steps)
lib.sw_predict_trjs(n_steps, dt)

# --- Get trajectory data ---
trj_data = np.zeros((n_steps, 3))
n_copied = lib.sw_get_trj_pos(1, ship_idx, trj_data.ctypes.data_as(POINTER(c_double)), n_steps)

print(f"Copied {n_copied} trajectory points for ship {ship_idx}.")
print(trj_data)

# --- Clean up ---
lib.sw_clear()
```