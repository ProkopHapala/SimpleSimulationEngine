# Problems to solve

## Ship Building 

* **Lua Scripting** - loading ship desings from Lua scripts.
* **Component system**  
   * **Workshop** to store all availabe material and components
* **Editor GUI**
   * Editor Gizmo 
* **Mass and material consumption**
  * Cut the shit into fintite elements made of particular materials. Calculate the density. 

## Truss Dynamics

  * Some space-ship layouts rely on rotation to keep afloat (centrifugal **tensegrity** structures to minimize number of heavy girders)
  * how the spaceship behave during manuevers (e.g. change of orientation)
  * how what effect have impulse (kick) from pulsed nuclear engines (e.g. [Orion project damper](https://en.wikipedia.org/wiki/Project_Orion_(nuclear_propulsion), [Medusa magnetic parachute](https://en.wikipedia.org/wiki/Nuclear_pulse_propulsion#Medusa)
  * Strengh - some components can break
  * Technical - [[Truss Simulation]]

## Damage Model
  * Probability of impact of various modules by projectiles
     * project projectiles on a plane. Store position of modules in a set of grid-points (to check co-location with the projecties)
     * sort modules (as 1D intervals) to check which will be hit first

## **Interaction with radiation**
  * **X-ray & Neutron Scattering** from nuclear engine (nuclear reactor, nuclear pulse capsule), or from other sources (e.g. enemy fire)
    * [[Monte Carlo Raytracing]]
  * **Thermal Radiation Scattering** from radiators, from sun or from enemy nuclear bombs
    * Diffuse reflection ([[Radiosity]]) 
    * Specular reflection

## Magnetic field and plasma dynamics

  * Coupling magnetic field with truss dynamics for tensegrity ship designs (Medusa-like)
  * Magnetic forces acting on conductors
     * Plasma nozzle for nuclear pulse propulsion
     * Magnetic confinement chamber for magnetic fusion propulsion 
  * Voltage generated on metallic components by movement in magnetic field
  * [[Magetic Pulse Plasma Nozzle Solver]]  - Pulse Plasma Nozzle to shield and direct expanding plasma from nuclear pulse (e.g. for inertial nuclear fusion propulsion)
     * we can assume that speed of plasma expansion is much faster than decay of magnetic field (i.e. decay of current) in the coils and wires. Therefore total magnetic field and total current is preserved. Current is iduced in wires and coils to compensate field pushed out by expanding plasma.   

## Orbit simulations

 * [[continuous thrust orbital transfer]] simulations using splines 

# Stand alone programs:

* [SolarSystemMap.cpp](https://github.com/ProkopHapala/SimpleSimulationEngine/blob/master/cpp/apps/OrbitalWar/SolarSystemMap.cpp) - **Map od solar system**s with many asteroides obtained from asteroid database (including realistic orbits and velocities)
* [spaceCraftEditor.cpp](https://github.com/ProkopHapala/SimpleSimulationEngine/blob/master/cpp/apps/OrbitalWar/spaceCraftEditor.cpp) - editor of an spacecraft which can be loaded from lua script.
* [spaceCraftDynamics.cpp](https://github.com/ProkopHapala/SimpleSimulationEngine/blob/master/cpp/apps/OrbitalWar/spaceCraftDynamics.cpp) - **Simulates dynamics of spacecraft** (e.g. when it is rotating or conducting other inner-manuevers)
* [spaceTactics.cpp](https://github.com/ProkopHapala/SimpleSimulationEngine/blob/master/cpp/apps/OrbitalWar/spaceTactics.cpp) - Example of using the engine to simulate a simple tactical situation where one spacecraft intercept other spacecraft around jupiter moons

# Technical tests / sketches:
* [test_OptContinuousThrust.cpp](https://github.com/ProkopHapala/SimpleSimulationEngine/blob/master/cpp/apps/OrbitalWar/test_OptContinuousThrust.cpp) - **Otimization of 
orbital trajectory** with propulsion engines with continuously variable thrust (e.e. ion engines). They are assumed to be limited by power (energy flow), the thrust can be increased at the cost of lower specific impulse (i.e. higher propelant consumption).
    * [TrajectoryVariation.h](https://github.com/ProkopHapala/SimpleSimulationEngine/blob/master/cpp/apps/OrbitalWar/TrajectoryVariation.h) - Solver for trajectroy optimization. 
* [test_SpaceFlightODE.cpp](https://github.com/ProkopHapala/SimpleSimulationEngine/blob/master/cpp/apps/OrbitalWar/test_SpaceFlightODE.cpp) - **simulates flight of space-craft** by integration of equations of motion using [ adaptive runge kuttaODE integrator](https://github.com/ProkopHapala/SimpleSimulationEngine/blob/master/cpp/common/dynamics/ODEintegrator.h). In particular it uses [Runge–Kutta–Fehlberg](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method).
* [orbitEditor.cpp](https://github.com/ProkopHapala/SimpleSimulationEngine/blob/master/cpp/apps/OrbitalWar/orbitEditor.cpp) - editor of trajectories (postion, velocity, acceleration) by mouse.

# Modules:

* [SpaceCraft.h](https://github.com/ProkopHapala/SimpleSimulationEngine/blob/master/cpp/apps/OrbitalWar/SpaceCraft.h) - Defines constrcution of spacecraft from components (girders, ropes, guns, sails, heat-shieds, fuel tanks etc. )
* [SpaceCraftDraw.h](https://github.com/ProkopHapala/SimpleSimulationEngine/blob/master/cpp/apps/OrbitalWar/SpaceCraftDraw.h) - Draws spacecraft made of individual components.
* [SpaceDraw.h](https://github.com/ProkopHapala/SimpleSimulationEngine/blob/master/cpp/apps/OrbitalWar/SpaceDraw.h) - Draw utilities for things like trajectries or planets usefull for space simulations.
* [SpaceWorld.h](https://github.com/ProkopHapala/SimpleSimulationEngine/blob/master/cpp/apps/OrbitalWar/SpaceWorld.h) - Space Simulation world with all relevant objects.
* [SpaceBodies.h](https://github.com/ProkopHapala/SimpleSimulationEngine/blob/master/cpp/apps/OrbitalWar/SpaceBodies.h) - Utilities for simulation of planets and other objects orbiting in solar system.
* [RublePile.h](https://github.com/ProkopHapala/SimpleSimulationEngine/blob/master/cpp/apps/OrbitalWar/RublePile.h) - Generate and store structure of an astoriede in space.
* [Asteroid.h](https://github.com/ProkopHapala/SimpleSimulationEngine/blob/master/cpp/apps/OrbitalWar/Asteroid.h) - Class for reading asteroid properties from [astorb database](ftp://ftp.lowell.edu/pub/elgb/astorb.html).
* [asteroidEngineering.h](https://github.com/ProkopHapala/SimpleSimulationEngine/blob/master/cpp/apps/OrbitalWar/asteroidEngineering.h) - Class to manipulate and store state of modified asteroides during colonization.
* [spaceCombat.h](https://github.com/ProkopHapala/SimpleSimulationEngine/blob/master/cpp/apps/OrbitalWar/spaceCombat.h) - class to simulate space combat. e.g. calculate performance of various guns, impact effects etc.




 

