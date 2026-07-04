# SailWar

Simulation of sailing ships — naval combat from the age of sail. 2D side-view with realistic sailing physics (wind, sail trim, buoyancy) and naval combat (cannon fire, damage).

## Physics

- Sailing aerodynamics: wind force on sails based on sail angle and wind direction
- Buoyancy and hydrodynamic drag
- Projectile ballistics for cannon fire
- Ship-to-ship collision

## Files

- **SailWar_main.cpp** — main application: naval combat, ship control, wind, cannon firing
- **GameWorld.cpp / .h** — world model: ships, wind, projectiles, combat
- **GameScreen.cpp / .h** — game screen/UI management
- **Frigate2D.cpp / .h** — frigate ship model: hull, sails, cannons, buoyancy
- **Yacht2D.cpp / .h** — yacht ship model: simpler sailing vessel
- **Gun.cpp / .h** — cannon: aiming, firing, reload
- **Projectile.cpp / .h** — cannonball: ballistics, impact, damage
- **test_SailPolar.cpp** — test: sail polar diagram (speed vs. wind angle)
- **test_buoyancy.cpp** — test: buoyancy calculation for ship hull
- **CMakeLists.txt** — build targets: `SailWar_main`, `test_SailPolar`, `test_buoyancy`
- **data/** — ship configuration
- **servis/** — service/debugging scripts
