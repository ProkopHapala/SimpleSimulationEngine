# SailWar_Multi

Multiplayer version of SailWar — age-of-sail naval combat with UDP networking. Client-server architecture using SDL2_net. Same sailing physics and combat as SailWar, with network synchronization for multiple players.

## Files

- **main.cpp** — entry point: launches client or server mode
- **SailWar_client.cpp** — network client: connects to server, sends input, receives state
- **SailWar_server.cpp** — network server: simulates world, broadcasts state to clients
- **SailWarWorld.cpp / .h** — networked world model: ships, wind, projectiles, state sync
- **Frigate2D.cpp / .h** — frigate ship model (same as SailWar)
- **Yacht2D.cpp / .h** — yacht ship model (same as SailWar)
- **Gun.cpp / .h** — cannon model (same as SailWar)
- **Projectile.cpp / .h** — projectile model with network sync
- **test_SailPolar.cpp** — sail polar diagram test (same as SailWar)
- **test_buoyancy.cpp** — buoyancy test (same as SailWar)
- **CMakeLists.txt** — build targets: `SailWar_client`, `SailWar_server` (requires SDL2_net, `WITH_NET=ON`)
- **data/** — ship configuration
- **servis/** — service/debugging scripts
