# TerrainSimplex Documentation

## Purpose
The TerrainSimplex class simulates terrain generation and erosion processes using simplex grids. It provides tools for creating noise-based terrain, modeling water flow, and applying erosion effects, making it suitable for applications like game development or procedural world generation.

## Algorithms
- **Simplex Noise Generation**: Uses fractional Brownian motion (fBm) with warp noise to create natural-looking heightmaps, scaled and normalized for terrain features.
- **Erosion Simulation**: Implements gradient-based flow erosion (with and without rain), contour-based water outflow, and droplet-based erosion to model terrain degradation and sediment transport, ensuring numerical stability through neighbor comparisons and boundary checks.
- **Indexing and Rasterization**: Employs hexagonal simplex grid indexing for efficient neighbor access and line rasterization, supporting fast queries and updates in the grid.

## Data Structures
- **Grid Arrays**: Stores terrain data in flat arrays (ground[], water[]) with dimensions nx×ny, using row-major indexing (iy*nx + ix)
- **Contour Tracking**: Uses parallel arrays (contour1[], contour2[], known[]) for efficient water flow pathfinding and erosion simulation
- **Neighbor Access**: Predefined neighbor offsets (neighs[]) and distances (neigh_dists[]) for hexagonal grid connectivity
- **Droplet State**: Tracks current droplet position (droplet_ix, droplet_iy) and height (droplet_h) for erosion simulation
- **Erosion Map**: Stores erosion values in a 2D array (erosion_map[]) to track terrain degradation over time
- **Sediment Map**: Tracks sediment transport and deposition in a 2D array (sediment_map[]) to model terrain changes
- **Water Table**: Simulates groundwater levels using a 2D array (water_table[]) to influence erosion and terrain formation
- **Heightmap**: Stores the generated terrain heightmap in a 2D array (heightmap[]) for visualization and further processing
- **Boundary Conditions**: Stores boundary conditions for the terrain simulation in a 2D array (boundary_conditions[]) to define the simulation domain

## Functions
- `raster_line()`: Rasterizes a line across the simplex grid to find intersection points and boundaries.
- `genTerrainNoise()`: Generates terrain height using simplex noise with parameters for scale, height scaling, frequency decay, and strength.
- `init_outflow()`: Initializes the water flow simulation by setting high water values and resetting known states.
- `outflow_step()`: Advances the water flow contour by extending paths to lower neighbors, simulating water movement.
- `extend_path()`: Evaluates and extends water flow paths if a lower elevation neighbor is found, updating water levels.
- `initErrosion()`: Sets up erosion parameters by clamping ground heights and initializing water levels.
- `flow_errosion_step()`: Performs a single erosion step with water flow, moving sediment based on gradients.
- `rain_and_evaporation()`: Simulates rain addition and water evaporation/mixing to update water levels across the grid.
- `flow_errosion_step_noRain()`: Conducts erosion without rain, counting pixels with water flow and updating terrain.
- `initDroplet()`: Initializes a droplet for erosion simulation at a random grid position.
- `droplet_step()`: Moves a droplet to the lowest neighboring height, eroding the terrain along the way.
- `errodeDroples()`: Runs multiple droplet erosion simulations for a specified number of droplets and steps.