# LTWorld.h - Tactical World Management System

## Purpose
Provides comprehensive world management for tactical battles, integrating terrain systems, unit management, and environmental simulation for large-scale military operations.

## Core Components

### 1. World Infrastructure

#### Terrain Systems
- **SimplexRuler**: Hexagonal grid for unit positioning
- **Ruler2DFast**: Rectangular grid for spatial queries
- **HydraulicGrid2D**: Water flow and river systems
- **PathFinder**: Optimal route planning between points
- **Terrain elevation**: `ground` array for height-based combat

#### Object Management
- **Static objects**: Buildings, trees, terrain features
- **Linear objects**: Walls, trenches, hedgerows
- **Unit types**: Comprehensive database of military units
- **Gun types**: Weapon system definitions

### 2. Unit Database System

#### Type Management
- **Unit types**: `std::vector<LTUnitType>` with full specifications
- **Gun types**: `std::vector<LTGunType>` for weapon systems
- **Object types**: `std::vector<LTObjectType>` for terrain features
- **Dictionary lookup**: Fast name-based access via unordered maps

#### Configuration Loading
- **File-based**: Text configuration for easy modding
- **Error handling**: Graceful failure on missing files
- **Validation**: Duplicate detection and warnings

### 3. Tactical Environment

#### Spatial Organization
- **Map squares**: Grid-based object placement
- **Linear objects**: Walls, rivers, roads affecting movement
- **Static objects**: Buildings and terrain features
- **Interaction radius**: Configurable engagement distances

#### Simulation Parameters
- **Time step**: Configurable simulation granularity
- **Damping**: Physics simulation stability
- **Frame rate**: Adjustable for performance/accuracy trade-offs

### 4. Tactical Queries

#### Spatial Analysis
- **Units in area**: Circle-based unit detection
- **Line objects**: Walls and barriers in range
- **Static objects**: Terrain features affecting combat
- **Surroundings**: Comprehensive environmental assessment

#### Deployment Optimization
- **Position finding**: Optimal unit placement
- **Formation**: Squad arrangement algorithms
- **Cover evaluation**: Protection assessment
- **Visibility**: Line-of-sight calculations

## Technical (from code)
Grounded in `cpp/apps/LandTactics/LTWorld.h`:

### Data Structures
- `class LTWorld` - main world manager
- `class LTMapSquare` - grid cell for object storage
- Vectors for unit types, gun types, and objects
- Unordered maps for fast dictionary lookup

### Integration Points
- **Terrain systems**: Shared with LandCraft for consistency
- **Path finding**: Leverages existing PathFinder implementation
- **Physics**: Rigid body simulation for units
- **Rendering**: OpenGL integration for visualization

### Configuration System
- `loadUnitTypes()` - unit database loading
- `loadGunTypes()` - weapon system definitions
- `initStaticObject()` - terrain feature initialization
- `initLinearObjects()` - wall/trench system setup

**Notes:**
- Shares terrain infrastructure with LandCraft for consistency
- Provides comprehensive tactical environment
- Supports large-scale battles with hundreds of units
- Configuration-driven for easy scenario creation
