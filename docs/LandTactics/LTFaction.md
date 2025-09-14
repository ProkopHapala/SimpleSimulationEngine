# LTFaction.h - Military Faction Management System

## Purpose
Manages military factions, their unit squads, and provides spatial queries for tactical operations and command structures.

## Core Components

### 1. Faction Organization

#### LTFaction Class
- **Purpose**: Represents a military force with distinct units and objectives
- **Key Attributes**:
  - `char* name` — faction identifier
  - `uint32_t color` — RGBA color used for rendering
  - `std::vector<LTSquad*> squads` — collection of controlled squads

#### Squad Management
- **Vector storage**: `std::vector<LTSquad*>` for efficient iteration
- **Ownership**: Factions manage squad lifecycle
- **Hierarchical structure**: Squads contain individual units

### 2. Spatial Query System

#### Unit Location Queries
- **getUnitAt(const Vec2d& p)**: Returns pointer to nearest `LTSquad` to position `p` (not individual `LTUnit`).
  - Implementation: Linear scan over `squads`, picks minimum `p.dist2(squad->pos)`.
  - Complexity: O(nSquads).
  - Purpose: Mouse selection / command of squads.

#### Area Queries
- **getUnitsInCircle(const Vec2d& pos, double R, std::vector<LTUnit*>& out)**: Appends pointers to all units within circle of radius R centered at `pos`.
  - Early rejection: skips squads if `pos.dist2(squad->pos) > sq(R + squad->radius)`.
  - Then iterates `squad->units` and checks `pos.dist2(u.pos) < R*R`.
  - Returns count of appended units.
  - Complexity: O(nSquads + totalUnitsInPassedSquads).

#### Tactical Applications
- **Selection**: Click-to-select unit mechanics
- **Targeting**: Identify units in weapon range
- **Formation**: Group selection and movement
- **Reconnaissance**: Detect enemy units in area

## Technical (from code)
Grounded in `cpp/apps/LandTactics/LTFaction.h`:

### Data Structures
- `char* name` - faction identifier
- `uint32_t color` - RGBA color for rendering
- `std::vector<LTSquad*> squads` - squad collection

### Query Methods
- `getUnitAt(const Vec2d& p)` - nearest squad selection
- `getUnitsInCircle(const Vec2d& pos, double R, std::vector<LTUnit*>& out)` - area search over squads/units

### Performance Characteristics
- **Nearest unit**: O(n) where n = total squads
- **Area search**: O(n*m) where n = squads, m = units per squad
- **Memory**: Minimal overhead with vector storage

**Notes:**
- Factions provide organizational structure for AI and player control
- Spatial queries enable intuitive unit selection and targeting
- Color system allows quick visual identification of forces
- `LTSquad` is referenced but not documented here; relies on fields `pos`, `radius`, and `units` [need-check: see `LTSquad.h` for exact API].
