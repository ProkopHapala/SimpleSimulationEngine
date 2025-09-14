# LTShelter.h - Cover and Terrain Object System

## Purpose
Provides tactical cover system for infantry and vehicles, including static terrain objects like buildings, trees, and fortifications that affect visibility and protection.

## Core Components

### 1. Static Object System

#### LTObjectType Base Class
- **Purpose**: Defines reusable static object templates
- **Key Features**:
  - `kind`: Object classification (tree, house, rocks)
  - `glo`: OpenGL display list for efficient rendering
  - Virtual render interface for extensibility

#### Object Types
- **LTRectHouseType**: Rectangular buildings with customizable dimensions
- **Rendering**: 2D top-down representation with orientation
- **Collision**: Precise geometric boundaries

### 2. Cover Mechanics

#### Linear Objects
- **LTLinearObject**: Walls, trenches, hedgerows
- **Cover value**: 0.0-1.0 scale for protection effectiveness
- **Width**: Physical thickness for collision detection
- **Intersection testing**: Line-of-sight blocking calculations

#### Shelter System
- **Protection zones**: Areas shielded from specific directions
- **Armor values**: Structural integrity for different materials
- **Entry points**: Units can occupy/leave shelters

### 3. Tactical Integration

#### Visibility System
- **Line-of-sight blocking**: Objects prevent target acquisition
- **Cover calculation**: Reduces incoming damage based on angle
- **Height consideration**: Multi-level terrain interaction

#### Movement Effects
- **Obstacle avoidance**: Units navigate around objects
- **Pathfinding**: Static objects as navigation nodes
- **Formation disruption**: Objects break unit formations

## Technical (from code)
Grounded in `cpp/apps/LandTactics/LTShelter.h`:

### Object Hierarchy
- `class LTObjectType` - base object type definition
- `class LTRectHouseType` - specific building implementation
- `class LTStaticObject` - individual object instances
- `class LTLinearObject` - line-based cover objects

### Geometric Operations
- `dp(const Vec2d& p)` - distance derivative for collision
- `intersection()` - line segment intersection testing
- `render()` - OpenGL display list generation

### Rendering System
- **Display lists**: Pre-compiled geometry for performance
- **Shape rendering**: Efficient 2D representation
- **Color coding**: Visual distinction of object types

**Notes:**
- Shelter system is WIP with TODO comments for future expansion
- Objects affect both visibility and movement
- Current implementation focuses on 2D top-down representation
