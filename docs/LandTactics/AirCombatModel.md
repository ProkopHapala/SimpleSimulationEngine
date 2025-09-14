# AirCombatModel.h - Aerial Combat Physics System

## Purpose
Provides realistic physics simulation for aerial combat, including aircraft performance, aerodynamics, and tactical maneuvering for air-to-air engagements.

## Core Components

### 1. Aircraft Performance Modeling

#### CombatAirCraft Structure
- **Physical parameters**:
  - `mass`: Total aircraft weight (kg)
  - `power`: Engine power output (W)
  - `area_wing`: Wing surface area (m²)
  - `area_hull`: Fuselage drag area (m²)
  - `aspect`: Wing aspect ratio for lift efficiency
  - `CD0`: Minimum drag coefficient
  - `CLmax`: Maximum lift coefficient before stall

#### Performance Calculations
- **Maximum speed**: Based on power-to-drag ratio
- **Turn rate**: Optimal turning performance at various speeds
- **Climb rate**: Vertical performance under different conditions
- **Energy management**: Kinetic vs potential energy trade-offs

### 2. Aerodynamic Physics

#### Force Calculations
- **Drag**: `Fd = Cd * (ρ/2) * v²` - air resistance
- **Lift**: `Fl = Cl * (ρ/2) * v²` - wing lift generation
- **Gravity**: `Fg = m * g` - constant downward force
- **Thrust**: Engine power converted to propulsive force

#### Performance Envelope
- **Speed optimization**: Maximum level flight speed
- **Turn optimization**: Best sustained turn rate
- **Climb optimization**: Maximum climb rate at stall boundary
- **Energy conservation**: Minimum energy loss maneuvers

### 3. Tactical Maneuvering

#### Turn Performance Analysis
- **Sustained turns**: Maximum turn rate without energy loss
- **Instantaneous turns**: Maximum turn rate at cost of energy
- **Turn radius**: Minimum turning radius at given speed
- **Energy bleeding**: Rate of energy loss in high-G maneuvers

#### Combat Maneuvers
- **Slip away**: Fast vs slow aircraft tactical analysis
- **Energy fighting**: Managing kinetic vs potential energy
- **Turn rate comparison**: Relative maneuvering advantage
- **Range optimization**: Optimal engagement distances

### 4. Physics Integration

#### Flight Dynamics
- **Climb curves**: Performance at various climb angles
- **Turn curves**: Performance envelope visualization
- **Energy management**: Balancing speed, altitude, and maneuverability
- **Propeller efficiency**: Power conversion to thrust

#### Environmental Factors
- **Air density**: Altitude effects on performance
- **Temperature**: Atmospheric condition impacts
- **Wind**: External force considerations

## Technical (from code)
Grounded in `cpp/common/CombatModels/AirCombatModel.h`:

### Physics Constants
- `const_airDensGround = 1.22` kg/m³ - sea level air density
- `const_GravAccel = 9.81` m/s² - gravitational acceleration

### Performance Methods
- `maxSpeed_simple()` - theoretical maximum speed
- `trunRate_simple()` - optimal sustained turn rate
- `climbRate_simple()` - maximum climb rate
- `climbRate_CLmax()` - climb rate at stall boundary

### Tactical Analysis
- **Slip away equations**: Mathematical analysis of turning combat
- **Energy calculations**: Kinetic and potential energy management
- **Thrust optimization**: Propeller efficiency vs speed

**Notes:**
- Models based on real aerodynamic principles
- Supports both propeller and jet aircraft
- Includes stall and energy limitations
- Provides framework for tactical AI decision making
