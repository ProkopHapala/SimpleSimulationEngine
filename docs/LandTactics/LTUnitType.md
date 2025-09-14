# LTUnitType.h - Military Unit Type Definitions

## Purpose
Defines comprehensive unit types for military vehicles and soldiers, including weapon systems, armor profiles, and mobility characteristics for realistic tactical simulation.

## Core Components

### 1. Weapon System Definitions

#### LTGunType Class
- **Purpose**: Defines individual weapon characteristics
- **Key Attributes**:
  - `rps`: Rate of fire (rounds per second)
  - `vMuzzle`: Muzzle velocity (1000 m/s default)
  - `caliber`: Barrel diameter (0.00762 m default)
  - `pMass`: Projectile mass (0.001 kg default)
  - `AP`: Armor penetration capability (100 mm default)
  - `ExDMG`: Explosive damage after penetration (4e+6 J default)
  - `dmgType`: Damage type (KE, HEAT, Frag, Fire)
  - `crossArea = caliber^2 * 0.25 * PI`
  - `balisicCoef = crossArea / pMass` (aerodynamic deceleration coefficient)
  - `spread`: stored as tangent-like factor; accuracy uses `sq(spread*dist)`
  - Methods: `fromString(const char*)`, `toStrCaptioned(char*)`, `updateAux()`

#### Ballistic Calculations
- **Velocity decay**: `exp(-balisicCoef*distance)`
- **Spread**: Quadratic increase with distance
- **Kinetic damage**: Based on remaining velocity and height difference
- **Ballistic coefficient**: `crossArea/pMass` for aerodynamic drag
 - Inline helpers (from header):
   - `getVelocityDecay(dist)`
   - `getSpread(dist)` returns `sq(spread*dist)`
   - `getKineticDamage(dist,dHeight)` returns `pMass*( v^2 - dHeight*GravityAcc )`

### 2. Unit Type Classification

#### UnitKind Enumeration
- **inf**: Infantry units
- **gun**: Artillery/anti-tank guns
- **tank**: Main battle tanks
- **stug**: Assault guns/tank destroyers
- **heli**: Helicopter units
 - String tables: `sUnitKind[5] = {"inf","gun","tank","stug","heli"}`
 - Damage types: `sDmgTypes[4] = {"KE","HEAT","Frag","Fire"}`

#### Armor Profiles
- **Directional armor**: Separate values for front, side, back, top, bottom
- **Realistic thickness**: Modeled in millimeters
- **Coverage**: Accounts for weak spots and angles
 - Fields: `armorFront, armorSide, armorBack, armorTop, armorBottom`

### 3. Mobility System

#### Performance Metrics
- **mass**: Total unit weight
- **enginePower**: Available propulsion power
- **maxSpeed**: Top speed on ideal terrain
- **Recovery rate**: How quickly units regain combat effectiveness
 - Size: `Vec3d sz` and `Vec3d szAreas` for dimensions/target areas [need-check: exact semantics]
 - Guns: `int nGun`, `LTGunType* guns[nGunMax]` (max 3)

### 4. Configuration System

#### File Format
```
UnitName UnitKind mass enginePower maxSpeed
armorFront armorSide armorBack armorTop armorBottom HP
nGuns gun1 gun2 ...
```

#### Loading System
- **Dictionary-based**: `UnitTypeDict` and `GunTypeDict` for fast lookup
- **Duplicate detection**: Warns about duplicate unit names
- **Cross-references**: Guns reference gun types by name
 - Methods: `LTUnitType::fromString(const char*, GunTypeDict&)`; `LTGunType::fromString(const char*)`
 - `vec2map()` helper populates maps from vectors; logs duplicates

## Technical (from code)
Grounded in `cpp/apps/LandTactics/LTUnitType.h`:

### Data Structures
- `enum UnitKind { inf, gun, tank, stug, heli }` - unit classification
- `enum DmgTypes { KE, HEAT, Frag, Fire }` - damage type system
- `class LTGunType` with ballistic physics calculations
- `class LTUnitType` with comprehensive unit attributes
 - `nGunMax = 3` (compile-time maximum)

### Configuration Loading
- `fromString()` methods for parsing text configuration
- `vec2map()` for building lookup dictionaries
- `printDict()` for debugging unit databases
 - `LTUnitType( const char* fname, GunTypeDict& )` convenience constructor

### Constants
- `static const int nGunMax = 3` - maximum weapons per unit
- String arrays for unit kind and damage type names

**Notes:**
- Configuration files use simple text format for easy modding
- Ballistic calculations include realistic physics for range and penetration
- Unit types support inheritance through composition rather than class hierarchy
 - Some fields like `szAreas` are used in hit/visibility computation [need-check: confirm in source]
