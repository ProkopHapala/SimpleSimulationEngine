# LTcommon.h - Common Tactical Definitions

## Purpose
Provides fundamental constants and utility functions shared across the LandTactics combat simulation, including damage calculation formulas and physical constants.

## Core Components

### 1. Physical Constants
- `static const double GravityAcc = 9.80665;` — standard gravitational acceleration (m/s²)
- `extern int default_font_texture;` — global font texture handle used by UI/plots

### 2. Damage Calculation System
- Signature: `inline double damage_ramp(double att, double def);`
- Implementation: `0.5*(1 + Treshold::r2( (att-def)/(att+def) ))`
- Purpose: Smooth damage multiplier from 0..1 based on attack vs defense ratio
- Behavior:
  - Balanced (att≈def) → ~0.5
  - Overmatch (att>>def) → →1.0
  - Undermatch (att<<def) → →0.0

## Technical (from code)
Grounded in `cpp/apps/LandTactics/LTcommon.h`:

- Defines `extern int default_font_texture` for global font access
- `static const double GravityAcc = 9.80665` - gravitational constant
- `inline double damage_ramp(double att, double def)` uses Treshold::r2 for smooth damage scaling

**Notes:**
- Damage ramp provides non-linear scaling for realistic combat outcomes
- Uses threshold functions for smooth transitions rather than hard breakpoints
- Actual use sites in combat code are in `LTUnit::getShot()` [need-check: confirm invocation path in .cpp]
