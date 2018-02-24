
#ifndef MinimalTacticsCommon_h
#define MinimalTacticsCommon_h

#include "fastmath.h"
#include "tresholdFunctions.h"

extern int   default_font_texture;

static const double GravityAcc = 9.80665;

class Unit;
class Faction;

inline double damage_ramp( double att, double def ){ return 0.5*(1+Treshold::r2( (att-def)/(att+def) )); }

#endif
