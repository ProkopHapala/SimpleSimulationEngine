
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <vector>
#include <unordered_map>
#include <string>

#include "fastmath.h"
#include "Vec3.h"
#include "Vec3.h"
#include "Mat3.h"
#include "VecN.h"

#include "testUtils.h"

#include "LTUnitType.h"
#include "LTUnit.h"

std::vector<LTGunType*>   gunTypes;
std::vector<LTUnitType*>  unitTypes;
std::vector<LTUnit*>      units;

GunTypeDict  gunTypeDict;
UnitTypeDict unitTypeDict;

// ============ Global Variables

extern "C"{
// ========= Grid initialization

int newWeaponType( const char* str ){
    LTGunType* o = new LTGunType(str);
    o->fromString(str);
    gunTypes.push_back(o);
    gunTypeDict[ o->name ] = o;
    return gunTypes.size()-1;
}

//int newUnitType( const char* str, int ngun, int* gunTs ){
int newUnitType( const char* str ){
    LTUnitType* o = new LTUnitType( str, gunTypeDict );
    //for(int i=i; i<ngun;  i++){ o->guns[i] = gunTypes[ gunTs[i]]; }
    unitTypes.push_back(o);
    unitTypeDict[ o->name ] = o;
    return unitTypes.size()-1;
}

int newUnit( int ityp ){
    LTUnit* o = new LTUnit( unitTypes[ityp] );
    units.push_back(o);
    return units.size()-1;
}

void setUnitPose( int i, double* pos, double* rot ){
    LTUnit* o = units[i];
    o->pos = *(Vec2d*)pos; 
    o->rot = *(Vec2d*)rot;
}

double shoot( int itarget, int iattacker, int igun ){
    LTUnit* attacker = units[iattacker];
    LTUnit* target   = units[itarget  ];
    return attacker->fireGun( igun, *target );
}

} // extern "C"
