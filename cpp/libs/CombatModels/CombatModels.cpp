
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


//#include "AirCombatModel.h"
//#include "NavalCombatModel.h"
//#include "SimpleBattleLine.h"
//#include "NavalCombatModel.h"

#include "MinimalDivisionLevel.h"

//std::unordered_map<std::string,WeaponType*> weaponTypes;
//std::unordered_map<std::string,UnitType*>   unitTypes;
//std::unordered_map<std::string,UnitState*>  units;


// Compare with units from:
// /apps/LandTactics/LTUnitType.h
// ultimately, there shoould be one unified system of units

std::vector<WeaponType*> weaponTypes;
std::vector<UnitType*>   unitTypes;
std::vector<UnitState*>  units;
Combat combat1;

// ============ Global Variables

extern "C"{
// ========= Grid initialization

int newWeaponType( const char* str ){
    WeaponType* o = new WeaponType();
    o->fromString(str);
    weaponTypes.push_back(o);
    return weaponTypes.size()-1;
}

int newUnitType( const char* str, int iprim, int isec ){
    UnitType* o = new UnitType();
    o->fromString(str);
    if(iprim>=0){ o->primary   = weaponTypes[iprim]; }else{ o->primary  =0; }
    if(isec >=0){ o->secondary = weaponTypes[isec ]; }else{ o->secondary=0; }
    unitTypes.push_back(o);
    return unitTypes.size()-1;
}

int newUnit( int ityp, int n ){
    UnitState* o = new UnitState();
    o->init( unitTypes[ityp] , n );
    units.push_back(o);
    return units.size()-1;
}

void prepareCombat( int natt, int ndef, int* atts, int* defs ){
    combat1.clear();
    //char stmp[8000];
    for(int i=0; i<natt; i++){
        int iu=atts[i];
        UnitState* u = units[ iu ];
        combat1.attacker.composition.units.push_back( u );
        //unit_tank.type->info(stmp,true); printf( "attacker  unit_tank  : \n %s\n", stmp );
    }
    for(int i=0; i<ndef; i++){
        int iu=atts[i];
        UnitState* u = units[ iu ];
        combat1.defender.composition.units.push_back( u );
        //unit_tank.type->info(stmp,true); printf( "attacker  unit_tank  : \n %s\n", stmp );
    }
}

void solveCombat( int nrounds, int ndist, double* dists, double dt_max, double advance_dist, double maxCammo ){
    iDEBUG = 1;
    combat1.conds.maxCamo=maxCammo;
    for(int i=0; i<ndist; i++){
        combat1.attacker.composition.reset();
        combat1.defender.composition.reset();
        combat1.start( 20.0 );   // dist
        for(int i=0; i<nrounds; i++){
            combat1.round( dt_max, advance_dist  );  // dt_max, advance_dist
        }
    }
}


} // extern "C"
