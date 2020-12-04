
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


#include "spaceCombat.h"

CombatAssembly battle;

// ============ Global Variables

extern "C"{
// ========= Grid initialization

    /*
    //    "%i %lf %lf %lf", &n, &layerDens, &spacing, &critEdens
    whippleShieldType ws1; ws1.fromString("2 10.0 0.5 150000");
    //                                  &area, &HPs, &HPexponent
    ProjectedTarget tg1;   tg1.fromString( "10.0  1.5e+8 2.5"); tg1.wshield = &ws1;

    // sscanf( s,              "%lf %lf"  , &mass, &caliber );
    ProjectileType pt1;    pt1.fromString( "0.15 0.12" );
    //sscanf( s, "%lf %lf %lf %lf %lf"  , &length, &maxForce, &maxPower, &scatter, &fireRate  );
    SpaceGunType   sg1;    sg1.fromString( "800 60000 1e+9 2e-4 10" );
    SpaceGun g1( 1, &sg1, &pt1 );

    DEBUG
    CombatAssembly battle;      DEBUG
    battle.addTarget( tg1 );    DEBUG
    battle.fireGun  ( g1  );    DEBUG
    //      dist[m]  accel[m/s^2]
    battle.colide( 1e+5,    0.1 );
    */

    void addTarget( const char* str_target, const char* str_shield ){
        whippleShieldType* ws = new whippleShieldType(); 
        ProjectedTarget*   tg = new ProjectedTarget();   
        tg->wshield = ws;
        ws->fromString( str_shield );
        tg->fromString( str_target );
        battle.addTarget(*tg);
    }

    void addGun( int n, const char* str_gun, const char* str_shot ){
        ProjectileType* pt = new ProjectileType();  pt->fromString( str_shot );
        SpaceGunType*   sg = new SpaceGunType();    sg->fromString( str_gun  ); 
        SpaceGun*       g  = new SpaceGun( n, sg, pt );
        battle.fireGun(*g);
    }

    void evaluateCombat( int n, double* dists,  double* accels, double* out ){
        int ntarget = battle.targets.size();
        for(int i=0; i<n; i++){
            battle.reset();
            battle.colide( dists[i],  accels[i] );
            for( int j=0; j<ntarget; j++ ){
                out[i*ntarget+j] = battle.targets[j].health;
            }
        }
    }

} // extern "C"
