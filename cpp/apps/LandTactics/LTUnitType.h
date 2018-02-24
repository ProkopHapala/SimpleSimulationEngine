
#ifndef LTUnitType_h
#define LTUnitType_h

#include <string>
#include <unordered_map>

#include <cstdlib>
#include <string.h>
#include <stdio.h>

//#include <math.h>
#include <fastmath.h>
#include "Vec2.h"
#include "Vec3.h"

#include "LTcommon.h"

static const int nGunMax = 3;

enum UnitKind { inf=1, gun=2, tank=3, stug=4, heli=5 };
enum DmgTypes { KE=1, HEAT=2, Frag=3, Fire=4 };

class LTGunType{ public:
    std::string name = "defaultGunType";
    double rps=0.1;         // [1] round per second
    int    nburst;
    //double coolDown;      // time
    double vMuzzle=1000;    // [m/s]
    double caliber=0.00762; // [m]
    double pMass   =0.001;  // [kg]  projectile mass
    double balisicCoef=caliber*caliber*0.25*M_PI*pMass;   // [m^2/kg] aerodynamic_decceleration a = balisicCoef * v^2 * airDensity
              // [kg]         projectile mass
    double AP=100;          // [mm]s
    double ExDMG=4e+6;      // [J] after penetration damage
    int    dmgType;         // e.g. KineticEnergy, Explosive, Flame ...

    //double period       = 20.0;
    //double range        = 0.0; // if this is <1.0 unit is not ranged;   We may modify this later based on physics
    //double damage       = 3.0;
    //double penetration  = 1.0;
    double spread       = 1.1; // tg(angular_spread),  spread = distance*fire_spread
    //double fStamina     = 0.8; // stamina cost of fire attack

    double getVelocityDecay(double dist)               { return exp(-balisicCoef*dist);  };
    double getSpread(double dist)                      { return sq(spread*dist); };
    double getKineticDamage(double dist,double dHeight){ return pMass * ( sq( vMuzzle * getVelocityDecay(dist) ) - dHeight*GravityAcc );   };

    void fromString( char * str_ ){
        //printf("1 \n");
        char *str    = strdup(str_);
        printf( "%s\n", str );
        char * token = strtok(str, ";"); //printf( ">>%s<<\n", token );
        name = strdup(token);
        //sscanf( token, "%[^\n]s", name );
        printf( "%s\n", name.c_str() );
        //printf( " basic \n" );
        token = strtok(NULL, ";"); //printf( ">>%s<<\n", token );
        sscanf( token, "%lf %li %lf %lf", &rps, &nburst, &spread, &vMuzzle, &pMass );
        //printf( "%lf %lf %lf %lf\n", mass, radius, max_speed, min_speed );
        //printf( " time \n" );
        token = strtok(NULL, ";"); //printf( ">>%s<<\n", token );
        sscanf( token, "%lf %lf", &AP, &ExDMG );
        //printf( "%lf %lf\n", stamina_regain, time_buffer );
        //printf( " melee \n" );

        std::string s_dmgType = token;
        dmgType = DmgTypes::KE;
        if      ( s_dmgType =="HEAT" ){ dmgType=DmgTypes::KE;   }
        else if ( s_dmgType =="Fire" ){ dmgType=DmgTypes::Fire; }
        else if ( s_dmgType =="Frag" ){ dmgType=DmgTypes::Frag; };

    };

    LTGunType(){};
    LTGunType( char * fname ){ fromString(fname); }

};

class LTUnitType{ public:

    UnitKind kind = UnitKind::inf;
    //char   * name = "default";
    std::string name = "defaultUnitType";

    //double radius_min      = 1.0, radius_max=3.0; // [m]
    //double width;
    //double length;
    //double height;
    //double Rturret;
    Vec3d sz;

    int nGun=1;
    LTGunType guns[nGunMax];

    // mobility
    double mass;
    double enginePower;
    double maxSpeed;
    //double speed           = 1.0;                 // [m/s]

    double armorFront;
    double armorSide;
    double armorBack;
    double armorTop;
    double armorBottom;

    // armor and defence
    double HP = 2.0;
    //double heal_prob        = 0.8;
    //double hit_area         = 1.0;
    double recovery_rate   = 10.0;

    void fromString( char * str_ ){
        //printf("1 \n");
        char *str    = strdup(str_);
        printf( "%s\n", str );
        char * token = strtok(str, ";"); //printf( ">>%s<<\n", token );
        name = strdup(token);

        token       = strtok(str, ";"); //printf( ">>%s<<\n", token );
        std::string s_kind = strdup(token);

        kind = UnitKind::inf;
        if      ( s_kind =="gun"  ){ kind=UnitKind::gun;  }
        else if ( s_kind =="tank" ){ kind=UnitKind::tank; }
        else if ( s_kind =="stug" ){ kind=UnitKind::stug; }
        else if ( s_kind =="heli" ){ kind=UnitKind::heli; };

        //sscanf( token, "%[^\n]s", name );
        printf( "%s\n", name.c_str() );
        //printf( " basic \n" );
        token = strtok(NULL, ";"); //printf( ">>%s<<\n", token );
        sscanf( token, "%lf %lf %lf", &sz.a, &sz.b, &sz.c );
        //printf( "%lf %lf %lf %lf\n", mass, radius, max_speed, min_speed );
        //printf( " time \n" );
        token = strtok(NULL, ";"); //printf( ">>%s<<\n", token );
        sscanf( token, "%lf %lf %lf", &mass, &maxSpeed, &enginePower );
        //printf( "%lf %lf\n", stamina_regain, time_buffer );
        //printf( " melee \n" );
        token = strtok(NULL, ";"); //printf( ">>%s<<\n", token );
        sscanf( token, "%lf %lf %lf %lf %lf", &armorFront, &armorSide, &armorBack, &armorTop, &armorBottom );
        //printf( "%lf %lf %lf %lf %lf %lf %lf %lf\n", melee_skill, melee_period, melee_range, melee_damage, melee_penetration, melee_fStamina, melee_push );
        //printf( " defence \n" );

        token = strtok(NULL, ";"); //printf( ">>%s<<\n", token );
        sscanf( token, "%lf", &HP );

        token = strtok(NULL, ";"); //printf( ">>%s<<\n", token );
        sscanf( token, "%li", &nGun );
    };

   LTUnitType(){};
   LTUnitType( char * fname ){ fromString(fname); }

 };

#endif
