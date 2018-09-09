
#ifndef LTUnitType_h
#define LTUnitType_h

#include <string>
#include <vector>
#include <unordered_map>

#include <cstdlib>
#include <string.h>
#include <stdio.h>

//#include <math.h>
#include <fastmath.h>
#include "Vec2.h"
#include "Vec3.h"


#include "IO_utils.h"
#include "LTcommon.h"

static const int nGunMax = 3;

enum UnitKind { inf, gun, tank, stug, heli };
enum DmgTypes { KE, HEAT, Frag, Fire };

static const char *sUnitKind[5] = { "inf", "gun", "tank", "stug", "heli" };
static const char *sDmgTypes[4] = { "KE","HEAT","Frag","Fire" };

class LTGunType;
class LTUnitType;

typedef std::unordered_map<std::string,LTGunType*>  GunTypeDict;
typedef std::unordered_map<std::string,LTUnitType*> UnitTypeDict;

/*
template <typename T>
int load2vector( char* fname, std::vector<T>& vec ){
    FILE * pFile;
    const int nbuff = 4096;
    char str[nbuff];
    pFile = fopen ( fname , "r");
    if (pFile == NULL){ printf("file not found: %s \n", fname ); return(-1); };
    while ( fgets( str , nbuff, pFile) != NULL ){
        printf("===str: >>%s<<\n", str);
        if (str[0]=='#') continue;
        T t = T(str);
        vec.push_back( t );
    }
    fclose(pFile);
    return vec.size();
}

template <typename T,typename T2>
int load2vectorDict( char* fname, std::vector<T>& vec, T2& dct ){
    FILE * pFile;
    const int nbuff = 4096;
    char str[nbuff];
    pFile = fopen ( fname , "r");
    if (pFile == NULL){ printf("file not found: %s \n", fname ); return(-1); };
    while ( fgets( str , nbuff, pFile) != NULL ){
        printf("===str: >>%s<<\n", str);
        if (str[0]=='#') continue;
        T t = T(str,dct);
        vec.push_back( t );
    }
    fclose(pFile);
    return vec.size();
}
*/

template <typename T,typename K>
int vec2map( std::vector<T>& vec, std::unordered_map<K,T*>& dct ){
    int ndup;
    printf( "vec2map: vec.size() %i \n", vec.size() );
    for( int i=0; i<vec.size(); i++ ){
        //printf(" vec2map vec[%i] \n", i );
        const T& t = vec[i];
        auto it = dct.find( t.name );
        if( it != dct.end() ){
            ndup++;
            printf( "WARRNING: ignoring duplicate: vec[%i].name=>>%s<< \n", i, t.name.c_str() );
        }else{
            //printf( "inserting >>%s<< \n", t.name.c_str() );
            //vec.push_back( t );
            dct.insert( (std::pair<K,T*>){ t.name, &vec[i] } );
        };
    }
    return ndup;
}

template <typename T,typename K>
void printDict( std::unordered_map<K,T*>& dct ){
    printf( "printDict dct.size() %i \n", dct.size() );
    for( std::pair<K,T*> p : dct ){
        printf( ">>%s<< >>%s<<\n", p.first.c_str(), p.second->name.c_str() );
    }
}

class LTGunType{ public:
    std::string name = "defaultGunType";
    double rps    = 0.1;           // [1] round per second
    int    nburst = 1;
    //double coolDown;             // time
    double vMuzzle     = 1000;     // [m/s]
    double caliber     = 0.00762;  // [m]
    double pMass       = 0.001;    // [kg]  projectile mass
    double crossArea   = caliber*caliber*0.25*M_PI;
    double balisicCoef = crossArea/pMass;   // [m^2/kg] aerodynamic_decceleration a = balisicCoef * v^2 * airDensity
              // [kg]         projectile mass
    double AP      = 100;   // [mm]s
    double ExDMG   = 4e+6;  // [J] after penetration damage
    int    dmgType = 0;     // e.g. KineticEnergy, Explosive, Flame ...

    //double period       = 20.0;
    //double range        = 0.0; // if this is <1.0 unit is not ranged;   We may modify this later based on physics
    //double damage       = 3.0;
    //double penetration  = 1.0;
    double spread       = 1.1; // tg(angular_spread),  spread = distance*fire_spread
    //double fStamina     = 0.8; // stamina cost of fire attack

    void fromString(const char * str_ );
    char* toStrCaptioned( char * sout );

    void updateAux();

    inline double getVelocityDecay(double dist)               const { return exp(-balisicCoef*dist);  };
    inline double getSpread(double dist)                      const { return sq(spread*dist); };
    inline double getKineticDamage(double dist,double dHeight)const { return pMass * ( sq( vMuzzle * getVelocityDecay(dist) ) - dHeight*GravityAcc );   };

    LTGunType(){};
    LTGunType( char * fname ){ fromString(fname); }

};

class LTUnitType{ public:

    int kind = UnitKind::inf;
    //char   * name = "default";
    std::string name = "defaultUnitType";

    //double radius_min      = 1.0, radius_max=3.0; // [m]
    //double width;
    //double length;
    //double height;
    //double Rturret;
    Vec3d sz;
    Vec3d szAreas;

    int nGun=1;
    LTGunType* guns[nGunMax];

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

    void  updateAux();
    void  fromString(const char * str_, GunTypeDict& gunTypeDict );
    char* toStrCaptioned( char * sout, bool bGunDetials );

    LTUnitType(){};
    LTUnitType( char * fname, GunTypeDict& gunTypeDict ){ fromString(fname, gunTypeDict); }

 };

#endif
