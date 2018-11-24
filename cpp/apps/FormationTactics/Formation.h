

#ifndef Formation_h
#define Formation_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "geom2D.h"
#include "Body2D.h"

#include "FormationTacticsCommon.h"
#include "Soldier.h"

const int VIEW_INJURY  = 1;
const int VIEW_STAMINA = 2;
const int VIEW_CHARGE  = 3;
const int VIEW_MORAL   = 4;

class Formation{ public:
    static constexpr double leaveBehindMargin = 20.0;
    static constexpr double maxDistFromCog    = 5.0;
    static constexpr double KDistFromCog      = 10.0;

    int          id;
    //char       * name;
    std::string  name;
    Faction    * faction;
    BattleLine * line;

    Vec2d fire_target;

    double maxWill    = 1.0;
    double bboxMargin = 0.0;
    double maxBbox2   = 4.0;

    bool   tAttack = 3.0;
    bool   melee   = true;
    double order   = 1.0;
    double moral   = 1.0;
    double stamina = 1.0; // stamina summary over all soldier
    bool shouldLeaveMenBehind = true;

    Vec2d center;
    Vec2d p00target,p01target;
    bool movingToTarget = false;
    //Vec2d pL,pR;
    Vec2d p00,p01,p10,p11;
    Vec2d dirFw,dirLf;
    double length,width;
    double kLength=3.0, kWidth=3.0;

    Rect2d bbox;
    Vec2d  cog;
    //Vec2d bbox_min, bbox_max;

    int nrows;
    int ncols;
    int nSoldiers;
    int nCapable,nAlive;
    Soldier * soldiers = NULL;

    // =========== function declaration

    void setTarget     ( const Vec2d& target );
    void setEnds       ( const Vec2d& pL_, const Vec2d& pR_, double width_ );
    void setCenterRot  ( const Vec2d& center_, const Vec2d& dirFw_ );
    void moveToTarget  ( );
    void moveBy        ( const Vec2d& dpos );
    void jumpToTarget  ( );

    void clean_temp    ( );
    void applyWillForce( Soldier& soldier );
    void applyWillForce( );

    int  leaveMenBehind   ( );
    bool eliminateInvalids( );
    void update_bbox      ( );
    int  interact         ( Formation * fb );
    int  interactInside   ( );
    void update           ( double dt );
    void render           ( const Vec3f& color, int view_type );

    char * reportSetup ( char * sout );
    char * reportStatus( char * sout );

    void setupSoldiers( SoldierType * type );
    void deploySoldiers( );

    Formation( int id_, int nrows_, int ncols_, SoldierType * type, Faction * faction_ );

    // =========== inline functions

    inline double willSaturation( double will ){
        if( will > maxWill ) return maxWill;
        return will;
    }

    inline bool checkMenBehind( ){
        double r2box = bbox.a.dist2( bbox.b );
        return r2box > ( maxBbox2*( width*width + length*length  ) ); // maxBbox2 is dimension less
    }

    inline void checkTarget( ){
        Vec2d d;
        d.set_sub( p00, p00target ); if( d.norm2() >  0.01 ) return;
        d.set_sub( p01, p01target ); if( d.norm2() >  0.01 ) return;
        movingToTarget = false;
    }

    Formation(){};
    ~Formation(){
        if( soldiers != NULL ) delete soldiers;
    }

};

struct SalvoTarget{
    double t;
    //int type;
    void * ptr;
};

struct SalvoBin{
    int nprj;
    double t0;
};

class Salvo{ public:
    Vec2d fw;
    Vec2d left;
    double xmin,xmax;
    double t_span = 2.0;
    //Vec2d p0,p1;
    double tgElev; // tangent of elevation
    double dbin;
    int    nbin=0, perBin=3;
    //double t;
    int nprj;

    SalvoBin*    bins;
    SalvoTarget* targets;
    // may have projectiles in bins

    // function declaration

    int  build        ( Formation& f, double d, char* buff );
    void clearTargets ( Formation& f );
    void gatherShots  ( Formation& f );
    void gatherTargets( Formation& f );

};

#endif
