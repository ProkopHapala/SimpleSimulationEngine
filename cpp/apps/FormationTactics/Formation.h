

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

class Formation{
	public:
    int          id;
    char       * name;
    Faction    * faction;
    BattleLine * line;

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

    void leaveMenBehind   ( );
    bool eliminateInvalids( );
    void update_bbox      ( );
    int  interact         ( Formation * fb );
    int  interactInside   ( );
    void update           ( double dt );
    void render           ( const Vec3f& color, int view_type );

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
        return r2box > ( maxBbox2*( width*width + length*length  ) );
    }

    inline void checkTarget( ){
        Vec2d d;
        d.set_sub( p00, p00target );
        if( d.norm2() >  0.01 ) return;
        d.set_sub( p01, p01target );
        if( d.norm2() >  0.01 ) return;
        movingToTarget = false;
    }

    Formation(){};
    ~Formation(){
        if( soldiers != NULL ) delete soldiers;
    }

};

#endif
