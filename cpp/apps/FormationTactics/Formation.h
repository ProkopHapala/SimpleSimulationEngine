

#ifndef Formation_h
#define Formation_h

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"
#include "Body2D.h"

#include "FormationTacticsCommon.h"
#include "Soldier.h"

class Formation{
	public:
    int          id;
    char       * name;
    Faction    * faction;
    BattleLine * line;

    double maxWill    = 1.0;
    double bboxMargin = 0.0;
    double maxBbox2   = 4.0;

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
    void interact         ( Formation * fb );
    void interactInside   ( );
    void update           ( double dt );
    void render           ( );

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
