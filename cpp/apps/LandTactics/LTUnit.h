
#ifndef LTUnit_h
#define LTUnit_h

#include "fastmath.h"
#include "tresholdFunctions.h"
#include "Vec2.h"
#include "Vec3.h"
#include "geom2D.h"
#include "Body2D.h"

#include "LTcommon.h"
#include "LTUnitType.h"
//#include "LTFaction.h"

class LTGun{ public:
    LTGunType * type = 0;
    int n=1;      // twin,quad ...
    int ammo=0;
    double timer   = 0;
    int attachedTo = 0; // attached to turret, hull or what ?
};

class LTUnit : public RigidBody2D { public:
    //UnitKind kind=UnitKind::inf;
    LTUnitType* type = NULL;
    LTGun guns[nGunMax]; // by default all guns attached to turrer

    bool alive=true;
    int  wound=0;

    double suppressed;

    Vec2d turret_dir = (Vec2d){1.0,0.0};

    // ===== inline functions

    //void move_to_goal ( double dt );
    //void fire_at_unit ( LTUnit * target );
    //void update       ( double dt );
    //double damage_ramp( double att, double def );
    //void setOpponent  ( LTUnit * opponent_ );
    //void setGoal      ( const Vec2d& goal_ );
    void setType( LTUnitType* type_ );
    void update ( double dt );
    void render ( uint32_t color, int iLOD);

    void getShot( const Vec3d& from, int nburst, double area, double damage_const, double damage_kinetic, double penetration );

    //void renderJob    ( uint32_t c );
    //void view();

    LTUnit(LTUnitType* type_){ type = type_; };

};

#endif
