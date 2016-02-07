

#ifndef Soldier_h
#define Soldier_h

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"
#include "Body2D.h"

#include "SoldierType.h"

class Soldier : public PointBody2D{
	public:
    bool alive   = true;
    bool capable = true;

    Vec2d willForce;

    SoldierType * type = NULL;
    double maxwf = 1.0;

    void clean_temp( ){ force.set(0.0); willForce.set(0.0); }

    void moveSoldier( double dt ){
        /*
        double  rwf2 = willForce.norm2();
        if( rwf2 > (maxwf*maxwf) ){
           willForce.mul( maxwf / sqrt( rwf2 ) );
        }
        */
        force.add( willForce );
        move_PointBody2D( dt );
    }

};

#endif
