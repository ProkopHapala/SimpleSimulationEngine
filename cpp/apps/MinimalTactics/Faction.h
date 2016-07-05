

#ifndef Faction_h
#define Faction_h

#include <vector>

#include "Vec2.h"
#include "Vec3.h"

#include "Unit.h"


class Faction{
	public:
    char  * name;
    Vec3f   color;

    std::vector<Unit*>  units;

    Unit* getUnitAt(const Vec2d& p ){
        int i=0;
        int imin=0;
        double r2min = 1e+300;
        for( Unit * u : units ){
            double r2 = p.dist2( u->pos );
            if( r2 < r2min ){ r2min=r2; imin=i; }
            i++;
        }
        printf( " imin %i r2min %f \n", imin, r2min );
        return units[imin];
    };

    Faction( char * name_, const Vec3f& color_ ){
        color.set( color_ );
        name = name_;
    };

};

#endif
