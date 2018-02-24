
#ifndef LTFaction_h
#define LTFaction_h

#include <vector>

#include "Vec2.h"
#include "Vec3.h"

#include "LTUnit.h"
#include "LTSquad.h"

class LTFaction{
	public:
    char  * name;
    //Vec3f   color;
    uint32_t color;

    std::vector<LTSquad*>  squads;

    LTSquad* getUnitAt(const Vec2d& p ){
        int i=0;
        int imin=0;
        double r2min = 1e+300;
        for( LTSquad * u : squads ){
            double r2 = p.dist2( u->pos );
            if( r2 < r2min ){ r2min=r2; imin=i; }
            i++;
        }
        printf( " imin %i r2min %f \n", imin, r2min );
        return squads[imin];
    };

    LTFaction( char * name_, uint32_t color_ ){
        //color.set( color_ );
        color = color_;
        name = name_;
    };

};

#endif
