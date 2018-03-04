
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

    int getUnitsInCircle( const Vec2d& pos, double R, std::vector<LTUnit*>& out ){
        double R2 = R*R;
        for( LTSquad* s : squads ){
            if( pos.dist2(s->pos) > sq( R + s->radius) ) continue;
            for( LTUnit& u : s->units ){
                //if( d.norm2() > sq( R + u->type->radius) )
                if( pos.dist2(u.pos) < R2 ) out.push_back(&u);
            }
        }
    }

    LTFaction( char * name_, uint32_t color_ ){
        //color.set( color_ );
        color = color_;
        name = name_;
    };

};

#endif
