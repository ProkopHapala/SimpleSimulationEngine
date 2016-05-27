

#ifndef BattleLine_h
#define BattleLine_h

#include <vector>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"
#include "Body2D.h"

#include "FormationTacticsCommon.h"
#include "Soldier.h"
#include "Formation.h"

class BattleLine{
	public:
    char    * name;
    Faction * faction;

    std::vector<Formation*> formations;

    void setTargetLine( const Vec2d& p1, const Vec2d& p2 ){
        Vec2d dp, p1i,p2i;
        dp.set_sub( p2, p1 );
        //double l = dp.norm();
        //double dl = l/
        dp.mul( 1.0/(formations.size()-1) );
        p1i.set    ( p1    );
        p2i.set_add( p1,dp );
        printf( " p1 (%3.3f,%3.3f) p2 (%3.3f,%3.3f)  dp (%3.3f,%3.3f) \n", p1.x,p1.y, p2.x,p2.y,  dp.x,dp.y );
        for( Formation* fm : formations ){
            fm->p00target.set( p1i );
            fm->p01target.set( p2i );
            printf( " (%3.3f,%3.3f) (%3.3f,%3.3f) \n", fm->p00target.x, fm->p00target.y, fm->p01target.x, fm->p01target.y  );
            p1i.set( p2i );
            p2i.add( dp  );
        }
    }

};

#endif
