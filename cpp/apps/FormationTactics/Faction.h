

#ifndef Faction_h
#define Faction_h

#include <vector>

//#include <SDL2/SDL.h>
//#include <SDL2/SDL_opengl.h>
//#include "Draw2D.h"

//#include "fastmath.h"
//#include "Vec2.h"
#include "Vec3.h"
//#include "geom2D.h"
//#include "Body2D.h"

#include "FormationTacticsCommon.h"
#include "Soldier.h"
#include "Formation.h"
#include "BattleLine.h"

class Faction{
	public:
    char  * name;
    Vec3f color;

    std::vector<BattleLine*> battleLines;
    std::vector<Formation*>  formations;

    void initFaction( int nFormations, int nrows_, int ncols_, std::vector<SoldierType>& soldierTypes, const Vec2d& p1, const Vec2d& p2, double width ){
        BattleLine* battleLine = new BattleLine();
        battleLine->formations.reserve( nFormations );
        formations.reserve( nFormations );
        battleLines.push_back( battleLine );
        for( int i=0; i<nFormations; i++ ){
            Formation * formation = new Formation( nrows_, ncols_, &(soldierTypes[0]), this );
            formation->width = width;
            formations .push_back( formation );
            battleLine->formations.push_back( formation );

        }
        battleLines[0]->setTargetLine( p1, p2 );
        int i=0;
        for( Formation * fm : formations ){
            fm->jumpToTarget();
            //printf( " \n", i,  );
            fm->deploySoldiers();
            i++;
        }
    }

    Faction( char * name_, const Vec3f& color_ ){
        color.set( color_ );
        name = name_;
    };

};

#endif
