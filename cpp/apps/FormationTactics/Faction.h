

#ifndef Faction_h
#define Faction_h

#include <vector>
#include <unordered_map>

//#include <SDL2/SDL.h>
//#include <SDL2/SDL_opengl.h>
//#include "Draw2D.h"

//#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
//#include "geom2D.h"
//#include "Body2D.h"

#include "FormationTacticsCommon.h"
#include "Soldier.h"
#include "Formation.h"
#include "BattleLine.h"

//#include "PolyLineFormation.h"

class Faction{
	public:
    char  * name;
    Vec3f   color;

    std::vector<BattleLine*> battleLines;
    std::vector<Formation*>  formations;

    //std::vector<PolyLineFormation*> polylines;

    Formation* getFormationAt(const Vec2d& p ){
        int i=0;
        int imin=0;
        double r2min = 1e+300;
        for(Formation * fm : formations ){
            double r2 = p.dist2( fm->center );
            if( r2 < r2min ){ r2min=r2; imin=i; }
            i++;
        }
        printf( " imin %i r2min %f \n", imin, r2min );
        return formations[imin];
    };

    void initFaction( int nFormations, int nrows_, int ncols_, std::vector<SoldierType>& soldierTypes, const Vec2d& p1, const Vec2d& p2, double width ){
        BattleLine* battleLine = new BattleLine();
        battleLine->formations.reserve( nFormations );
        formations.reserve( nFormations );
        battleLines.push_back( battleLine );
        for( int i=0; i<nFormations; i++ ){
            //char * name; asprintf(&name,"%3iof%s", i, name );
            Formation * formation = new Formation( i, nrows_, ncols_, &(soldierTypes[0]), this );
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

    int initFaction( char * fname, Vec2d pos0, Vec2d dir0, const SoldierTypeDict& name2soldierType ){
        printf("Faction::initFaction(%s)\n", fname );
        BattleLine* battleLine = new BattleLine();
        battleLine->formations.reserve( 20 );
        //printf( "bl.fms.size %i \n", battleLine->formations.size() );
        //formations.reserve( nFormations );
        battleLines.push_back( battleLine );

        FILE * pFile;
        const int nbuff = 1024;
        char str[nbuff];
        pFile = fopen ( fname , "r");
        if (pFile == NULL){ printf("ERROR in Faction::initFaction(%s) file not found \n", fname ); exit(0); };
        //fgets ( str, nbuff, pFile);
        //printf(">>%s<<\n", str);

        int i=0;
        while ( fgets ( str , nbuff, pFile) != NULL ){
            if (str[0]=='#') continue;
            //printf( "2 bl.fms.size %i \n", battleLine->formations.size() );
            //printf(">>%s<<\n", str);
            //Vec2d pos, dir;

            int nrow, ncol;
            char* token;
            std::string name = strtok( str, ";");
            //printf(">>%s<<\n", name.c_str() );
            //printf( "3 bl.fms.size %i \n", battleLine->formations.size() );
            SoldierType * typ;// =  &(soldierTypes[0]);
            auto found = name2soldierType.find( name );
            if( found->second ){ typ = found->second; }else{ printf("cannot found SoldierType : %s\n", name ); exit(0); };
            //typ->toStrCaptioned(str); puts(str);

            token = strtok( NULL, ";");
            //printf(">>%s<<\n", token);
            //Vec2d p1, p2;
            //sscanf( token, "%i %i  %lf %lf  %lf %lf", &nrow, &ncol, &p1.x, &p1.y, &p2.x, &p2.y );
            //printf(" %i %i  %lf %lf  %lf %lf \n",  nrow, ncol, p1.x, p1.y, p2.x, p2.y );
            //p1.mul_cmplx(dir0); p1.add(pos0);
            //p2.mul_cmplx(dir0); p2.add(pos0);
            Vec2d p,d;
            sscanf( token, "%i %i  %lf %lf  %lf %lf", &nrow, &ncol, &p.x, &p.y, &d.x, &d.y );
            printf(" %i %i  %lf %lf  %lf %lf \n",  nrow, ncol, p.x, p.y, d.x, d.y );
            d.mul_cmplx(dir0); p.add(pos0);

            Formation * fm = new Formation( i, nrow, ncol, typ, this );
            formations .push_back( fm );
            battleLine->formations.push_back( fm );
            //char srtBuf[256]; sprintf(srtBuf, "\n", , );
            fm->name = typ->name + std::to_string( i );
            //printf( " fms.size() %i %i\n",battleLine->formations.size(), formations.size());

            fm->width  = (nrow+1)*typ->radius;
            fm->length =  ncol*typ->radius;
            //printf("DEBUG 1\n");
            fm->setCenterRot( p, d ); //printf("DEBUG 2\n");
            //fm->p00target.set( p1 );
            //fm->p01target.set( p2 );
            //fm->setTarget( pos );
            //fm->jumpToTarget();     //printf("DEBUG 3\n");
            fm->deploySoldiers();     //printf("DEBUG 4\n");
            //printf("DEBUG 4\n");
            i++;
        }
        fclose(pFile);
        //exit(0);

        return i;
    }

    Faction( char * name_, const Vec3f& color_ ){
        color.set( color_ );
        name = name_;
    };

};

#endif
