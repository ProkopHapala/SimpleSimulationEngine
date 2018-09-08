
#ifndef LTSurroundings_h
#define LTSurroundings_h

#include <vector>

#include "LTUnit.h"
#include "LTShelter.h"

class LTsurrounding{ public:
    // ## Considerations:
    // - maybe we should rather use unordered_set to efficiently check if one object is inserted multiple times ?
    std::vector<LTStaticObject*> objs;
    std::vector<LTLinearObject*> lobjs;
    std::vector<LTUnit*>         enemies;
    std::vector<LTUnit*>         coleagues;

    bool   bConstr=false;
    Vec2d  ConstrPos;
    double ConstrRad =  1.0;
    double ConstrE   = -1.0;

    double RSearchLocal  = 5.0;
    double RSearchGlobal = 50.0;

    double Rdeploy = 1.0; // desired average deployment distanct

    double cover2E   =  1.0;
    double overlap2E = -1.0;

    void clear(){ objs.clear(); lobjs.clear(); enemies.clear(); coleagues.clear(); }

    double unitPosFittness( const LTUnit* u, const Vec2d& p )const{
        //printf(" === unitPosFittness \n");
        double cover = 0;
        for( LTLinearObject* lo : lobjs ){
            double r2    = lo->dp(p).norm2();
            double r2max = sq(lo->width);
            //printf( " %f %f \n", r2, r2max  );
            if( r2 < r2max ){
                //E+= getCover( *u, *lo );
                cover = fmax( cover, lo->cover * (r2max-r2)/r2max );
            }
        }
        double overlap = 0;
        double invR2 = 1/(Rdeploy*Rdeploy);
        for( LTUnit* uu : coleagues ){
            if(uu!=u){
                //overlap += 1/( 1.0 + invR2*uu->pos.dist2( p ) ); // lorenz repulsion of neighbors
                overlap += 1/( 1.0 + invR2*uu->goal_pos.dist2( p ) );
            }
        }
        //printf( "overlap %f (%i) cover %f (%i)  \n", overlap, coleagues.size(), cover, lobjs.size() );
        double E = overlap2E*overlap + cover2E*cover;
        if( bConstr ){
            E += ConstrE * sq( p.dist2(ConstrPos)/sq(ConstrRad) );
        }
        return E;
    };

    void tryFindBetterPlace( LTUnit* u, int n, bool bJumpToGoal )const{

        // rather use goal_pos instead of pos

        Vec2d   p0 = u->goal_pos;
        double  E0 = unitPosFittness(u,p0);
        for(int i=0; i<n; i++){
            Vec2d p = p0;
            double fstrat = randf();
            if( fstrat<0.8 ){ p.add( randf(-RSearchLocal ,RSearchLocal),  randf(-RSearchLocal ,RSearchLocal)  ); }
            else{             p.add( randf(-RSearchGlobal,RSearchGlobal), randf(-RSearchGlobal,RSearchGlobal) ); }
            double E = unitPosFittness(u,p);
            if( E>E0 ){ E0 = E; p0 = p; };
        }
        u->goal_pos = p0;
        //if( bJumpToGoal ) u->pos = u->goal_pos;
    }

};

#endif
