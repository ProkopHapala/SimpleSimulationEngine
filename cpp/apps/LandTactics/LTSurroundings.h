
#ifndef LTSurroundings_h
#define LTSurroundings_h

#include <vector>

#include "LTUnit.h"
#include "LTShelter.h"

class LTsurrounding{
    std::vector<LTStaticObject*> objs;
    std::vector<LTLinearObject*> lobjs;
    std::vector<LTUnit*>         enemies;
    std::vector<LTUnit*>         coleagues;

    double RSearchLocal  = 5.0;
    double RSearchGlobal = 50.0;

    double unitPosFittness( const LTUnit* u, const Vec2d& p ){
        double E = 0.0;
        for( LTLinearObject* lo : lobjs ){
            if( lo->dp(p).norm2() < sq(lo->width) ){
                //E+= getCover( *u, *lo );
                E += 1.0;
            }
        }
        for( LTUnit* uu : coleagues ){
            if(uu!=u){
                E -= 1/( 1.0+uu->pos.dist2( u->pos ) ); // lorenz repulsion of neighbors
            }
        }
    };

    void tryFindBetterPlace( LTUnit* u, int n ){
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
    }

};

#endif
