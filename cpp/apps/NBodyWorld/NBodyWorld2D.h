

#ifndef NBodyWorld2D_h
#define NBodyWorld2D_h

#include <unordered_set>

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"
#include "HashMap2D.h"
#include "Body2D.h"

inline void pairwiseForce( const Vec2d& pa, const Vec2d& pb, double qq, Vec2d& Fout ){
    const double r2max = 0.5d;
    Vec2d d;
    //d.set_sub( pb, pa );
    d.set_sub( pa, pb );
    double r2 = d.norm2();
    if( r2 < r2max ){
        double mr2 = r2max-r2;
        double fr  = ( r2max/r2 + qq - 2.0 )*mr2*mr2;
        Fout.set( d.x*fr, d.y*fr );
    }else{
        Fout.set( 0.0, 0.0 );
    }
};

class Particle2D: public PointBody2D{
    public:
    double charge;

};


class NBodyWorld{
	public:
    double  dt_frame = 1.0;
    int    per_frame = 1;

	int nParticles;
	Particle2D* particles;

	HashMap2D<Particle2D> map;

    std::unordered_set<ULONG> activeCells;
    std::unordered_set<ULONG> activeCellsNeighbors;

    int nActiveParticles;
    Particle2D** activeParticles;

    void init();
    void update();
    void simulationStep_BruteForce( double dt );
    void simulationStep( double dt );
    void activateCell  ( ULONG i );
    void activateAroundParticle( Particle2D* pi, ULONG& icell_old );
    void assembleForces( ULONG i );
    void assembleForces_offside( ULONG i, ULONG j, UINT ni, Particle2D** buf_i );

};

#endif
