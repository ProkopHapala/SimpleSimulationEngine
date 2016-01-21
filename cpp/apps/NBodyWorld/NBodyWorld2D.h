

#ifndef NBodyWorld2D_h
#define NBodyWorld2D_h

#include <unordered_set>

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"
#include "HashMap2D.h"
#include "Body2D.h"

inline void shortRangeForce( const Vec2d& pa, const Vec2d& pb, Vec2d& Fout ){
    Vec2d d; d.set_sub( pb, pa );
    double r2 = d.norm2();
    double fr = 1-r2;
           fr = fr*fr;
    Fout.set( d.x*fr, d.y*fr );
};

class Particle2D: public PointBody2D{
    public:
    double charge;

};


class NBodyWorld{
	public:

	int nParticles;
	Particle2D* particles;

	HashMap2D<Particle2D> map;

    std::unordered_set<ULONG> activeCells;
    std::unordered_set<ULONG> activeCellsNighbors;

    int nActiveParticles;
    Particle2D** activeParticles;

    void init();
    void simulationStep( double dt );
    void activateParticles( ULONG i );
    void assembleForces   ( ULONG i );
    void assembleForces_offside( ULONG i, ULONG j, UINT ni, Particle2D** buf_i );

};

#endif
