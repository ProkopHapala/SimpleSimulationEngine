

#ifndef NBodyWorld2D_h
#define NBodyWorld2D_h

#include <unordered_set>

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"
#include "HashMap2D.h"
#include "Body2D.h"


inline void stringForce( const Vec2d& pa, const Vec2d& pb, double k, Vec2d& Fout ){
    Vec2d d;
    d.set_sub( pb, pa );
    Fout.set( d.x*k, d.y*k );
};


inline bool pairwiseForce( const Vec2d& pa, const Vec2d& pb, double qq, Vec2d& Fout ){
    const double r2max = 4.0d;
    Vec2d d;
    //d.set_sub( pb, pa );
    d.set_sub( pa, pb );
    double r2 = d.norm2();
    if( r2 < r2max ){
        double mr2 = r2max-r2;
        double fr  = ( 1/(r2+0.01) + qq - 0.1 )*mr2*mr2;
        Fout.set( d.x*fr, d.y*fr );
        return true;
    }else{
        Fout.set( 0.0, 0.0 );
        return false;
    }
};


const double f2conv = 4.0;
const double v2conv = 0.01;

class Particle2D: public PointBody2D{
    public:
    double charge;
    int    stepsConverged;

/*
    inline bool converged( ){
        if ( force.norm2() > force2conv ) return false;
        if ( vel  .norm2() > vel2conv   ) return false;
        return true;
    }
*/

};


class NBodyWorld{
	public:
    double dt_frame  = 0.2;
    int    per_frame = 10;
    double damping   = 0.2;

    double anchorStiffness = 0.5;

    double damp;
    double dt;

    double v2max;
    double f2max;

	int nParticles;
	Particle2D* particles;

	HashMap2D<Particle2D> map;

    std::unordered_set<ULONG> activeCells;
    std::unordered_set<ULONG> activeCellsNeighbors;

    int nActiveParticles;
    Particle2D** activeParticles;

    Vec2d anchor;
    Particle2D* picked = NULL;

    int n_moves, n_interactions;

    void init();
    void update();
    void simulationStep_BruteForce( double dt );
    void simulationStep_semiBruteForce( double dt );
    void simulationStep( double dt );
    void activateCell  ( ULONG i );
    void activateAroundParticle( Particle2D* pi, ULONG& icell_old );
    void assembleForces( ULONG i );
    void assembleForces_offside( ULONG i, ULONG j, UINT ni, Particle2D** buf_i );
    bool moveParticle     ( Particle2D* pi );
    void moveParticleDebug( Particle2D* pi, int i );
    void checkHashMapConsistency( );

    inline void setSimParams( double dt_frame_, double per_frame_, double damping_ ){
        dt_frame  = dt_frame_;
        per_frame = per_frame_;
        damping   = damping_;
        evalAuxSimParams();
        //printf( " dt_frame, per_frame,  dt, damp " );
    }

    inline void evalAuxSimParams(){
        const double dampMin = 0.5;
        dt   = dt_frame / per_frame;
        damp = ( 1 - damping * dt );
        if( damp < dampMin ){ damp = dampMin; }
    };

};

#endif
