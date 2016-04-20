

#ifndef BlockHouseWorld_h
#define BlockHouseWorld_h

#include <unordered_set>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

class BlockHouse{
	public:
    std::vector<Formation*> formations;
    TerrainCubic   terrain;

	//HashMap2D<Particle2D> map;

    double dt_frame  = 0.1;
    int    per_frame = 1;
    double damping   = 0.2;
    double damp;
    double dt;

    void init();
    void update();
    void simulationStep( double dt );

    //void assembleForces( ULONG i );
    //void assembleForces_offside( ULONG i, ULONG j, UINT ni, Particle2D** buf_i );
    //bool moveParticle     ( Particle2D* pi );

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
