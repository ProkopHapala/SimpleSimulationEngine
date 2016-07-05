

#ifndef TacWorld_h
#define TacWorld_h

#include <unordered_set>

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"
#include "HashMap2D.h"
#include "Body2D.h"

#include "TileBuffer2D.h"
#include "Ruler2DFast.h"
#include "TerrainCubic.h"

#include "Unit.h"

class TacWorld{
	public:
    //std::vector<BattleLine*> battleLines;
    std::vector<Unit*>    units;
    std::vector<Faction*> factions;

    static constexpr int nTypesMax = 16;
    UnitType      unitTypes[nTypesMax];
    int           nTypes = 0;

    double RmaxInteract = 1.5;

    Ruler2DFast colruler;

    TerrainCubic   terrain;

    int nSoldierInteractions = 0;
    int nSoldiers = 0;



	//HashMap2D<Particle2D> map;

    double dt_frame  = 0.1;
    int    per_frame = 1;
    double damping   = 0.2;
    double damp;
    double dt;

    void init();
    void update();
    void simulationStep( double dt );
    int  getUnitAt( const Vec2d& p, Faction * faction );

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
