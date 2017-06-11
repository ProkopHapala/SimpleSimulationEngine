

#ifndef FormationWorld_h
#define FormationWorld_h

#include <unordered_map>

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"
#include "HashMap2D.h"
#include "Body2D.h"

#include "TileBuffer2D.h"
#include "Ruler2DFast.h"
#include "TerrainCubic.h"

#include "FormationTacticsCommon.h"
#include "Formation.h"
#include "BattleLine.h"
#include "Faction.h"



class FormationWorld{
	public:
    //std::vector<BattleLine*> battleLines;
    std::vector<Formation*>  formations;
    std::vector<Faction*>    factions;
    std::vector<SoldierType> soldierTypes;

    std::unordered_map<std::string,SoldierType*> name2soldierType;

    double RmaxInteract = 1.5;

    Ruler2DFast colruler;
    TileBuffer2D<Soldier*,16,16,64> colbuf;

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
    int  formationInteractions( );
    int  formationInteractions_buff( );

    void refreshFormations( );

    int  loadSoldierTypes( char * fname );
    void updateSoldierTypeMap();

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
