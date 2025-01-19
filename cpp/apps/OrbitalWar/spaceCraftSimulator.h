#ifndef SpaceCraftSimulator_h
#define SpaceCraftSimulator_h

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "Draw3D.h"
#include "SDL_utils.h"
#include "AppSDL2OGL_3D.h"
#include "GUI.h"
#include "SpaceCraft.h"
#include "SpaceCraft2Mesh2.h"
#include "SpaceCraftDraw.h"
#include "TrussDynamics_d.h"
#include "TrussDynamics_f.h"
//#include "OCL_Orb.h"

namespace SpaceCrafting {

int readlist(const char* input, std::vector<int>& lst) {
    int num;
    const char* p = input;
    int i=0;
    while (*p) {
        if (sscanf(p, "%d", &num) == 1) {
            //printf("readlist() [%i]: %d\n", i, num);
            lst.push_back(num);
            while (*p && *p != ',') p++;
            if (*p == ',') p++;
            i++;
        } else {
            // Handle invalid input
            printf("Invalid input at: %s\n", p);
            break;
        }
    }
    return i;
}

void init_workshop(){
    //StickMaterial *o = new StickMaterial();
    workshop.add_Material     ( "Steel", 7.89e+3, 1.2e+9, 1.2e+9, 200.0e+9, 200.0e+9, 0.85, 800 );
    workshop.add_StickMaterial( "GS1_long", "Steel", 0.1, 0.005, 0.0 );
}

void to_TrussDynamics_d( TrussDynamics_d& sim, Mesh::Builder2& mesh, int nfix=0, int* fixPoints=0 ){
    printf("#### ==== SpaceCraftSimulator::to_TrussDynamics_d() sim.linSolveMethod=%i \n", sim.linSolveMethod );
    exportSim( sim, mesh, workshop );
    if( sim.nPoint==0 ){ printf( "ERROR in SpaceCraftSimulator::to_TrussDynamics_d() sim.nPoint=%i => exit() \n", sim.nPoint ); exit(0); };
    //if(fixPoints.size()>0) sim.setFixPoints( fixPoints.size(), fixPoints.data() );
    if(nfix>0) sim.setFixPoints( nfix, fixPoints );
    if( ( sim.linSolveMethod == (int)TrussDynamics_d::LinSolveMethod::Cholesky       ) ||
        ( sim.linSolveMethod == (int)TrussDynamics_d::LinSolveMethod::CholeskySparse ) ){
        sim.prepare_LinearSystem( true, true, true, 256 );
    }
    sim.cleanVel();
    sim.cleanForce();
    sim.addAngularVelocity( sim.pos0, sim.omega.f*sim.omega.w );   
    //if( omega.norm2()>1e-16 )sim.addAngularVelocity( p0, omega );
    //sim.addAngularVelocity( p0, Vec3dZ*0.01 );   
    //sim.addAngularVelocity2( p0, Vec3dZ*0.000001 ); // accelerate more the end
};

void to_TrussDynamics_f(TrussDynamics_f& sim, Mesh::Builder2& mesh, int nfix=0, int* fixPoints=0 ){
    printf("#### ==== SpaceCraftSimulator::to_TrussDynamics_f() sim.linSolveMethod=%i \n", sim.linSolveMethod);
    exportSim( sim, mesh, workshop );
    if( sim.nPoint==0 ){ printf( "ERROR in SpaceCraftSimulator::to_TrussDynamics_f() sim.nPoint=%i => exit() \n", sim.nPoint ); exit(0); };
    if(nfix>0) sim.setFixPoints( nfix, fixPoints );
    //if(fixPoints.size()>0) sim.setFixPoints( fixPoints.size(), fixPoints.data() );
    if( ( sim.linSolveMethod == (int)TrussDynamics_d::LinSolveMethod::Cholesky       ) ||
        ( sim.linSolveMethod == (int)TrussDynamics_d::LinSolveMethod::CholeskySparse ) ){
        sim.prepare_LinearSystem( true, true, true, 256 );
    }
    //sim.linSolveMethod = (int)TrussDynamics_d::LinSolveMethod::CG;
    sim.linSolveMethod = (int)TrussDynamics_d::LinSolveMethod::Cholesky;
    //sim.linSolveMethod = (int)TrussDynamics_d::LinSolveMethod::CholeskySparse;
    sim.cleanVel();
    sim.cleanForce();
    //if( omega.norm2()>1e-16 )sim.addAngularVelocity( p0, omega );
    //sim.addAngularVelocity( p0, Vec3fZ*0.01 );
};


// =======================================================================
class SpaceCraftSimulator { public:

    Mesh::Builder2 mesh;
    TrussDynamics_d sim;
    TrussDynamics_f       sim_f;

    std::vector<int> fixPoints;

    bool   bDouble = false;
    bool   bRun    = false;
    double time=0;

    // ---- Functions
    void         reloadShip ( const char* fname );
    virtual void initSimulators( double dt=0.1, Vec3d p0=Vec3dZero, Vec3d omega=Vec3dZero );
    virtual void initSimDefault();

    virtual TrussDynamics_d* getTrussSim  (){ return &sim;   };
    virtual TrussDynamics_f* getTrussSim_f(){ return &sim_f; };
    virtual Mesh::Builder2*  getMesh(){ return &mesh; };

    //virtual void mouseHandling( );
    //void drawSim  ( TrussDynamics_d&   sim );
    //void drawSim_f( TrussDynamics_f& sim );
	//SpaceCraftSimulator();

};

void SpaceCraftSimulator::reloadShip( const char* fname ){
    printf("#### START reloadShip('%s')\n", fname );
    theSpaceCraft->clear();                       DEBUG
    //luaL_dofile(theLua, "data/spaceshil1.lua"); DEBUG
    if( Lua::dofile(theLua,fname) ){ printf( "ERROR in reloadShip() Lua::dofile(%s) \n", fname ); exit(0); }
    printf( "Lua::dofile(%s) DONE \n", fname );
    theSpaceCraft->checkIntegrity();
    mesh.clear();
    BuildCraft_truss( mesh, *theSpaceCraft, 30.0 );
    mesh.printSizes();
    printf("#### END reloadShip() \n");
};

void SpaceCraftSimulator::initSimulators( double dt, Vec3d p0, Vec3d omega ){
    to_TrussDynamics_d  ( sim,   mesh, fixPoints.size(), fixPoints.data() );
    to_TrussDynamics_f( sim_f, mesh, fixPoints.size(), fixPoints.data() );
    //sim.cleanVel( Quat4d{0.0,1.0,0.0,0.0} );
    //sim.addAngularVelocity( p0, Vec3fZ*0.01 );

}

void SpaceCraftSimulator::initSimDefault(){
    init_workshop ();
    //makeTrussShape( mesh, 3, 10, 100.0, 10.0 );
    makeTrussShape( mesh, 2, 100, 100.0, 10.0 );
    initSimulators();
}

}// namespace SpaceCrafting

#endif // SpaceCraftSimulator_h
