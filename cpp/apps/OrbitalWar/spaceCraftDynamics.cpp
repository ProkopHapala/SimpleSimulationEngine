
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "globals.h"
#include "Draw2D.h"
#include "Draw3D.h"
#include "SDL_utils.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"

//int verbosity = 0;

//#include "Truss.h"
#include "SpaceCraft.h"
#include "SpaceCraft2Mesh2.h"
#include "SpaceCraftDraw.h"
//#include "SoftBody.h"

//#include "SphereSampling.h"
//#include "DrawSphereMap.h"
//#include "Draw3D_Surf.h"

#include "AppSDL2OGL_3D.h"
#include "GUI.h"
#include "IO_utils.h"
#include "argparse.h"
#include "testUtils.h"

//#include "EditSpaceCraft.h"
//#include "TriangleRayTracer.h"
//#include "Radiosity.h"

#include "Tree.h"

#include "spaceCraftEditorUtils.h"
#include "OrbSim_d.h"
#include "OrbSim_f.h"
#include "SpaceCraftDynamicsApp.h"

using namespace SpaceCrafting;

enum class EDIT_MODE:int{ vertex=0, edge=1, component=2, size }; // http://www.cprogramming.com/c++11/c++11-nullptr-strongly-typed-enum-class.html

double elementSize  = 5.;

// Mesh::Builder2 mesh;
// OrbSim         sim;
// OrbSim_f       sim_f;
SpaceCraftSimulator W;

// ===================== MAIN

int main(int argc, char *argv[]){
    // example: use like : ./spaceCraftEditor -s data/ship_ICF_interceptor_1.lua
    printf( "argc %i \n", argc );
    SDL_DisplayMode dm = initSDLOGL( 8 );
	int junk;
	SpaceCraftDynamicsApp * app = new SpaceCraftDynamicsApp( junk, dm.w-150, dm.h-100 );
    app->bindSimulators( &W ); 

    LambdaDict funcs;
    funcs["-float"   ]={0,[&](const char** ss){ app->bDouble=false; printf( "COMMAND LINE: -float  app->bDouble=%i\n", app->bDouble );  } };
    funcs["-double"  ]={0,[&](const char** ss){ app->bDouble=true;  printf( "COMMAND LINE: -double app->bDouble=%i\n", app->bDouble );  } };
    funcs["-fix"     ]={1,[&](const char** ss){ int n =  readlist( ss[0], W.fixPoints); printf("COMMAND LINE: -fix[%i]{%s}\n", n, ss[0] );  } };
    funcs["-perframe"]={1,[&](const char** ss){            sscanf( ss[0], "%i", &app->perFrame ); printf( "COMMAND LINE: -perframe(%i) \n", app->perFrame ); } };
    funcs["-nsolve"  ]={1,[&](const char** ss){ int nsolv; sscanf( ss[0], "%i", &nsolv ); printf( "COMMAND LINE: -nsolve(%i) \n", nsolv ); W.sim_f.nSolverIters=nsolv; W.sim.nSolverIters=nsolv;  } };
    funcs["-method"  ]={1,[&](const char** ss){ int im;    sscanf( ss[0], "%i", &im    ); printf( "COMMAND LINE: -method(%i) \n", im    ); W.sim_f.linSolveMethod=im;  W.sim.linSolveMethod=im;   } };
    funcs["-bmix"    ]={1,[&](const char** ss){ int istart; float bmix;  sscanf( ss[0], "%i,%f", &istart, &bmix ); W.sim.mixer.b_end=bmix; W.sim.mixer.istart=istart; printf( "COMMAND LINE: -bmix( istart:%i bmix: %f ) \n", W.sim.mixer.istart, W.sim.mixer.b_end );    } };
    funcs["-dt"      ]={1,[&](const char** ss){ float dt;  sscanf( ss[0], "%f", &dt ); W.sim.dt=dt; W.sim_f.dt=dt; printf( "COMMAND LINE: -dt( dt: %f ) \n", W.sim.dt );    } };
    
    funcs["-shape"   ]={1,[&](const char** ss){ 
        int ishape, nseg; int nret=sscanf( ss[0],"%i,%i", &ishape, &nseg ); if(nret<2){ printf( "ERROR reading argument -shape nret(%i)<2 %s \n", nret, ss[0] ); exit(0);}
        printf( "COMMAND LINE: -shape( ishape: %i nseg: %i ) \n", ishape, nseg );
        init_workshop ();
        makeTrussShape( W.mesh, ishape, nseg, 100.0, 10.0 );
        W.initSimulators();
    }};
    funcs["-s"]={1,[&](const char** ss){ 
        W.reloadShip( ss[0] );
        W.initSimulators( );
    }};
    process_args( argc, argv, funcs );
    if( W.sim.nPoint == 0 ){ W.initSimDefault(); }

	//thisApp = new SpaceCraftDynamicsApp( junk , 800, 600 );
	app->loop( 1000000 );
	return 0;
}
















