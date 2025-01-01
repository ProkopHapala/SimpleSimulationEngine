
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "Draw3D.h"
#include "SDL_utils.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"

int verbosity = 0;

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

#include "spaceCraftDynamicsShared.h"

using namespace SpaceCrafting;

enum class EDIT_MODE:int{ vertex=0, edge=1, component=2, size }; // http://www.cprogramming.com/c++11/c++11-nullptr-strongly-typed-enum-class.html

double elementSize  = 5.;

Mesh::Builder2 mesh;
OrbSim         sim;

// ===================== MAIN

int main(int argc, char *argv[]){
    // example: use like : ./spaceCraftEditor -s data/ship_ICF_interceptor_1.lua
    printf( "argc %i \n", argc );
    LambdaDict funcs;
    funcs["-s"]={1,[&](const char** ss){ 
        reloadShip( ss[0], mesh );
        to_OrbSim( sim, mesh );
    }}; 
    SDL_DisplayMode dm = initSDLOGL( 8 );
	int junk;
	SpaceCraftDynamicsApp * app = new SpaceCraftDynamicsApp( junk, dm.w-150, dm.h-100 );
    app->_mesh = &mesh;
    app->_sim  = &sim;
    process_args( argc, argv, funcs );
    if( sim.nPoint == 0 ){ app->initSimDefault(); }

	//thisApp = new SpaceCraftDynamicsApp( junk , 800, 600 );
	app->loop( 1000000 );
	return 0;
}
















