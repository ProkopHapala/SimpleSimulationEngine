
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
#include "TrussDynamics_d.h"
#include "TrussDynamics_f.h"
#include "SpaceCraftDynamicsApp.h"

using namespace Mesh;
using namespace SpaceCrafting;

enum class EDIT_MODE:int{ vertex=0, edge=1, component=2, size }; // http://www.cprogramming.com/c++11/c++11-nullptr-strongly-typed-enum-class.html

double elementSize  = 5.;

SpaceCraftSimulator W;

Vec3d wheel_speed      {0.0,0.0,0.0};
Vec3d wheel_speed_setup{ 5.0, 5.0, 5.0 };

void SpaceCraftControl(double dt){ applySliders2sim( *theSpaceCraft, W.sim, (double*)&wheel_speed ); }

void solve_float(){
    //printf( "solve_float() \n" );
    Quat4f* psa =  W.sim_f.ps_pred;
    Quat4f* psb =  W.sim_f.ps_cor;
    for(int itr=0; itr<W.sim_f.nSolverIters; itr++){ 
        W.sim_f.updateJacobi_lin( psa, psb, W.sim_f.bvec );
        for(int i=0; i<W.sim_f.nPoint; i++){  psa[i] = psb[i]; }
    }
}

int make_Skeleton( SpaceCraft* theSpaceCraft, int nnode, Vec3d* node_pos, double* node_size, int ngirdes, Vec2i* girdes, int* girer_nsegs,  int nropes, Vec2i* ropes ){
    printf("### make_Skeleton() creating %d nodes, %d girders, %d ropes\n", nnode, ngirdes, nropes);
    
    // Initialize workshop if needed (you might want to move this outside)
    init_workshop();
    
    int nodeEdgeType = 0;
    for(int i=0; i<nnode; i++){
        theSpaceCraft->add_Node(node_pos[i], node_size[i], nodeEdgeType );
        //printf("Added node %d at position (%g,%g,%g)\n", i, node_pos[i].x, node_pos[i].y, node_pos[i].z);
    }
    
    // Create girders
    Vec3d up = {0.0, 0.0, 1.0}; // Default up direction
    int nseg = 5;               // Number of segments along girder
    int mseg = 3;               // Number of segments around girder
    Vec2d wh = {1.0, 1.0};      // Width and height
    //Quat4i stickTypes = {1, 2, 3, 4}; // Default stick material types
    Quat4i stickTypes = {0,0,0,0}; // Default stick material types
    
    for(int i=0; i<ngirdes; i++){
        Vec2i e = girdes[i];
        // Calculate vector between nodes to determine up direction
        Vec3d dir = node_pos[e.y] - node_pos[e.x];
        dir.normalize();
        // Get a proper up vector (not parallel to dir)
        if(fabs(dir.z) < 0.9){ up = {0.0, 0.0, 1.0}; } 
        else                 { up = {1.0, 0.0, 0.0}; }
        int nseg_i = nseg;
        if(girer_nsegs){ nseg_i = girer_nsegs[i]; }
        theSpaceCraft->add_Girder(e.x, e.y, up, nseg_i, mseg, wh, stickTypes);
        //printf("  Added girder between nodes %d and %d\n", e.x, e.y);
    }
    
    // // Create ropes
    // double thick = 0.1; // Default thickness
    // for(int i=0; i<nropes; i++){
    //     int node1 = ropes[i].a;
    //     int node2 = ropes[i].b;
    //     theSpaceCraft->add_Rope(node1, node2, thick);
    //     printf("  Added rope between nodes %d and %d\n", node1, node2);
    // }
    
    return 1; // Success
}

void keyStateHandling_local( const Uint8 *keys ){
    // This is a bit of a hack, assuming wheel_speed is a global in the main app file
    // A better solution would be a proper control system.
    wheel_speed = Vec3dZero;
    if( keys[ SDL_SCANCODE_KP_5 ] ){ wheel_speed.y=-wheel_speed_setup.y; }
    if( keys[ SDL_SCANCODE_KP_8 ] ){ wheel_speed.y= wheel_speed_setup.y; }
    if( keys[ SDL_SCANCODE_KP_4 ] ){ wheel_speed.x=-wheel_speed_setup.x; }
	if( keys[ SDL_SCANCODE_KP_6 ] ){ wheel_speed.x= wheel_speed_setup.x; }
	if( keys[ SDL_SCANCODE_KP_7 ] ){ wheel_speed.z=-wheel_speed_setup.z; }
	if( keys[ SDL_SCANCODE_KP_9 ] ){ wheel_speed.z= wheel_speed_setup.z; }
    if( keys[ SDL_SCANCODE_KP_0 ] ){ wheel_speed.x=0; wheel_speed.y=0; wheel_speed.z=0; }
}



// ===================== MAIN

int main(int argc, char *argv[]){

    // disable stdout buffering
    setbuf(stdout, NULL);
    //Or use the more flexible setvbuf:
    //setvbuf(stdout, NULL, _IONBF, 0); 

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
    funcs["-G"       ]={1,[&](const char** ss){ Vec3d  G;  sscanf( ss[0], "%lf,%lf,%lf", &G.x, &G.y, &G.z ); W.sim.accel.f=G; W.sim_f.accel.f=(Vec3f)G; printf( "COMMAND LINE: -G( sim.accel: %f %f %f ) \n", W.sim.accel.x, W.sim.accel.y, W.sim.accel.z );    } };
    funcs["-omega"   ]={1,[&](const char** ss){ Quat4d o;  sscanf( ss[0], "%lf,%lf,%lf", &o.x, &o.y, &o.z ); o.w=o.f.normalize(); W.sim.omega=o; W.sim_f.omega=(Quat4f)o; printf( "COMMAND LINE: -omega( sim.omega: %f %f %f %f ) \n", W.sim.omega.x, W.sim.omega.y, W.sim.omega.z, W.sim.omega.w );    } };
    funcs["-drag"    ]={1,[&](const char** ss){ sscanf( ss[0], "%lf", &W.sim.Cdrag );  printf( "COMMAND LINE: -drag( sim.Cdrag: %f ) \n", W.sim.Cdrag );    } };
   
    funcs["-debug_orig"]={0,[&](const char** ss){ bUseOriginalLuaWrappers=true; printf( "COMMAND LINE: -debug_orig \n" ); } };

    funcs["-shape"   ]={1,[&](const char** ss){ 
        int ishape, nseg; int nret=sscanf( ss[0],"%i,%i", &ishape, &nseg ); if(nret<2){ printf( "ERROR reading argument -shape nret(%i)<2 %s \n", nret, ss[0] ); exit(0);}
        printf( "COMMAND LINE: -shape( ishape: %i nseg: %i ) \n", ishape, nseg );
        init_workshop ();
        makeTrussShape( W.mesh, ishape, nseg, 100.0, 10.0 );
        W.initSimulators();
    }};

    funcs["-oct_nodes"] = {0, [&](const char**){
        printf("funcs[-oct_nodes]: Manual Construction Blocks Test:\n");
        //testConstructionBlocks(bb, truss);
        const int nnodes = 5;
        Vec3d nodes[nnodes] = {
            {0.0,    0.0,    0.0}, // 0
            {0.0,    0.0,  200.0}, // 1
            {0.0,    0.0, -150.0}, // 2
            {0.0, -100.0,    0.0}, // 3
            {0.0,  100.0,    0.0}  // 4
        };
        double node_sizes[nnodes] = {10.0, 5., 5.0, 5.0, 5.0 };
        const int nGirders = 4;
        Vec2i girdes[nGirders] = {{0,1}, {0,2}, {0,3}, {0,4}};
        int girer_nsegs[nGirders] = {8, 6, 4, 4};
        const int nRopes = 4;
        Vec2i ropes[nRopes] = { {1,3}, {1,4}, {2,3}, {2,4} };

        // ----- here we should create Nodes, girders and ropes in theSpaceCraft ( see SpaceCrafting::SpaceCraft in SpaceCraft.h )
        make_Skeleton( theSpaceCraft, nnodes, nodes, node_sizes, nGirders, girdes, girer_nsegs, nRopes, ropes );

        // Attach a wheel to the girders
        const int   wheel_girders[4] = {2,3,1,0};
        // We must provide at least 3 points to define the circle. The 4th can be calculated by intersection.
        // We set the first three to the midpoint (0.5) of their respective girders.
        const float wheel_pos[4]     = {0.5f, 0.5f, 0.3f, -1.0f};
        Vec3d p0_ring_center_guess = {0.0, 0.0, -50.0}; // A point "near" the desired ring center to help resolve which side of the girder to attach to
        theSpaceCraft->make_Ring2(wheel_girders, wheel_pos, p0_ring_center_guess, 32, {2.0, 2.0}, "GS1_long", {0,0,0,0}, 0);

        // ----- Build mesh from SpaceCraft and TrussDynamics_d to prepare simulation

        CMesh oct=(CMesh){Solids::Octahedron_nverts,Solids::Octahedron_nedges,Solids::Octahedron_ntris,Solids::Octahedron_nplanes, Solids::Octahedron_verts, Solids::Octahedron_edges, Solids::Octahedron_tris, Solids::Octahedron_planes, Solids::Octahedron_planeVs};
        theSpaceCraft->nodeMeshes.push_back(oct);
        
        W.mesh.clear();
        BuildCraft_blocks( W.mesh, *theSpaceCraft, 30.0 );
        W.mesh.printSizes();
        W.initSimulators();

        // --- Initialize slider paths now that the mesh and simulation points exist.
        printf("--- Initializing Slider Paths ---\n");
        for( Ring* o : theSpaceCraft->rings ){
            // bSelf=true: The path is on the ring itself. bShared=true: All sliders on this ring share the same path.
            o->updateSlidersPaths( true, true, W.sim.points );
        }
        sliders2edgeverts( *theSpaceCraft, W.sim );

        W.sim.user_update = SpaceCraftControl;

        // ---------- Backup from constructionBlockApp.cpp, this should perhaps go to BuildCraft_truss in SpaceCraft2Mesh2.h
        //bool bUseSpecialPlanes=false;
        // CMesh oct=(CMesh){Solids::Octahedron_nverts,Solids::Octahedron_nedges,Solids::Octahedron_ntris,Solids::Octahedron_nplanes, Solids::Octahedron_verts, Solids::Octahedron_edges, Solids::Octahedron_tris, Solids::Octahedron_planes, Solids::Octahedron_planeVs};
        // truss.facingNodes( oct, nnodes, node_positions, chs );
        // truss.bridgeFacingPolygons( nedges, edges, node_positions, 4, chs );
        // printf("  truss.printSizes():  "); truss.printSizes();
        // truss.write_obj("truss.obj", ObjMask::Verts | ObjMask::Edges | ObjMask::Polygons );
        //printf("Extruding chunk %i by 2.0 units...\n", chs.a);
        //truss.extrudeFace(chs.a, 2.0);
    }};

    funcs["-s"]={1,[&](const char** ss){ 
        W.reloadShip( ss[0] );
        W.initSimulators( );
    }};
    process_args( argc, argv, funcs );
    if( W.sim.nPoint == 0 ){ W.initSimDefault(); }

    W.sim.extern_b = W.sim_f.bvec;
    W.sim.extern_x = W.sim_f.ps_pred;
    W.sim.extern_solve = &solve_float;   // function pointer extern_solve

	//thisApp = new SpaceCraftDynamicsApp( junk , 800, 600 );
	app->loop( 1000000 );
	return 0;
}
