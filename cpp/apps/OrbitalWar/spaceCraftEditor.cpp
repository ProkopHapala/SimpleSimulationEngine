
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

//int verbosity = 0;


#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"

#include "testUtils.h"

#include "SpaceCraft.h"
//#include "SpaceCraft2Mesh.h"   // deprecated
#include "MeshBuilder2.h"
#include "SpaceCraft2Mesh2.h"
#include "SoftBody.h"

#include "Truss.h"
#include "SpaceCraft2Truss.h" // deprecated

#include "SpaceCraftDraw.h"

#include "SphereSampling.h"
#include "DrawSphereMap.h"
#include "Draw3D_Surf.h"

#include "IO_utils.h"

#include "EditSpaceCraft.h"

#include "OrbSim_f.h"
#include "OrbSim_d.h"
#include "TriangleRayTracer.h"
#include "Radiosity.h"

#include "Tree.h"

#include "spaceCraftEditorUtils.h"

#include "AppSDL2OGL_3D.h"
#include "GUI.h"
#include "SpaceCraftGUI.h"

#include "argparse.h"


#include "Lingebra.h"

// ======================  Global Variables & Declarations

using namespace SpaceCrafting;

bool bShipReady = false;
Mesh::Builder2 mesh2;
OrbSim sim;
int glo_truss=0, glo_capsula=0, glo_ship=0;
char str_tmp[8096];
double elementSize  = 5.;
bool bRun = false;

int iFrame = 0;

Vec3d wheel_speed       = {0.0,0.0,0.0};
Vec3d wheel_speed_setup = { 5.0, 5.0, 5.0 };


class LinSolverOrbSim : public LinSolver{ public:
    OrbSim* sim=0;
    virtual void dotFunc( int n, double * x, double * Ax ) override {
        //printf( "LinSolverOrbSim::dotFunc(n=%i) \n", n );
        //sim->dot_Linearized_neighs2(n, x, Ax);
        sim->dot_Linearized_bonds(n, x, Ax);
    }
};
LinSolverOrbSim linSolver;

void SpaceCraftControl(double dt){ applySliders2sim( *theSpaceCraft, sim, (double*)&wheel_speed ); }

void runSim( OrbSim& sim, int niter=100 ){
    long t0 = getCPUticks();
    if(bRun){
        //printf( "runSim() linSolveMethod=%i\n", sim.linSolveMethod  );
        if( sim.linSolveMethod == (int)OrbSim::LinSolveMethod::Force ){
            //sim.run( niter, 1e-3, 1e-8 );
            sim.run_omp( niter, false, 1e-3, 1e-4 );
            //sim.run_omp( niter, false, 1e-6, 1e-4 );
        }else{
            sim.run_LinSolve( niter );
        }
        //sim.run_omp( 1, true, 1e-3, 1e-4 );
        double T = (getCPUticks()-t0)*1e-6;
        //printf( "runSim() DONE T=%g[ms] %g[ms/iter] niter=%i,nP=%i,nE=%i \n", T, T/niter, niter, sim.nPoint, sim.nNeighMax );
        //printf( "runSim() DONE T=%g[ms] %g[ms/iter] niter=%i,nP=%i,nE=%i \n", T, T/niter, niter, sim.nPoint, sim.nNeighMax );
        sim.updateInveriants(false);
    }
    
    sim.evalBondTension();
    renderPointForces( sim.nPoint, sim.points, sim.forces, 1e-3 );

    if(sim.pointBBs.ncell>0) updatePointBBs( sim.pointBBs, sim.BBs, sim.points,            true );  // It crashes here because of the wrong obj2cell mapping
    if(sim.edgeBBs .ncell>0) updateEdgeBBs ( sim.edgeBBs,  sim.BBs, sim.bonds, sim.points, false );
    
}

// ======================  Free Functions

void renderShip(){
    //printf( "SpaceCraftEditorApp.cpp::renderShip() \n" );
    if(glo_ship){ glDeleteLists(glo_ship,1); };
    glo_ship = glGenLists(1);
    glNewList( glo_ship, GL_COMPILE );
    //glColor3f(0.2,0.2,0.2);

    //  If we use drawSpaceCraft() there are no problems with backface lighting
    //drawSpaceCraft( *theSpaceCraft, 1, false, true );

    /*
    // If we use drawSpaceCraft_Mesh() there are problems with backface lighting
    Mesh::Builder mesh;
    glColor3f(1.0,1.0,1.0);
    drawSpaceCraft_Mesh( *theSpaceCraft, mesh, 1, false, true, (Vec3d){0.5,0.5,0.5} );
    mesh.newSub();
    mesh.printSizes();
    Draw3D::drawMeshBuilder( mesh );
    mesh.write_obj( "ship.obj" );
    */

    //mesh2.max_size = 30.0;
    //sim.printAllNeighs();
    //theSpaceCraft->printAll_girders();
    //theSpaceCraft->updateSliderPaths();
    
    glEnable(GL_LIGHTING);
    glEnable     ( GL_LIGHT0           );
    glEnable     ( GL_NORMALIZE        );
    Draw3D::drawMeshBuilder2( mesh2, 0b110, 1, true, true );

    /*
    Draw3D::color(Vec3d{1.0f,0.0f,1.0f});
    for(int i=0; i<theSpaceCraft->sliders.size(); i++){
        const Slider& o = theSpaceCraft->sliders[i];
        Draw3D::drawLineStrip( o.path.n, o.path.ps, sim.points, o.path.closed );
    }
    */
    
    //Draw3D::color(Vec3d{1.0,0.,1.});
    //for( const Node& o : theSpaceCraft->nodes ){ Draw3D::drawPointCross( o.pos, 5 ); }

    /*
    radiositySolver.clearTriangles();
    //theSpaceCraft->plate2raytracer( theSpaceCraft->shields[0], radiositySolver, elementSize, true );
    //theSpaceCraft->plate2raytracer( theSpaceCraft->shields[1], radiositySolver, elementSize, true );
    //theSpaceCraft->plate2raytracer( theSpaceCraft->shields[2], radiositySolver, elementSize, true );
    theSpaceCraft->toRayTracer( radiositySolver, elementSize );

    // radiosity samples
    glColor3f(.0,.0,.0);
    glBegin(GL_LINES);
    TriangleRayTracer& rt = radiositySolver;
    for( const SurfElement& s: rt.elements ){
        Draw3D::vertex( s.pos );
        Draw3D::vertex( s.pos+s.normal*elementSize );
    }
    glEnd();
    */

    glEndList();
}

void solveTrussCG(){
    LinSolverOrbSim& ls = linSolver;
    int nitr  = 3;
    double dt=0.05;
    sim.cleanVel();
    for(int itr=0; itr<nitr; itr++){ 
        printf( "--------- move[%i]\n", itr );
        printf("pre-move ps:"); print_vector(sim.nPoint, (double*)sim.points, 4,0,3 );
        //printf("pre-move fs:"); print_vector(sim.nPoint, (double*)sim.forces, 4,0,3 );
        for(int i=0; i<sim.nPoint; i++){ 
            sim.vel   [i].f.add_mul( sim.forces[i].f, dt); 
            sim.points[i].f.add_mul( sim.vel   [i].f, dt); 
        }
        printf("post-move ps:"); print_vector(sim.nPoint, (double*)sim.points, 4,0,3 );
        //printf("post-move vs:"); print_vector(sim.nPoint, (double*)sim.vel   , 4,0,3 );
        //printf("linearize ps:"); print_vector(sim.nPoint, (double*)sim.points,  4,0,3 );
        sim.prepareLinearizedTruss_ling(linSolver.b); 
        printf("f(Ax=f):     "); VecN::print_vector(ls.n, ls.b);
        VecN::set( ls.n, 0.0, ls.x );
        int ncg = ls.solve_CG_( 10, 1.e-6 );
        printf("dp: "); VecN::print_vector( ls.n, ls.x );
        double dd = VecN::sum2( ls.n, ls.x ); printf( "move[%i,%i] |d|=%g \n", itr, ncg, sqrt(dd) );
        Vec3d * dx = (Vec3d*)ls.x;
        for(int i=0; i<sim.nPoint; i++){ // correct point positions (ignoring fixed points)
            if( ! (sim.kFix[i]>1.0) ){ 
                //dx[i].set(0.0);
                sim.points[i].f.add( dx[i] ); 
            }
        }        
    }
}

void makeTestTruss(){
    printf( "#### makeTestTruss() \n" );
    mesh2.clear();
    int stickType = 0;
    //workshop.stickMaterials.vec[stickType]->k = 0.01;
    mesh2.stick  (    {-2.0, 0.0, 0.0}, 
                      {-1.0,-0.1, 0.0}, stickType );
    mesh2.stickTo( 1, { 0.0,-0.2, 0.0}, stickType );
    mesh2.stickTo( 2, { 1.0,-0.1, 0.0}, stickType );
    mesh2.stickTo( 3, { 2.0, 0.0, 0.0}, stickType );
    //for(int i=0; i<mesh2.verts.size(); i++){ Vec3d& p = mesh2.verts[i].pos; printf( "vert[%i] i %i p(%g,%g,%g)\n", i, p.x, p.y, p.z );}
    exportSim( sim, mesh2, workshop );
    _realloc ( sim.dpos, sim.nBonds );
    _realloc( sim.points_bak, sim.nPoint );
    sim.points[0].w = 1e+300; // fix first point
    sim.points[4].w = 1e+300; // fix last point
    sim.cleanForce();
    sim.cleanVel();
    sim.forces[2].f.y = -5.0;
    printf( "#### makeTestTruss() DONE \n" );
}

void reloadShip( const char* fname  ){
    printf("#### START reloadShip('%s')\n", fname );
    theSpaceCraft->clear();                  // clear all components
    if( Lua::dofile(theLua,fname) ){ printf( "ERROR in reloadShip() Lua::dofile(%s) \n", fname ); exit(0); }
    theSpaceCraft->checkIntegrity();

    mesh2.clear();
    BuildCraft_truss( mesh2, *theSpaceCraft, 30.0 );

    exportSim( sim, mesh2, workshop );
    if( ( sim.linSolveMethod == (int)OrbSim::LinSolveMethod::Cholesky      ) ||
        ( sim.linSolveMethod == (int)OrbSim::LinSolveMethod::CholeskySparse ) ){
        //sim.dt = 0.01;
        //sim.accel = Quat4d{0.0,0.0,0.0,0.0};
        sim.prepare_LinearSystem( true, true, true, 256 );
    }
    sim.cleanVel();
    sim.cleanForce();

    linSolver.realloc( sim.nPoint*3, true );
    linSolver.sim = &sim;

    makeBBoxes( *theSpaceCraft, sim );
    makePointCunks( sim.edgeBBs, sim.bonds, sim.pointChunks );

    sim.user_update = SpaceCraftControl;
    sliders2edgeverts( *theSpaceCraft, sim );
    renderShip();
    //sim.updateInveriants(true);
    sim.updateInveriants(false);
    bShipReady = true;
    printf("#### END reloadShip('%s')\n", fname );
};

// ====================== Class Definitions

class SpaceCraftEditorApp : public AppSDL2OGL_3D { public:

	class OnSelectLuaShipScript : public GUIEventCallback{ public:
        virtual int GUIcallback(GUIAbstractPanel* caller) override {
            ((DropDownList*)caller)->selectedToStr(str+sprintf(str,"data/"));
            reloadShip(str);
            return 0;
        }
    };

    int perFrame = 10;

	bool bSmoothLines = 1;
	bool bWireframe   = 1;

	//DropDownList lstLuaFiles;
    GUI gui;
    BoundGUI*  compGui   =0;
    GirderGUI* girderGui =0;
    PlateGUI*  plateGui  =0;
    PickerUI   picker;

    DropDownList* lstLuaFiles=0;
    OnSelectLuaShipScript onSelectLuaShipScript;

    EDIT_MODE edit_mode = EDIT_MODE::component;
    //EDIT_MODE edit_mode = EDIT_MODE::vertex;
    //EDIT_MODE edit_mode = EDIT_MODE::edge;
    int picked = -1;
    Vec3d mouse_ray0;

    int picked_block = -1;

    //void selectCompGui();

    // https://stackoverflow.com/questions/29145476/requiring-virtual-function-overrides-to-use-override-keyword

	virtual void draw   () override;
	virtual void drawHUD() override;
	//virtual void mouseHandling( );
	//virtual void camera();
	virtual void eventHandling   ( const SDL_Event& event  ) override;
	virtual void keyStateHandling( const Uint8 *keys ) override;
    //virtual void mouseHandling( );

	SpaceCraftEditorApp( int& id, int WIDTH_, int HEIGHT_ , int argc, char *argv[]);

};

SpaceCraftEditorApp::SpaceCraftEditorApp( int& id, int WIDTH_, int HEIGHT_, int argc, char *argv[] ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    //Lua1.init();
    fontTex       = makeTexture    ( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    GUI_fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );
    Draw::fontTex = fontTex;

    plateGui  = (PlateGUI* )gui.addPanel( new PlateGUI ( WIDTH-105, 5, WIDTH-5, fontSizeDef*2+2) );
    girderGui = (GirderGUI*)gui.addPanel( new GirderGUI( WIDTH-105, 5, WIDTH-5, fontSizeDef*2+2) );

    //truss.loadXYZ(  "data/octShip.xyz" );
    //DropDownList* lstLuaFiles = new DropDownList( "lua files",20,HEIGHT_-100,200,5); gui.addPanel(lstLuaFiles);
    lstLuaFiles = new DropDownList( "lua files",20,HEIGHT_-100,200,5); gui.addPanel(lstLuaFiles);
    lstLuaFiles->setCommand([this](GUIAbstractPanel* panel){ onSelectLuaShipScript.GUIcallback(panel); });

    TreeView* tvDir      = new TreeView    ( "DirView",20,HEIGHT_-400,200,20 ); gui.addPanel(tvDir);
    dir2tree(tvDir->root, "data" );
    tvDir->updateLines();

    ogl_asteroide=glGenLists(1);
	glNewList( ogl_asteroide, GL_COMPILE );
    drawAsteroide( 32, 50, 0.1, false );
    //drawAsteroide( 32, 50, 0.1, true );
    glEndList();

    ogl_geoSphere=glGenLists(1);
	glNewList( ogl_geoSphere, GL_COMPILE );
    drawAsteroide( 16, 0, 0.0, true );
    //drawAsteroide( 32, 50, 0.1, true );
    glEndList();

    //std::vector<std::string> luaFiles;
    listDirContaining( "data", ".lua", lstLuaFiles->labels );

    //sim.linSolveMethod = (int)OrbSim::LinSolveMethod::Cholesky;
    sim.linSolveMethod = (int)OrbSim::LinSolveMethod::Force;

    VIEW_DEPTH = 10000.0;
    zoom = 1000.0;

    picker.picker = &sim;   picker.Rmax=10.0;
    theSpaceCraft = new SpaceCraft();
    initSpaceCraftingLua();
    //if(argc<=1)reloadShip( "data/ship_ICF_marksman_2.lua" );
    if(argc<=1){
        //reloadShip( "data/test_materials.lua" ); 
        const char* fname = "data/test_materials.lua";
        if( Lua::dofile(theLua,fname) ){ printf( "ERROR in reloadShip() Lua::dofile(%s) \n", fname ); exit(0); }
        makeTestTruss();
        printf("zoom %g \n", zoom );
        zoom = 10.0;
    }

    //onSelectLuaShipScript.GUIcallback(lstLuaFiles);
    printf( "### SpaceCraftEditorApp() DONE\n" );
}

void SpaceCraftEditorApp::draw(){
    //printf( " ==== frame %i \n", frameCount );
    //glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    //glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
    //glClearColor( 0.0f, 0.0f, 0.0f, 1.0f );
    glClearColor( 0.8f, 0.8f, 0.8f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	//glDisable(GL_DEPTH_TEST);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	glColor3f(1.0,0.5,1.0);
	//if(glo_capsula) glCallList(glo_capsula);

    //Vec3d cam.rot.c;
    //float lightPos   []{ 1.0f, -1.0f, 1.0f, 0.0f  };
    glLightfv( GL_LIGHT0, GL_POSITION,  (float*)&cam.rot.c  );
    //glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    iFrame = frameCount;

	/*
    // to make nice antialiased lines without supersampling buffer
    // see  https://www.opengl.org/discussion_boards/showthread.php/176559-GL_LINE_SMOOTH-produces-bold-lines-%28big-width%29
    // https://www.opengl.org/discussion_boards/showthread.php/170072-GL_POLYGON_SMOOTH-is-it-that-bad
    if(bSmoothLines){
        glEnable (GL_BLEND);
        //glColor4f(0.0,0.0,0.0,0.1);
        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable (GL_LINE_SMOOTH);
        glHint (GL_LINE_SMOOTH_HINT, GL_FASTEST);
        glLineWidth(0.5);
        //glLineWidth(1.5);
	}
	//glHint (GL_LINE_SMOOTH_HINT, GL_NICEST);
    if(bWireframe){
        glEnable(GL_POLYGON_SMOOTH); // THIS WILL MAKE WIREFRAME AS SIDE EFFECT
        //glHint(GL_POLYGON_SMOOTH_HINT, GL_FASTEST);
        //glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
    }
    */

    glDisable(GL_LIGHTING);

    //Draw3D::drawMatInPos( sim.I, sim.cog, (Vec3d){sqrt(1/sim.I.xx),sqrt(1/sim.I.yy),sqrt(1/sim.I.zz)} );

    if(bShipReady)Draw3D::drawMatInPos( Mat3dIdentity, sim.cog, Vec3d{100.,100.,100.} );

    glDisable(GL_CULL_FACE);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    /*
    glLineWidth(1.0); glColor3f(1.0,0.0,0.0);
    sim.cleanForce();
    SpaceCraftControl(0.1);
    sim.evalEdgeVerts();
    renderPointForces( sim.nPoint, sim.points, sim.forces, 1e-3 );
    */

    // Render simulation
    glLineWidth(1.0); 
    runSim( sim, bShipReady?perFrame:1 );
    renderTruss( sim.nBonds, sim.bonds, sim.points, sim.strain, 1000.0 );

    glLineWidth(3.0); glColor3f(0.0,0.5,0.0);
    //printf( "sim.damped_bonds.size() %i \n", sim.damped_bonds.size() );
    for(int i: sim.damped_bonds){
        int2 b = sim.bonds[i];
        Draw3D::drawLine( sim.points[b.x].f, sim.points[b.y].f );
    }
    glLineWidth(1.0);

    glPointSize(5.0);
    glBegin(GL_POINTS);
    for(int i: sim.damped_points){ Draw3D::vertex( sim.points[i].f ); }
    glEnd();
    glPointSize(1.0);

    
    glColor3f(0.0,0.0,0.0);
    if(bShipReady==false)renderPoinSizes( sim.nPoint, sim.points, 1.0 );

    glDisable(GL_DEPTH_TEST);

    picked_block = sim.pick_BBox( picker.ray0, picker.hray, 10000.0, 1 );

    if(picked_block>=0){
        glPointSize(5);
        Quat8d bb = sim.BBs[picked_block];
        glColor3f( 0.0,1.0,0.0 ); Draw3D::drawBBox( bb.lo.f, bb.hi.f );
        glColor3f( 0.0,1.0,1.0 ); renderEdgeBox ( picked_block, sim.edgeBBs,  sim.bonds, sim.points );
        glColor3f( 1.0,0.0,1.0 ); renderPointBox( picked_block, sim.pointChunks,         sim.points );
        glColor3f( 0.0,0.0,1.0 ); renderPointBox( picked_block, sim.pointBBs,            sim.points );
    }

    //pointLabels( mesh.verts.size(), &mesh.verts[0].pos, 0.1, 0.0, fontTex, 10.0, 0 );
    
    glColor3f(0.0,0.5,1.0);
    //drawSliderBonds( sim );  // draw the lines connected to the two endpoints of the actual active edge of she slider path
    drawSliderPaths( *theSpaceCraft, sim  ); // Draw Line connected to interpolated position of slider on its path 

    Draw3D::color( Vec3d{0.0,0.0,1.0} );
    //for(int i=0; i<mesh2.verts.size(); i++){  Draw3D::drawInt( mesh2.verts[i].pos, i, Draw::fontTex, 0.02 );}


    // ----------- Picking
    picker.hray = (Vec3d)(cam.rot.c);
    picker.ray0 = (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y);
    glLineWidth(5.0);
    if     (picker.edit_mode == EDIT_MODE::vertex){ if( picker.picked>=0 ){ Vec3d p = *(Vec3d*)picker.getPickedObject(); glColor3f(0.0,1.0,0.0); Draw3D::drawPointCross( p, 10.0 );                              } }
    else if(picker.edit_mode == EDIT_MODE::edge  ){ if( picker.picked>=0 ){ Vec2i b = *(Vec2i*)picker.getPickedObject(); glColor3f(0.0,1.0,0.0); Draw3D::drawLine      ( sim.points[b.x].f, sim.points[b.y].f ); } }
    //mouse_ray0 = (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y);
    //printf( "%i\n", EDIT_MODE::vertex );
    // if(picked>=0){
    //     switch(edit_mode){
    //         case EDIT_MODE::vertex:
    //             glColor3f(1.0,1.0,1.0); Draw3D::drawPointCross( truss.points[picked], 0.3 );
    //             if ( SDL_GetMouseState(NULL, NULL) & SDL_BUTTON(SDL_BUTTON_LEFT) ){ Draw3D::drawLine(truss.points[picked], mouse_ray0); }
    //             break;
    //         case EDIT_MODE::edge  : glColor3f(1.0,1.0,1.0); auto ed = truss.edges[picked]; Draw3D::drawLine( truss.points[ed.a], truss.points[ed.b] ); break;
    //     }
    // }

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
    if( (picker.picked>=0) && (edit_mode==EDIT_MODE::component) ){
        glColor3f(0,1.0,0);
        drawPicked( *theSpaceCraft, picked );
    }
};

void SpaceCraftEditorApp::drawHUD(){
    glDisable( GL_LIGHTING );
    glDisable(GL_DEPTH_TEST);

    // void Draw::drawText( const char * str, int itex, float sz, Vec2i block_size ){

    sprintf(str_tmp, "time=%10.5f[s] mass=%g cog(%g,%g,%g) vcog(%g,%g,%g) L(%g,%g,%g) torq(%g,%g,%g) |F|=%g \n", sim.time, sim.mass, sim.cog.x,sim.cog.y,sim.cog.z, sim.vcog.x,sim.vcog.y,sim.vcog.z, sim.L.x,sim.L.y,sim.L.z, sim.torq.x,sim.torq.y,sim.torq.z, sim.F_residual );
    //sprintf( str_tmp, "time= %10.5f[s] \n ", sim.time );
    Draw::drawText( str_tmp, fontTex, fontSizeDef,  {WIDTH,HEIGHT-20}  );

    gui.draw();
    //glColor3f(1.0f,1.0f,1.0f);   txtStatic.view3D( {5,5}, fontTex, 8 );
    //glPopMatrix();
}

void SpaceCraftEditorApp::keyStateHandling( const Uint8 *keys ){
    //Mat3d camMat;
    //qCamera.toMatrix_T(camMat);
    //Draw3D::drawMatInPos(  (Mat3f)camMat, (Vec3d){0.0,0.0,0.0} );
    //qCamera.toMatrix(camMat);
    if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.dyaw  (  keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.dyaw  ( -keyRotSpeed ); }
	//if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.dpitch(  keyRotSpeed ); }
	//if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.dpitch( -keyRotSpeed ); }
    if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.roll(  keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.roll( -keyRotSpeed ); }

    if( keys[ SDL_SCANCODE_W ] ){ cam.pos.add_mul( cam.rot.b, +0.05*zoom ); }
	if( keys[ SDL_SCANCODE_S ] ){ cam.pos.add_mul( cam.rot.b, -0.05*zoom );  }
	if( keys[ SDL_SCANCODE_A ] ){ cam.pos.add_mul( cam.rot.a, -0.05*zoom );  }
	if( keys[ SDL_SCANCODE_D ] ){ cam.pos.add_mul( cam.rot.a, +0.05*zoom );  }
    if( keys[ SDL_SCANCODE_LEFTBRACKET  ] ){ theSpaceCraft->nodes[7]->calong-=0.001; theSpaceCraft->nodes[7]->updateBound(); }
    if( keys[ SDL_SCANCODE_RIGHTBRACKET ] ){ theSpaceCraft->nodes[7]->calong+=0.001; theSpaceCraft->nodes[7]->updateBound(); }

    wheel_speed = Vec3dZero;
	if( keys[ SDL_SCANCODE_KP_5 ] ){ wheel_speed.y=-wheel_speed_setup.y; }
    if( keys[ SDL_SCANCODE_KP_8 ] ){ wheel_speed.y= wheel_speed_setup.y; }
    if( keys[ SDL_SCANCODE_KP_4 ] ){ wheel_speed.x=-wheel_speed_setup.x; }
	if( keys[ SDL_SCANCODE_KP_6 ] ){ wheel_speed.x= wheel_speed_setup.x; }
	if( keys[ SDL_SCANCODE_KP_7 ] ){ wheel_speed.z=-wheel_speed_setup.z; }
	if( keys[ SDL_SCANCODE_KP_9 ] ){ wheel_speed.z= wheel_speed_setup.z; }

};

void SpaceCraftEditorApp::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    if(event.type == SDL_MOUSEWHEEL){
        if     (event.wheel.y > 0){ zoom/=VIEW_ZOOM_STEP; }
        else if(event.wheel.y < 0){ zoom*=VIEW_ZOOM_STEP; }
    }
    gui.onEvent(mouseX,mouseY,event);
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_m:  edit_mode = (EDIT_MODE)((((int)edit_mode)+1)%((int)EDIT_MODE::size)); printf("edit_mode %i\n", (int)edit_mode); break;
                case SDLK_LEFTBRACKET  : circ_inc(picked_block, sim.edgeBBs.ncell ); break;
                case SDLK_RIGHTBRACKET : circ_dec(picked_block, sim.edgeBBs.ncell ); break;
                case SDLK_m: picker.switch_mode(); break;
                //case SDLK_h:  warrior1->tryJump(); break;
                case SDLK_l:
                    //reloadShip( );
                    onSelectLuaShipScript.GUIcallback(lstLuaFiles);
                    break;
                case SDLK_SPACE: bRun = !bRun; break;
                case SDLK_KP_0:  sim.cleanVel(); break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    /*
                    switch(edit_mode){
                        case EDIT_MODE::vertex    : picked = truss.pickVertex( mouse_ray0, (Vec3d)cam.rot.c, 0.5  ); printf("picked %i\n", picked); break;
                        case EDIT_MODE::edge      : picked = truss.pickEdge  ( mouse_ray0, (Vec3d)cam.rot.c, 0.25 ); printf("picked %i\n", picked); break;
                        case EDIT_MODE::component :
                            picked = theSpaceCraft->pick( mouse_ray0, (Vec3d)cam.rot.c, 3.0 );
                            //printf("picked %i %i\n", picked, theSpaceCraft->pickedTyp );
                            selectCompGui();
                            //printf("onEvent %li %li \n", girderGui, compGui );
                            if( compGui ){
                                if(picked>=0){ compGui->bindLoad( theSpaceCraft->getPicked(picked) ); }
                                //else         { compGui->unbind( ); };
                            }
                            break;
                    }; break;
                    */
                   picker.pick();
                case SDL_BUTTON_RIGHT: break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    /*
                    switch(edit_mode){
                        case EDIT_MODE::vertex: int ip2 = truss.pickVertex( mouse_ray0, (Vec3d)cam.rot.c, 0.5  ); if((picked>=0)&(ip2!=picked)); truss.edges.push_back((TrussEdge){picked,ip2,0}); break;
                        //case EDIT_MODE::edge  : picked = truss.pickEdge  ( mouse_ray0, camMat.c, 0.25 ); printf("picked %i\n", picked); break;
                    }; break;
                    */
                case SDL_BUTTON_RIGHT:break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );

    //printf( "compGui %li \n", compGui );
    if(compGui) if( compGui->check() ){ renderShip(); }

}

// ===================== MAIN

int main(int argc, char *argv[]){    
    printf( "argc %i \n", argc );
    SDL_DisplayMode dm = initSDLOGL( 8 );
	int junk;
	SpaceCraftEditorApp * app = new SpaceCraftEditorApp( junk , dm.w-150, dm.h-100, argc, argv );
    //app->bindSimulators( &W ); 

    LambdaDict funcs;
    funcs["-s"]={1,[&](const char** ss){ reloadShip( ss[0] ); }}; 
    funcs["-perframe"]={1,[&](const char** ss){            sscanf( ss[0], "%i", &app->perFrame );              printf( "COMMAND LINE: -perframe(%i) \n", app->perFrame ); } };
    funcs["-method"  ]={1,[&](const char** ss){ int im;    sscanf( ss[0], "%i", &im );  sim.linSolveMethod=im; printf( "COMMAND LINE: -method(%i)   \n", im            ); } };
    funcs["-dt"      ]={1,[&](const char** ss){ float dt;  sscanf( ss[0], "%f", &dt );  sim.dt=dt;             printf( "COMMAND LINE: -dt( dt: %f ) \n", sim.dt        ); } };
    
    //funcs["-bmix"    ]={1,[&](const char** ss){ int istart; float bmix;  sscanf( ss[0], "%i,%f", &istart, &bmix ); W.sim.mixer.b_end=bmix; W.sim.mixer.istart=istart; printf( "COMMAND LINE: -bmix( istart:%i bmix: %f ) \n", W.sim.mixer.istart, W.sim.mixer.b_end );    } };
    //funcs["-fix"     ]={1,[&](const char** ss){ int n =  readlist( ss[0], W.fixPoints); printf("COMMAND LINE: -fix[%i]{%s}\n", n, ss[0] );  } };
    //funcs["-nsolve"  ]={1,[&](const char** ss){ int nsolv; sscanf( ss[0], "%i", &nsolv ); printf( "COMMAND LINE: -nsolve(%i) \n", nsolv ); W.sim_f.nSolverIters=nsolv; W.sim.nSolverIters=nsolv;  } };    

    process_args( argc, argv, funcs );

	//app = new SpaceCraftEditorApp( junk , 800, 600 );
	app->loop( 1000000 );
	return 0;
}
