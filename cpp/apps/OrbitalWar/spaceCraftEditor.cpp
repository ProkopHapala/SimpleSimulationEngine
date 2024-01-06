
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

int verbosity = 0;


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

#include "OrbSim.h"
#include "TriangleRayTracer.h"
#include "Radiosity.h"

#include "Tree.h"

#include "spaceCraftEditorUtils.h"

#include "AppSDL2OGL_3D.h"
#include "GUI.h"
#include "SpaceCraftGUI.h"

#include "argparse.h"

// TODO:
// Collision Detection:
//   Try Collisions using [kBoxes.h]    tested in    [test_BoxAndSweep.cpp]
// 

// ======================  Global Variables & Declarations

using namespace SpaceCrafting;

//SpaceCraft craft;
//Truss      truss;

bool bShipReady = false;
Mesh::Builder2 mesh2;
OrbSim_f sim;
int glo_truss=0, glo_capsula=0, glo_ship=0;
char str_tmp[8096];
double elementSize  = 5.;
bool bRun = false;

Vec3d wheel_speed       = {0.0,0.0,0.0};
//Vec3d wheel_speed_setup = { 0.1, 0.1, 0.1 };
//Vec3d wheel_speed_setup = { 0.5, 0.5, 0.5 };
Vec3d wheel_speed_setup = { 5.0, 5.0, 5.0 };

void SpaceCraftControl(double dt){
    // wheel_speed
    //printf( "SpaceCraftControl() wheel_speed(%g,%g,%g)\n", wheel_speed.x, wheel_speed.y, wheel_speed.z );
    // --- move the sliders & update corresponding EdgeVerts
    for( int i=0; i<theSpaceCraft->sliders.size(); i++ ){
        Slider* o = theSpaceCraft->sliders[i];

        int icon = o->icontrol;
        if(icon>=0){ 
            if(icon<3){ o->speed = wheel_speed.array[icon]; }
        }

        EdgeVertBond& ev = sim.edgeVertBonds[i];

        Vec3f  d  = sim.points[ev.verts.y].f - sim.points[ev.verts.x].f;
        Vec3f  dv = sim.vel[ev.verts.z].f - sim.vel[ev.verts.y].f*ev.c + sim.vel[ev.verts.x].f*(1-ev.c);
        double l = d.norm();
        float  f = d.dot( ev.f )/l; // force along the slider path
        float  v = d.dot( dv   )/l; // velocity along the slider path

        o->move( dt, l, v, f );
        
        // -- update corresponding EdgeVerts
        
        ev.c = o->path.fromCur( ev.verts.x, ev.verts.y );
        ev.verts.z = o->ivert;
        //ev.K = 1000.0;
        //o->updateEdgeVerts( sim.points );

    }
}


void runSim( OrbSim_f& sim, int niter=100 ){
    long t0 = getCPUticks();
    if(bRun){
        /*
        for(int itr=0; itr<niter; itr++){
            sim.cleanForce();   
            //sim.evalTrussForces_neighs();
            sim.evalTrussForces_neighs2();
            //sim.evalTrussForces_bonds();
            //sim.evalTrussCollisionImpulses_bonds(0.5);
            //sim.evalTrussCollisionImpulses_bonds(1.0);
            //sim.evalTrussCollisionImpulses_bonds(0.1);
            sim.applyForceCentrifug( {0.,0.,0.}, {0.0,0.0,1.0}, 5e-2 );
            //sim.applyForceRotatingFrame( {0.,0.,0.}, {0.0,0.0,1.0}, 0.1 );
            //sim.move_GD( 0.00001 );
            //sim.move_MD( 1e-3, 1e-5 );
            //sim.move_MD( 1e-3, 1e-8 );
            sim.move_MD( 1e-3, 1e-8 );
            //sim.move_GD( 1e-7 );
        }
        */
        //sim.run( niter, 1e-3, 1e-8 );
        //sim.run( niter, 1e-3, 1e-4 );
        //sim.run_omp( niter, false, 1e-3, 1e-4 );
        //sim.run_omp( niter, true, 1e-3, 1e-5 );
        //sim.run_omp( niter, true, 1e-3, 1e-4 );
        sim.run_omp( niter, true, 1e-3, 1e-6 );
        //sim.run_omp( 1, true, 1e-3, 1e-4 );
        double T = (getCPUticks()-t0)*1e-6;
        //printf( "runSim() DONE T=%g[ms] %g[ms/iter] niter=%i,nP=%i,nE=%i \n", T, T/niter, niter, sim.nPoint, sim.nNeighMax );
        //printf( "runSim() DONE T=%g[ms] %g[ms/iter] niter=%i,nP=%i,nE=%i \n", T, T/niter, niter, sim.nPoint, sim.nNeighMax );
        sim.updateInveriants(false);
    }
    
    sim.evalBondTension();
    //renderPoinSizes( sim.nPoint, sim.points, 0.001 );
    //renderPointForces( sim.nPoint, sim.points, sim.forces, 1e-6 );
    renderPointForces( sim.nPoint, sim.points, sim.forces, 1e-3 );
    //renderPointForces( sim.nPoint, sim.points, sim.forces, 1e-4 );
    //renderPointForces( sim.nPoint, sim.points, sim.forces, 1.0 );

    if(sim.pointBBs.ncell>0) updatePointBBs( sim.pointBBs, sim.BBs, sim.points,            true );  // It crashes here because of the wrong obj2cell mapping
    if(sim.edgeBBs .ncell>0) updateEdgeBBs ( sim.edgeBBs,  sim.BBs, sim.bonds, sim.points, false );
    
}

// ======================  Free Functions

void renderEdgeBox( int ib, Buckets& buckets, int2* edges, Quat4f* points){
    //const Quat8f& bb = BBs[ib];
    //Draw3D::drawBBox( bb.lo.f, bb.hi.f );
    int i0 = buckets.cellI0s[ib];
    int n  = buckets.cellNs [ib];
    //printf( "---- renderEdgeBox[%i] n=%i i0=%i \n", ib, n, i0 );
    for(int i=0; i<n; i++){
        int ie = buckets.cell2obj[i+i0];
        int2& e = edges[ie];
        //printf( "renderEdgeBox[%i] ie %i e(%i,%i) \n", ib, ie, e.x, e.y );
        Draw3D::drawLine( points[e.x].f, points[e.y].f );
    }
}

void renderPointBox( int ib, Buckets& buckets, Quat4f* points){
    //const Quat8f& bb = BBs[ib];
    //Draw3D::drawBBox( bb.lo.f, bb.hi.f );
    int i0 = buckets.cellI0s[ib];
    int n  = buckets.cellNs [ib];
    //printf( "---- renderEdgeBox[%i] n=%i i0=%i \n", ib, n, i0 );
    glBegin(GL_POINTS);
    for(int i=0; i<n; i++){
        int ip = buckets.cell2obj[i+i0];
        //printf( "renderEdgeBox[%i] ie %i e(%i,%i) \n", ib, ie, e.x, e.y );
        Draw3D::vertex( points[ip].f );   
    }
    glEnd();
}

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
    drawSpaceCraft_Mesh( *theSpaceCraft, mesh, 1, false, true, (Vec3f){0.5,0.5,0.5} );
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
    Draw3D::color(Vec3f{1.0f,0.0f,1.0f});
    for(int i=0; i<theSpaceCraft->sliders.size(); i++){
        const Slider& o = theSpaceCraft->sliders[i];
        Draw3D::drawLineStrip( o.path.n, o.path.ps, sim.points, o.path.closed );
    }
    */
    
    //Draw3D::color(Vec3f{1.0,0.,1.});
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


void makeBBoxes( const SpaceCraft& craft, OrbSim_f& sim ){
    // // --- Bounding boxes
    int nbuck  = exportBuckets( craft,             0, 16, true );
    //printf( "nbuck %i \n", nbuck );
    sim.recallocBBs( nbuck );
    sim.pointBBs.cleanO2C(0);  // by default we assign all points to cell 0
    int nbuck_ = exportBuckets( craft, &sim.pointBBs, 16, true );
    //sim.pointBBs.printCells();
    //sim.pointBBs.printObjCellMaping(0,100);   // see if obj2cell is incorrect at the beginning
    sim.pointBBs.updateCells();
    //printf("### pointBBs.printCells() \n"); sim.pointBBs.printCells();
    sim.edgesToBBs();
    sim.edgeBBs.updateCells();
    //printf("### edgeBBs.printCells() \n"); sim.edgeBBs.printCells();
    updatePointBBs( sim.pointBBs, sim.BBs, sim.points,            true );  // It crashes here because of the wrong obj2cell mapping
    updateEdgeBBs ( sim.edgeBBs,  sim.BBs, sim.bonds, sim.points, false );
    //updateEdgeBBs ( sim.edgeBBs, sim.BBs, sim.bonds, sim.points, true );
    //sim.printBBs();
}

void makeTestTruss(){
    mesh2.clear();
    int stickType = 0;
    mesh2.stick  ( {0.0,0.0,0.0}, {10.0,0.0,0.0}, stickType );
    mesh2.stickTo( 0,             {5.0,10.0,0.0}, stickType );
    mesh2.edge   ( 1, 2, stickType );
    workshop.stickMaterials.vec[stickType]->preStrain = 0.01;

    exportSim( sim, mesh2, workshop );
    sim.updateInveriants(true);

    bShipReady = false;
    renderShip();
}

void reloadShip( const char* fname  ){
    theSpaceCraft->clear();                  // clear all components
    //luaL_dofile(theLua, "data/spaceshil1.lua");
    printf("#### START reloadShip('%s')\n", fname );
    //Lua::dofile(theLua,fname);
    if( Lua::dofile(theLua,fname) ){ printf( "ERROR in reloadShip() Lua::dofile(%s) \n", fname ); exit(0); }
    //printf( "Lua::dofile(%s) DONE \n", fname );
    theSpaceCraft->checkIntegrity();
    //luaL_dostring(theLua, "print('LuaDEBUG 1'); n1 = Node( {-100.0,0.0,0.0} ); print('LuaDEBUG 2'); print(n1)");
    //theSpaceCraft->toTruss();
    //truss.clear();
    //toTruss(*theSpaceCraft, truss);
    long t0;

    t0 = getCPUticks();
    mesh2.clear();
    BuildCraft_truss( mesh2, *theSpaceCraft, 30.0 );
    printf( "BuildCraft_truss() DONE T=%g[ms] \n", (getCPUticks()-t0)*1e-6 );
    mesh2.printSizes();

    t0 = getCPUticks();
    exportSim( sim, mesh2, workshop );

    
    printf( "exportSim() DONE T=%g[ms] \n", (getCPUticks()-t0)*1e-6 );

    t0 = getCPUticks();
    makeBBoxes( *theSpaceCraft, sim );
    makePointCunks( sim.edgeBBs, sim.bonds, sim.pointChunks );
    sim.pointChunks.printCells();
    printf( "makeBBoxes() DONE T=%g[ms] \n", (getCPUticks()-t0)*1e-6 );

    sim.user_update = SpaceCraftControl;

    printf( "reloadShip().updateSlidersPaths \n" );
    // update ring slider paths
    for( Ring* o : theSpaceCraft->rings ){
        o->updateSlidersPaths( true, true, sim.points );
    }
    printf( "reloadShip().SlidersToEdgeVerts \n" );
    // ----- Sliders to EdgeVerts
    sim.nEdgeVertBonds =  theSpaceCraft->sliders.size();
    sim.edgeVertBonds  = new EdgeVertBond[ sim.nEdgeVertBonds ];
    for( int i=0; i<theSpaceCraft->sliders.size(); i++ ){
        Slider* o = theSpaceCraft->sliders[i];
        EdgeVertBond& ev = sim.edgeVertBonds[i];
        ev.c = o->path.fromCur( ev.verts.x, ev.verts.y );
        ev.verts.z = o->ivert;
        //ev.K = 1000.0;
        ev.K = 100000.0;
        //ev.K = 0.0;

        o->speed = 1.0;
        //o->speed = 0.1;
        //o->updateEdgeVerts( sim.points );
    }

    // --- kick to rings
    // Vec3f v0 = {1.0,0.0,0.0};
    // for( Ring* o : theSpaceCraft->rings ){
    //     for(int iv=o->pointRange.a; iv<o->pointRange.b; iv++){
    //         sim.vel[iv].f = v0; 
    //     }
    // }

    renderShip();

    sim.updateInveriants(true);
    bShipReady = true;
    printf("#### END reloadShip('%s')\n", fname );
};

/*
template<typename T>
class Driver{ public:
    T* master;
    T* slave;
    void bind(T* master_, T* slave_ ){ master=master_; slave=slave_; };
    void apply(){ slave=master; }
};
*/

// ====================== Class Definitions




class SpaceCraftEditorApp : public AppSDL2OGL_3D { public:

	class OnSelectLuaShipScript : public GUIEventCallback{ public:
        virtual int GUIcallback(GUIAbstractPanel* caller) override {
            ((DropDownList*)caller)->selectedToStr(str+sprintf(str,"data/"));
            reloadShip(str);
            return 0;
        }
    };

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
    //lstLuaFiles->onSelect = new OnSelectLuaShipScript();
    lstLuaFiles->onSelect = &onSelectLuaShipScript;

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

    /*
    Vec3f pf;
    Vec3d pd;
    //pd.set(1.5,1.800001,1.9);
    pd.set(1.5,1.800000001,1.9);
    //pd = pf;
    //pd = (Vec3d)pf;
    pf = (Vec3f)pd;
    print(pf);printf(" // pf\n");
    print(pd);printf(" // pd\n");
    //exit(0);
    */

    //camera();

    picker.picker = &sim;   picker.Rmax=10.0;
    theSpaceCraft = new SpaceCraft();
    initSpaceCraftingLua();
    //if(argc<=1)reloadShip( "data/ship_ICF_marksman_2.lua" );
    if(argc<=1){
        reloadShip( "data/test_materials.lua" ); 
        makeTestTruss();
    }

    //onSelectLuaShipScript.GUIcallback(lstLuaFiles);

    /*
    glo_truss = makeTruss(truss);

    glo_capsula = glGenLists(1);
    glNewList( glo_capsula, GL_COMPILE );
    //Draw3D::drawCylinderStrip  ( 16, 10, 10, {0.0,0.0,-16.0,}, {0.0,0.0,16} );
    Draw3D::drawCapsula( (Vec3f){0.0,0.0,-1.0}, (Vec3f){0.0,0.0,1.0}, 2.0, 1.0, 0.7, 0.7, 0.2, 32, false );
    glEndList();
    //delete [] ups;
    */


    VIEW_DEPTH = 10000.0;
    zoom = 1000.0;
    printf( "### SpaceCraftEditorApp() DONE\n" );
}

//void SpaceCraftEditorApp::camera(){
//    camera_FreeLook( camPos );
//}

/*
void SpaceCraftEditorApp::selectCompGui(){
    if(compGui)compGui->close();
    switch( (ComponetKind)theSpaceCraft->pickedTyp ){
        case ComponetKind::Radiator:
        case ComponetKind::Shield:
            compGui=plateGui;
            break;
        case ComponetKind::Girder:
            compGui=girderGui;
            //printf("setGirderGui %li %li \n", girderGui, compGui );
            break;
        //case ComponetKind::Rope:
        default:
            compGui=0;
            break;
    }
    if(compGui)compGui->open();
}
*/

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

    //Vec3f cam.rot.c;
    //float lightPos   []{ 1.0f, -1.0f, 1.0f, 0.0f  };
    glLightfv( GL_LIGHT0, GL_POSITION,  (float*)&cam.rot.c  );
    //glLightfv(GL_LIGHT0, GL_POSITION, light_position);


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

    //Draw3D::drawMatInPos( sim.I, sim.cog, (Vec3f){sqrt(1/sim.I.xx),sqrt(1/sim.I.yy),sqrt(1/sim.I.zz)} );

    if(bShipReady)Draw3D::drawMatInPos( Mat3fIdentity, sim.cog, Vec3f{100.,100.,100.} );

    /*
    glLineWidth(3.0);
    Draw3D::color( Vec3f{1.0,0.0,1.0} );
    //drawSpaceCraft_sliderPaths( *theSpaceCraft, sim.points );
    //drawSpaceCraft_nodes( theSpaceCraft, const Quat4f* ps=0, float sz=10, int i0=0, int i1=-1 );
    Draw3D::color(Vec3f{0.f,0.5f,0.f});  drawSpaceCraft_nodes( *theSpaceCraft, 0, 10 );
    glLineWidth(5.0);
    Draw3D::color(Vec3f{0.5f,0.5f,0.f}); drawSpaceCraft_nodes( *theSpaceCraft, sim.points, 5 );
    //printf( "nodes.size() = %i \n", theSpaceCraft->nodes.size() );
    
    for(SpaceCrafting::Node* nd: theSpaceCraft->nodes ){
        //printf( "node[%i] %li \n", nd->id, (long)(nd->boundTo) );
        //printf( "node[%i] \n", nd->id );
        if(nd->boundTo==0) continue;
        //printf( "node[%i] %li bt.id=%i\n", nd->id, (long)(nd->boundTo), nd->boundTo->id );
        if((nd->boundTo->id>1000)||(nd->boundTo->id<0)){ printf("ERROR node[%id]->boundTo->id==%i\n", nd->id, nd->boundTo->id ); exit(0); }
        //Node* nd = theSpaceCraft->nodes[7];
        Draw3D::color(Vec3f{0.f,0.5f,0.5f});  Draw3D::drawPointCross( nd->boundTo->nodes.x->pos*(1-nd->calong) + nd->boundTo->nodes.y->pos*nd->calong, 5 );
        Mat3d rot;
        nd->boundTo->rotMat( rot );
        Vec3d pos = nd->boundTo->nodes.x->pos*(1-nd->calong) + nd->boundTo->nodes.y->pos*nd->calong;
        Draw3D::drawMatInPos( rot, pos, Vec3dOne*10.0 );
        glLineWidth(3.0);
        Draw3D::color(Vec3f{1.f,1.0f,1.0f});
        Draw3D::drawVecInPos( ((Girder*)nd->boundTo)->up*30.0, pos );
    }
    */

    // { // Draw Radiators anchor points
    //     glPointSize(10);
    //     const Radiator& o =  *theSpaceCraft->radiators[0];
    //     const Girder& g1  =  *theSpaceCraft->girders[o.g1];
    //     const Girder& g2  =  *theSpaceCraft->girders[o.g2];
    //     Draw3D::color(Vec3f{1.f,0.f,0.f}); drawPointRange( 10, g1.pointRange, 4, 0, o.g1span, sim.points );
    //     Draw3D::color(Vec3f{0.f,0.f,1.f}); drawPointRange( 10, g2.pointRange, 4, 1, o.g2span, sim.points );
    // }
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
    runSim( sim, bShipReady?100:1 );
    renderTruss( sim.nBonds, sim.bonds, sim.points, sim.strain, 1000.0 );
    glColor3f(0.0,0.0,0.0);
    if(bShipReady==false)renderPoinSizes( sim.nPoint, sim.points, 1.0 );

    // --- render bounding boxes
    // glColor3f(0.0,0.5,0.0);
    // for(int i=0; i<sim.nBBs; i++){
    //     const Quat8T<float> bb = sim.BBs[i];
    //     Draw3D::drawBBox( bb.lo.f, bb.hi.f );
    // }

    glDisable(GL_DEPTH_TEST);

    picked_block = sim.pick_BBox( picker.ray0, picker.hray, 10000.0, 1 );

    if(picked_block>=0){
        glPointSize(5);
        Quat8f bb = sim.BBs[picked_block];
        glColor3f( 0.0,1.0,0.0 ); Draw3D::drawBBox( bb.lo.f, bb.hi.f );
        glColor3f( 0.0,1.0,1.0 ); renderEdgeBox ( picked_block, sim.edgeBBs,  sim.bonds, sim.points );
        glColor3f( 1.0,0.0,1.0 ); renderPointBox( picked_block, sim.pointChunks,         sim.points );
        glColor3f( 0.0,0.0,1.0 ); renderPointBox( picked_block, sim.pointBBs,            sim.points );
    }

    // --- draw nodes
    // glPointSize(10);
    // glColor3f(1.0,0.0,1.0);
    // glBegin(GL_POINTS);
    // for(int i=0; i<theSpaceCraft->nodes.size(); i++){ 
    //     theSpaceCraft->nodes[i]->pos = (Vec3d)sim.points[ theSpaceCraft->nodes[i]->ivert ].f;
    //     Draw3D::vertex( theSpaceCraft->nodes[i]->pos );
    //     //Draw3D::vertex( sim.points[i].f ); 
    // }
    // for(int i=0; i<theSpaceCraft->nodes.size(); i++){   // labels
    //     Draw3D::drawText( std::to_string(i).c_str(), (Vec3f)theSpaceCraft->nodes[i]->pos, fontTex, 0.04, 0 );
    //     //Draw3D::vertex( sim.points[i].f ); 
    // }
    // glEnd();

    // --- draw slider bonds
    glLineWidth(3.0);
    // glColor3f(0.0,0.5,0.0);
    // // render EdgeVertBonds
    // glBegin(GL_LINES);
    // for(int i=0; i<sim.nEdgeVertBonds; i++){
    //     EdgeVertBond& ev = sim.edgeVertBonds[i];
    //     Vec3f a = sim.points[ ev.verts.x ].f;
    //     Vec3f b = sim.points[ ev.verts.y ].f;
    //     Vec3f c = sim.points[ ev.verts.z ].f;
    //     Draw3D::vertex( a );  Draw3D::vertex( c );
    //     Draw3D::vertex( b );  Draw3D::vertex( c );
    // }
    // glEnd();
    // --- draw slider paths
    glLineWidth(3.0);
    
    //for( const Slider* o: theSpaceCraft->sliders){ 
    for( int i=0; i<theSpaceCraft->sliders.size(); i++ ){
        const Slider* o = theSpaceCraft->sliders[i];
        // glColor3f(0.0,0.5,1.0);
        // drawSliderPath( o, sim.points, 10 ); 
        
        glColor3f(1.0,0.0,1.0);
        Vec3f p  = sim.getEdgeVertPos( i );
        Vec3f p0 = sim.points[ o->ivert ].f;
        Draw3D::drawLine( p0, p );
        
    }
    //glColor3f(0.0,0.5,1.0); glPointSize(20); glBegin(GL_POINTS); for( const Slider* o: theSpaceCraft->sliders){ Draw3D::vertex( sim.points[o->ivert].f ); } glEnd();
    //glPointSize(3000/zoom); glBegin(GL_POINTS); for( const Slider* o: theSpaceCraft->sliders){ Draw3D::vertex( sim.points[o->ivert].f ); } glEnd();

    //if(glo_truss) glCallList(glo_truss);
    //if(glo_ship ) glCallList(glo_ship);

    //pointLabels( mesh.verts.size(), &mesh.verts[0].pos, 0.1, 0.0, fontTex, 10.0, 0 );
    

    Draw3D::color( Vec3f{0.0,0.0,1.0} );
    //for(int i=0; i<mesh2.verts.size(); i++){  Draw3D::drawInt( mesh2.verts[i].pos, i, Draw::fontTex, 0.02 );}

    /*
    if(ogl_asteroide){
        glPushMatrix();
        glScalef(100,100.0,100.0);
        glCallList(ogl_asteroide);
        glPopMatrix();
    }
    */

    //Mat3d camMat;
    //qCamera.toMatrix_T(camMat);
    //Draw3D::drawMatInPos( cam.rot, (Vec3f){0.0,0.0,0.0} );

    picker.hray = (Vec3d)(cam.rot.c);
    picker.ray0 = (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y);

    glLineWidth(5.0);
    if     (picker.edit_mode == EDIT_MODE::vertex){ if( picker.picked>=0 ){ Vec3f p = *(Vec3f*)picker.getPickedObject(); glColor3f(0.0,1.0,0.0); Draw3D::drawPointCross( p, 10.0 );                              } }
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

    //glColor3f(0.0f,0.0f,0.0f); drawTruss( truss.edges.size(), &truss.edges[0], &truss.points[0] );
    //glColor3f(1.0f,1.0f,1.0f); Draw3D::drawPoints( truss.points.size(), &truss.points[0], 0.1 );

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

//void SpaceCraftEditorApp::keyStateHandling( const Uint8 *keys ){ };
/*
void SpaceCraftEditorApp::mouseHandling( ){
    SDL_GetMouseState( &mouseX, &mouseY ); mouseY=HEIGHT-mouseY;
    mouse_t   = (mouseX-20 )/timeScale + tstart;
    mouse_val = (mouseY-200)/valScale;
    sprintf( curCaption, "%f %f\0", mouse_t, mouse_val );
    int ipoint_ = binSearchFrom<double>(mouse_t,splines.n,splines.ts);
    if( (splines.ts[ipoint_+1]-mouse_t)<(mouse_t-splines.ts[ipoint_]) ) ipoint_++;
    if(ipoint_!=ipoint){
        ipoint=ipoint_;
        char buff[100];
        Vec3d r,v;
        r.set( splines.CPs[0][ipoint],splines.CPs[1][ipoint],splines.CPs[2][ipoint] );
        v.set( splines.getPointDeriv(ipoint,0), splines.getPointDeriv(ipoint,1), splines.getPointDeriv(ipoint,2) );
        sprintf(buff, "%i %f r(%3.3f,%3.3f,%3.3f) v(%3.3f,%3.3f,%3.3f)", ipoint, splines.ts[ipoint], r.x,r.y,r.z,   v.x,v.y,v.z );
        txtStatic.inputText = buff;
    }
    //printf("%i %i %3.3f %3.3f \n", mouseX,mouseY, mouse_t, mouse_val );
}
*/



void SpaceCraftEditorApp::keyStateHandling( const Uint8 *keys ){
    //Mat3d camMat;
    //qCamera.toMatrix_T(camMat);
    //Draw3D::drawMatInPos(  (Mat3f)camMat, (Vec3f){0.0,0.0,0.0} );
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

LambdaDict funcs;
SpaceCraftEditorApp * app;

int main(int argc, char *argv[]){

    printf( "argc %i \n", argc );

    // example: use like : ./SpaceCraftEditorApp -s data/ship_ICF_interceptor_1.lua
    //funcs["-s"]={1,[&](const char** ss){ app->reloadShip( ss[0] ); }}; 
    funcs["-s"]={1,[&](const char** ss){ reloadShip( ss[0] ); }}; 

	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	// https://www.opengl.org/discussion_boards/showthread.php/163904-MultiSampling-in-SDL
	//https://wiki.libsdl.org/SDL_GLattr
	//SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 8);
    //SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 8);
    //SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 8);
    //SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE, 8);
    //SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    //SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    //SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 1);
    //SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 2);
    //SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 4);
    SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 8);
    //SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 16);
    glEnable(GL_MULTISAMPLE);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
    SDL_DisplayMode dm;
    SDL_GetDesktopDisplayMode(0, &dm);
	app = new SpaceCraftEditorApp( junk , dm.w-150, dm.h-100, argc, argv );
    

    process_args( argc, argv, funcs );

	//app = new SpaceCraftEditorApp( junk , 800, 600 );
	app->loop( 1000000 );
	return 0;
}
















