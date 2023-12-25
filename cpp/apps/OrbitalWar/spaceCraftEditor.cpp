
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

#include "Truss.h"
#include "SpaceCraft.h"
#include "SpaceCraftDraw.h"
#include "SpaceCraft2Mesh.h"
#include "SpaceCraft2Mesh2.h"
#include "SoftBody.h"

#include "SphereSampling.h"
#include "DrawSphereMap.h"
#include "Draw3D_Surf.h"

#include "AppSDL2OGL_3D.h"
#include "GUI.h"

#include "IO_utils.h"

#include "EditSpaceCraft.h"

#include "OrbSim.h"
#include "TriangleRayTracer.h"
#include "Radiosity.h"

#include "Tree.h"

#include "spaceCraftEditorUtils.h"

#include "argparse.h"



// TODO:
// Collision Detection:
//   Try Collisions using [kBoxes.h]    tested in    [test_BoxAndSweep.cpp]
// 



// ======================  Global Variables & Declarations

using namespace SpaceCrafting;
enum class EDIT_MODE:int{ vertex=0, edge=1, component=2, size }; // http://www.cprogramming.com/c++11/c++11-nullptr-strongly-typed-enum-class.html

//SpaceCraft craft;
Truss      truss;

Mesh::Builder2 mesh2;
OrbSim_f sim;
int glo_truss=0, glo_capsula=0, glo_ship=0;
//char str[8096];
double elementSize  = 5.;



// Render 

void renderTruss(int nb, int2* bonds, Quat4f* ps, float* strain=0, float sc=1.0 ){
    glBegin(GL_LINES);
    for(int i=0; i<nb; i++ ){
        //printf( "renderTruss()[%i] \n", i );
        int2 b =  bonds[i];
        if(strain){
            float f=strain[i]*sc;
            if(f>0){  Draw3D::color(Vec3f{f,0,0}); }else{ Draw3D::color(Vec3f{0,0,f}); };
        } 
        Draw3D::vertex( ps[b.x].f );
        Draw3D::vertex( ps[b.y].f );
        //printf( "renderTruss[%i](%i,%i) p(%g,%g,%g) p(%g,%g,%g)\n", i, b.x, b.y,  ps[b.x].f.x,ps[b.x].f.y,ps[b.x].f.z,   ps[b.y].f.x,ps[b.y].f.y,ps[b.y].f.z );
    }
    glEnd();
}

// ======================  Free Functions

void renderShip(){
    printf( "spaceCraftEditor.cpp::renderShip() \n" );
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
    mesh2.clear();
    BuildCraft_truss( mesh2, *theSpaceCraft, 30.0 );
    mesh2.printSizes();
    exportSim( sim, mesh2, workshop );
    //sim.printAllNeighs();
    //theSpaceCraft->printAll_girders();
    theSpaceCraft->updateSliderPaths();
    
    glEnable(GL_LIGHTING);
    glEnable     ( GL_LIGHT0           );
    glEnable     ( GL_NORMALIZE        );
    Draw3D::drawMeshBuilder2( mesh2, 0b110, 1, true, true );
    
    /*
    Draw3D::color(Vec3f{1.0f,0.0f,1.0f});
    for(int i=0; i<theSpaceCraft->sliders.size(); i++){
        const Slider2& o = theSpaceCraft->sliders[i];
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

void reloadShip( const char* fname  ){
    theSpaceCraft->clear();                  // clear all components
    //luaL_dofile(theLua, "data/spaceshil1.lua");
    printf("#### START reloadShip('%s')\n", fname );
    Lua::dofile(theLua,fname);
    printf( "Lua::dofile(%s) DONE \n", fname );
    //luaL_dostring(theLua, "print('LuaDEBUG 1'); n1 = Node( {-100.0,0.0,0.0} ); print('LuaDEBUG 2'); print(n1)");
    theSpaceCraft->toTruss();
    renderShip();
    printf("#### END reloadShip('%s')\n", fname );
};

void drawPicked( const SpaceCraft& craft, int ipick ){
    switch( (ComponetKind)craft.pickedTyp){
        case ComponetKind::Radiator :{
            //printf("drawPicked radiator[%i] \n", ipick );
            drawPlateContour( craft.radiators[ipick], &craft.nodes[0], &craft.girders[0], false );
            }break;
        case ComponetKind::Shield:{
            //printf("drawPicked shield[%i] \n", ipick );
            drawPlateContour( craft.shields  [ipick], &craft.nodes[0], &craft.girders[0], false );
            }break;
        case ComponetKind::Girder:{
            const Girder& o = craft.girders[ipick];
            Draw3D::drawLine(craft.nodes[o.p0].pos,craft.nodes[o.p1].pos);
            } break;
        case ComponetKind::Rope:{
            const Rope& o = craft.ropes[ipick];
            Draw3D::drawLine(craft.nodes[o.p0].pos,craft.nodes[o.p1].pos);
            } break;
    }
}

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

class GUIPanelWatcher{ public:
    GUIPanel* master;
    void*   slave;
    bool bInt=false;
    void bind    ( GUIPanel* master_, void* slave_, bool bInt_ ){ master=master_; slave=slave_; bInt=bInt_; };
    void bindLoad( GUIPanel* master_, void* slave_, bool bInt_ ){ bind(master_,slave_,bInt_); load(); master_->redraw=true; };
    void apply(){ if(bInt){ *(int*)slave  =  (int )master->value;   }else{ *(double*)slave = master->value;   };                     };
    void load (){ if(bInt){ master->value = *(int*)slave;           }else{ master->value   = *(double*)slave; }; master->val2text(); };
    bool check(){ if(master&&slave) if( master->redraw ){ apply(); return true; }; return false; }
};

class ComponetGUI:public MultiPanel{ public:
    //static const int np=4;
    bool binded=false;
    GUIPanelWatcher* drivers=0;
    ComponetGUI(const std::string& caption, int xmin, int ymin, int xmax, int dy,int nsub):MultiPanel(caption,xmin,ymin,xmax,dy,nsub){
        opened=false;
    }
    virtual void view()override{ if(opened)MultiPanel::view(); };
    virtual void bindLoad(ShipComponent* o)=0;
    void unbind(){
        binded=false;
        close();
        //opened=false;
        //redraw=true; tryRender();
    }
    bool check(){
        if(!binded) return false;
        bool bChanged=false;
        for(int i=0; i<nsubs; i++){
            bChanged |= drivers[i].check();
        }
        return bChanged;
    }
};

class PlateGUI:public ComponetGUI{ public:
    //static const int np=4;
    PlateGUI(int xmin, int ymin, int xmax, int dy):ComponetGUI("Plate",xmin,ymin,xmax,dy,4){
        opened=false;
        subs[0]->caption="g1.a ";
        subs[1]->caption="g1.b ";
        subs[2]->caption="g2.a ";
        subs[3]->caption="g2.b ";
    }
    void bindLoad(ShipComponent* o_)override{
        Plate* o = (Plate*) o_;
        //printf( "PlateGUI %li \n", o );
        if(drivers==0)drivers=new GUIPanelWatcher[nsubs];
        drivers[0].bindLoad(subs[0],&o->g1span.a,false);
        drivers[1].bindLoad(subs[1],&o->g1span.b,false);
        drivers[2].bindLoad(subs[2],&o->g2span.a,false);
        drivers[3].bindLoad(subs[3],&o->g2span.b,false);
        binded=true;
        open();
        //redraw=true; tryRender();
    }
};

class GirderGUI:public ComponetGUI{ public:
    //static const int np=4;
    GirderGUI(int xmin, int ymin, int xmax, int dy):ComponetGUI("Girder",xmin,ymin,xmax,dy,2){
        opened=false;
        subs[0]->caption="nseg ";
        subs[1]->caption="mseg ";
    }
    void bindLoad(ShipComponent* o_)override{
        Girder* o = (Girder*) o_;
        //printf( "GirderGUI %li \n", o );
        if(drivers==0)drivers=new GUIPanelWatcher[nsubs];
        drivers[0].bindLoad(subs[0],&o->nseg,true);
        drivers[1].bindLoad(subs[1],&o->mseg,true);
        binded=true;
        open();
        //redraw=true; tryRender();
    }
};


class SpaceCraftEditGUI : public AppSDL2OGL_3D { public:

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
    ComponetGUI*  compGui   =0;
    GirderGUI*  girderGui =0;
    PlateGUI*   plateGui  =0;

    DropDownList* lstLuaFiles=0;
    OnSelectLuaShipScript onSelectLuaShipScript;

    EDIT_MODE edit_mode = EDIT_MODE::component;
    //EDIT_MODE edit_mode = EDIT_MODE::vertex;
    //EDIT_MODE edit_mode = EDIT_MODE::edge;
    int picked = -1;
    Vec3d mouse_ray0;

    void selectCompGui();

    // https://stackoverflow.com/questions/29145476/requiring-virtual-function-overrides-to-use-override-keyword

	virtual void draw   () override;
	virtual void drawHUD() override;
	//virtual void mouseHandling( );
	//virtual void camera();
	virtual void eventHandling   ( const SDL_Event& event  ) override;
	virtual void keyStateHandling( const Uint8 *keys ) override;
    //virtual void mouseHandling( );

	SpaceCraftEditGUI( int& id, int WIDTH_, int HEIGHT_ , int argc, char *argv[]);

};

SpaceCraftEditGUI::SpaceCraftEditGUI( int& id, int WIDTH_, int HEIGHT_, int argc, char *argv[] ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

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

    theSpaceCraft = new SpaceCraft();
    initSpaceCraftingLua();
    if(argc<=1)reloadShip( "data/ship_ICF_interceptor_1.lua" );
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
}

//void SpaceCraftEditGUI::camera(){
//    camera_FreeLook( camPos );
//}

void SpaceCraftEditGUI::selectCompGui(){
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


void SpaceCraftEditGUI::draw(){
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
    glLineWidth(0.5);
    renderTruss( sim.nBonds, sim.bonds, sim.points, sim.strain, 1000.0 );
    glLineWidth(3.0);
    Draw3D::color( Vec3f{1.0,0.0,1.0} );
    drawSpaceCraft_sliderPaths( *theSpaceCraft, sim.points );

    
    {
        glPointSize(10);
        const Radiator& o =  theSpaceCraft->radiators[0];
        const Girder& g1  =  theSpaceCraft->girders[o.g1];
        const Girder& g2  =  theSpaceCraft->girders[o.g2];
        Draw3D::color(Vec3f{1.f,0.f,0.f}); drawPointRange( 10, g1.poitRange, 4, 0, o.g1span, sim.points );
        Draw3D::color(Vec3f{0.f,0.f,1.f}); drawPointRange( 10, g2.poitRange, 4, 1, o.g2span, sim.points );
    }

    glDisable(GL_CULL_FACE);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

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

    //printf( "%i\n", EDIT_MODE::vertex );
    if(picked>=0){
        switch(edit_mode){
            case EDIT_MODE::vertex:
                glColor3f(1.0,1.0,1.0); Draw3D::drawPointCross( truss.points[picked], 0.3 );
                if ( SDL_GetMouseState(NULL, NULL) & SDL_BUTTON(SDL_BUTTON_LEFT) ){ Draw3D::drawLine(truss.points[picked], mouse_ray0); }
                break;
            case EDIT_MODE::edge  : glColor3f(1.0,1.0,1.0); auto ed = truss.edges[picked]; Draw3D::drawLine( truss.points[ed.a], truss.points[ed.b] ); break;
        }

    }

    mouse_ray0 = (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y);
    //glColor3f(0.0f,0.0f,0.0f); drawTruss( truss.edges.size(), &truss.edges[0], &truss.points[0] );
    //glColor3f(1.0f,1.0f,1.0f); Draw3D::drawPoints( truss.points.size(), &truss.points[0], 0.1 );


    

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
    if( (picked>=0) && (edit_mode==EDIT_MODE::component) ){
        glColor3f(0,1.0,0);
        drawPicked( *theSpaceCraft, picked );
    }

};

void SpaceCraftEditGUI::drawHUD(){
    glDisable( GL_LIGHTING );
    glDisable(GL_DEPTH_TEST);

    gui.draw();
    //glColor3f(1.0f,1.0f,1.0f);   txtStatic.view3D( {5,5}, fontTex, 8 );
    //glPopMatrix();
}

//void SpaceCraftEditGUI::keyStateHandling( const Uint8 *keys ){ };
/*
void SpaceCraftEditGUI::mouseHandling( ){
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



void SpaceCraftEditGUI::keyStateHandling( const Uint8 *keys ){
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
};

void SpaceCraftEditGUI::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    //if(event.type == SDL_MOUSEWHEEL){
    //    if     (event.wheel.y > 0){ zoom*=VIEW_ZOOM_STEP; }
    //    else if(event.wheel.y < 0){ zoom/=VIEW_ZOOM_STEP; }
    //}

    if(event.type == SDL_MOUSEWHEEL){
        if     (event.wheel.y > 0){ zoom/=VIEW_ZOOM_STEP; }
        else if(event.wheel.y < 0){ zoom*=VIEW_ZOOM_STEP; }
    }
    gui.onEvent(mouseX,mouseY,event);
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_m:  edit_mode = (EDIT_MODE)((((int)edit_mode)+1)%((int)EDIT_MODE::size)); printf("edit_mode %i\n", (int)edit_mode); break;
                //case SDLK_h:  warrior1->tryJump(); break;
                case SDLK_l:
                    //reloadShip( );
                    onSelectLuaShipScript.GUIcallback(lstLuaFiles);
                    break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
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
                case SDL_BUTTON_RIGHT: break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    switch(edit_mode){
                        case EDIT_MODE::vertex: int ip2 = truss.pickVertex( mouse_ray0, (Vec3d)cam.rot.c, 0.5  ); if((picked>=0)&(ip2!=picked)); truss.edges.push_back((TrussEdge){picked,ip2,0}); break;
                        //case EDIT_MODE::edge  : picked = truss.pickEdge  ( mouse_ray0, camMat.c, 0.25 ); printf("picked %i\n", picked); break;
                    }; break;
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
SpaceCraftEditGUI * app;

int main(int argc, char *argv[]){

    printf( "argc %i \n", argc );

    // example: use like : ./spaceCraftEditor -s data/ship_ICF_interceptor_1.lua
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
	app = new SpaceCraftEditGUI( junk , dm.w-150, dm.h-100, argc, argv );
    

    process_args( argc, argv, funcs );

	//app = new SpaceCraftEditGUI( junk , 800, 600 );
	app->loop( 1000000 );
	return 0;
}
















