
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

#include "Truss.h"
#include "SpaceCraft.h"
#include "SpaceCraft2Mesh2.h"
#include "SpaceCraftDraw.h"
#include "SoftBody.h"

//#include "SphereSampling.h"
//#include "DrawSphereMap.h"
//#include "Draw3D_Surf.h"

#include "AppSDL2OGL_3D.h"
#include "GUI.h"
#include "IO_utils.h"
#include "testUtils.h"

//#include "EditSpaceCraft.h"

//#include "TriangleRayTracer.h"
//#include "Radiosity.h"

#include "Tree.h"

#include "spaceCraftEditorUtils.h"
#include "OrbSim_d.h"

using namespace SpaceCrafting;

enum class EDIT_MODE:int{ vertex=0, edge=1, component=2, size }; // http://www.cprogramming.com/c++11/c++11-nullptr-strongly-typed-enum-class.html

//SpaceCraft craft;
Truss      truss;
SoftBody   body;
std::vector<BondType> bondTypes;

int glo_truss=0;

//char str[8096];
double elementSize  = 5.;

Mesh::Builder2 mesh;
OrbSim         sim;

/**
 * Creates a ship truss structure assuming that the ship is composed of a central spine and peripheral maneuvering pendulums all connected by ropes.
 * 
 * @param truss Truss object to store the truss structure.
 * @param n number of central spine points. 
 * @param m number of peripheral maneuvering pendulum points. (legs)
 * @param L length of each truss segment.
 * @param perRope The number of ropes per connection.
 */
void makeShipTruss1( Truss& truss, int n, int m, double L, int perRope ){
    printf( "makeShipTruss1() n=%i m=%i perRope=%i L=%g\n", n, m, perRope, L );
    // central spine
    for(int i=0; i<(n+1); i++){
        truss.points.push_back( (Vec3d){0.,0.,i*L} ); // 1
        if(i>0) truss.edges.push_back( (TrussEdge){ i-1,i,   0 }   );
    }
    // peripheral maneuvering pendulum points (legs)
    int ip0=truss.points.size();
    for(int i=0; i<n; i++){
        for(int j=0; j<m; j++){ 
            Vec2d rot; rot.fromAngle( (j  + 0.5*i )*2*M_PI/m );
            truss.points.push_back( (Vec3d){rot.x*L,rot.y*L,L*(i+.5)} );
        }
    }
    // ropes
    for(int i=0; i<n; i++){
        int ip=ip0+i*m;
        for(int j=0; j<m; j++){
            int ie = ip+j;
            // conection legs to spine
                     truss.addRope( i  , ie, 1, perRope );
            if(i<n){ truss.addRope( i+1, ie, 1, perRope ); }
            // connection between legs
            int ie2;
            if(j==0){ ie2=ie+m-1; }else{ ie2=ie-1; }
            truss.addRope( ie, ie2, 2, perRope );
        }
    }
    // ToDo : radiators between lengs ?
}

/**
 * Converts a Truss object (which define abstract structure of truss composed of points and edges, and faces) to a SoftBody object (which is a physical model of the truss used in simulations - evaulates forces and moves points). 
 * 
 * @param truss     Truss object to convert.
 * @param bondTypes array of BondType objects.
 * @param body      The SoftBody object to store the converted data and do the simulation later.
 */
void truss2SoftBody( const Truss& truss, BondType* bondTypes, SoftBody& body ){
    printf( "truss2SoftBody() npoints %i nedge %i \n", truss.points.size(), truss.edges.size() );
    body.allocate( truss.points.size(), truss.edges.size() );
    delete [] body.points;
    body.points = (Vec3d*)&(truss.points[0]);      // this will help to visualizatio
    for(int i=0; i<truss.points.size(); i++){
        body.points[i] = truss.points[i];
    }
    for(int i=0; i<truss.edges.size(); i++){
        //printf( "truss2SoftBody : bond[%i]] \n", i );
        const TrussEdge& ed = truss.edges[i];
        Bond&             b = body.bonds [i];
        b.id     = i;
        b.i      = ed.a;
        b.j      = ed.b;
        b.type   = bondTypes + ed.type;
        b.broken = false;
    }
    body.prepareBonds ( true );
    body.preparePoints( true, -1, -1 );
    body.gravity=Vec3dZero;
}

void makeShip1( int n, int m, int perRope, double L, Truss& truss, SoftBody& body, BondType* bondTypes ){
    printf( "makeShip1()\n" );
    makeShipTruss1( truss, n, m, L, perRope );
    truss2SoftBody( truss, &bondTypes[0], body );
}

void makeShip_Wheel( int nseg=8){
    printf("makeShip_Wheel()\n");
    StickMaterial *o = new StickMaterial();
    //Material{ name="Kevlar", density=1.44e+3, Spull=3.6e+9, Spush=0.0,    Kpull=154.0e+9, Kpush=0.0,      reflectivity=0.6,  Tmelt=350 }
    //Material{ name="Steel" , density=7.89e+3, Spull=1.2e+9, Spush=1.2e+9, Kpull=200.0e+9, Kpush=200.0e+9, reflectivity=0.85, Tmelt=800 }
    //st1  = StickMaterial( "GS1_long", "Steel", 0.1,  0.005 )
    //st2  = StickMaterial( "GS1_perp", "Steel", 0.05, 0.003 )
    //st3  = StickMaterial( "GS1_in",   "Steel", 0.04, 0.002 )
    //st4  = StickMaterial( "GS1_out",  "Steel", 0.04, 0.002 )

    workshop.add_Material     ( "Steel", 7.89e+3, 1.2e+9, 1.2e+9, 200.0e+9, 200.0e+9, 0.85, 800 );
    workshop.add_StickMaterial( "GS1_long", "Steel", 0.1, 0.005, 0.0 );

    Vec3d p0{0.0,0.0,0.0};
    Vec3d p1{1.0,0.0,0.0};
    Vec3d ax{0.0,0.0,1.0};
    
    //BuildCraft_truss( mesh, *theSpaceCraft, 30.0 );
    mesh.clear();
    //mesh.block();
    //mesh.wheel( p0, p1, ax, nseg, 0.2 );
    //wheel( mesh, p0, p1, ax, nseg, Vec2d{0.2,0.2}, Quat4i{0,0,0,0} );
    wheel( mesh, p0, p1, ax, 3, Vec2d{0.2,0.2}, Quat4i{0,0,0,0} );
    //wheel( mesh, o->pose.pos, o->pose.pos+o->pose.rot.b*o->R, o->pose.rot.c, o->nseg, o->wh, o->st );
    //Quat4i& b = mesh.blocks.back();
    //o->pointRange = {b.x,(int)mesh.verts.size()};
    //o->stickRange = {b.y,(int)mesh.edges.size()};

    mesh.printSizes();

    //int nneighmax_min = 16;
    int nneighmax_min = 8;
    exportSim( sim, mesh, workshop,  nneighmax_min );
    for(int i=0; i<sim.nPoint; i++) sim.points[i].w=1.0;
    for(int i=0; i<sim.nBonds; i++) sim.bparams[i].y=10000.0;
    //sim2.printAllNeighs();

    double omega = 1.0;
    double dt    = 0.05;

    mat2file<double>( "bond_params.log",  sim.nBonds,4, (double*)sim.bparams  );
    

    sim.dt = dt;
    int n = sim.nPoint;
    mat2file<int>( "neighs_before.log",  n, sim.nNeighMax,      sim.neighs,     "%5i " );
    sim.prepare_Cholesky( 0.05, 32 );
    mat2file<int>( "neighs_after.log",   n, sim.nNeighMaxLDLT,  sim.neighsLDLT, "%5i " );

    mat2file<double>( "PDmat.log",  n,n, sim.PDmat  );
    mat2file<double>( "LDLT_L.log", n,n, sim.LDLT_L );
    mat2file<double>( "LDLT_D.log", n,1, sim.LDLT_D );

    sim.cleanVel();
    sim.addAngularVelocity(  p0, ax*omega );
    //apply_torq( sim2.nPoint, p0, ax*omega, sim2.points, sim2.vel );  
    mat2file<double>( "vel_1.log", n,4, (double*)sim.vel );
    
    //exportSim( sim, mesh, workshop );
    //sim.printAllNeighs();

    printf("#### END makeShip_Whee()\n" );
    //exit(0);
};

/**
 * Draws a SoftBody object (which is a physical model of the truss used in simulations - evaulates forces and moves points).
 * 
 * @param body SoftBody object to draw.
 * @param vsc  Scale factor for drawing velocities.
 * @param fsc  Scale factor for drawing forces.
*/
void drawSoftBody( SoftBody& body, float vsc, float fsc ){
    glBegin(GL_LINES);
    glColor3f(.0,.0,.0);
    for(int i=0; i<body.nbonds; i++){
        Bond& b = body.bonds[i];
        //Draw::color_of_hash( b.type->id*1545 + 456 );
        float c =  b.type->kTens*0.5 - b.type->kPress;
        glColor3f(c,c,c);
        //if(c<0){ glLineWidth(2); }else{ glLineWidth(1); }
        Draw3D::vertex( body.points[b.i] );
        Draw3D::vertex( body.points[b.j] );
    }
    for(int i=0; i<body.npoints; i++){
        if(fsc>0){ glColor3f(.9,.0,.0); Draw3D::vertex( body.points[i] ); Draw3D::vertex( body.points[i]+body.forces[i]    *fsc ); }
        if(vsc>0){ glColor3f(.0,.0,.9); Draw3D::vertex( body.points[i] ); Draw3D::vertex( body.points[i]+body.velocities[i]*vsc );}
    }
    if(body.kinks){
        glColor3f(1.0,.0,1.0);
        float f=0.1;
        for(int i=0; i<body.nkink; i++){
            const Kink& k=body.kinks[i];
            Draw3D::vertex( body.points[k.c]*(1-f) + body.points[k.a]*f );
            Draw3D::vertex( body.points[k.c]*(1-f) + body.points[k.b]*f );
        }
    }
    glEnd();
}

// =======================================================================
class SpaceCraftDynamicsApp : public AppSDL2OGL_3D { public:

    bool bDrawTrj=false;
    double time=0;
    int perFrame = 10;
    // https://stackoverflow.com/questions/29145476/requiring-virtual-function-overrides-to-use-override-keyword

	virtual void draw   () override;
	virtual void drawHUD() override;
	//virtual void mouseHandling( );
	//virtual void camera();
	virtual void eventHandling   ( const SDL_Event& event  ) override;
	virtual void keyStateHandling( const Uint8 *keys ) override;
    //virtual void mouseHandling( );

    void makeSoftBody();
    void drawBody();
    void drawSim();

	SpaceCraftDynamicsApp( int& id, int WIDTH_, int HEIGHT_, int argc, char *argv[] );

};

void SpaceCraftDynamicsApp::makeSoftBody(){
    printf("SpaceCraftDynamicsApp::makeSoftBody()\n");
    double modul = 100.0e+9;
    bondTypes.push_back( BondType::stick( 0, 0.1 ,  modul, 2.0e+9 , 2e+3 ) );
    bondTypes.push_back( BondType::rope ( 1, 0.01,  modul, 20.0e+9, 4e+3 ) );
    bondTypes.push_back( BondType::rope ( 2, 0.01,  modul, 20.0e+9, 4e+3 ) );
    //makeShipTruss1( truss, 3, 3, 100.0 );
    int n=1;
    int m=3;
    //makeShipTruss1( truss, n, m, 100.0, 1 );
    makeShipTruss1( truss, n, m, 100.0, 2 ); // multiple segments per rope - seems unstable
    truss2SoftBody( truss, &bondTypes[0], body );
    //body.findKinks( 0.5, 100000000.0 );
    body.findKinks( 100000.0, 100000000.0 );
    /*
    for(int i=0; i<truss.edges.size(); i++){
        TrussEdge& ed = truss.edges[i];
        Vec3d& p1 = truss.points[ed.a];
        Vec3d& p2 = truss.points[ed.b];
        printf( "edge[%i] (%i,%i) (%g,%g,%g) (%g,%g,%g) \n", i, ed.a, ed.b, p1.x,p1.y,p1.z, p2.x,p2.y,p2.z );
    }
    */
    double pendulumWeight = 300.0; // [kg]
    float speed = 10.0;
    //printf( "body.npoints %i \n" , body.npoints );
    for(int i=n+1; i<n+1+n*m; i++){
        //Vec3d& p = body.points[i];
        //double r2 = sq(p.x) + sq(p.y);
        //body.velocities[i].set_cross( p, Vec3dZ*speed );
        body.mass[i] += pendulumWeight;
    }
    body.preparePoints( true, -1.0, -1.0 );
    for(int i=0; i<body.npoints; i++){
        //printf( "body.npoints[%i]  \n" , i );
        Vec3d& p  = body.points[i];
        double r2 = sq(p.x) + sq(p.y);
        //body.points[i].addRandomCube( 5.0 ); // DEBUG RANDOMNESS
        body.velocities[i].set_cross( p, Vec3dZ*speed );
        printf( "point[%i] mass %g[kg] \n" , i, body.mass[i] );
    }
    //body.prepareBonds ( true );
    //body.preparePoints( true, -1, -1 );
    body.updateInvariants();

    perFrame = 100;
    body.dt        = 0.00001;
    body.damp        = 0.0;
    // body.damp_stick  = 0.5;
    body.viscosity = 0.0;
    printf("### SpaceCraftDynamicsApp::makeSoftBody() DONE\n");
}

void SpaceCraftDynamicsApp::drawSim(){
    sim.run_Cholesky(1);
    renderTruss( sim.nBonds, sim.bonds, sim.points, sim.strain, 1000.0 );
};

SpaceCraftDynamicsApp::SpaceCraftDynamicsApp( int& id, int WIDTH_, int HEIGHT_, int argc, char *argv[] ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {
    //Lua1.init();
    fontTex     = makeTexture    ( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    GUI_fontTex = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );
    if(argc<=1){
        //makeSoftBody();
        //reloadShip( "data/ship_ICF_interceptor_1.lua" );
        makeShip_Wheel();
    }
    //exit(0);
    VIEW_DEPTH = 10000.0;
    zoom = 10.0;
}

void SpaceCraftDynamicsApp::drawBody(){
//int npull = 2; int ipulls[]{2,3};
	//int npull = 4; int ipulls[]{2,3,4,5};
	int npull = 4; int ipulls[]{1,2,3,4};
	int iref  = 4;

	//body.fmax = 100000.0;

	double ang0 = atan2(body.points[iref].x,body.points[iref].y);
    for(int itr=0; itr<perFrame; itr++){
        time+=body.dt;
        int pullDir = (((int)(time))%2)*2-1;
        for(int i=0; i<npull; i++){
            body.bonds[ipulls[i]].l0*=(1 + 0.4*body.dt*pullDir );
        }
        //body.step( );
        body.cleanForces();
        body.evalBondForces();
        body.evalKinkForces();   // currently crashes
        body.move_LeapFrog();
    }
    double ang1  = atan2(body.points[iref].x,body.points[iref].y);
    if(ang1<ang0)ang1+=M_PI*2;
    double omega = (ang1-ang0)/(body.dt*perFrame);

    // TODO:
    // ERROR: there is some minor assymetry in construction of the spaceCraft
    //        some is in mass distribution, but even for constant mass distribution angular momentum is not completely along Z-axis

    Vec3d cog   = body.evalCOG();
    Vec3d vcog  = body.evalCOGspeed();
    Vec3d L     = body.evalAngularMomentum();
    double Ekin = body.evalEkin();

    //printf( "cog (%g,%g,%g) \n", cog.x,cog.y,cog.z );
    //printf( "p   (%g,%g,%g) \n", ptot.x,ptot.y,ptot.z );
    //printf( "L   (%g,%g,%g) \n", L.x,L.y,L.z );
    //if(frameCount<10)
    if(frameCount%100==0)printf( "E %g cog (%g,%g,%g) p (%g,%g,%g) L (%g,%g,%g) \n", Ekin, cog.x,cog.y,cog.z, vcog.x,vcog.y,vcog.z, L.x,L.y,L.z );

    double junk;
    //glRotatef( modf(time*1.5,&junk)*360.0,  0.,0.,1.0 );

    //for(int i=0; i<body.npoints; i++){ Vec3d f = body.forces[i]; printf("force[i] (%g,%g,%g) \n", i, f.x,f.y,f.z ); };
	//drawTrussDirect( truss );
	if(frameCount%200==0)bDrawTrj=!bDrawTrj;
	if(bDrawTrj){
        //glClearColor( 0.8f, 0.8f, 0.8f, 0.1f );
        //glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
        glBegin(GL_POINTS);
        glColor3f(0,0,0);
        for(int i=0; i<body.npoints; i++){ Draw3D::vertex(body.points[i]);  };
        glEnd();
        Draw3D::drawPointCross( body.points[0], 3.0 );
	}else{
        drawSoftBody( body, 0.02, 0.00001  );
        //drawSoftBody( body, 0.02, 1.0  );
        glLineWidth(2);
        for(int i=0; i<npull; i++){
            Bond& b = body.bonds[ipulls[i]];
            glColor3f(0.,1.,0.); Draw3D::drawLine( body.points[b.i], body.points[b.j] );
        }
        glLineWidth(1);
        Draw3D::drawAxis(10.0);
    }
}

void SpaceCraftDynamicsApp::draw(){
    //printf( " ==== frame %i \n", frameCount );
    //glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    //glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
    //glClearColor( 0.0f, 0.0f, 0.0f, 1.0f );
    glClearColor( 0.8f, 0.8f, 0.8f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    //drawBody();
    drawSim();

	//if(!bDrawTrj)glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	//glDisable(GL_DEPTH_TEST);
	//glEnable(GL_DEPTH_TEST);
};

void SpaceCraftDynamicsApp::drawHUD(){
    glDisable( GL_LIGHTING );
    glDisable(GL_DEPTH_TEST);
    //gui.draw();
    //glColor3f(1.0f,1.0f,1.0f);   txtStatic.view3D( {5,5}, fontTex, 8 );
    //glPopMatrix();
}


void SpaceCraftDynamicsApp::keyStateHandling( const Uint8 *keys ){

    //Mat3d camMat;
    //qCamera.toMatrix_T(camMat);
    //Draw3D::drawMatInPos(  (Mat3f)camMat, (Vec3f){0.0,0.0,0.0} );
    //qCamera.toMatrix(camMat);

    if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.dyaw  (  keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.dyaw  ( -keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.dpitch(  keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.dpitch( -keyRotSpeed ); }

    if( keys[ SDL_SCANCODE_W ] ){ cam.pos.add_mul( cam.rot.b, +0.05*zoom ); }
	if( keys[ SDL_SCANCODE_S ] ){ cam.pos.add_mul( cam.rot.b, -0.05*zoom );  }
	if( keys[ SDL_SCANCODE_A ] ){ cam.pos.add_mul( cam.rot.a, -0.05*zoom );  }
	if( keys[ SDL_SCANCODE_D ] ){ cam.pos.add_mul( cam.rot.a, +0.05*zoom );  }

};

void SpaceCraftDynamicsApp::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    //if(event.type == SDL_MOUSEWHEEL){
    //    if     (event.wheel.y > 0){ zoom*=VIEW_ZOOM_STEP; }
    //    else if(event.wheel.y < 0){ zoom/=VIEW_ZOOM_STEP; }
    //}

    if(event.type == SDL_MOUSEWHEEL){
        if     (event.wheel.y > 0){ zoom/=VIEW_ZOOM_STEP; }
        else if(event.wheel.y < 0){ zoom*=VIEW_ZOOM_STEP; }
    }
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_t: bDrawTrj=!bDrawTrj; glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT ); break;
                case SDLK_m:  break;
                //case SDLK_h:  warrior1->tryJump(); break;
                case SDLK_l:break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:break;
                case SDL_BUTTON_RIGHT: break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:  break;
                case SDL_BUTTON_RIGHT: break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );

}

// ===================== MAIN

SpaceCraftDynamicsApp * thisApp;

int main(int argc, char *argv[]){
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
	thisApp = new SpaceCraftDynamicsApp( junk , dm.w-150, dm.h-100, argc, argv );
	//thisApp = new SpaceCraftDynamicsApp( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















