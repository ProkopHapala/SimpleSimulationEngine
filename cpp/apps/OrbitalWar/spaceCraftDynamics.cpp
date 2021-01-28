
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

#include "Truss.h"
#include "SpaceCraft.h"
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

using namespace SpaceCrafting;

enum class EDIT_MODE:int{ vertex=0, edge=1, component=2, size }; // http://www.cprogramming.com/c++11/c++11-nullptr-strongly-typed-enum-class.html

//SpaceCraft craft;
Truss      truss;
SoftBody   body;
std::vector<BondType> bondTypes;

int glo_truss=0;

//char str[8096];
double elementSize  = 5.;


void makeShipTruss1_simple( Truss& truss, int n, int m, double L ){
    // central spine
    for(int i=0; i<(n+1); i++){
        truss.points.push_back( (Vec3d){0.,0.,i*L} ); // 1
        if(i>0) truss.edges.push_back( (TrussEdge){ i-1,i,   0 }   );
    }
    // peripheral maneuvering pendulums
    for(int i=0; i<n; i++){
        for(int j=0; j<m; j++){
            Vec2d rot; rot.fromAngle( (j  + 0.5*i )*2*M_PI/m );
            truss.points.push_back( (Vec3d){rot.x*L,rot.y*L,L*(i+.5)} );
            int ie = truss.points.size()-1;
            // conection legs to spine
                     truss.edges.push_back( (TrussEdge){ i  ,ie, 1 }   );
            if(i<n){ truss.edges.push_back( (TrussEdge){ i+1,ie, 1 }   ); }
            // connection between legs
            int ie2;
            if(j==0){ ie2=ie+m-1; }else{ ie2=ie-1; }
            truss.edges.push_back( (TrussEdge){ ie2, ie, 2 }   );
        }
    }
    // ToDo : radiators between lengs ?
}

void makeShipTruss1( Truss& truss, int n, int m, double L, int perRope ){
    // central spine
    for(int i=0; i<(n+1); i++){
        truss.points.push_back( (Vec3d){0.,0.,i*L} ); // 1
        if(i>0) truss.edges.push_back( (TrussEdge){ i-1,i,   0 }   );
    }
    // peripheral maneuvering pendulum points
    int ip0=truss.points.size();
    for(int i=0; i<n; i++){
        for(int j=0; j<m; j++){
            Vec2d rot; rot.fromAngle( (j  + 0.5*i )*2*M_PI/m );
            truss.points.push_back( (Vec3d){rot.x*L,rot.y*L,L*(i+.5)} );
        }
    }
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

void truss2SoftBody( const Truss& truss, BondType* bondTypes, SoftBody& body ){
    printf( "npoints %i nedge %i \n", truss.points.size(), truss.edges.size() );
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
    makeShipTruss1( truss, n, m, L, perRope );
    truss2SoftBody( truss, &bondTypes[0], body );
}

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
        if(fsc>0){ glColor3f(.0,.0,.9); Draw3D::vertex( body.points[i] ); Draw3D::vertex( body.points[i]+body.forces[i]    *fsc ); }
        if(vsc>0){ glColor3f(.9,.0,.0); Draw3D::vertex( body.points[i] ); Draw3D::vertex( body.points[i]+body.velocities[i]*vsc );}
    }
    glEnd();
}

class SpaceCraftDynamicsApp : public AppSDL2OGL_3D { public:

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

	SpaceCraftDynamicsApp( int& id, int WIDTH_, int HEIGHT_ );

};

SpaceCraftDynamicsApp::SpaceCraftDynamicsApp( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    //Lua1.init();
    fontTex     = makeTexture    ( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    GUI_fontTex = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );

    //	            id, linearDensity, kPress,kTens  sPress,sTens;  // strength
    //bondTypes.push_back( BondType::stick( 0, 0.1 ,  200e+9, 2.0e+9  ) );
    //bondTypes.push_back( BondType::rope ( 0, 0.01,  200e+9, 20.0e+9 ) );
    //bondTypes.push_back( BondType::rope ( 0, 0.01,  200e+9, 20.0e+9 ) );

    double modul = 100.0e+9;
    bondTypes.push_back( BondType::stick( 0, 0.1 ,  modul, 2.0e+9 , 2e+3 ) );
    bondTypes.push_back( BondType::rope ( 1, 0.01,  modul, 20.0e+9, 4e+3 ) );
    bondTypes.push_back( BondType::rope ( 2, 0.01,  modul, 20.0e+9, 4e+3 ) );

    //makeShipTruss1( truss, 3, 3, 100.0 );

    int n=2;
    int m=3;
    //makeShipTruss1( truss, n, m, 100.0, 1 );
    makeShipTruss1( truss, n, m, 100.0, 2 ); // multiple segments per rope - seems unstable
    truss2SoftBody( truss, &bondTypes[0], body );

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
        //body.mass[i] += pendulumWeight;
    }
    body.preparePoints( true, 1.0, 1.0 );
    for(int i=0; i<body.npoints; i++){
        //printf( "body.npoints[%i]  \n" , i );
        Vec3d& p  = body.points[i];
        double r2 = sq(p.x) + sq(p.y);
        //if( r2>1.0 ){
        body.velocities[i].set_cross( p, Vec3dZ*speed );
            //body.mass[i] += pendulumWeight;
        //}
        printf( "point[%i] mass %g[kg] \n" , i, body.mass[i] );
    }
    body.updateInvariants();

    perFrame = 100;
    body.dt        = 0.00001;
    body.damp      = 0.0;
    body.viscosity = 0.0;


    VIEW_DEPTH = 10000.0;
    zoom = 1000.0;
}



void SpaceCraftDynamicsApp::draw(){
    //printf( " ==== frame %i \n", frameCount );
    //glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    //glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
    //glClearColor( 0.0f, 0.0f, 0.0f, 1.0f );
    glClearColor( 0.8f, 0.8f, 0.8f, 1.0f );

	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable(GL_DEPTH_TEST);
	//glEnable(GL_DEPTH_TEST);

	//int npull = 2; int ipulls[]{2,3};
	int npull = 4; int ipulls[]{2,3,4,5};
	int iref  = 4;

	double ang0 = atan2(body.points[iref].x,body.points[iref].y);
    for(int itr=0; itr<perFrame; itr++){

        time+=body.dt;

        /*
        int pullDir = (((int)(time))%2)*2-1;
        for(int i=0; i<npull; i++){
            body.bonds[ipulls[i]].l0*=(1 + 0.2*body.dt*pullDir );
        }
        */
        //body.step( );
        body.cleanForces();
        body.evalBondForces();
        //body.move_LeapFrog();
    }
    double ang1  = atan2(body.points[iref].x,body.points[iref].y);
    if(ang1<ang0)ang1+=M_PI*2;
    double omega = (ang1-ang0)/(body.dt*perFrame);

    Vec3d cog  = body.evalCOG();
    Vec3d vcog = body.evalCOGspeed();
    Vec3d L    = body.evalAngularMomentum();
    //printf( "cog (%g,%g,%g) \n", cog.x,cog.y,cog.z );
    //printf( "p   (%g,%g,%g) \n", ptot.x,ptot.y,ptot.z );
    //printf( "L   (%g,%g,%g) \n", L.x,L.y,L.z );
    if(frameCount<10)printf( "cog (%g,%g,%g) p (%g,%g,%g) L (%g,%g,%g) \n", cog.x,cog.y,cog.z, vcog.x,vcog.y,vcog.z, L.x,L.y,L.z );

    double junk;
    //glRotatef( modf(time*1.5,&junk)*360.0,  0.,0.,1.0 );

	//drawTrussDirect( truss );
	drawSoftBody( body, 0.02, 0.00001  );

	glLineWidth(2);
	for(int i=0; i<npull; i++){
        Bond& b = body.bonds[ipulls[i]];
        glColor3f(0.,1.,0.); Draw3D::drawLine( body.points[b.i], body.points[b.j] );
    }
    glLineWidth(1);

    Draw3D::drawAxis(10.0);

	glDisable(GL_LIGHTING);
	//glEnable(GL_LIGHTING);


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
	thisApp = new SpaceCraftDynamicsApp( junk , dm.w-150, dm.h-100 );
	//thisApp = new SpaceCraftDynamicsApp( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















