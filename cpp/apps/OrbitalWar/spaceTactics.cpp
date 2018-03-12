
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

#include "spline_hermite.h"
#include "SplineManager.h"


#include "SpaceBodies.h"
#include "SpaceWorld.h"

#include "AppSDL2OGL_3D.h"
#include "GUI.h"
#include "testUtils.h"

/*


#### Ideas:
 -

 - Optimal placement of units is given by
    - effective range of weapons and
    - narrow angle of shields

    *       *       *        *   //army 1
    ..
    . .
    .  .
    .   .
    .    .
    .     .
    .      .
    *       *       *        *   //army 1

*/

/*
void drawPolyLineRef( int n, Vec3d * ps, Vec3d * ps_ref ){   // closed=false
    //printf("%i %i\n", n, closed );
    glBegin(GL_LINE_STRIP);
    for(int i=0; i<n; i++){
        glVertex3d( ps[i].x-ps_ref[i].x, ps[i].y-ps_ref[i].y, ps[i].z-ps_ref[i].z );
    };
    glEnd();
};
*/

const int nThrusts = 6;
const double ThrustTimes[nThrusts  ] = {-10.0,35000.0,40000.0,60000.0,65000.0,100000.0};
const double Thrusts    [nThrusts*3] = {
    0.000, 0.0,0.0,
    0.001, 0.0,0.0,
    -0.200, 0.0,0.0,
    -0.100, 0.0,0.0,
    -0.001, 0.0,0.0,
    0.000, 0.0,0.0
};

class BodyInteraction{ public:
    int i,j,kind;
};

class SpaceTactics : public AppSDL2OGL_3D { public:
    SplineManager splines;

    SpaceWorld world;

    bool dragging=false;
    int     nsamp=0;
    double tstart, tend, dt;


    std::vector<BodyInteraction> interactions;

    float scF = 1.0;
    double view_scale = 1/1e+9;
    Vec3d  view_shift = (Vec3d){0.0,0.0,0.0};


    int iTrjMin=0,iTrjMax=0;
    bool bRefTrj = false;
    int trjStep = 10;
    SpaceBody* referenceBody = 0;

    int          fontTex;
    GUITextInput txtStatic;
    //GUITextInput txtPop;
    char curCaption[100];

    //DEGUB
    //Vec3d *  dpos;
    //Vec3d * ddpos;

    void drawTrj    ( int n, Vec3d* ps );
    void drawBodyTrj( SpaceBody& b );
    void drawPlanet ( SpaceBody& b );
    //void drawInteraction( BodyInteraction& bi );
    void drawInteractionTrj( BodyInteraction& bi );

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	virtual void keyStateHandling( const Uint8 *keys );
    //virtual void mouseHandling( );

	SpaceTactics( int& id, int WIDTH_, int HEIGHT_ );

	inline Vec3d transformTrjPos( Vec3d* ps, int i){
        //printf("%f %f %f \n", ps[i].x,ps[i].y,ps[i].z );
	    if( referenceBody ){
            if( bRefTrj ){ return (ps[i]-view_shift-referenceBody->trjPos[i])*view_scale; }
            else{          return (ps[i]-view_shift-referenceBody->pos      )*view_scale; }
        }else{             return (ps[i]-view_shift                         )*view_scale; }
	};

};

SpaceTactics::SpaceTactics( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    //world.addPlanet( "Sun"  , 100000000000.0, 1.0, { 0.0, 0.0, 0.0}, {0.0,0.0,0.0} );
    //world.addPlanet( "Earth",   1000000000.0, 0.1, {10.0, 0.0, 0.0}, {0.0,1.0,0.0} );

    double vJup = 13.1e+3;
    double rJup = 7.784120E+011;

    world.addPlanet( "Sun"      ,   1.99E+030,  696.0000e+6,  { 0.0, 0.0, 0.0},          {0.0,0.0,0.0}    );
    world.addPlanet( "Jupiter"  ,   1.90E+27,     71.4920e+6, { 0.0,         rJup, 0.0}, {vJup,0.0  ,0.0} );
    world.addPlanet( "Io"       ,   8.93190E+23, 1.83000E+06, { 4.21700E+08, rJup, 0.0}, {vJup,17330,0.0} );
    world.addPlanet( "Europa"   ,   4.80000E+23, 1.56080E+06, { 6.71034E+08, rJup, 0.0}, {vJup,13740,0.0} );
    world.addPlanet( "Ganymede" ,   1.48190E+24, 2.63120E+06, { 1.07041E+09, rJup, 0.0}, {vJup,10880,0.0} );
    world.addPlanet( "Calisto"  ,   1.07590E+24, 2.41030E+06, { 1.88271E+09, rJup, 0.0}, {vJup,8204 ,0.0} );

    referenceBody = &world.planets[4];

    world.addShip( "ShipA1",   10.00E+3, 2.5E+03, (Vec3d){ -4e+6, 0, 0.0} +referenceBody->pos, (Vec3d){ 0.0, -5.6e+3,0.0}+referenceBody->vel );
    world.addShip( "ShipA2",   10.00E+3, 2.5E+03, (Vec3d){ -3e+6, 0, 0.0} +referenceBody->pos, (Vec3d){ 0.0, -5.6e+3,0.0}+referenceBody->vel );
    world.addShip( "ShipA3",   10.00E+3, 2.5E+03, (Vec3d){ -5e+6, 0, 0.0} +referenceBody->pos, (Vec3d){ 0.0, -5.6e+3,0.0}+referenceBody->vel );

    //world.addShip( "ShipB1",   10.00E+3, 2.5E+03, (Vec3d){ 3e+8, rJup, 0.0}, (Vec3d){ vJup, 26.6e+3,0.0} );
    //world.addShip( "ShipB2",   10.00E+3, 2.5E+03, (Vec3d){ 4e+8, rJup, 0.0}, (Vec3d){ vJup, 22.6e+3,0.0} );
    //world.addShip( "ShipB3",   10.00E+3, 2.5E+03, (Vec3d){ 5e+8, rJup, 0.0}, (Vec3d){ vJup, 20.6e+3,0.0} );

    world.addShip( "ShipB1",   10.00E+3, 2.5E+03, (Vec3d){ 0, rJup-4.5e+8, 0.0}, (Vec3d){ vJup+ 20.6e+3,0.0,0.0} );
    world.addShip( "ShipB2",   10.00E+3, 2.5E+03, (Vec3d){ 0, rJup-4.6e+8, 0.0}, (Vec3d){ vJup+ 20.6e+3,0.0,0.0} );
    world.addShip( "ShipB3",   10.00E+3, 2.5E+03, (Vec3d){ 0, rJup-4.7e+8, 0.0}, (Vec3d){ vJup+ 20.6e+3,0.0,0.0} );

    SpaceBody& jup = world.planets[1];
    world.intertialTansform( jup.pos*-1.0, jup.vel*-1.0 );
    bRefTrj = false;
    //referenceBody = &jup;

    world.allocateODE ();
    world.allocateTrjs( 500 );
    //world.allocateTrjs( 20 );

    interactions.push_back( {0,3, 0} );


    double dt = 200.0;
    nonUni2spline( 0.0, dt, nThrusts, ThrustTimes, (Vec3d*)Thrusts, world.trj_n, world.ships[3].trjThrust );

    world.predictTrjs( 100000000, dt );
    //exit(0);

    //for( SpaceBody& b : world.ships ){
    //    for(int i=0; i<world.trj_n; i++){ printf("%s[%i] T= %f %f %f \n", b.name.c_str(), i, b.trjThrust[i].x, b.trjThrust[i].y, b.trjThrust[i].z ); }
    //}
    //exit(0);

    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    //txtPop   .inputText = "txtPop";
    txtStatic.inputText = "txtStatic";

    zoom = 5;

}

void SpaceTactics::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	//glEnable(GL_DEPTH_TEST);

    //bool bAnimate = false;
    bool bAnimate = true;
	if(bAnimate){
        int nPhases = 10;
        int nChunk = (world.trj_n/nPhases);
        iTrjMin = ( ((frameCount/10)%nPhases)*nChunk );
        //iTrjMax=_min(iTrjMin+nChunk,world.trj_n);
        iTrjMax=_min(iTrjMin+trjStep,world.trj_n);
	}else{
        iTrjMin=0; iTrjMax=world.trj_n;
	}


	glColor3f(0.0f,0.0f,0.0f); for(SpaceBody& b : world.planets ){ drawBodyTrj( b ); drawPlanet( b ); }
    glColor3f(0.0f,0.0f,1.0f); for(SpaceBody& b : world.ships   ){ drawBodyTrj( b ); }

    glColor3f(0.0f,0.5f,0.0f); SpaceBody& b = world.planets[2];
    //printf("%s.trj \n", b.name.c_str() );
    //for(int i=0; i<world.trj_n; i++){ printf("%i %f %f %f \n ", i, b.trjPos[i].x, b.trjPos[i].y, b.trjPos[i].z ); };
    //exit(0);

    for( BodyInteraction& bi : interactions ){ drawInteractionTrj( bi ); }

};

void SpaceTactics::drawHUD(){
    glDisable ( GL_LIGHTING );
    glColor3f(1.0f,1.0f,1.0f);   txtStatic.view3D( {5,5}, fontTex, 8 );
}

void SpaceTactics::eventHandling( const SDL_Event& event ){
    switch( event.type ){
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:  LMB = true;  break;
                case SDL_BUTTON_RIGHT: RMB = true;  break;
            };  break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:  LMB = false; break;
                case SDL_BUTTON_RIGHT: RMB = false; break;
            }; break;
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_ESCAPE:   quit(); break;
                //case SDLK_SPACE:    STOP = !STOP; printf( STOP ? " STOPED\n" : " UNSTOPED\n"); break;
                case SDLK_KP_MINUS: view_scale/=VIEW_ZOOM_STEP; break;
                case SDLK_KP_PLUS:  view_scale*=VIEW_ZOOM_STEP; break;
            } break;
        case SDL_QUIT: quit(); break;
    };
};

void SpaceTactics::keyStateHandling( const Uint8 *keys ){

    if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.dyaw  (  keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.dyaw  ( -keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.dpitch(  keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.dpitch( -keyRotSpeed ); }

	double step = 0.05/view_scale;
    if( keys[ SDL_SCANCODE_W  ] ){ view_shift.y +=step; }
	if( keys[ SDL_SCANCODE_S  ] ){ view_shift.y -=step; }
	if( keys[ SDL_SCANCODE_A  ] ){ view_shift.x -=step; }
	if( keys[ SDL_SCANCODE_D  ] ){ view_shift.x +=step; }
	//printf(" view_shift %g %g %g \n", view_shift.x, view_shift.y, view_shift.z);
};

//void SpaceTactics::keyStateHandling( const Uint8 *keys ){ };

void SpaceTactics::drawPlanet( SpaceBody& b ){
    Vec3d p = b.pos;
    p.sub(view_shift);
    if( referenceBody ) p.sub(referenceBody->pos);
    //Draw3D::drawPointCross( p*view_scale, b.radius*view_scale );
    p.mul(view_scale);
    Draw3D::drawPointCross( p, 0.1 );
    //printf( "%s R %f (%f,%f,%f) \n", b.name.c_str(), b.radius*view_scale, p.x, p.y, p.z );
    //Draw3D::drawSphereOctLines(2, b.radius*view_scale, p );
    Draw3D::drawSphere_oct(16, b.radius*view_scale, p );
    //Draw3D::drawSphere_oct(3, 0.2, p );

    //Draw3D::drawSphere_oct(3, 0.2, {0.0,0.0,0.0} );
    //Draw3D::drawSphere_oct(3, 0.2, {1.0,2.0,1.0} );
    //Draw3D::drawSphere_oct(3, 0.2, p);
}

void SpaceTactics::drawTrj( int n, Vec3d * ps ){
    glBegin(GL_LINE_STRIP);
    Vec3f p;
    if( referenceBody ){
        if( bRefTrj ){ for(int i=0; i<n; i++){
            convert( (ps[i] - referenceBody->trjPos[i] - view_shift)*view_scale, p );
            //printf( "%i : %f %f %f \n", i, p.x, p.y, p.z );
            glVertex3d( p.x, p.y, p.z );
        }}else{ for(int i=0; i<n; i++){
            convert( (ps[i] - referenceBody->pos - view_shift)*view_scale, p );
            glVertex3d( p.x, p.y, p.z );
        }}
    }else{ for(int i=0; i<n; i++){
            convert( (ps[i] - view_shift)*view_scale, p );
            glVertex3d( p.x, p.y, p.z );
    }}
    glEnd();

};

void SpaceTactics::drawBodyTrj( SpaceBody& b ){
    int n = world.trj_n;
    glBegin(GL_LINE_STRIP);
    Vec3d p,f;
    for(int i=0; i<n; i++){ p=transformTrjPos( b.trjPos, i ); glVertex3d( p.x, p.y, p.z ); };
    glEnd();
    if( b.trjThrust ){
        double sc = scF * view_scale * 1e+9;
        glBegin(GL_LINES );
        for(int i=iTrjMin; i<iTrjMax; i+=trjStep ){ p=transformTrjPos( b.trjPos, i ); f=b.trjThrust[i]; glVertex3d( p.x, p.y, p.z ); glVertex3d( p.x+f.x*sc, p.y+f.y*sc, p.z+f.z*sc );
            //if( &b == &world.ships[0] ) printf( "thrust %s[%i] (%f,%f,%f) \n", b.name.c_str(), i, f.x, f.y, f.z );
        };
        glEnd();
    }
};


void SpaceTactics::drawInteractionTrj( BodyInteraction& bi ){
    int n = world.trj_n;
    glBegin(GL_LINE_STRIP);
    Vec3d p;
    glBegin(GL_LINES );
    for(int i=iTrjMin; i<iTrjMax; i+=trjStep ){
        p=transformTrjPos( world.ships[bi.i].trjPos, i ); glVertex3d( p.x, p.y, p.z );
        p=transformTrjPos( world.ships[bi.j].trjPos, i ); glVertex3d( p.x, p.y, p.z );
        //if( &b == &world.ships[0] ) printf( "thrust %s[%i] (%f,%f,%f) \n", b.name.c_str(), i, f.x, f.y, f.z );
    };
    glEnd();
};


// ===================== MAIN

SpaceTactics * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new SpaceTactics( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















