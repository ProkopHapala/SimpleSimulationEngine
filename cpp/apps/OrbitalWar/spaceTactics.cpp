
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

char strBuf[0x10000];

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
    typedef AppSDL2OGL_3D Super;
    SplineManager splines;

    SpaceWorld world;

    bool dragging=false;
    int     nsamp=0;
    double tstart, tend, dt;


    std::vector<BodyInteraction> interactions;

    float scF = 1.0; float scSz = 0.2;
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

    double timeCur = 0.0;
    double tmin    = 0;
    double tmax    = 100;

    //DEGUB
    //Vec3d *  dpos;
    //Vec3d * ddpos;

    int iShip=0;

    void drawTrj    ( int n, Vec3d* ps );
    void drawBodyTrj( SpaceBody& b );
    void drawPlanet ( SpaceBody& b, int iTrj, double du);
    void drawShip   ( SpaceBody& b, int iTrj, double du);
    //void drawInteraction( BodyInteraction& bi );
    void drawInteractionTrj( BodyInteraction& bi );

	virtual void draw   ();
	virtual void drawHUD();
	virtual void mouseHandling( );
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


	void setTimeRange(double tmin_, double tmax_){ tmin=tmin_; tmax=tmax_; }

	inline double x2t(double x){ return x*(tmax-tmin)/WIDTH; };
	inline double t2x(double t){ return WIDTH*t/(tmax-tmin); };

};

SpaceTactics::SpaceTactics( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    //world.addPlanet( "Sun"  , 100000000000.0, 1.0, { 0.0, 0.0, 0.0}, {0.0,0.0,0.0} );
    //world.addPlanet( "Earth",   1000000000.0, 0.1, {10.0, 0.0, 0.0}, {0.0,1.0,0.0} );

    printf( "DEBUG 0 \n" );

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

    referenceBody = &world.planets[1];

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
    iShip = 3;

    double dt = 200.0;
    nonUni2spline( 0.0, dt, nThrusts, ThrustTimes, (Vec3d*)Thrusts, world.trj_n, world.ships[3].trjThrust );

    world.predictTrjs( 100000000, dt );

    setTimeRange(world.trj_t0, world.trj_t0+world.trj_dt*world.trj_n );
    //exit(0);

    //for( SpaceBody& b : world.ships ){
    //    for(int i=0; i<world.trj_n; i++){ printf("%s[%i] T= %f %f %f \n", b.name.c_str(), i, b.trjThrust[i].x, b.trjThrust[i].y, b.trjThrust[i].z ); }
    //}
    //exit(0);

    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    //txtPop   .inputText = "txtPop";
    txtStatic.inputText = "txtStatic";

    /*
    double nsec=0.00015;
    for(int i=0; i<20; i++){
        printf( "nsec %e %e \n", nsec, Time::niceUnitAbove(nsec) );
        nsec*=10;
    }
    exit(0);
    */
    zoom = 5;

}

void SpaceTactics::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	//glEnable(GL_DEPTH_TEST);

    //bool bAnimate = false;

    /*
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
	*/

	world.trjTime2index(timeCur,iTrjMin); iTrjMax=iTrjMin+1;

	glColor3f(0.0f,0.0f,0.0f); for(SpaceBody& b : world.planets ){ drawBodyTrj( b ); drawPlanet( b, iTrjMin, 0 ); }
    glColor3f(0.0f,0.0f,1.0f); for(SpaceBody& b : world.ships   ){ drawBodyTrj( b ); drawShip  ( b, iTrjMin, 0 );  }

    glColor3f(0.0f,0.5f,0.0f); SpaceBody& b = world.planets[2];
    //printf("%s.trj \n", b.name.c_str() );
    //for(int i=0; i<world.trj_n; i++){ printf("%i %f %f %f \n ", i, b.trjPos[i].x, b.trjPos[i].y, b.trjPos[i].z ); };
    //exit(0);

    for( BodyInteraction& bi : interactions ){ drawInteractionTrj( bi ); }

    for(SpaceBody& b : world.planets ){ b.getTrjPos(iTrjMin,0); };
    glColor3f(0.0f,0.0f,1.0f); for(SpaceBody& b : world.ships   ){ drawBodyTrj( b ); }

};

void SpaceTactics::drawHUD(){
    glDisable ( GL_LIGHTING );
    //glColor3f(1.0f,1.0f,1.0f);   txtStatic.view3D( {5,5}, fontTex, 8 );

    //int   nt   =(tmax-tmin)/Time::niceUnitBelow( (tmax-tmin)/40 );
    //float tstep=(tmax-tmin)/nt;
    float tstep=Time::niceUnitAbove( (tmax-tmin)/40.0 );
    int   nt = (int)((tmax-tmin)/tstep);
    //printf( "tstep %f nt %i   | %f %f \n", tstep, nt, tmax, tmin );

    //float xstep=WIDTH/(float)nt;
    //glPushMatrix();
    glColor3f(1.0f,1.0f,1.0f);
    glBegin(GL_LINES);
    float t,x;
    for( int i=0; i<nt; i++ ){
        t = i*tstep+tmin;
        x = t2x(t);
        glVertex3f(x,20,0); glVertex3f(x,0,0);
        for(int j=0; j<9; j++){ t+=0.1*tstep; x=t2x(t); glVertex3f(x,10,0); glVertex3f(x,0,0); }
    };
    glEnd();
    for( int i=0; i<nt; i++ ){
        t = i*tstep+tmin;
        x = t2x(t);
        Time::toStr( t, strBuf);
        Draw2D::drawText( strBuf, {x,10}, {80,15}, fontTex, 8 );
    }
    glColor3f(0.0,1.0,0.0);
    x=t2x(timeCur); Draw2D::drawLine({x,30},{x,0});
    Time::toStr( timeCur, strBuf);
    Draw2D::drawText( strBuf, {x,20}, {80,15}, fontTex, 8 );


    SpaceBody& thisShip = world.ships[iShip];
    if( thisShip.trjThrust ){
        //float scT = scF*100000;
        float scT = scF*100;
        //printf( "thisShip %i \n", thisShip.trjThrust );
        glBegin(GL_LINES);
        Vec3d oT;
        float ox;
        float T0=50.0;
        for(int i=0; i<world.trj_n; i++){
            double t = world.ind2time(i);
            double x = t2x(t);
            Vec3d T  = thisShip.getThrust(i,0);
            //printf( "T = (%f,%f,%f)\n", T.x,T.y,T.z);
            glColor3f(1.0,0.0,0.0); glVertex3f(ox,T0+oT.x*scT,0); glVertex3f(x,T0+T.x*scT,0);
            glColor3f(0.0,1.0,0.0); glVertex3f(ox,T0+oT.y*scT,0); glVertex3f(x,T0+T.y*scT,0);
            glColor3f(0.0,0.0,1.0); glVertex3f(ox,T0+oT.z*scT,0); glVertex3f(x,T0+T.z*scT,0);
            oT = T;
            ox=x;
        }
        glEnd();
    }


    //Draw::drawText( "abcdefghijklmnopqrstuvwxyz \n0123456789 \nABCDEFGHIJKLMNOPQRTSTUVWXYZ \nxvfgfgdfgdfgdfgdfgdfg", fontTex, 8, {10,5} );
    //glTranslatef( 10.0,HEIGHT-20,0.0  ); glColor3f(1.0,0.0,1.0); currentFormation->reportSetup (strBuf); Draw::drawText( strBuf, fontTex, 8, {80,15} );
    //glTranslatef( 300.0,      0.0,0.0 ); glColor3f(0.0,0.5,0.0); currentFormation->reportStatus(strBuf); Draw::drawText( strBuf, fontTex, 8, {80,15} );
    //glTranslatef( 300.0,      0.0,0.0 ); glColor3f(0.0,0.5,0.8); currentFormation->soldiers[0].type->toStrCaptioned(strBuf); Draw::drawText( strBuf, fontTex, 8, {80,15} );
    //glPopMatrix();


}

void SpaceTactics::mouseHandling( ){
    SDL_GetMouseState( &mouseX, &mouseY );   mouseY=HEIGHT-mouseY;
    mouse_begin_x = (2*mouseX-WIDTH )*zoom/HEIGHT;
    mouse_begin_y = (2*mouseY-HEIGHT)*zoom/HEIGHT;
    int mx,my; Uint32 buttons = SDL_GetRelativeMouseState( &mx, &my);
    //printf( " %i %i \n", mx,my );
    if ( buttons & SDL_BUTTON(SDL_BUTTON_RIGHT)) {
        Quat4d q; q.fromTrackball( 0, 0, -mx*mouseRotSpeed, my*mouseRotSpeed );
        qCamera.qmul_T( q );
    }

    if ( buttons & SDL_BUTTON(SDL_BUTTON_LEFT)) {
        if( mouseY<20 ){ timeCur = x2t(mouseX); Time::toStr(timeCur,strBuf); printf( "%s\n", strBuf );  };
    }
    //Super::mouseHandling();
};

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

void SpaceTactics::drawPlanet( SpaceBody& b, int iTrj, double du ){
    Vec3d p;
    if(iTrj>0){
        p = b.getTrjPos(iTrj,du);
        if( referenceBody ) p.sub( referenceBody->getTrjPos(iTrj,du) );
    }else{
        p = b.pos;
        if( referenceBody ) p.sub(referenceBody->pos);
    }
    p.sub(view_shift);
    p.mul(view_scale);
    Draw3D::drawPointCross( p, 0.1 );
    Draw3D::drawSphere_oct(16, b.radius*view_scale, p );
}


void SpaceTactics::drawShip( SpaceBody& b, int iTrj, double du ){
    Vec3d p;
    if(iTrj>0){
        p = b.getTrjPos(iTrj,du);
        if( referenceBody ) p.sub( referenceBody->getTrjPos(iTrj,du) );
    }else{
        p = b.pos;
        if( referenceBody ) p.sub(referenceBody->pos);
    }
    p.sub(view_shift);
    p.mul(view_scale);
    Mat3d mat = b.rotMat; mat.mul(b.sizes); mat.mul(scSz);
    Draw3D::drawMatInPos ( mat, p );
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
	SDL_DisplayMode dm;
    SDL_GetDesktopDisplayMode(0, &dm);
	thisApp = new SpaceTactics( junk , dm.w-150, dm.h-100 );
	//thisApp = new SpaceTactics( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















