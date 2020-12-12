
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


//#include "Asteroid.h"


#include "SpaceBodies.h"
#include "SpaceWorld.h"
#include "RublePile.h"

#include "asteroidEngineering.h"

#include "SpaceDraw.h"

#include "AppSDL2OGL_3D.h"
#include "ScreenSDL2OGL_3D.h"
#include "GUI.h"
#include "testUtils.h"

char strBuf[0x10000];


/*

================   GAME NAME  ==================
------------------------------------------------
!!!!!!!!!!!!!!!!  SPACE ROCK  !!!!!!!!!!!!!!!!!!
-------------------------------------------------

*/

void drawEifelTower(){
    float h4=325;
    float w3=10.0,   h3=276.13;
    float w2=20.0,   h2=115.73;
    float w1=60.24/2,h1=57.63;
    float w0=124.9/2,h0=0.0;
    glBegin(GL_LINES);
    glVertex3f( w3, w3,h3); glVertex3f(0,0,h4);
    glVertex3f(-w3, w3,h3); glVertex3f(0,0,h4);
    glVertex3f( w3,-w3,h3); glVertex3f(0,0,h4);
    glVertex3f(-w3,-w3,h3); glVertex3f(0,0,h4);

    glVertex3f( w2, w2,h2); glVertex3f( w3, w3,h3);
    glVertex3f(-w2, w2,h2); glVertex3f(-w3, w3,h3);
    glVertex3f( w2,-w2,h2); glVertex3f( w3,-w3,h3);
    glVertex3f(-w2,-w2,h2); glVertex3f(-w3,-w3,h3);

    glVertex3f( w1, w1,h1); glVertex3f( w2, w2,h2);
    glVertex3f(-w1, w1,h1); glVertex3f(-w2, w2,h2);
    glVertex3f( w1,-w1,h1); glVertex3f( w2,-w2,h2);
    glVertex3f(-w1,-w1,h1); glVertex3f(-w2,-w2,h2);

    glVertex3f( w0, w0,h0); glVertex3f( w1, w1,h0);
    glVertex3f(-w0, w0,h0); glVertex3f(-w1, w1,h0);
    glVertex3f( w0,-w0,h0); glVertex3f( w1,-w1,h0);
    glVertex3f(-w0,-w0,h0); glVertex3f(-w1,-w1,h0);

    glEnd();
}




class PlanetView : public ScreenSDL2OGL_3D { public:
    RublePile* asteroid = 0;
    int fontTex;

    PlanetView( int& id, int WIDTH_, int HEIGHT_ ): ScreenSDL2OGL_3D( id, WIDTH_, HEIGHT_ ){}

    ~PlanetView(){
        printf( "~PlanetView \n" );
        if(asteroid)delete asteroid;
        //ScreenSDL2OGL_3D::~ScreenSDL2OGL_3D();
        printf( "~PlanetView DONE \n" );
    }

    virtual void draw(){
        glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

        glEnable(GL_DEPTH_TEST);

        if(frameCount==1){
            cam.zmax = asteroid->radius * 4.0;
            zoom = 25000.0f;
        }
        //printf("cam.zoom %g \n", cam.zoom  );

        double r = asteroid->radius;
        Draw3D::drawBBox( {0.0f,0.0f,0.0f}, asteroid->radius );

        //asteroid->relaxBoulders(0.25,1);

        //printf( "PlanetView.draw() \n" );
        SpaceDraw::rublePile( *asteroid );
        //printf( "PlanetView::draw \n"  );
        //Draw3D::drawAxis( 1.0 );

        glColor3f(0.48f,0.48f,0.48f); Draw3D::drawRectGridLines( {20,20}, {-10000.0 ,-10000.0},  {1000.0,0.0,0.0},  {0.0,1000.0,0.0} );
        glColor3f(0.52f,0.52f,0.52f); Draw3D::drawRectGridLines( {20,20}, {-100000.0,-100000.0}, {10000.0,0.0,0.0}, {0.0,10000.0,0.0} );

        glColor3f(0.0f,0.0f,0.0f);
        //glTranslatef( asteroid->radius , asteroid->radius , asteroid->radius );
        glTranslatef( 0,0, asteroid->radius*1.5 );
        drawEifelTower();

    }

    virtual void drawHUD(){
        glPushMatrix();
        glDisable(GL_DEPTH_TEST);
        glColor3f( 0.0,0.0,0.0 );
        glTranslatef( 10.0,HEIGHT-20,0.0  ); asteroid->infoToBuff(strBuf); Draw::drawText( strBuf, fontTex, fontSizeDef, {80,15} );
        glPopMatrix();
    }


};

class SolarSystemMap : public AppSDL2OGL_3D { public:
    typedef AppSDL2OGL_3D Super;
    SplineManager splines;

    SpaceWorld world;
    double epoch;

    bool dragging=false;
    int     nsamp=0;
    double tstart, tend, dt;

    int ipick = 1;

    GLint olg_orbits = 0;

    PlanetView* planetView = 0;

    //float scF = 1.0; float scSz = 0.2;
    //double view_scale = 1/1e+9;
    //Vec3d  view_shift = (Vec3d){0.0,0.0,0.0};


    //int iTrjMin=0,iTrjMax=0;
    //bool bRefTrj = false;
    //int trjStep = 10;
    SpaceBody* referenceBody = 0;

    int          fontTex;
    GUITextInput txtStatic;
    //GUITextInput txtPop;
    //char curCaption[100];

    //double timeCur = 0.0;
    //double tmin    = 0;
    //double tmax    = 100;

    //DEGUB
    //Vec3d *  dpos;
    //Vec3d * ddpos;

    //int iShip=0;

	virtual void draw   ();
	virtual void drawHUD();
	virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );
    //virtual void mouseHandling( );

    virtual void removeChild(ScreenSDL2OGL* child);

	SolarSystemMap( int& id, int WIDTH_, int HEIGHT_ );

	/*
	inline Vec3d transformTrjPos( Vec3d* ps, int i){
        //printf("%f %f %f \n", ps[i].x,ps[i].y,ps[i].z );
	    if( referenceBody ){
            if( bRefTrj ){ return (ps[i]-view_shift-referenceBody->trjPos[i])*view_scale; }
            else{          return (ps[i]-view_shift-referenceBody->pos      )*view_scale; }
        }else{             return (ps[i]-view_shift                         )*view_scale; }
	};
    */

	//void setTimeRange(double tmin_, double tmax_){ tmin=tmin_; tmax=tmax_; }

	//inline double x2t(double x){ return x*(tmax-tmin)/WIDTH; };
	//inline double t2x(double t){ return WIDTH*t/(tmax-tmin); };

};

void SolarSystemMap::removeChild(ScreenSDL2OGL* child){
    if(child==planetView){ planetView=0; };
    AppSDL2OGL_3D::removeChild(child);
}

SolarSystemMap::SolarSystemMap( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    srand(45454);

    //world.addPlanet( "Sun"  , 100000000000.0, 1.0, { 0.0, 0.0, 0.0}, {0.0,0.0,0.0} );
    //world.addPlanet( "Earth",   1000000000.0, 0.1, {10.0, 0.0, 0.0}, {0.0,1.0,0.0} );

    world.addPlanet( "Sun"      ,   1.99E+030,  696.0000e+6,  { 0.0, 0.0, 0.0},          {0.0,0.0,0.0}    );

    world.planets.reserve(3000);

    SpaceBody* sun = &world.planets.back();
    //referenceBody =

    //                                    name           a                e       inc         lan         apa         ma      R           mass
    world.planets.push_back( SpaceBody( "Mercury", sun,	0.39 *const_AU,	0.21,	0.1222580,	0.8435468,	1.3518701,	4.4026077,2439.7e+3, 3.30E+23) );
    world.planets.push_back( SpaceBody( "Venus",   sun, 0.72 *const_AU,	0.01,	0.0592489,	1.3383305,	2.2956836,	3.1761455,6051.8e+3, 4.87E+24) );
    world.planets.push_back( SpaceBody( "Earth",   sun, 1.00 *const_AU,	0.02,	0.0000009,	-0.1965352,	1.7967674,	1.7534337,6378.14e+3,5.97E+24) );
    world.planets.push_back( SpaceBody( "Mars",    sun, 1.52 *const_AU,	0.09,	0.0322992,	0.8653088,	5.8650191,	6.2038308,3396.19e+3,6.42E+23) );
    world.planets.push_back( SpaceBody( "Jupiter", sun, 5.20 *const_AU,	0.05,	0.0227818,	1.7550359,	0.2575033,	0.6004697,71492e+3,  1.90E+27) );
    world.planets.push_back( SpaceBody( "Saturn",  sun, 9.54 *const_AU,	0.05,	0.0433620,	1.9847019,	1.6132417,	0.8716928,60268e+3,  5.68E+26) );
    world.planets.push_back( SpaceBody( "Uranus",  sun, 19.19*const_AU,	0.05,	0.0134366,	1.2955558,	2.9838889,	5.4669329,25559e+3,  8.68E+25) );
    world.planets.push_back( SpaceBody( "Neptune", sun, 30.07*const_AU,	0.01,	0.0308778,	2.2989772,	0.7848981,	5.3211603,24764e+3,  1.02E+26) );
    world.planets.push_back( SpaceBody( "Pluto",   sun, 39.48*const_AU,	0.25,	0.2991800,	1.9251587,	3.9107027,	4.1700944,1151e+3,   1.31E+22) );

    world.load_astorb( "data/astorb_with_diameter.dat", 3000 );

    /*

    OrbitalElements el;
    //Orbit* orb;
    for(int i=0; i<5; i++){
        el = OrbitalElements( 2.0, 0.5, 0.1, 0.0 + i*0.5, 0.0 );
        char str[40];
        sprintf(str, "body_%0i", i );
        world.planets.push_back( SpaceBody( str, referenceBody, new Orbit(el), 1.0) ); //world.planets.back().orbit=orb;
    }
    */

    //world.intertialTansform( sun.pos*-1.0, sun.vel*-1.0 );

    //bRefTrj = false;


    //fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    fontTex = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );

    //txtPop   .inputText = "txtPop";
    txtStatic.inputText = "txtStatic";


    zoom = 5;

    std::vector<int> view_orbit_sel{1,2,3,4,5,6,7,8,9};
    olg_orbits = glGenLists(1);
    glNewList(olg_orbits, GL_COMPILE);
    for(int i : view_orbit_sel){
        if(world.planets[i].orbit)
        SpaceDraw::orbit( 64, *world.planets[i].orbit, 0, 2*M_PI );
    }
    glEndList();



}

void SolarSystemMap::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	//glEnable(GL_DEPTH_TEST);

	//world.trjTime2index(timeCur,iTrjMin); iTrjMax=iTrjMin+1;

	double epochPerFrame = 0.1;

    glColor3f(0.0f,0.0f,0.0f);
    epoch = frameCount*epochPerFrame;
	SpaceDraw::asPoints( world.planets.size(), &world.planets[0], epoch );
	SpaceDraw::asCrosses( 10, &world.planets[0], epoch, 0.1 );

	glColor3f(0.52f,0.52f,0.52f);
	Draw3D::drawRectGridLines( {40,40}, {-5.0,-5.0}, {0.25,0.0,0.0}, {0.0,0.25,0.0} );

	glColor3f(0.5f,0.7f,0.5f);
	if(ipick>=0){
        SpaceDraw::asCrosses( 1, &world.planets[ipick], epoch, 0.1 );
        SpaceDraw::orbit    ( 64, *world.planets[ipick].orbit, 0, 2*M_PI );
	}

	glColor3f(0.52f,0.52f,0.8f);
	glCallList( olg_orbits );

	/*
	SpaceBody& picked_body = world.planets[ipick];
	if( picked_body.orbit ){
        printf( "picked body : %s \n", picked_body.name.c_str() );
        SpaceDraw::orbit( 100, *picked_body.orbit, 0, 2*M_PI );
	}
	*/

	//glColor3f(0.0f,0.0f,0.0f); for(SpaceBody& b : world.planets ){ drawBodyTrj( b ); drawPlanet( b, iTrjMin, 0 ); }
    //glColor3f(0.0f,0.0f,1.0f); for(SpaceBody& b : world.ships   ){ drawBodyTrj( b ); drawShip  ( b, iTrjMin, 0 );  }

    //glColor3f(0.0f,0.5f,0.0f); SpaceBody& b = world.planets[2];
    //printf("%s.trj \n", b.name.c_str() );
    //for(int i=0; i<world.trj_n; i++){ printf("%i %f %f %f \n ", i, b.trjPos[i].x, b.trjPos[i].y, b.trjPos[i].z ); };
    //exit(0);

    //for( BodyInteraction& bi : interactions ){ drawInteractionTrj( bi ); }

    //for(SpaceBody& b : world.planets ){ b.getTrjPos(iTrjMin,0); };
    //glColor3f(0.0f,0.0f,1.0f); for(SpaceBody& b : world.ships   ){ drawBodyTrj( b ); }

    Draw3D::drawAxis(1.0);


};

void SolarSystemMap::drawHUD(){
    glDisable ( GL_LIGHTING );
    //glColor3f(1.0f,1.0f,1.0f);   txtStatic.view3D( {5,5}, fontTex, 8 );

    //int   nt   =(tmax-tmin)/Time::niceUnitBelow( (tmax-tmin)/40 );
    //float tstep=(tmax-tmin)/nt;
    //float tstep=Time::niceUnitAbove( (tmax-tmin)/40.0 );
    //int   nt = (int)((tmax-tmin)/tstep);
    //printf( "tstep %f nt %i   | %f %f \n", tstep, nt, tmax, tmin );



}

void SolarSystemMap::mouseHandling( ){
    //printf( "window[%i] HasFocus \n", id );
    SDL_GetMouseState( &mouseX, &mouseY );   mouseY=HEIGHT-mouseY;
    mouse_begin_x = (2*mouseX-WIDTH )*zoom/HEIGHT;
    mouse_begin_y = (2*mouseY-HEIGHT)*zoom/HEIGHT;
    int mx,my; Uint32 buttons = SDL_GetRelativeMouseState( &mx, &my);
    //printf( " %i %i \n", mx,my );
    if ( buttons & SDL_BUTTON(SDL_BUTTON_RIGHT)) {
        Quat4f q; q.fromTrackball( 0, 0, -mx*mouseRotSpeed, my*mouseRotSpeed );
        qCamera.qmul_T( q );
    }
    qCamera.toMatrix(cam.rot);
    if ( buttons & SDL_BUTTON(SDL_BUTTON_LEFT)) {
        //if( mouseY<20 ){ timeCur = x2t(mouseX); Time::toStr(timeCur,strBuf); printf( "%s\n", strBuf );  };
        Vec3d mpos = (Vec3d)( cam.pos + cam.rot.a*mouse_begin_x +   cam.rot.b*mouse_begin_y );
        mpos.mul( 1./SpaceDraw::zoom );
        printf( " zoom %g mpos (%g,%g,%g)\n", SpaceDraw::zoom, mpos.x,mpos.y,mpos.z  );
        ipick = world.pickPlanet( mpos, (Vec3d)cam.rot.c, epoch );

        if(ipick>=0){
            int junk;
            if(!planetView){
                printf( " new PlanetView \n" );
                planetView = new PlanetView( junk, 800, 600 );
                planetView->fontTex = fontTex;
                planetView->parent    = this;
                planetView->iinparent = child_windows.size();
                child_windows.push_back(planetView);
                printf( " child_windows.push_back(planetView) \n" );
            }else{
                if(planetView->asteroid) delete planetView->asteroid;
            }

            planetView->asteroid       =  new RublePile();
            planetView->asteroid->body = &world.planets[ipick];
            double mass = world.planets[ipick].mass;
            printf( "mass %g \n", mass );
            srand( ipick + 4654646 );
            int nfrag = 1<<(rand()%4);
            planetView->asteroid->makeRublePile(nfrag,8,20,mass,mass*0.01,mass,2.0);
            //for(Boulder& b : planetView->asteroid->boulders ){ printf( "pos(%g,%g,%g,) rot(%g,%g,%g,)(%g,%g,%g,)(%g,%g,%g,)\n", b.pos.x, b.pos.y, b.pos.z,   b.rotMat.a.x,b.rotMat.a.y,b.rotMat.a.z,  b.rotMat.b.x,b.rotMat.b.y,b.rotMat.b.z,  b.rotMat.c.x,b.rotMat.c.y,b.rotMat.c.z   ); }
            //exit(0);
            //planetView->asteroid->relaxBoulders(0.5,10);
            planetView->asteroid->relaxBoulders(0.25,50);

        }
    }
    //Super::mouseHandling();
};

void SolarSystemMap::eventHandling( const SDL_Event& event ){
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
                //case SDLK_KP_MINUS: view_scale/=VIEW_ZOOM_STEP; break;
                //case SDLK_KP_PLUS:  view_scale*=VIEW_ZOOM_STEP; break;
                case SDLK_KP_MINUS: zoom*=VIEW_ZOOM_STEP; break;
                case SDLK_KP_PLUS:  zoom/=VIEW_ZOOM_STEP; break;
            } break;
        case SDL_QUIT: quit(); break;
    };
};

/*
void SolarSystemMap::keyStateHandling( const Uint8 *keys ){

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
*/


// ===================== MAIN

SolarSystemMap * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	SDL_DisplayMode dm;
    SDL_GetDesktopDisplayMode(0, &dm);
	thisApp = new SolarSystemMap( junk , dm.w-150, dm.h-100 );
	//thisApp = new SolarSystemMap( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















