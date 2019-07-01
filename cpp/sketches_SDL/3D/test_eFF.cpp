
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw3D.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "VecN.h"


#include "AppSDL2OGL_3D.h"
#include "testUtils.h"
#include "SDL_utils.h"
#include "Plot2D.h"

//#include "MMFF.h"

//#include "eFF.h"
#include "e2FF.h"

#define R2SAFE  1.0e-8f

/*
int pickParticle( int n, Vec3d * ps, const Mat3d& cam, double R ){
    double tmin =  1e+300;
    int imin    = -1;
    for(int i=0; i<n; i++){
        double x = cam.a.dot(ps[i]);
        doyble y = cam.b.dot(ps[i]);
        double ti = raySphere( ray0, hRay, R, ps[i] );
        if(ti<tmin){ imin=i; tmin=ti; }
    }
    return imin;
}
*/


// ======= THE CLASS

class TestAppRARFF: public AppSDL2OGL_3D { public:

    //RigidAtom     atom1;
    //RigidAtomType type1,type2;

    bool bRun = false;

//    EFF ff;
    E2FF ff;
    int ipicked  = -1, ibpicked = -1;

    //Plot2D plot1;

    //double Emin,Emax;
    //int     npoints;
    //Vec3d*  points  =0;
    //double* Energies=0;
    //Vec3d * Forces  =0;

    int      fontTex;

    virtual void draw   ();
    virtual void drawHUD();
    //virtual void mouseHandling( );
    virtual void eventHandling   ( const SDL_Event& event  );
    //virtual void keyStateHandling( const Uint8 *keys );

    TestAppRARFF( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppRARFF::TestAppRARFF( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );

    //ff.loadFromFile_bas( "data/CH4.bas" );
    //ff.loadFromFile_bas( "data/C2H6.bas" );
    ff.loadFromFile_bas( "data/C2H6_e2FF.bas" );
    //ff.loadFromFile_bas( "data/C2.bas" );
    //ff.loadFromFile_bas( "data/H2.bas" );
    //ff.loadFromFile_bas( "data/C2e2.bas" );
    //ff.loadFromFile_bas( "data/H-e.bas" );
    printf( " ff.na, ff.ne %i %i \n", ff.na, ff.ne );

    ff.clearForce();
    ff.clearVel();
    ff.forceEE();
    ff.forceAE();
    ff.forceAA();
    
    //ff.move( 0.01, 0.9 );

}

void TestAppRARFF::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);

    ff.clearForce();
    ff.forceEE();
    ff.forceAE();
    ff.forceAA();
    //if(bRun) ff.move( 0.1, 0.9 );
    if(bRun) ff.run( 1, 0.1, 0.5 );

    Vec3d d = ff.apos[0]-ff.apos[1];

    //printf("C1-C2 %g C1-e %g C2-e %g \n", (ff.apos[0]-ff.apos[1]).norm(), 
    //                                      (ff.apos[0]-ff.epos[0]).norm(), 
    //                                      (ff.apos[1]-ff.epos[0]).norm() );

    double fsc = 1.0;
    glColor3f(0.0,0.0,0.0);
    for(int i=0; i<ff.na; i++){
        Draw3D::drawPointCross( ff.apos  [i]    , ff.aQ  [i]*0.1 );
        Draw3D::drawVecInPos(   ff.aforce[i]*fsc, ff.apos[i] );
        //printf( " %i %f %f %f %f  \n", i, ff.aQ[i], ff.apos[i].x,ff.apos[i].y,ff.apos[i].z );
        //printf( " %i %f %f %f %f  \n", i, ff.aQ[i], ff.aforce[i].x, ff.aforce[i].y, ff.aforce[i].z );
    }

    glColor3f(1.0,1.0,1.0);
    //for(int i=0; i<ff.ne; i++){
    //    Draw3D::drawPointCross( ff.epos  [i],     0.1 );
    //    Draw3D::drawVecInPos(   ff.eforce[i]*fsc, ff.epos[i] );
    //}
    for(int i=0; i<ff.ne; i+=2){
        Draw3D::drawLine(ff.epos[i],ff.epos[i+1] );
    }

    //exit(0);

};


void TestAppRARFF::drawHUD(){


	glTranslatef( 100.0,100.0,0.0 );
	glScalef    ( 10.0,10.00,1.0  );
	//plot1.view();

}

/*
void TestAppRARFF::mouseHandling( ){
    int mx,my; Uint32 buttons = SDL_GetRelativeMouseState( &mx, &my);
    if ( buttons & SDL_BUTTON(SDL_BUTTON_LEFT)) {
        
    }
    AppSDL2OGL_3D::mouseHandling( );
};
*/

void TestAppRARFF::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
                case SDLK_SPACE: bRun = !bRun;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //ipicked = pickParticle( ff.natoms, ff.apos, ray0, (Vec3d)cam.rot.c , 0.5 );
                break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //ibpicked = pickParticle( ff.natoms, ff.apos, ray0, (Vec3d)cam.rot.c , 0.5 );
                    //printf( "dist %i %i = ", ipicked, ibpicked );
                    break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}

// ===================== MAIN

TestAppRARFF* thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppRARFF( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















