
/// @file @brief  This is a physics demo showcasing a method for resolving collisions and packing objects together. It uses an iterative Gauss-Seidel solver (`SphereGaussSeidel.h`) to enforce non-penetration constraints between a collection of spheres or boxes. The user can watch as the objects, which may start in an overlapping state, quickly settle into a valid, tightly-packed configuration.
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <unordered_map>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "fastmath.h"
#include "Vec3.h"
//#include "Mat3.h"
#include "quaternion.h"
//#include "raytrace.h"
#include "Solids.h"

#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"

#include "SphereGaussSeidel.h"

#include "testUtils.h"

// ======================  TestApp


void tetraStar(){
    uint32_t cpoly = 0xFFFFFFFF;
    uint32_t cwire = 0xFF000000;
    Draw3D::drawMesh( Solids::Tetrahedron, cpoly, cwire );
    glPushMatrix(); glRotated(60, 1, 1, 1);  Draw3D::drawMesh( Solids::Tetrahedron, cpoly, cwire ); glPopMatrix();
    glPushMatrix(); glRotated(60,-1, 1, 1);  Draw3D::drawMesh( Solids::Tetrahedron, cpoly, cwire ); glPopMatrix();
    glPushMatrix(); glRotated(60, 1,-1, 1);  Draw3D::drawMesh( Solids::Tetrahedron, cpoly, cwire ); glPopMatrix();
    glPushMatrix(); glRotated(60, 1, 1,-1);  Draw3D::drawMesh( Solids::Tetrahedron, cpoly, cwire ); glPopMatrix();
    glPushMatrix(); glRotated(45, 1,0,0);    Draw3D::drawMesh( Solids::Tetrahedron, cpoly, cwire ); glPopMatrix();
    glPushMatrix(); glRotated(45, 0,1,0);    Draw3D::drawMesh( Solids::Tetrahedron, cpoly, cwire ); glPopMatrix();
    glPushMatrix(); glRotated(45, 0,0,1);    Draw3D::drawMesh( Solids::Tetrahedron, cpoly, cwire ); glPopMatrix();
    glScalef(-1,-1,-1);
    glColor3f(1,1,1);
    Draw3D::drawMesh( Solids::Tetrahedron, cpoly, cwire );
    glPushMatrix(); glRotated(60, 1, 1, 1);  Draw3D::drawMesh( Solids::Tetrahedron, cpoly, cwire ); glPopMatrix();
    glPushMatrix(); glRotated(60,-1, 1, 1);  Draw3D::drawMesh( Solids::Tetrahedron, cpoly, cwire ); glPopMatrix();
    glPushMatrix(); glRotated(60, 1,-1, 1);  Draw3D::drawMesh( Solids::Tetrahedron, cpoly, cwire ); glPopMatrix();
    glPushMatrix(); glRotated(60, 1, 1,-1);  Draw3D::drawMesh( Solids::Tetrahedron, cpoly, cwire ); glPopMatrix();
    glPushMatrix(); glRotated(45, 1,0,0);    Draw3D::drawMesh( Solids::Tetrahedron, cpoly, cwire ); glPopMatrix();
    glPushMatrix(); glRotated(45, 0,1,0);    Draw3D::drawMesh( Solids::Tetrahedron, cpoly, cwire ); glPopMatrix();
    glPushMatrix(); glRotated(45, 0,0,1);    Draw3D::drawMesh( Solids::Tetrahedron, cpoly, cwire ); glPopMatrix();
}

class TestAppSphereGaussSeidel : public AppSDL2OGL_3D { public:

    SphereGaussSeidel gs;

    int ipick=0;

    bool running = false;
    int perFrame = 10;

	virtual void draw   ();
	virtual void eventHandling   ( const SDL_Event& event  );

	TestAppSphereGaussSeidel ( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppSphereGaussSeidel::TestAppSphereGaussSeidel( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {
    //glEnable(GL_AUTO_NORMALS);
    //fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );

    int n = 100;
    gs.realloc_all(n);
    Vec3d p = Vec3dZero;

    for( int i=0; i<n; i++ ){
        gs.Rs [i] = randf(0.5,2.0);
        //gs.Rs [i] = randf(1.0,1.0);
        gs.pos[i] = Vec3dZero;
        gs.pos[i].addRandomCube( 10.0 );
    }

    //delay = 250;
}

void TestAppSphereGaussSeidel::draw   (){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glEnable(GL_DEPTH_TEST);

	glMatrixMode(GL_MODELVIEW);
    //glColor3f(0,0,1); Draw3D::drawLines( msh.nedge, (int*)msh.edges, msh.verts );
    setLightingNormal();
    //setLightingRGB();

    if(running){
        double eConv = 1e-3;
        double err = gs.eval_errors();
        //int nmove = gs.update_pos( err*0.1, 0.5, 0.5 );
        double maxMove = gs.update_pos( fmax( err*0.1, eConv ), 0.3, 5.0 );
        //if(nmove>0)
        printf( "==== frame %i err %g maxMove %g \n", frameCount, err, maxMove );
        //if( err<1e-6 ) running=false;
        if( maxMove<eConv ) running=false;
    }

    glColor3f(0,0,0);
    for(int i=0; i<gs.n; i++){
        glColor3f(0,0,0);
        if(i==ipick)glColor3f(0,1,0);
        Draw3D::drawBBox( (Vec3f)gs.pos[i], gs.Rs[i] );
        //glColor3f(1,1,1);Draw3D::drawPointCross( (Vec3f)gs.pos[i], gs.Rs[i] );
        int ioff=i*gs.N_NEIGH;
        for(int k=0; k<gs.N_NEIGH; k++){
            double dk=gs.neighDist[ioff+k];
            if(dk>0){ glColor3f(0,0,1); }else{ glColor3f(1,0,0); };
            if(k>=3)dk*=-1;
            Vec3d d=Vec3dZero;
            int kk = k%3;
            //printf( "i %i k %i kk %i  d %g \n", i, k, kk, dk );
            d.array[kk]=dk;
            Draw3D::drawArrow( gs.pos[i], gs.pos[i]+d, 0.1 );
        }
    }


    /*
    Vec3d dir = (Vec3d)cam.rot.c;
    int i = assignCubicFace( dir );

    printf( "i %i \n", i );
    Draw3D::drawMesh( Solids::Cube, 0xFFFFFFFF, 0xFF000000 );

    Vec3d arrow = Vec3dZero;
    double dk=-10.0;
    if(i>=3)dk=-dk;
    arrow.array[i%3]=dk;
    glColor3f(1,0,0);
    Draw3D::drawVecInPos( arrow, {0,0,0} );
    Draw3D::drawArrow( {0,0,0}, arrow, 0.2 );
    Draw3D::drawAxis(4);
    //tetraStar();
    */

    Draw3D::drawAxis(10);

};

void TestAppSphereGaussSeidel::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    int rot;
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_x:  qCamera.set( 0.0, 0.0, 0.0 ,1.0 ); break;
                case SDLK_y:  qCamera.set( 0.7071, 0.0, 0.0 ,0.7071 );; break;
                case SDLK_z:  qCamera.set( 0.0, 0.7071, 0.0 ,0.7071 ); break;
                case SDLK_LEFTBRACKET:  ipick++; if(ipick>=gs.n)ipick=0; printf("ipick %i \n", ipick ); break;
                case SDLK_RIGHTBRACKET: ipick--; if(ipick<=0)ipick=gs.n; printf("ipick %i \n", ipick ); break;
                case SDLK_SPACE: running = !running; printf( "running %i \n", running ); break;
            }
            break;
        /*
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    pickParticle( world.picked );
                break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    world.picked = NULL;
                    break;
            }
            break;
        */
    };
    AppSDL2OGL_3D::eventHandling( event );
}



// ===================== MAIN

TestAppSphereGaussSeidel * testApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestAppSphereGaussSeidel( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}
