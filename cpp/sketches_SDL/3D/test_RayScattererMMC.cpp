
/// @file @brief  This program demonstrates volumetric light scattering using a Monte Carlo method, implemented in `RayScatter.h`. It simulates the transport of light or particles through a medium by tracing the random paths of many individual rays as they scatter. This is useful for rendering effects like fog, smoke, or murky water.
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw3D.h"
#include "Solids.h"

#include "VecN.h"
#include "NumT.h"
//#include "Opt3d.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"
#include "geom3D.h"
//#include "Body.h"

#include "Mesh.h"
//#include "CMesh.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

#include "RayScatter.h"

#include "VecNm.h"



void drawTrinagles( const RayScatterWorld& scat  ){
    int n=scat.triangles.size();
    glColor3f( 0.,1.,0.); Draw3D::drawTriangles( n, (int*)&scat.triangles[0], &scat.points[0],  1 );
    glColor3f( 0.,0.,1.); Draw3D::drawTriangles( n, (int*)&scat.triangles[0], &scat.points[0],  2  );
};

// ============= Application

class TestAppRayScattererMMC : public AppSDL2OGL_3D { public:

    //int defaultObjectShape;
    RayScatterWorld scene;
    int ogl;
    bool bIntegrate=false;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppRayScattererMMC( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppRayScattererMMC::TestAppRayScattererMMC( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    ScatterMatrial* mat;
    scene.materials.resize(2);
    //mat = &scene.materials[0]; mat->scatterProb=0.0; //mat->setLorenz(double width);
    //mat = &scene.materials[0]; mat->scatterProb=2.0; mat->setLorenz(0.25);
    mat = &scene.materials[0]; mat->mfp=1.0; mat->setLorenz(0.25);

    CMesh msh;
    ogl = Draw::list( 0 );
    msh = Solids::Tetrahedron;
    scene.addMesh( msh, {-1,0} );
    scene.addMesh( msh, {-1,0} );

    VecNm<double> vnm( msh.nvert, 3, (double*)&scene.points[scene.points.size()-msh.nvert] );
    Vec3d v= (Vec3d){-2.0,0.5,0.0};
    vnm.shift( (double*)&v );

    for(int i=0; i<scene.points.size(); i++){ Vec3d p = scene.points[i]; printf("p[%i] (%g,%g,%g) \n", p.x,p.y,p.z ); }

    //Draw3D::drawMesh( msh, 0, 0xFF00FF00 );
    //glColor3f( 0.,0.,1.); Draw3D::drawTriangles( msh.ntri, (int*)msh.tris, msh.verts,  true );
    drawTrinagles( scene );

    glColor3f( 1.,1.,1.);
    //msh = Solids::Tetrahedron;
    double tg=0.1;
    scene.bDraw = true;
    int nray=5;
    for(int i=0; i<nray; i++){
        printf( "============ ray[%i] \n", i );
        Vec3d ray0 = (Vec3d){10.0,0.0,0.0};
        Vec3d hRay = (Vec3d){-1.0,randf(-tg,tg),randf(-tg,tg)}; hRay.normalize();
        Draw3D::drawPointCross( ray0, 0.1  );
        Draw3D::drawVecInPos  ( hRay, ray0 );
        scene.scatter( 10, ray0, hRay );
    }

    glEndList();


}

void TestAppRayScattererMMC::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	//glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    //glEnable( GL_LIGHTING );
    //glEnable(GL_DEPTH_TEST);

    glDisable( GL_LIGHTING   );
    glDisable( GL_DEPTH_TEST );

    int    nray   = 50;
    int    perRay = 10;
    double tg   = 0.1;
    if(bIntegrate){
        glDisable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
        glColor4f(1.0,1.0,1.0,0.05);
        //msh = Solids::Tetrahedron;
        scene.bDraw = true;
        for(int i=0; i<nray; i++){
            Vec3d ray0 = (Vec3d){10.0,0.0,0.0};
            Vec3d hRay = (Vec3d){-1.0,randf(-tg,tg),randf(-tg,tg)}; hRay.normalize();
            Draw3D::drawPointCross( ray0, 0.1  );
            Draw3D::drawVecInPos  ( hRay, ray0 );
            scene.scatter( perRay, ray0, hRay );
        }
    }else{
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
        glCallList(ogl);
    }

};


void TestAppRayScattererMMC::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
                case SDLK_m:  bIntegrate=!bIntegrate; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}


void TestAppRayScattererMMC::drawHUD(){
    glDisable ( GL_LIGHTING );
    /*
    glColor3f( 0.0f, 1.0f, 0.0f );
    glBegin( GL_LINES );
        float whalf = WIDTH *0.5;
        float hhalf = HEIGHT*0.5;
        glVertex3f( whalf-10,hhalf, 0 ); glVertex3f( whalf+10,hhalf, 0 );
        glVertex3f( whalf,hhalf-10, 0 ); glVertex3f( whalf,hhalf+10, 0 );
    glEnd();
    */
}

// ===================== MAIN

TestAppRayScattererMMC * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppRayScattererMMC( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
