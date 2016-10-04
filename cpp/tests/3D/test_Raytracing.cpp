
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw3D.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"
//#include "Body.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"



class TriangleMesh{
    public:
    int ntris;
    int npoints;
    Vec3d * points;
    Vec3i * tris;

    void init(int npoints_, int ntris_){ npoints=npoints_; ntris=ntris_; points = new Vec3d[npoints]; tris = new Vec3i[ntris]; }

    double ray( const Vec3d &ray0, const Vec3d &hRay, bool& inside, Vec3d& hitpos, Vec3d& normal ){
        //int    imin  = 0;
        double t_min = 1e+300;
        //Vec3d hitpos_min,normal_min;
        for(int i=0; i<ntris; i++ ){
            Vec3i itri = tris[i];
            Vec3d A = points[itri.x];
            Vec3d B = points[itri.y];
            Vec3d C = points[itri.z];
            Vec3d hitpos_, normal_;
            bool inside_;
            double t = rayTriangle( ray0, hRay, A, B, C, inside_, hitpos_, normal_ );
            if( inside & ( t<t_min ) ){
                t_min = t;
                inside = true;
                hitpos = hitpos_;
                normal = normal_;
            }
        };
    };

};


// ============= Application

class TestAppRaytracing : public AppSDL2OGL_3D {
	public:
    //MultiFight3DWorld world;
    double dvel = 10.0;

    //std::vector<KinematicBody*> objects;
    TriangleMesh mesh;

    int nRays = 5;
    Vec3d   raySource;
    Vec3d * hRays;

    int defaultObjectShape;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppRaytracing( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppRaytracing::TestAppRaytracing( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    double range   = 5.0;
    //mesh.init( 30, 10 );
    //for( int i=0; i<mesh.npoints; i++ ){ mesh.points[i].set( randf(-range,range), randf(-range,range), randf(-range,range) ); }
    //for( int i=0; i<mesh.ntris;   i++ ){ mesh.tris  [i].set( i%mesh.npoints, (i+1)%mesh.npoints, (i+2)%mesh.npoints ); }
    double scatter = 5.0;
    int ntris      = 10;
    mesh.init( 3*ntris, ntris );
    for( int i=0; i<mesh.ntris; i++ ){
        int i3 = i*3;
        Vec3d A; A.set(randf(-range,range), randf(-range,range), randf(-range,range));
        mesh.points[i3  ].set( A );
        mesh.points[i3+1].set( A ); mesh.points[i3+1].add( randf(-scatter,scatter), randf(-scatter,scatter), randf(-scatter,scatter) );
        mesh.points[i3+2].set( A ); mesh.points[i3+2].add( randf(-scatter,scatter), randf(-scatter,scatter), randf(-scatter,scatter) );
        mesh.tris  [i].set( i3, i3+1, i3+2 );
    }

    double dray = 0.3;
    hRays = new Vec3d[0];
    raySource    .set(3.0, 3.0, 3.0);
    Vec3d rayDir0; rayDir0.set(raySource*-1);
    for( int i=0; i<nRays; i++ ){
        hRays[i].set( rayDir0 );
        hRays[i].add( randf(-dray,dray), randf(-dray,dray), randf(-dray,dray) );
        hRays[i].normalize();
    }

    defaultObjectShape = glGenLists(1);
    glNewList( defaultObjectShape , GL_COMPILE );
        glEnable( GL_LIGHTING );
        glColor3f( 0.8f, 0.8f, 0.8f );
        Draw3D::drawTriangles( mesh.ntris, (int*)mesh.tris, mesh.points );
    glEndList();

    zoom = 16.0;

    printf("initialization DONE !");

}

void TestAppRaytracing::draw(){
    printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    glEnable( GL_LIGHTING );
    glCallList( defaultObjectShape );

    glDisable ( GL_LIGHTING );
    for( int i=0; i<nRays; i++ ) {
        bool  inside;
        Vec3d hitpos, normal;
        glColor3f( 0.8f, 0.0f, 0.0f ); Draw3D::drawLine( raySource, raySource + hRays[i]*100  );
        double t = mesh.ray( raySource, hRays[i], inside, hitpos, normal );
        if (inside){
            glColor3f( 0.0f, 0.0f, 0.8f );
            Draw3D::drawPointCross( hitpos, 0.2 );
            Draw3D::drawVecInPos  ( normal, hitpos );
        }
    }


};


void TestAppRaytracing::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}


void TestAppRaytracing::drawHUD(){
    glDisable ( GL_LIGHTING );
    glColor3f( 0.0f, 1.0f, 0.0f );
    glBegin( GL_LINES );
    float whalf = WIDTH *0.5;
    float hhalf = HEIGHT*0.5;
    glVertex3f( whalf-10,hhalf, 0 ); glVertex3f( whalf+10,hhalf, 0 );
    glVertex3f( whalf,hhalf-10, 0 ); glVertex3f( whalf,hhalf+10, 0 );
    glEnd();
}

// ===================== MAIN

TestAppRaytracing * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppRaytracing( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















