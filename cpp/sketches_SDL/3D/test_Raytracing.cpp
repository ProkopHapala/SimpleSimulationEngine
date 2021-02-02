
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

    double ray( const Vec3d &ray0, const Vec3d &hRay, Vec3d& normal ){
        /*
        //int    imin  = 0;
        Vec3d hX,hY;
        hRay.getSomeOrtho( hX, hY );

        double t_min = 1e+300;
        //Vec3d hitpos_min,normal_min;
        for(int i=0; i<ntris; i++ ){
            Vec3i itri = tris[i];
            Vec3d A = points[itri.x];
            Vec3d B = points[itri.y];
            Vec3d C = points[itri.z];
            Vec3d normal_;
            bool inside_;
            //double t = rayTriangle( ray0, hRay, A, B, C, inside_, hitpos_, normal_ );
            double t = rayTriangle2( ray0, hRay, hX, hY, A, B, C, normal_ );
            //printf( "t=%f\n", t );
            inside_ = (t<0.9e+300 )&&(t>0);
            if( inside_ && ( t<t_min ) ){
                t_min = t;
                normal = normal_;
            }
        }
        return t_min;
        */
        int imin;
        return rayTriangles( ray0, hRay, ntris, tris, points, normal, imin );
    }

};


// ============= Application

class TestAppRaytracing : public AppSDL2OGL_3D {
	public:
    //MultiFight3DWorld world;
    double dvel = 10.0;

    //std::vector<KinematicBody*> objects;
    TriangleMesh mesh;

    const static int ngons   = 5;
    int              vertinds[ngons]   = {0,1,2,3,4};
    double           verts   [ngons*3] = { -1.0,-1.0,5.0,    2.0,-1.0,5.0,    2.0,1.0,5.0,   1.0,2.0,5.0,   -1.0,1.0,5.0     };

    int nRays = 3;
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

    double dray = 3.8;
    raySource    .set(3.0, 3.0, 3.0);
    Vec3d rayDir0; rayDir0.set(raySource*-1);
    hRays = new Vec3d[nRays];
    for( int i=0; i<nRays; i++ ){
        hRays[i].set( rayDir0 );
        hRays[i].add( randf(-dray,dray), randf(-dray,dray), randf(-dray,dray) );
        hRays[i].normalize();
    }

/*
    mesh.init( 3, 1 );
    mesh.points[0].set( 1.0,0.0,0.0 );
    mesh.points[1].set( 0.0,1.0,0.0 );
    mesh.points[2].set( 0.0,0.0,1.0 );
    mesh.tris  [0].set( 0, 1, 2 );

    nRays = 1;
    hRays = new Vec3d[nRays];
    raySource.set(-3.0, -3.0, -3.0);
    hRays[0].set( 1.0, 1.0, 1.0 );
    hRays[0].normalize();
*/
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
    glEnable(GL_DEPTH_TEST);
    double t;
    Vec3d hitpos, normal;

    /*
    glCallList( defaultObjectShape );
    glDisable ( GL_LIGHTING );
    for( int i=0; i<nRays; i++ ) {
        glColor3f( 0.8f, 0.0f, 0.0f ); Draw3D::drawLine( raySource, raySource + hRays[i]*100  );

        bool  inside=false;
        Vec3d hitpos, normal;
        double t = mesh.ray( raySource, hRays[i], inside, hitpos, normal );

        if (inside){
            glColor3f( 0.0f, 0.0f, 0.8f );
            Draw3D::drawPointCross( hitpos, 0.2 );
            Draw3D::drawVecInPos  ( normal, hitpos );
        }
    }
    */

    Vec3d ray0; ray0.set_add( (Vec3d)cam.pos,{5.3,-5.8,-5.9});

    glCallList( defaultObjectShape );
    t = mesh.ray( ray0, (Vec3d)cam.rot.c, normal );   //raySphere( camPos, camMat.c, 2.0, * );
    hitpos.set_add_mul( ray0, (Vec3d)cam.rot.c, t);
    if( t<t_inf ){
        glColor3f( 0.0f, 0.0f, 0.8f );
        Draw3D::drawPointCross( hitpos, 0.2 );
        Draw3D::drawVecInPos  ( normal*10, hitpos );
    };

    t = rayPolygon( ray0, (Vec3d)cam.rot.c, (Vec3d)cam.rot.a, (Vec3d)cam.rot.b, ngons, vertinds, (Vec3d*)verts, normal );   //raySphere( camPos, camMat.c, 2.0, * );
    hitpos.set_add_mul( ray0, (Vec3d)cam.rot.c, t);
    if( t<t_inf ){
        glColor3f( 0.0f, 0.0f, 0.8f );
        Draw3D::drawPointCross( hitpos, 0.2 );
        Draw3D::drawVecInPos  ( normal, hitpos );
        glColor3f(1.0f,0.0f,0.0f);
    }else{ glColor3f(1.0f,1.0f,1.0f); };
    Draw3D::drawPlanarPolygon( ngons, vertinds, (Vec3d*)verts );

    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
    glColor3f( 1.0f, 1.0f, 1.0f );
    Draw3D::drawPointCross( ray0, 1.0 );
    Draw3D::drawVecInPos  ( (Vec3d)cam.rot.c, ray0 );

    //STOP = true;

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
















