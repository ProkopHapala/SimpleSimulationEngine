
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "SDL_utils.h"
#include "Draw.h"
#include "Draw3D.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"
#include "Mesh.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

// ============= Application

class TestAppMesh : public AppSDL2OGL_3D { public:

    Mesh mesh;
    Mesh mesh2;

    Plane3D plane1,plane2,plane3,plane4;
    LineInterval3d l12;

    bool dragging;
    Vec2f mouse0;
    int ipicked;

    int      fontTex;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppMesh( int& id, int WIDTH_, int HEIGHT_ );

};

#include <unistd.h>
TestAppMesh::TestAppMesh( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    getcwd( str, 1024);
    printf("WORKING_DIRECTORY = >>%s<<", str);

    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );

    mesh.fromFileOBJ( "common_resources/turret.obj" );
    //mesh.fromFileOBJ( "common_resources/tank_hull.obj" );

    mesh.findEdges(); //exit(0);

    mesh.polygons[2]->printPoints();
    mesh.polygons[8]->printPoints();
    mesh.insertEdgeVertex( 14 );
    mesh.points[mesh.points.size()-1].add(0.5,0.5,0.0);
    mesh.polygons[2]->printPoints();
    mesh.polygons[8]->printPoints();

    int ip = mesh.colapseEdge(12);
    mesh.cleanRemovedPoints();

    mesh.polygonsToTriangles( false );
    mesh.tris2normals(true);
    printf("initialization DONE ! \n");

    mesh.rendered_shape = glGenLists(1);
    glNewList( mesh.rendered_shape , GL_COMPILE );
        glEnable( GL_LIGHTING );
        glColor3f( 0.8f, 0.8f, 0.8f );

        Draw3D::drawMesh( mesh );
        glColor3f(0.0f,0.0f,0.9f);
        for(int i=0; i<mesh.points.size(); i++){
            Draw3D::drawVecInPos( mesh.normals[i], mesh.points[i] );
        }
        /*
        glBegin(GL_TRIANGLES);
        for( Vec3i tri : mesh.triangles ){
            Vec3f p,n;
            convert( mesh.points[tri.a], p ); convert( mesh.normals[tri.a], n ); glNormal3f( n.x, n.y, n.z ); glVertex3f( p.x, p.y, p.z );
            convert( mesh.points[tri.b], p ); convert( mesh.normals[tri.b], n ); glNormal3f( n.x, n.y, n.z ); glVertex3f( p.x, p.y, p.z );
            convert( mesh.points[tri.c], p ); convert( mesh.normals[tri.c], n ); glNormal3f( n.x, n.y, n.z ); glVertex3f( p.x, p.y, p.z );
        };
        */
        for(int i=0; i<mesh.points.size(); i++){
            Draw3D::drawVecInPos( mesh.normals[i], mesh.points[i] );
        }
        glEnd();
    glEndList();

    // TO DO : convex hull
    /*
    int npoints = 20;
    Vec3d * points = new Vec3d[npoints];
    for(int i=0; i<npoints; i++){
        points[i].set( randf(-1.0,1.0), randf(-1.0,1.0), randf(-1.0,1.0) );
    }

    mesh.rendered_shape = glGenLists(1);
    glNewList( mesh.rendered_shape , GL_COMPILE );
    glDisable(GL_LIGHTING);
    glColor3f(1.0,1.0,1.0);
    glBegin( GL_POINTS );
        for(int i=0; i<npoints; i++){ glVertex3f( points[i].x, points[i].y, points[i].z ); }
    glEnd();
    glEndList();
    */

    /*
    plane1 = {1.0,0.0,0.0, -1.0};
    plane2 = {0.0,1.0,0.0, -1.0};
    l12.fromPlanes( plane1.normal, plane1.iso, plane2.normal, plane2.iso );

    plane3 = {0.0,0.0,+1.0, -3.0};
    plane4 = {0.0,0.0,-1.0, -2.0};

    l12.trim(plane3.normal, plane3.C);
    l12.trim(plane4.normal, plane4.C);
    printf( "l12 t (%g,%g) \n", l12.t0, l12.t1 );
    */

    /*
    for(int i=0; i<nplanes; i++){
        planes[i].normal.fromRandomSphereSample();
        planes[i].normal.normalize();
    }
    */

    for(int i=0; i<18; i++){
        float d=0.3;
        planes[i].normal.add(randf(-d,d),randf(-d,d),randf(-d,d));
        planes[i].normal.normalize();
    }

    mesh2.fromPlanes( nplanes, planes );
    printf( "mesh2: npoints %i nedges %i \n", mesh2.points.size(), mesh2.edges.size() );

    for( MeshEdge& edge: mesh2.edges ){
        Vec3d& a =  mesh.points[ edge.verts.a ];
        Vec3d& b =  mesh.points[ edge.verts.b ];
        printf( "edge %i(%g,%g,%g) %i(%g,%g,%g) \n", edge.verts.a, a.x,a.y,a.z,   edge.verts.b, b.x,b.y,b.z );
    }

    //exit(0);

}

void TestAppMesh::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );


    glColor3f(0.0,0.0,0.0);
    for( MeshEdge& edge: mesh2.edges ){
        Vec3d& a =  mesh2.points[ edge.verts.a ];
        Vec3d& b =  mesh2.points[ edge.verts.b ];
        //Draw3D::drawPointCross( l12.p0, 0.1 ); Draw3D::drawVecInPos( l12.hdir, l12.p0 );
        //printf( "edge %i(%g,%g,%g) %i(%g,%g,%g) \n", edge.verts.a, a.x,a.y,a.z,   edge.verts.b, b.x,b.y,b.z );
        Draw3D::drawLine( a, b );
    }

    for( Vec3d& p: mesh2.points ){
        glEnable(GL_DEPTH_TEST);
        Draw3D::drawPointCross( p, 0.1 );
    }

    int i=0;
    for( Polygon* pl: mesh2.polygons ){
        Draw::color_of_hash(i*4456464+54844); i++;
        glBegin(GL_TRIANGLE_FAN);
        for( int i : pl->ipoints ){
            Vec3d& v = mesh2.points[i];
            glVertex3f( v.x,v.y,v.z  );
        }
        glEnd();
    }

};


void TestAppMesh::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
                case SDLK_t:
                    plane1.normal.fromRandomSphereSample();
                    plane2.normal.fromRandomSphereSample();
                    float perturb = 0.3;
                    plane2.normal = plane1.normal*(perturb-1) + plane2.normal*perturb;
                    plane2.normal.normalize();
                    plane1.iso = randf(-0.7,-1.8);
                    plane1.iso = randf(-0.7,-1.8);
                    perturb = 0.5;
                    l12.p0.fromRandomCube(0.5);
                    l12.hdir.fromRandomSphereSample();
                    l12.hdir = plane1.normal*(perturb-1) + l12.hdir*perturb;
                    l12.hdir.normalize();
                    l12.infiniteSpan();
                    l12.trim(plane1.normal, plane1.C);
                    l12.trim(plane2.normal, plane2.C);
                    printf( "%g %g (%g,%g%g) (%g,%g%g)\n", l12.t0, l12.t1,    l12.p0.x,l12.p0.y,l12.p0.z,    l12.hdir.x,l12.hdir.y,l12.hdir.z );
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    Vec3d ray0;
                    ray0.set_lincomb( 1, mouse_begin_x, mouse_begin_y, (Vec3d)cam.pos, (Vec3d)cam.rot.a, (Vec3d)cam.rot.b );
                    //glColor3f(1,1,1); Draw3D::drawPointCross( ray0, 0.2 );
                    ipicked= mesh.pickVertex( ray0, (Vec3d)cam.rot.c );
                    mouse0.set(mouse_begin_x, mouse_begin_y);
                    dragging=true;
                    break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    if(dragging){
                        Vec3d d;
                        d.set_lincomb( 1, mouse_begin_x, mouse_begin_y, (Vec3d)cam.pos, (Vec3d)cam.rot.a, (Vec3d)cam.rot.b );
                        d.sub(mesh.points[ipicked]);
                        double c = d.dot((Vec3d)cam.rot.c);  d.add_mul((Vec3d)cam.rot.c, -c);
                        mesh.points[ipicked].add(d);

                        glDeleteLists(mesh.rendered_shape, 1);
                        mesh.rendered_shape = glGenLists(1);
                        glNewList( mesh.rendered_shape , GL_COMPILE );
                            glEnable( GL_LIGHTING );
                            glColor3f( 0.8f, 0.8f, 0.8f );
                            Draw3D::drawMesh( mesh );
                            mesh.tris2normals(true);
                            glColor3f(0.0f,0.0f,0.9f);
                            for(int i=0; i<mesh.points.size(); i++){
                                Draw3D::drawVecInPos( mesh.normals[i], mesh.points[i] );
                            }
                        glEndList();
                    }
                    dragging=false;
                    break;
            }
    };
    AppSDL2OGL::eventHandling( event );
}

void TestAppMesh::drawHUD(){
    glDisable ( GL_LIGHTING );

}

// ===================== MAIN

TestAppMesh * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppMesh( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















