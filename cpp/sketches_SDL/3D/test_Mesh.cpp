
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "SDL_utils.h"
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

class TestAppMesh : public AppSDL2OGL_3D {
	public:

    Mesh mesh;

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

TestAppMesh::TestAppMesh( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

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
    printf("initialization DONE !");

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

    //exit(0);

}

void TestAppMesh::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    glEnable( GL_LIGHTING );
    glEnable(GL_DEPTH_TEST);
    glCallList(mesh.rendered_shape);

    glDisable( GL_LIGHTING );
    Draw3D::drawAxis(1.0);

    glDisable( GL_DEPTH_TEST );
    /*
    // this should work for perspective
    Vec3d hRay; //hRay.set_lincomb(  camMat.a, camMat.b, camMat.b, );
    camMat.dot_to( {mouseX-WIDTH*0.5f, mouseY-HEIGHT*0.5f, zoom }, hRay);
    Draw3D::drawLine      ( camPos, camPos+hRay );
    Draw3D::drawPointCross( camPos+hRay, 0.2 );
    int ip = mesh.pickVertex( camPos, camMat.c );
    */

    glColor3f(0.0,0.0,1.0);
    char str[256];
    for(int i=0; i<mesh.points.size(); i++){
        sprintf(str,"%i\0",i);
        Vec3d& p = mesh.points[i];
        Draw3D::drawText(str, p, fontTex, 0.03, 0);
    }
    glColor3f(1.0,0.0,0.0);
    for(int i=0; i<mesh.polygons.size(); i++){
        Vec3d c = mesh.faceCog( i );
        sprintf(str,"%i\0",i);
        Vec3d& p = mesh.points[i];
        Draw3D::drawText(str, c, fontTex, 0.03, 0);
    }

    glColor3f(0.0,0.7,0.0);
    for(int i=0; i<mesh.edges.size(); i++){
        MeshEdge& ed = mesh.edges[i];
        Draw3D::drawLine( mesh.points[ed.verts.a], mesh.points[ed.verts.b] );
        Vec3d c = (mesh.points[ed.verts.a]+mesh.points[ed.verts.b])*0.5;
        sprintf(str,"%i\0",i);
        Vec3d& p = mesh.points[i];
        Draw3D::drawText(str, c, fontTex, 0.03, 0);
    }

    //glColor3f(0,0,0); Draw3D::drawPointCross( mesh.points[ipicked], 0.2 );

};


void TestAppMesh::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    Vec3d ray0;
                    ray0.set_lincomb( 1, mouse_begin_x, mouse_begin_y, camPos, camMat.a, camMat.b );
                    //glColor3f(1,1,1); Draw3D::drawPointCross( ray0, 0.2 );
                    ipicked= mesh.pickVertex( ray0, camMat.c );
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
                        d.set_lincomb( 1, mouse_begin_x, mouse_begin_y, camPos, camMat.a, camMat.b );
                        d.sub(mesh.points[ipicked]);
                        double c = d.dot(camMat.c);  d.add_mul(camMat.c, -c);
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
















