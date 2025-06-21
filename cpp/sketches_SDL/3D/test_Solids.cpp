
/// @file @brief  A simple demo program for displaying and testing the generation of primitive shapes. It uses `Solids.h` to create various polyhedra, including the Platonic solids (cube, icosahedron, etc.), and renders them using a basic mesh class (`CMesh.h`). It's a good starting point for verifying that the 3D rendering pipeline is working correctly.
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "fastmath.h"
#include "Vec2.h"
#include "CMesh.h"
#include "Solids.h"

#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"

#include "GLUtils.h"

// ======================  TestApp

class TestAppSolids : public AppSDL2OGL_3D {
	public:

	int point_cloud;
	int shape;

	// ---- function declarations

	virtual void draw   ();

	TestAppSolids( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppSolids::TestAppSolids( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

/*
    //normals =
	shape=glGenLists(1);
	glNewList( shape, GL_COMPILE );
        printf( " Solids::nTetrahedron_tris %i \n", Solids::tetrahedron.nTris );
        glEnable    ( GL_LIGHTING );
        glShadeModel( GL_FLAT     );
        glColor3f( 0.8, 0.8, 0.8 );
        Draw3D::drawTriangles( Solids::tetrahedron.nTris, (int*)(&Solids::tetrahedron.tris[0][0]), Solids::tetrahedron.verts );
	glEndList();
*/

    //CMesh msh = Solids::RhombicDodecahedron;

    //const CMesh& msh = Solids::Tetrahedron;
    //const CMesh& msh = Solids::Octahedron;
    //const CMesh& msh = Solids::Cube;
    //const CMesh& msh = Solids::Icosahedron;
    const CMesh& msh = Solids::RhombicDodecahedron;


	shape=glGenLists(1);
	glNewList( shape, GL_COMPILE );
        //printf( " Solids::nTetrahedron_tris %i \n", Solids::nTetrahedron_tris );

        glPushMatrix();
        glDisable ( GL_LIGHTING );
        glColor3f( 1.0f, 0.0f, 1.0f );
        /*
        Draw3D::drawLines    ( Solids::Tetrahedron_nedges, Solids::Tetrahedron_edges, Solids::Tetrahedron_verts );
        glTranslatef(  2.0f, 0.0f, 0.0f );
        Draw3D::drawLines    ( Solids::Octahedron_nedges, Solids::Octahedron_edges, Solids::Octahedron_verts );
        glTranslatef( -4.0f, 0.0f, 0.0f );
        Draw3D::drawLines    ( Solids::Cube_nedges, Solids::Cube_edges, Solids::Cube_verts );
        */
        //Draw3D::drawLines    ( Solids::RhombicDodecahedron_nedges, Solids::RhombicDodecahedron_edges, Solids::RhombicDodecahedron_verts );

        //Draw3D::drawLines    ( Solids::Icosahedron_nedges, (int*)Solids::Icosahedron_edges, Solids::Icosahedron_verts );

        Draw3D::drawLines( msh.nedge, (int*)msh.edges, msh.verts );

        glPopMatrix();

        glPushMatrix();
        glEnable    ( GL_LIGHTING );
        glShadeModel( GL_FLAT     );
        glColor3f( 0.8f, 0.8f, 0.8f );
        //Draw3D::drawTriangles( Solids::Tetrahedron_ntris,  Solids::Tetrahedron_tris,  Solids::Tetrahedron_verts );
        /*
        Draw3D::drawPolygons( Solids::Tetrahedron_nfaces, Solids::Tetrahedron_ngons, Solids::Tetrahedron_faces, Solids::Tetrahedron_verts );
        glTranslatef(  2.0f, 0.0f, 0.0f );
        Draw3D::drawPolygons( Solids::Octahedron_nfaces,  Solids::Octahedron_ngons,  Solids::Octahedron_faces,  Solids::Octahedron_verts  );
        glTranslatef(  -4.0f, 0.0f, 0.0f );
        Draw3D::drawPolygons( Solids::Cube_nfaces,        Solids::Cube_ngons,        Solids::Cube_faces,        Solids::Cube_verts        );
        */
        //Draw3D::drawPolygons( Solids::RhombicDodecahedron_nfaces,        Solids::RhombicDodecahedron_ngons,        Solids::RhombicDodecahedron_faces,        Solids::RhombicDodecahedron_verts        );

        Draw3D::drawPolygons( msh.nfaces,  msh.ngons, msh.faces, msh.verts );

        //Draw3D::drawPolygons( Solids::Icosahedron_nfaces,        Solids::Icosahedron_ngons,        Solids::Icosahedron_faces,        Solids::Icosahedron_verts        );

        /*
        int nVert = countVerts( Solids::Icosahedron_nfaces, Solids::Icosahedron_ngons );
        GLfloat * verts   = new GLfloat[nVert*3];
        GLfloat * normals = new GLfloat[nVert*3];
        hardFace( Solids::Icosahedron_nfaces, Solids::Icosahedron_ngons, Solids::Icosahedron_faces, Solids::Icosahedron_verts, verts, normals );
        */

        int nVert = countVerts( msh.nfaces, msh.ngons );
        GLfloat * verts   = new GLfloat[nVert*3];
        GLfloat * normals = new GLfloat[nVert*3];
        hardFace( msh.nfaces, msh.ngons, msh.faces, msh.verts, verts, normals );


        Vec3f * verts_   = (Vec3f*)verts;
        Vec3f * normals_ = (Vec3f*)normals;
        for(int i=0; i<nVert; i++){
            //Draw3D::drawVecInPos( ((Vec3f*)normals)[i], ((Vec3f*)verts)[i] );
            Draw3D::drawVecInPos( normals_[i], verts_[i] );
            printf("%g %g %g   %g %g %g\n", verts_[i].x, verts_[i].y, verts_[i].z,   normals_[i].x,normals_[i].y,normals_[i].z );
        }

        glPopMatrix();

	glEndList();


}

void TestAppSolids::draw   (){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );


    //GLfloat lightPos[] = {cam.pos.x*-100.0,cam.pos.y*-100.0,cam.pos.z*-100.0};
    //GLfloat lightPos[] = {cam.pos.x*100.0,cam.pos.y*100.0,cam.pos.z*100.0};
    //printf( "(%f,%f,%f)\n", cam.pos.x, cam.pos.y, cam.pos.z );
    //printf( "(%f,%f,%f)\n", cam.rot.b.x,cam.rot.b.y,cam.rot.b.z );
    //GLfloat lightPos[] = { -cam.rot.b.x,-cam.rot.b.y,-cam.rot.b.z, 0.0 };
    GLfloat lightPos[] = { -cam.rot.c.x,-cam.rot.c.y,-cam.rot.c.z, 0.0 };
    //GLfloat lightPos[] = { 1.0, 0.0, 0.0, 0.0 };
	glLightfv( GL_LIGHT0, GL_POSITION, lightPos );


	glCallList( shape );

	glDisable ( GL_LIGHTING );

	Draw3D::drawAxis ( 3.0f );

};

// ===================== MAIN

TestAppSolids * testApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestAppSolids( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}
