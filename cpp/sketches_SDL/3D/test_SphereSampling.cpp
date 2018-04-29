
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "fastmath.h"
#include "Vec2.h"
#include "Solids.h"

#include "Draw.h"
#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"

#include "SDL_utils.h"
#include "GLUtils.h"

// ======================  TestApp
const int nsamp = 32;
const int npix  = 8*nsamp*nsamp;
uint32_t pix   [npix];


void sampleOctahedron( Vec3d p, uint8_t& iface, double& a, double& b ){
    //iface = (p.x>0) | ( ((uint8_t)(p.y>0))<<1) | ( ((uint8_t)(p.z>0) )<<2);
    iface = ( ((uint8_t)(p.x>0))<<2) | ( ((uint8_t)(p.y>0))<<1) | ( ((uint8_t)(p.z>0) ));
    double renorm = 1.0d/( fabs(p.x) + fabs(p.y) + fabs(p.z) );
    a = fabs(p.x);
    b = fabs(p.y);
}

void sampleTri( Vec3d p, Vec3i tri, Vec3d* verts, Vec3d& c ){
    c.a = p.dot( verts[tri.a] );
    c.b = p.dot( verts[tri.b] );
    c.c = p.dot( verts[tri.c] );
    c.mul( 1/(c.a + c.b + c.c) );
}

class TestAppSolids : public AppSDL2OGL_3D {
	public:

	int point_cloud;
	int shape;

	char str[2048];
	int  fontTex;

	// ---- function declarations

	virtual void draw   ();

	TestAppSolids( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppSolids::TestAppSolids( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
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

    for( int i=0; i<npix; i++ ){ pix[i]=0; }

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
        glPopMatrix();

        glPushMatrix();
        glEnable    ( GL_LIGHTING );
        glShadeModel( GL_FLAT     );
        glColor3f( 0.8f, 0.8f, 0.8f );
        //Draw3D::drawTriangles( Solids::Tetrahedron_ntris,  Solids::Tetrahedron_tris,  Solids::Tetrahedron_verts );

        //Draw3D::drawPolygons( Solids::Tetrahedron_nfaces, Solids::Tetrahedron_ngons, Solids::Tetrahedron_faces, Solids::Tetrahedron_verts );
        //glTranslatef(  2.0f, 0.0f, 0.0f );
        Draw3D::drawPolygons( Solids::Octahedron_nfaces,  Solids::Octahedron_ngons,  Solids::Octahedron_faces,  Solids::Octahedron_verts  );
        //glTranslatef(  -4.0f, 0.0f, 0.0f );
        //Draw3D::drawPolygons( Solids::Cube_nfaces,        Solids::Cube_ngons,        Solids::Cube_faces,        Solids::Cube_verts        );

        //Draw3D::drawPolygons( Solids::RhombicDodecahedron_nfaces,        Solids::RhombicDodecahedron_ngons,        Solids::RhombicDodecahedron_faces,        Solids::RhombicDodecahedron_verts        );
        //Draw3D::drawPolygons( Solids::Icosahedron_nfaces,        Solids::Icosahedron_ngons,        Solids::Icosahedron_faces,        Solids::Icosahedron_verts        );

        /*
        int nVert = countVerts( Solids::Icosahedron_nfaces, Solids::Icosahedron_ngons );
        GLfloat * verts   = new GLfloat[nVert*3];
        GLfloat * normals = new GLfloat[nVert*3];
        hardFace( Solids::Icosahedron_nfaces, Solids::Icosahedron_ngons, Solids::Icosahedron_faces, Solids::Icosahedron_verts, verts, normals );
        Vec3f * verts_ = (Vec3f*)verts;
        Vec3f * normals_ = (Vec3f*)normals;
        for(int i=0; i<nVert; i++){
            //Draw3D::drawVecInPos( ((Vec3f*)normals)[i], ((Vec3f*)verts)[i] );
            Draw3D::drawVecInPos( normals_[i], verts_[i] );
            printf("%g %g %g   %g %g %g\n", verts_[i].x, verts_[i].y, verts_[i].z,   normals_[i].x,normals_[i].y,normals_[i].z );
        }
        */

        glPopMatrix();

	glEndList();

	zoom = 3.0;

}

void TestAppSolids::draw   (){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glEnable(GL_DEPTH_TEST);
	glCallList( shape );


	glColor3f( 0.0f,1.0f,1.0f );
	Vec3d p = camMat.c*-3.0;
	Draw3D::drawPointCross( p, 0.1 );

	uint8_t iface; double a,b;
	sampleOctahedron( p, iface, a, b );

	//printf( "iface %i | (%f,%f,%f) \n", iface, p.x, p.y, p.z );

	Vec3i ivs = Solids::Octahedron_tris[iface];
	Vec3d c;
	sampleTri( p, ivs, Solids::Octahedron_verts, c );

    if ( SDL_GetMouseState(NULL, NULL) & SDL_BUTTON(SDL_BUTTON_LEFT) ) {
        //SDL_Log("Mouse Button 1 (left) is pressed.");
        uint32_t * pixF = pix + nsamp*nsamp*iface;
        //float step = 1.0/nsamp;
        int ia = (int)(c.a*nsamp);
        int ib = (int)(c.b*nsamp);
        int i = ia*nsamp + ib;
        pixF[i] = 0xffffffff;
    }


	double fsc = 1.1;
	Vec3d p_;
	p_.set_lincomb( c.a,c.b,c.c,  Solids::Octahedron_verts[ivs.a], Solids::Octahedron_verts[ivs.b], Solids::Octahedron_verts[ivs.c] );
	p_.mul(fsc);
	Draw3D::drawLine( Solids::Octahedron_verts[ivs.a]*fsc, p_ );
	Draw3D::drawLine( Solids::Octahedron_verts[ivs.b]*fsc, p_ );
	Draw3D::drawLine( Solids::Octahedron_verts[ivs.c]*fsc, p_ );





	glColor3f( 1.0f,0.0f,1.0f );
	fsc = 1.1;
    //Draw3D::drawTriangle( Solids::Octahedron_verts[ivs.a]*fsc, Solids::Octahedron_verts[ivs.b]*fsc, Solids::Octahedron_verts[ivs.c]*fsc );

    glDisable ( GL_LIGHTING );

    glColor3f( 1.0f,0.0f,0.0f );

    for( int i=0; i<Solids::Octahedron_ntris; i++ ){
        ivs = Solids::Octahedron_tris[i];
        //p   = ( Solids::Octahedron_verts[ivs.a] + Solids::Octahedron_verts[ivs.b] + Solids::Octahedron_verts[ivs.c] ) * 0.3;
        //sprintf( str, "%i (%i,%i,%i)", i, (int)(p.x>0), (int)(p.y>0), (int)(p.z>0) );
        //sprintf( str, "%i_%i%i%i", i, (int)(p.x>0), (int)(p.y>0), (int)(p.z>0) );
        //Draw3D::drawText(str, p, fontTex, 0.03, 0);

        uint32_t * pixF = pix + nsamp*nsamp*i;

        glBegin(GL_POINTS);
        Vec3d& a = Solids::Octahedron_verts[ivs.a];
        Vec3d& b = Solids::Octahedron_verts[ivs.b];
        Vec3d& c = Solids::Octahedron_verts[ivs.c];
        float step = 1.0/nsamp;
        for( int ia = 0; ia<nsamp; ia++ ){
            float ca = ia*step;
            for( int ib = 0; ib<(nsamp-ia); ib++ ){
                float cb = ib*step;
                float cc = 1-ca-cb;
                p = (a*ca + b*cb + c*cc);
                p.normalize()*2.0;
                Draw::setRGB ( pixF[ ia*nsamp + ib ] );
                glVertex3f( (float)p.x, (float)p.y, (float)p.z );
                //printf( "iab %i %i %i |  (%f,%f,%f) \n", i, ia, ib, p.x, p.y, p.z );
            }
        }
        glEnd();

    }
    //exit(0);

    glColor3f( 0.0f,1.0f,0.0f );
    for( int i=0; i<Solids::Octahedron_nverts; i++ ){
        sprintf( str, "%i", i );
        Draw3D::drawText(str, Solids::Octahedron_verts[i], fontTex, 0.03, 0);
    }


    //if( SDL_MOUSEBUTTONDOWN ){ printf("mouse down %i \n", SDL_MOUSEBUTTONDOWN  ); };



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
















