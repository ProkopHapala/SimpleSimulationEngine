
/// @file @brief This is a math and geometry visualization that deals with mapping orientations on a sphere. It uses an icosahedron as a base and, with the help of Quaternions (`Quat4f.h`), visualizes a consistent set of local coordinate systems (u,v directions) on the surface. This is a key technique for applying seamless textures or vector fields (like wind patterns) to a sphere without singularities at the poles.
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
#include "RotationMesh.h"

#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"

// ==== functions

void quats2matrices( int n ){
    float  d = 1.0f/n;
    glBegin   ( GL_POINTS   );
    glColor3f ( 0.0f, 0.0f, 0.0f );
    for( int ix=-n; ix<n; ix++ ){
        float x = (ix+0.5f) * d;
        for( int iy=-n; iy<n; iy++ ){
            float y = (iy+0.5f) * d;
            for( int iz=-n; iz<n; iz++ ){
                float z = (iz+0.5f) * d;
                float r2 = x*x + y*y + z*z;
                if( r2 < 1.0f ){
                    Quat4f q;
                    Mat3f  M;
                    q.set( x, y, z, sqrt(1.0f - r2) );
                    q.toMatrix( M );
                    glColor3f ( M.a.x, M.a.y, M.a.z );
                    glVertex3f( (float)M.a.x, (float)M.a.y, (float)M.a.z );
                    //glVertex3f( (float)M.b.x, (float)M.b.y, (float)M.b.z );
                    //glVertex3f( (float)M.c.x, (float)M.c.y, (float)M.c.z );
                }
            }
        }
    }
    glEnd();
}

void matrices2quats( int n ){
    glBegin   ( GL_POINTS  );
    for( int i=0; i<n; i++ ){
        Mat3f  M;   M.fromRand( {randf(),randf(),randf()} );
        Quat4f q;   q.fromMatrix(M);
        glVertex3f( (float)q.x, (float)q.y, (float)q.z );
    }
    glEnd();
}


// ======================  TestApp

class TestAppQuatRotSample : public AppSDL2OGL_3D {
	public:

	int point_cloud;
	int sphere;
	int vobL;
	RotationMesh rmesh;

	int iselected = 0;

	// ---- function declarations

	virtual void draw   ();
    virtual void eventHandling( const SDL_Event& event  );

	TestAppQuatRotSample( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppQuatRotSample::TestAppQuatRotSample( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    /*
	point_cloud=glGenLists(1);
	glNewList( point_cloud, GL_COMPILE );
        glDisable ( GL_LIGHTING );
        //quats2matrices( 10 );
        glColor3f(0.0f,0.0f,0.0f); matrices2quats( 100000 );
	glEndList();
    */

    /*
	sphere=glGenLists(1);
	glNewList( sphere, GL_COMPILE );
		//glEnable    ( GL_LIGHTING );
		//glShadeModel( GL_SMOOTH   );
		glDisable ( GL_LIGHTING );
		glColor3f ( 0.5f, 0.5f, 0.5f );
		Draw3D::drawSphere_oct( 8, 0.95f, {0.0d,0.0d,0.0d} );
	glEndList();
	*/

	vobL=glGenLists(1);
    glNewList( vobL, GL_COMPILE );
		glDisable ( GL_LIGHTING );
		glBegin(GL_LINE_STRIP);
		glVertex3d( 0.0,0.0,0.0 );
		glVertex3d( 0.0,0.0,1.0 );
		glVertex3d( 0.0,0.1,1.0 );
		glEnd();
	glEndList();


	int nroll=4;
	int ndir=Solids::Icosahedron_nverts;
	rmesh.allocate(ndir*nroll,60);
	rmesh.fromDirsNroll(ndir,nroll, (Vec3d*)Solids::Icosahedron_verts, {1.0,0.0,0.0});

	int nntot = rmesh.findNeighs( 0.3 );
	printf("nntot %i \n", nntot );

    /*
    int nroll=4;
	int ndir=Solids::Cube_nverts;
	rmesh.allocate(ndir*nroll,60);
	rmesh.fromDirsNroll(ndir,nroll, (Vec3d*)Solids::Cube_verts, {1.0,0.0,0.0});
	*/
}

void TestAppQuatRotSample::draw   (){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	// draw rotation states
	glColor3f(0.9f,0.9f,0.9f);
	for(int i=0; i<rmesh.n; i++){
        Draw3D::drawShape( vobL, {0.0,0.0,0.0}, rmesh.rots[i] );
        //Draw3D::drawShape( vobL, {0.0,0.0,0.0}, rmesh.mrots[i] );
	}

    // draw neighbor connections
    glColor3f(0.0f,0.0f,0.9f);
    glBegin(GL_LINES);
    Quat4f q;
    Mat3f  mat;
    Vec3f  p;
    for(int i=0; i<rmesh.n; i++){
        //Draw3D::drawShape( vobL, {0.0,0.0,0.0}, rmesh.rots[i] );
        Vec3f pi;
        //convert(rmesh.mrots[i].c, pi );
        //printf("%i %i \n", i, rmesh.nneighs[i] );
        convert(rmesh.rots[i],q); q.toMatrix_T( mat );
        mat.dot_to( {0.0,0.1,1.0}, pi );
        //printf( "%i (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f,%3.3f)\n", i, pi.x,pi.y,pi.z,  q.x,q.y,q.z,q.w );

        //if(i==iselected){ glColor3f(0.9f,0.0f,0.9f); }else{ glColor3f(0.0f,0.0f,0.9f); };
        if(i!=iselected) continue;
        for(int j=0; j<rmesh.nneighs[i]; j++){
            int neighij = rmesh.neighs[i][j];
            Vec3f pj;
            //convert(rmesh.mrots[neighij].c, pj );
            convert(rmesh.rots[neighij],q); q.toMatrix_T( mat );
            mat.dot_to( {0.0,0.1,1.0}, pj );
            glVertex3f(pi.x,pi.y,pi.z);
            glVertex3f(pj.x,pj.y,pj.z);
        }
        //exit(0);
	}
	glEnd();
	//exit(0);

	/*
	glColor3f(0.1f,0.1f,0.1f);
	glPushMatrix();
	glScalef(0.5,0.5,0.5);
	Draw3D::drawLines    ( Solids::Icosahedron_nedges, Solids::Icosahedron_edges, Solids::Icosahedron_verts );
	glPopMatrix();
	//Draw3D::drawLines    ( Solids::Cube_nedges, Solids::Cube_edges, Solids::Cube_verts );
	*/


	//glCallList( sphere );

	//glDisable ( GL_LIGHTING );
	//glCallList( point_cloud );
	//Draw3D::drawAxis ( 3.0f );

};


void TestAppQuatRotSample::eventHandling( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_LEFTBRACKET:  iselected++; if(iselected>=rmesh.n)iselected=rmesh.n-1; break;
                case SDLK_RIGHTBRACKET: iselected--; if(iselected<0)iselected=0; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}

// ===================== MAIN

TestAppQuatRotSample * testApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestAppQuatRotSample( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}
