
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "fastmath.h"
#include "Vec2.h"

#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"

// ======================  TestApp

class TestAppQuatRotSample : public AppSDL2OGL_3D {
	public:

	int point_cloud;
	int sphere;

	// ---- function declarations

	virtual void draw   ();

	TestAppQuatRotSample( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppQuatRotSample::TestAppQuatRotSample( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

/*
	point_cloud=glGenLists(1);
	glNewList( point_cloud, GL_COMPILE );
		glDisable ( GL_LIGHTING );
		glBegin( GL_POINTS );
		for( int i=0; i<100; i++ ){
			glColor3f( randf(-1.0,1.0), randf(-1.0,1.0), randf(-1.0,1.0) );
			glVertex3f( randf(-1.0,1.0), randf(-1.0,1.0), randf(-1.0,1.0) );
		}
		glEnd();
	glEndList();
*/

	int    n = 10;
	float  d = 1.0f/n; 
	
	point_cloud=glGenLists(1);
	glNewList( point_cloud, GL_COMPILE );
		glDisable ( GL_LIGHTING );
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
	glEndList();

	sphere=glGenLists(1);
	glNewList( sphere, GL_COMPILE );
		//glEnable    ( GL_LIGHTING );
		//glShadeModel( GL_SMOOTH   );
		glDisable ( GL_LIGHTING );
		glColor3f ( 0.5f, 0.5f, 0.5f );
		Draw3D::drawSphere_oct( 8, 0.95f, {0.0d,0.0d,0.0d} );
	glEndList();
}

void TestAppQuatRotSample::draw   (){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glCallList( sphere );

	glDisable ( GL_LIGHTING );
	glCallList( point_cloud );
	Draw3D::drawAxis ( 3.0f );

};

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
















