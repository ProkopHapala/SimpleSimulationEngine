
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
#include "SceneGraph.h"

#include "AppSDL2OGL_3D.h"

#include "GLUtils.h"

using namespace Scene;

// ======================  TestApp

class TestAppSolids : public AppSDL2OGL_3D {
	public:

    Scene::Group scene_root;

	int point_cloud;
	int shape;

	// ---- function declarations

	virtual void draw   ();

	TestAppSolids( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppSolids::TestAppSolids( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {


    //Line* l1 = new Line  ({1.0,0.0,0.0},{-1.0,0.0,0.0});
    Group* g1 = new Group(  );
    g1->objs.push_back( new Line  ({1.0,0.0,0.0},{-1.0,0.0,0.0}) );
    //g1->objs.push_back( new Sphere({0.0,0.0,2.0},0.5) );
    //g1->objs.push_back( new Cone( 6, {0.0,0.0,-1.0},{0.0,0.0,+2.0},1.0,0.2) );
    g1->objs.push_back( new Capsula( 12, {0.0,0.0,-1.0},{0.0,0.0,+2.0}, 1.0,0.5, M_PI_2, M_PI_2 ));

    /*
    //Node* n1 = new Node( g1, {00.0,0.0,0.0} );
    //Node* n2 = new Node( g1, {10.0,0.0,0.0} );
    //Node* n3 = new Node( g1, {20.0,0.0,0.0} );

    scene_root.objs.push_back( new Node( g1, {00.0,0.0,0.0} ) );
    scene_root.objs.push_back( new Node( g1, {10.0,0.0,0.0} ) );
    scene_root.objs.push_back( new Node( g1, {20.0,0.0,0.0} ) );

    //scene_root.objs.push_back( new Line  ({1.0,0.0,0.0},{-1.0,0.0,0.0}) );
    //scene_root.objs.push_back( new Sphere({0.0,2.0,0.0},0.5) );
    */

    make_circle( &scene_root, 8, M_PI, g1, {0.0,0.0,3.0}, {5.0,0.0,0.0}, {0.0,0.0,1.0} );

    shape = Draw::list();
    glColor3f(1.0f,0.0f,0.0f);
    //l1->render();
    scene_root.render();
    glEndList();

}

void TestAppSolids::draw   (){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    //GLfloat lightPos[] = { -cam.rot.c.x,-cam.rot.c.y,-cam.rot.c.z, 0.0 };
    //GLfloat lightPos[] = { 1.0, 0.0, 0.0, 0.0 };
	//glLightfv( GL_LIGHT0, GL_POSITION, lightPos );

	glCallList( shape );

	//glDisable ( GL_LIGHTING );
	//Draw3D::drawAxis ( 3.0f );

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
















