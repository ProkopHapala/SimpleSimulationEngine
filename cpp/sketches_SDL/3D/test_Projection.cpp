
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "Mat4.h"
#include "quaternion.h"


#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"

// ==== functions

// ======================  TestApp

const int nps = 100000;
Vec3f ps[nps];

class TestAppProjection : public AppSDL2OGL_3D {
	public:
	int ogl1=0,ogl2=0;

	bool view_projected = true;


	// ---- function declarations

	virtual void draw   ();
    virtual void eventHandling( const SDL_Event& event  );

	TestAppProjection( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppProjection::TestAppProjection( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    Mat4f m; m.getPerspectiveMatrix( -10.0f, 10.0f, -7.0f, 7.0f, 40.0f, 50.0 );

    //m.print();

    for(int i=0; i<nps; i++){
        ps[i] = (Vec3f){ randf(-30.0f,30.0f), randf(-30.0f,30.0f), randf(0.0f,100.0f) };
    }

    ogl1 = glGenLists(1);
    glNewList(ogl1, GL_COMPILE);
    glBegin(GL_POINTS);
    for(int i=0; i<nps; i++){
        Quat4f p = (Quat4f){ps[i].x,ps[i].y,ps[i].z,1.0f};
        glColor3f(0.5+(p.x/60.0), 0.5+(p.y/60.0), (p.z)/100.0);
        glVertex3f( p.x, p.y, p.z );
    }
    glEnd();
    glEndList();

    ogl2 = glGenLists(1);
    glNewList(ogl2, GL_COMPILE);
    glBegin(GL_POINTS);
    for(int i=0; i<nps; i++){
        Quat4f p = (Quat4f){ps[i].x,ps[i].y,ps[i].z,1.0f};
        Quat4f p_;
        m.dot_to( p, p_ );
        p_.mul(1/p_.w);
        glColor3f(0.5+(p.x/60.0), 0.5+(p.y/60.0), (p.z)/100.0);
        glVertex3f( p_.x, p_.y, p_.z );
        /*
        //printf( " %i (%g,%g,%g) (%g,%g,%g) \n",  i, p.x, p.y, p.z,   p_.x, p_.y, p_.z );
        //if( (p_.x>-1.0) && (p_.x<1.0) && (p_.y>-1.0) && (p_.y<1.0) && (p_.z>-1.0) && (p_.x<1.0) ){
            //glColor3f(1.0f,1.0f,1.0f); glVertex3f( p.x, p.y, -p.z );
            glColor3f(1.0f,0.0f,0.0f); glVertex3f( p_.x, p_.y, -p_.z );
            //printf( " %i (%g,%g,%g) (%g,%g,%g) TRUE  \n",  i, p.x, p.y, p.z,   p_.x, p_.y, p_.z );
        //}else{
            //printf( " %i (%g,%g,%g) (%g,%g,%g) FALSE \n",  i, p.x, p.y, p.z,   p_.x, p_.y, p_.z );
        //}
        */
    }
    glEnd();
    glEndList();

    m.print();

}

void TestAppProjection::draw   (){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glDisable(GL_LIGHTING);

	if(view_projected){glCallList( ogl2 );}else{glCallList( ogl1 );}


	Draw3D::drawAxis(10.0);

};


void TestAppProjection::eventHandling( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_m: view_projected=!view_projected; printf("view_projected %i\n", view_projected ); break;
                //case SDLK_LEFTBRACKET:  iselected++; if(iselected>=rmesh.n)iselected=rmesh.n-1; break;
                //case SDLK_RIGHTBRACKET: iselected--; if(iselected<0)iselected=0; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}

// ===================== MAIN

TestAppProjection * testApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestAppProjection( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}
















