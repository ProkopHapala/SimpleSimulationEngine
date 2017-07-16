
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

const int nps = 500000;
Vec3f ps [nps];
Vec3f ps_[nps];


inline bool clip_cond( Vec3f p ){
    return (p.x>-1.0) && (p.x<1.0) && (p.y>-1.0) && (p.y<1.0) && (p.z>-1.0) && (p.z<1.0);
    //return (p.x>-1.0) && (p.x<1.0) && (p.y>-1.0) && (p.y<1.0) && (p.z>-1.0) && (p.z<1.0);
    //return (p.x>-0.5) && (p.x<0.5) && (p.y>-0.5) && (p.y<0.5) && (p.z>0.0) && (p.z<0.5);
}


class TestAppProjection : public AppSDL2OGL_3D {
	public:
	int oAll =0,oSel =0;
	int oAll_=0,oSel_=0;

	Vec3f corner = (Vec3f){2.0,0.5,2.0f};

	bool view_projected = false;

	// ---- function declarations

	virtual void draw   ();
    virtual void eventHandling( const SDL_Event& event  );

	TestAppProjection( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppProjection::TestAppProjection( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    float far = 30.0;

    Mat4f m;
    m.getPerspectiveMatrix( -corner.x, corner.x, -corner.y, corner.y, -corner.z, -far );

    //m.print();

    for(int i=0; i<nps; i++){
        Quat4f p,p_;
        p.set( randf(-10.0f,10.0f), randf(-10.0f,10.0f), randf(-50.0f,50.0f), 1.0f );
        m.dot_to( p, p_ );
        p_.mul(1.0f/p_.w);
        ps [i].set(p.x, p.y, p.z );
        ps_[i].set(p_.x,p_.y,p_.z);
    }

    oAll = glGenLists(1);
    glNewList(oAll, GL_COMPILE);
    glBegin(GL_POINTS);
    glColor3f(0.0f,0.0f,1.0f);
    for(int i=0; i<nps; i++){
        glVertex3f( ps[i].x, ps[i].y, ps[i].z );
    }
    glEnd();
    glEndList();

    oSel = glGenLists(1);
    glNewList(oSel, GL_COMPILE);
    glBegin(GL_POINTS);
    glColor3f(1.0f,0.0f,0.0f);
    for(int i=0; i<nps; i++){
        if( clip_cond( ps_[i] ) ){
            glVertex3f( ps[i].x, ps[i].y, ps[i].z );
        }
    }
    glEnd();
    glEndList();

    oAll_ = glGenLists(1);
    glNewList(oAll_, GL_COMPILE);
    glBegin(GL_POINTS);
    glColor3f(0.0f,0.0f,1.0f);
    for(int i=0; i<nps; i++){
        glVertex3f( ps_[i].x, ps_[i].y, ps_[i].z );
    }
    glEnd();
    glEndList();

    oSel_ = glGenLists(1);
    glNewList(oSel_, GL_COMPILE);
    glBegin(GL_POINTS);
    glColor3f(1.0f,0.0f,0.0f);
    for(int i=0; i<nps; i++){
        if( clip_cond( ps_[i] ) ){
            glVertex3f( ps_[i].x, ps_[i].y, ps_[i].z );
        }
    }
    glEnd();
    glEndList();

    m.print();

}

void TestAppProjection::draw   (){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glDisable(GL_LIGHTING);
	glDisable(GL_DEPTH_TEST);

	if(view_projected){
        glCallList( oAll_ ); glCallList( oSel_ );
    }else{
        glCallList( oAll );glCallList( oSel );

        glColor3f(0.0f,0.0f,0.0f);
        Draw3D::drawPointCross( (Vec3f){-corner.x,-corner.y,corner.z}, 0.1 );
        Draw3D::drawPointCross( (Vec3f){-corner.x,+corner.y,corner.z}, 0.1 );
        Draw3D::drawPointCross( (Vec3f){+corner.x,-corner.y,corner.z}, 0.1 );
        Draw3D::drawPointCross( (Vec3f){+corner.x,+corner.y,corner.z}, 0.1 );

        /*
        glColor3f(1.0f,1.0f,1.0f);
        Draw3D::drawPointCross( (Vec3f){-corner.x,-corner.y,-corner.z}*2, 0.1 );
        Draw3D::drawPointCross( (Vec3f){-corner.x,+corner.y,-corner.z}*2, 0.1 );
        Draw3D::drawPointCross( (Vec3f){+corner.x,-corner.y,-corner.z}*2, 0.1 );
        Draw3D::drawPointCross( (Vec3f){+corner.x,+corner.y,-corner.z}*2, 0.1 );
        */
    }

	Draw3D::drawAxis(10.0);

    glColor3f(0.5f,0.0f,0.0f); Draw3D::drawScale( {0.0,0.0,0.0}, {10.0,0.0,0.0}, {0.0,1.0,1.0}, 1.0, 0.1, 0.1 );
    glColor3f(0.0f,0.5f,0.0f); Draw3D::drawScale( {0.0,0.0,0.0}, {0.0,10.0,0.0}, {1.0,0.0,1.0}, 1.0, 0.1, 0.1 );
    glColor3f(0.0f,0.0f,0.5f); Draw3D::drawScale( {0.0,0.0,0.0}, {0.0,0.0,10.0}, {1.0,1.0,0.0}, 1.0, 0.1, 0.1 );

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
















