
// References:
//  - https://stackoverflow.com/questions/8923174/opengl-vao-best-practices
//  - http://www.swiftless.com/tutorials/opengl4/4-opengl-4-vao.html
//  - http://www.lastrayofhope.co.uk/2011/07/30/using-vertex-array-objects/
//  - http://www.openglsuperbible.com/2013/12/09/vertex-array-performance/
//  - http://antongerdelan.net/opengl/vertexbuffers.html
//  - https://en.wikipedia.org/wiki/Vertex_buffer_object


#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>
#include <GL/glew.h>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

#include "Solids.h"

//#include "SceneOGL3.h"
#include "ScreenSDL2OGL3.h"
#include "AppSDL2OGL3.h"

/*
#include "Mesh.h"
#include "Solids.h"
#include "GLfunctions.h"
#include "GLobjects.h"
#include "GLObject.h"
#include "GL3Utils.h"
*/

#include "IO_utils.h"

// ====================================
//      VAOsTestApp
// ====================================

GLuint makeVAO( int nVert, Vec3f* verts, Vec3f* colors ){

    GLuint vao;
    GLuint vbo_verts;
    GLuint vbo_colors;

    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    glGenBuffers(1, &vbo_verts );
    glBindBuffer(GL_ARRAY_BUFFER, vbo_verts);
    glBufferData(GL_ARRAY_BUFFER, nVert *3* sizeof(GLfloat), verts, GL_STATIC_DRAW);
    //glBindBuffer(GL_ARRAY_BUFFER, vbo_verts);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);

    glGenBuffers(1, &vbo_colors );
    glBindBuffer(GL_ARRAY_BUFFER, vbo_colors );
    glBufferData(GL_ARRAY_BUFFER, nVert *3* sizeof(GLfloat), colors, GL_STATIC_DRAW);
    //glBindBuffer(GL_ARRAY_BUFFER, bo_colors);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, NULL);

    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);

    return vao;
}

class VAOsTestApp : public AppSDL2OGL3, public SceneOGL3 { public:

    const Uint8 *scanKeys;
    Uint32 mouseButtons;

    float camDist = 100.0;

    Shader *sh1=0;

    GLuint vao1=0,vao2=0;
    //GLMesh *glmesh=0,*glScreenMesh=0,*glCamMesh=0;

	//virtual void update();
    virtual void eventHandling   ( const SDL_Event& event  );

    virtual void draw( Camera& cam );

	VAOsTestApp(int W, int H);

};

VAOsTestApp::VAOsTestApp(int W, int H):AppSDL2OGL3(W,H),SceneOGL3(){

    for( ScreenSDL2OGL3* screen: screens ) screen->scenes.push_back( this );

    sh1=new Shader();
    //sh1->init( "common_resources/shaders/shade3D.glslv",   "common_resources/shaders/shade3D.glslf"   );
    sh1->init( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/color3D.glslf"   );
    sh1->getDefaultUniformLocation();

    sh1->use();
    glUniform3f(sh1->getUloc("lightColor"   ), 0.5,0.45,0.4 );
    glUniform3f(sh1->getUloc("diffuseColor" ), 1.0,1.0,1.0  );
    glUniform3f(sh1->getUloc("ambientColor" ), 0.2,0.25,0.3 );
    glUniform3f(sh1->getUloc("specularColor"), 2.0,2.0,2.0  );
    glUniform3f(sh1->getUloc("lightPos"     ), 10.0,-10.0,-10.0 );
    //glUniform3f(sh1->getUloc("lightPos"     ), 10.0,10.0,10.0 );

    const int nVerts = 100;
    Vec3f verts [nVerts];
    Vec3f colors[nVerts];

    for(int i=0; i<nVerts; i++){
        float x = 0.1*i;
        verts [i].set( x, sin(x*5),cos(x*2) );
        colors[i].set( sin(x*5),0.0, 1.0 );
    }

    vao1 = makeVAO( nVerts, verts, colors );

    for(int i=0; i<nVerts; i++){
        float x = 0.1*i;
        verts [i].set( sin(x*5),x ,cos(x*2) );
        colors[i].set( sin(x*5),0.0, 1.0 );
    }

    vao2 = makeVAO( nVerts, verts, colors );

};

void VAOsTestApp::draw( Camera& cam ){

    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);

    sh1->use();
    setCamera( *sh1, cam );
    sh1->setModelPoseT( Vec3dZero, Mat3dIdentity );
    GLuint ucolor = sh1->getUloc("baseColor");
    glUniform4f( ucolor, 0.0, 0.0, 0.0, 1.0 );

    //drawElements( draw_mode, inds, nInds );

    glBindVertexArray( vao1 );
    glDrawArrays( GL_LINE_STRIP, 0, 100);

    glBindVertexArray( vao2 );
    glDrawArrays( GL_LINE_STRIP, 0, 100);

};

void VAOsTestApp::eventHandling( const SDL_Event& event  ){
    switch( event.type ){
        case SDL_KEYDOWN :
        switch( event.key.keysym.sym ){
            case SDLK_ESCAPE   : SDL_Quit(); exit(1); break;
            case SDLK_SPACE    : STOP = !STOP; printf( STOP ? " STOPED\n" : " UNSTOPED\n"); break;
        }; break;
        case SDL_QUIT: SDL_Quit(); exit(1); break;

        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_RIGHT:
                   // mouseSteer = true;
                    break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_RIGHT:
                  //  mouseSteer = false;
                 //   pilot->resetSteer();
                    break;
            }
            break;
    }

};

VAOsTestApp * app;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	SDL_ShowCursor(SDL_DISABLE);
    SDL_DisplayMode dm;
    SDL_GetDesktopDisplayMode(0, &dm);
    app = new VAOsTestApp( dm.w-150, dm.h-100 );
    app->loop( 1000000 );
    app->quit();
    return 0;
}



