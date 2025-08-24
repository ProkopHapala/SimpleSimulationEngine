

//#define SPEED_TEST

#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>
//#include <SDL2/SDL.h>
//#include <SDL2/SDL_opengl.h>
#include <GL/glew.h>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "Mat4.h"

#include "Solids.h"

#include "Shader.h"
#include "GLObject.h"
#include "DrawOGL3.h"
#include "SceneOGL3.h"
#include "ScreenSDL2OGL3.h"
#include "AppSDL2OGL3.h"

#include "Noise.h"
#include "Mesh.h"
#include "Solids.h"

#include "GLfunctions.h"
#include "Shader.h"
#include "GLobjects.h"
#include "GLObject.h"
#include "GL3Utils.h"

#include "TerrainOGL3.h"
#include "testUtils.h"
#include "IO_utils.h"

#include "DebugView.h"

DEBUG_VIEW_DEFINE();

class LandscapeTestApp : public AppSDL2OGL3, public SceneOGL3 { public:

    const Uint8 *scanKeys;
    Uint32 mouseButtons;

    // Bring base-class overload draw(Camera&) into scope to avoid hiding warnings
    using SceneOGL3::draw;

    TerrainOGL3       terrain1;
    TerrainOGL3_patch terrain2;

    TerrainOGL3_normals terrNormals;

    float camDist = 100.0;

    //Shader *sh1=0,*shTx=0;
    //GLMesh *glmesh=0,*glScreenMesh=0,*glCamMesh=0;

	//virtual void update();
    virtual void eventHandling   ( const SDL_Event& event  );

    //virtual void draw( Camera& cam );
    virtual void draw( );

	LandscapeTestApp(int W, int H);

};

LandscapeTestApp::LandscapeTestApp(int W, int H):AppSDL2OGL3(W,H),SceneOGL3(){

    for( ScreenSDL2OGL3* screen: screens ) screen->scenes.push_back( this );

    DEBUG_VIEW_INIT();

    int imgH = 100;
    int imgW = 100;
    float * height_map = new float[imgH*imgW];
    for( int iy=0; iy<imgH; iy++ ){
        for( int ix=0; ix<imgW; ix++ ){
            float x = ix*2*(M_PI/imgW);
            float y = iy*2*(M_PI/imgH);
            float f = cos(x*10.0)*cos(y*10.0)*0.5 + 0.5;
            //height_map[ iy*imgW + ix ] = sin(x)*sin(y)*0.5 + 0.5;
            //height_map[ iy*imgW + ix ] = cos(x*20.0)*cos(y*20.0)*0.5 + 0.5;
            //height_map[ iy*imgW + ix ] = cos(x*20.0)*0.5 + 0.5;
            //height_map[ iy*imgW + ix ] = cos(y*40.0)*0.4 + 0.5 + x*0.01;
            //height_map[ iy*imgW + ix ] = 0;
            //height_map[ iy*imgW + ix ] =randf();
            //height_map[ iy*imgW + ix ] = f*f*f*f;

            height_map[ iy*imgW + ix ] = 1/(1+f*f*10000.0);
        }
    }

    terrNormals.init( {40,40}, 100.0,  {imgW, imgH},  height_map , 1.0 );
    terrNormals.mapScale.z = 50.0;
    terrNormals.mapScale.x = 0.0002;
    terrNormals.mapScale.y = 0.0002;
    terrNormals.derivScale = 0.02;
    terrNormals.txStep = (Vec2f){ 1.0/imgW, 1.0/imgH };

    //newTexture2D( txHeight, imgW, imgH, height_map, GL_RED, GL_FLOAT );
    terrain1.init( {50,100}, 100.0,  {imgW, imgH},  height_map   );
    //terrain2.init( {50,100}, 100.0,  {imgW, imgH},  height_map   );

    terrain2.init( {40,120}, 100.0,  {imgW, imgH},  height_map , 1.0, false );
    terrain2.mapScale.z = 50.0;
    terrain2.mapScale.x = 0.0002;
    terrain2.mapScale.y = 0.0002;
    terrain2.derivScale = 0.02;
    terrain2.txStep = (Vec2f){ 1.0/imgW, 1.0/imgH };

    delete [] height_map;

    camDist = 2.0;
    screen->cam.aspect = screen->HEIGHT / (float)screen->WIDTH;

};

//void LandscapeTestApp::update(){
//    AppSDL2OGL3::update();
//}

void LandscapeTestApp::draw(  ){
    Camera cam;

    //long time_start = getCPUticks();

    glClearColor(0.8, 0.8, 0.8, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );

    glEnable( GL_DEPTH_TEST );
    glDepthFunc   ( GL_LESS );

    /*
    terrain1.pos.x = camPos.x;
    terrain1.pos.z = camPos.z;
    terrain1.setViewRange( {mRot.c.x, mRot.c.z}, 0.3 );
    terrain1.sh.use();
    terrain1.sh.set_camPos( (float*)&camPos );
    terrain1.sh.set_camMat( (float*)&camMat );
    terrain1.draw();
    */

    //terrain2.mesh->draw_mode = GL_LINE_STRIP;
    //terrain2.pos.x = cam.pos.x;
    //terrain2.pos.z = cam.pos.z;
    //terrain1.setViewRange( {mRot.c.x, mRot.c.z}, 0.3 );
    //terrain2.sh.use();
    //terrain2.sh.set_camPos( (float*)&camPos );
    //terrain2.sh.set_camMat( (float*)&camMat );
    //terrain2.draw();

    terrain2.flexibility = 0.0;
    terrain2.mesh->draw_mode = GL_TRIANGLES;
    terrain2.draw(cam);


    terrNormals.mesh->draw_mode = GL_LINES;
    terrNormals.flexibility = 20.0;
    terrNormals.draw(cam);

    glPointSize(3.0);
    terrNormals.mesh->draw_mode = GL_POINTS;
    terrNormals.flexibility = 0.0;
    terrNormals.draw(cam);


    /*
    glPointSize(3.0);
    terrain2.flexibility = 30.0;
    //terrain2.mesh->draw_mode = GL_LINE_STRIP;
    terrain2.mesh->draw_mode = GL_POINTS;
    terrain2.draw(cam);
    //SDL_GL_SwapWindow(window);
    //long time_end = getCPUticks();
    */

    DEBUG_mesh->addLine( (Vec3f){0.0,0.0,0.0},  (Vec3f){terrain2.viewMin.x*1000.0, 0.0, terrain2.viewMin.y*1000.0}, {1.0,0.0,0.0} );
    DEBUG_mesh->addLine( (Vec3f){0.0,0.0,0.0},  (Vec3f){terrain2.viewMax.x*1000.0, 0.0, terrain2.viewMax.y*1000.0}, {0.0,0.0,1.0} );

    DEBUG_drawCamera( cam, {1.0,0.0,1.0} );

    DEBUG_draw( cam, (Vec3d)terrain2.pos, Mat3dIdentity );

};

void LandscapeTestApp::eventHandling( const SDL_Event& event  ){
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

LandscapeTestApp * app;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	SDL_ShowCursor(SDL_DISABLE);
    SDL_DisplayMode dm;
    SDL_GetDesktopDisplayMode(0, &dm);
    app = new LandscapeTestApp( dm.w-150, dm.h-100 );
    app->loop( 1000000 );
    app->quit();
    return 0;
}
