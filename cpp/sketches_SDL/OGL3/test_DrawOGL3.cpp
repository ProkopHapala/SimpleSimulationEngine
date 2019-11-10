
// Copied tutorial from
//  http://www.opengl-tutorial.org/intermediate-tutorials/billboards-particles/particles-instancing/
//  https://github.com/opengl-tutorials/ogl/tree/master/tutorial18_billboards_and_particles

#include <stdlib.h>
#include <stdio.h>

#include <vector>
#include <algorithm>

#include <GL/glew.h>
#include <SDL2/SDL.h>

#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "Mat4.h"

#include "GL3Utils.h"
#include "GLObject.h"
#include "GLfunctions.h"
#include "GLInstances.h"
#include "Shader.h"

#include "Solids.h"
#include "CMesh.h"

#include "Mesh.h"
#include "Solids.h"
#include "GLfunctions.h"
#include "GLobjects.h"
#include "GLObject.h"
#include "Shader.h"

#include "Shader.h"
#include "GLObject.h"
#include "DrawOGL3.h"
#include "SceneOGL3.h"
#include "ScreenSDL2OGL3.h"
#include "AppSDL2OGL3.h"

//#include "Draw3D_Surf.h"
//#include "SpaceCraftDraw.h"
//#include "DrawSphereMap.h"



#define _DEBUG_VIEW_ 1
#include "DebugView.h"  //do we need it ?
DEBUG_VIEW_DEFINE()


/*

Relative Cost of State Changes
https://computergraphics.stackexchange.com/questions/37/what-is-the-cost-of-changing-state

Render Tarrget (FBO)         60 000  /s
Shader Program              300 000  /s
Render Output Unit ROP    1 000 000  /s     (Change Blend Options : https://stackoverflow.com/questions/25505996/opengl-state-redundancy-elimination-tree-render-state-priorities )
Texture Bindings          1 500 000  /s
Vertex Format             2 000 000  /s
UBO Bindings              3 000 000  /s
Vertex Bindings           5 000 000  /s
Uniform Updates          10 000 0000 /s

*/














class HorizontTestApp : public AppSDL2OGL3, public SceneOGL3 { public:

    const Uint8 *scanKeys;
    Uint32 mouseButtons;

    FrameBuffer  frameBuff1;

    float camDist = 100.0;
    float screenDist = 100.0;

    //Shader *sh1=0,*shTx=0;
    //GLMesh *glmesh=0,*glScreenMesh=0,*glCamMesh=0;


    Shader *sh1,*shDebug,*shTx;
    GLMesh *glmesh,*gledges,*msh_normals, *glDebug, *glTxDebug;

	//virtual void update();
    virtual void eventHandling   ( const SDL_Event& event  );

    virtual void draw( Camera& cam );

	HorizontTestApp(int W, int H);

};

HorizontTestApp::HorizontTestApp(int W, int H):AppSDL2OGL3(W,H),SceneOGL3(){

    for( ScreenSDL2OGL3* screen: screens ) screen->scenes.push_back( this );

    DEBUG_VIEW_INIT()

    //for( ScreenSDL2OGL3* screen: screens ) screen->scenes.push_back( this );
    printf("DEBUG 3 \n");

    shDebug=new Shader();
    //shDebug->init( "common_resources/shaders/const3D.glslv",   "common_resources/shaders/const3D.glslf"   );
    shDebug->init( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/color3D.glslf"   );
    shDebug->getDefaultUniformLocation();

    sh1=new Shader();
    //sh1->init( "common_resources/shaders/const3D.glslv",   "common_resources/shaders/const3D.glslf"   );
    sh1->init( "common_resources/shaders/shade3D.glslv",   "common_resources/shaders/shade3D.glslf"   );
    sh1->getDefaultUniformLocation();

    shTx=new Shader();
    //sh1->init( "common_resources/shaders/const3D.glslv",   "common_resources/shaders/const3D.glslf"   );
    shTx->init( "common_resources/shaders/texture3D.glslv",   "common_resources/shaders/texture.glslf"   );
    shTx->getDefaultUniformLocation();

    sh1->use();
    glUniform3f(sh1->getUloc("lightColor"   ), 0.5,0.45,0.4 );
    glUniform3f(sh1->getUloc("diffuseColor" ), 1.0,1.0,1.0  );
    glUniform3f(sh1->getUloc("ambientColor" ), 0.2,0.25,0.3 );
    glUniform3f(sh1->getUloc("specularColor"), 2.0,2.0,2.0  );
    glUniform3f(sh1->getUloc("lightPos"     ), 10.0,-10.0,-10.0 );
    //glUniform3f(sh1->getUloc("lightPos"     ), 10.0,10.0,10.0 );

    GLMeshBuilder mshDebug;
    mshDebug.addLine      ( (Vec3f){0.0,0.0,0.0}, {10.0,10.0,10.0}, {1.0,0.0,0.0} );
    mshDebug.addPointCross( {0.0,0.0,0.0}, 1.0, {0.0,0.0,1.0} );
    glDebug = mshDebug.makeLineMesh();

    // ========== Mesh Builder

    GLMeshBuilder mshbuild;

    //UVFunc2smooth( {10,10}, {0.0,-M_PI*0.0}, {M_PI*1.0,M_PI*0.5}, uvfunc , mshbuild );
    //UVFunc2wire( {10,10}, {0.0,-M_PI*0.5}, {M_PI*2.0,M_PI*0.5}, uvfunc , mshbuild );
    //UVFunc2wire( {10,10}, {0.0,-M_PI*0.5}, {M_PI*1.8,M_PI*0.4}, uvfunc , mshbuild );
    //UVFunc2smooth( {10,10}, {0.0,-M_PI*0.5}, {M_PI*1.8,M_PI*0.4}, uvfunc, mshbuild );
    //UVFunc2smooth( {20,20}, {0.0,-M_PI*0.499}, {M_PI*2.0,M_PI*0.499}, uvfunc, mshbuild );
    //UVFunc2smooth( {32,16}, {0.0,0.0}, {1.0,M_PI*2.00}, uvfunc, mshbuild );
    //UVFunc2smooth( {32,16}, {0.01,0.0}, {0.99,M_PI}, uvfunc, mshbuild );

    //Cone2Mesh( {20,20}, {0.01,0.0}, {0.99,M_PI*2.0},     1.0,0.2,3.0, false, mshbuild );
    //Sphere2Mesh  ( {20,20}, {-M_PI*0.4,0.0}, {M_PI*0.4,M_PI},     1.0        , false, mshbuild );
    //Torus2Mesh   ( {20,20}, {0.0,0.0},       {M_PI*2.0,M_PI*2.0}, 0.5,1.5,     false, mshbuild );
    Teardrop2Mesh( {20,16}, {0.01,0.0}, {0.99,M_PI*2.0},   0.8,0.2,6.0, 0.5, false, mshbuild );    mshbuild.moveSub ( 0, {0.0,0.0,1.0} );
    //Sphere2Mesh  ( {20,20}, {-M_PI*0.49,0.0}, {M_PI*0.49,M_PI},     1.0 , false, mshbuild ); mshbuild.scaleSub( 1, {0.5,0.5,1.0} );
    Teardrop2Mesh( {10,8}, {0.01,0.0}, {0.99,M_PI}, 0.5,0.1,1.0, 0.0, false, mshbuild ); mshbuild.moveSub( 1, {0.0,0.45,-1.5} );
    float naca1[4]={2.0,0.15,0.0,0.0};
    float naca2[4]={1.0,0.15,0.0,0.0};
    NACASegment2Mesh( {20,2}, {-1.0,0.0}, {1.0,1.0}, naca1, naca2, 5.0, 0.0, false, mshbuild );

    //NACASegment2Mesh( {20,2}, {-1.0,0.0}, {1.0,1.0}, naca1,naca2, 5.0, false, mshbuild );
    mshbuild.duplicateSub( 2 ); mshbuild.scaleSub( 3, {-1.0,1.0,1.0} );
    mshbuild.duplicateSub( 2 ); mshbuild.scaleSub( 4, { 0.4,0.5,0.5} ); mshbuild.moveSub( 4, {0.0,0.0,-4.0} );
    mshbuild.duplicateSub( 2 ); mshbuild.scaleSub( 5, {-0.4,0.5,0.5} ); mshbuild.moveSub( 5, {0.0,0.0,-4.0} );
    mshbuild.duplicateSub( 2 ); mshbuild.scaleSub( 6, {-0.4,0.5,0.5} ); mshbuild.moveSub( 6, {0.0,0.0,-4.0} );  mshbuild.rotateSub( 6,{0.0,0.0,0.0},{0.0,0.0,1.0}, -M_PI*0.5 );


    Parabola2Mesh( {20,10}, {-1.5,-0.5*M_PI}, {1.0,0}, 2.0, 4.0, 0.0, false, mshbuild );

    Hyperbola2Mesh( {20,20}, {-1.5,0.0}, {1.0,M_PI}, 1.5, 2.0, 4.0, 0.0, false, mshbuild );



    //mshbuild.moveSub( 0, {1.0,2.0,3.0} );
    //mshbuild.rotateSub( 0, {0.0,0.0,0.0}, {1.0,0.0,0.0}, M_PI*0.5 );
    //mshbuild.scaleSub( 0, {1.0,0.5,0.25} );

    glmesh = mshbuild.makeGLmesh();
    //glmesh = mshbuild.normals2GLmesh(0.1);
    //glmesh->draw_mode = GL_LINES;
    msh_normals = mshbuild.normals2GLmesh(0.1);
    //glmesh = new GLMesh();
    //glmesh->init_wireframe( Solids::Octahedron );

    /*
    Camera& cam = screens[0]->cam;
    cam.zmin = 1.0; cam.zmax = 1000.0; cam.zoom = 5.00f;
    cam.aspect = screens[0]->HEIGHT/(float)screens[0]->WIDTH;
    screen->camLookAt = new Vec3f(); (*screen->camLookAt)={0.0,0.0,0.0};
    */


    Camera& cam = screens[0]->cam;
    cam.zmin = 1.0;
    cam.zmax = 1000.0;
    cam.zoom = 5.00f;
    cam.aspect = screens[0]->HEIGHT/(float)screens[0]->WIDTH;


    /*
    glTxDebug = new GLMesh();
    //DEFAULT_Bilboard_verts, DEFAULT_Bilboard_verts[]
    //glTxDebug->init( 6, 0,  NULL, DEFAULT_Bilboard_verts, NULL, NULL, DEFAULT_Bilboard_UVs);
    glTxDebug->init( 6, 0,  NULL, DEFAULT_Bilboard_verts, NULL, NULL, DEFAULT_Bilboard_UVs_2x2);
    */

    //return 0;

};



void HorizontTestApp::draw( Camera& cam ){

    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);

    //DEBUG_mesh->clear();
    /*
    for(int i=0; i<12; i++){ DEBUG_mesh->addLine(
        ((Vec3d){randf(),randf(),randf()})+myCraft->pos,
        ((Vec3d){randf(),randf(),randf()})+myCraft->pos,
        {randf(),randf(),randf()} );
    }
    */

    sh1->use();
    //Mat3f mrot; mrot.setOne();

    //cam.lookAt( myCraft->pos, 20.0 );
    cam.lookAt( (Vec3f){0.0,0.0,0.0}, 20.0 );

    //cam.lookAt( (Vec3d){0.0,0.0,0.0}, 20.0 );
    //setCamera( *sh1, cam );

    setCamera( *sh1, cam );
    sh1->setModelPoseT( Vec3dZero, Mat3dIdentity );


    //cam.rot = (Mat3f)myCraft->rotMat;
    //cam.lookAt( myCraft->pos, 20.0 );
    //sh1->setModelPose( myCraft->pos, myCraft->rotMat );

    int narg;
    GLuint ucolor = sh1->getUloc("baseColor");
    glUniform4f( ucolor, 0.0, 0.0, 0.0, 1.0 );
    glmesh->draw();
    /*
    narg = glmesh->preDraw();
    glmesh->drawRaw(GL_TRIANGLES, (frameCount*9)%glmesh->nInds , glmesh->nInds );
    glmesh->postDraw( narg );
    */

    glUniform4f( ucolor, 1.0, 0.0, 0.0, 1.0 );
    msh_normals->draw();
    /*
    narg = msh_normals->preDraw();
    msh_normals->drawRaw(GL_LINES, (frameCount*8)%msh_normals->nVerts , msh_normals->nVerts );
    msh_normals->postDraw( narg );
    */

    /*
    shDebug->use();
    setCamera( *shDebug, cam );
    //shDebug->setModelPoseT( myCraft->pos, myCraft->rotMat );
    glDebug->draw();
    */

    /*
    shTx->use();
    setCamera(*shTx, cam);
    //shTx->setModelPoseT( myCraft->pos, myCraft->rotMat );
    /glTxDebug->draw();
    */

    //shTx->setModelPoseT( myCraft->pos, Mat3dIdentity*10.0 );

    //DEBUG_draw(cam,myCraft->pos,myCraft->rotMat);
    //DEBUG_draw(cam,Vec3dZero,Mat3dIdentity);
};

void HorizontTestApp::eventHandling( const SDL_Event& event  ){
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




HorizontTestApp * app;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	SDL_ShowCursor(SDL_DISABLE);
    SDL_DisplayMode dm;
    SDL_GetDesktopDisplayMode(0, &dm);
    app = new HorizontTestApp( dm.w-150, dm.h-100 );
    app->loop( 1000000 );
    app->quit();
    return 0;
}
