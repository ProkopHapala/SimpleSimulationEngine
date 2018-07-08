

//#define SPEED_TEST

#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>
//#include <SDL2/SDL.h>
//#include <SDL2/SDL_opengl.h>
#include <GL/glew.h>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

#include "Solids.h"

#include "Shader.h"
#include "GLObject.h"
#include "DrawOGL3.h"
#include "SceneOGL3.h"
#include "ScreenSDL2OGL3.h"
#include "AppSDL2OGL3.h"

#include "Mesh.h"
#include "Solids.h"
#include "GLfunctions.h"
#include "GLobjects.h"
#include "GLObject.h"
#include "GL3Utils.h"


//#include "SimplexRuler.h"
//#include "Ruler2DFast.h"
//#include "TerrainHydraulics.h"
//#include "Terrain25D.h"
//#include "TerrainOGL3.h"

#include "IO_utils.h"

// ====================================
//      AeroCraftGUI
// ====================================


void drawMeshArray( Vec2i ns, Vec2d sz, const GLMesh& glmsh, const Shader& sh,  float zCut , bool side ){
    for(int ix=-ns.x; ix<ns.x; ix++){
        for(int iz=0; iz<ns.y; iz++){
            Vec3d pos = Vec3d{ix*sz.x,0.0,iz*sz.y};
            if( (pos.z<zCut) != side ){
                //shTx->setModelPoseT( Vec3dZero, Mat3dIdentity );
                pos.y = sin(iz*0.5)*cos(ix*0.5)*sin(iz*0.1)*cos(ix*0.1)*5.0   - 5.0;
                sh.set_modelPos( pos );
                glmsh.draw();
            }
        }
    }
}

class AeroCraftGUI : public AppSDL2OGL3, public SceneOGL3 { public:

    const Uint8 *scanKeys;
    Uint32 mouseButtons;

    FrameBuffer  frameBuff1;

    // Terrain - maybe move to Shooter
    //SimplexRuler       ruler;
    //Ruler2DFast        square_ruler;
    //TerrainHydraulics  hydraulics;
    //HydraulicGrid2D      hydraulics;
    //double * ground    = NULL;
    //TerrainOGL3 terrain1;

    float camDist = 100.0;

    Shader *sh1=0,*shTx=0;
    GLMesh *glmesh=0,*glScreenMesh=0,*glCamMesh=0;

	//virtual void update();
    virtual void eventHandling   ( const SDL_Event& event  );

    virtual void draw( Camera& cam );

	AeroCraftGUI(int W, int H);

};

AeroCraftGUI::AeroCraftGUI(int W, int H):AppSDL2OGL3(W,H),SceneOGL3(){

    //DEBUG_mesh = new GLMeshBuilder();
    //DEBUG_VIEW_INIT()

    for( ScreenSDL2OGL3* screen: screens ) screen->scenes.push_back( this );

    //shDebug=new Shader();
    //shDebug->init( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/color3D.glslf"   );
    //shDebug->getDefaultUniformLocation();

    sh1=new Shader();
    sh1->init( "common_resources/shaders/shade3D.glslv",   "common_resources/shaders/shade3D.glslf"   );
    sh1->getDefaultUniformLocation();

    shTx=new Shader();
    shTx->init( "common_resources/shaders/texture3D.glslv",   "common_resources/shaders/texture.glslf"   );
    shTx->getDefaultUniformLocation();

    sh1->use();
    glUniform3f(sh1->getUloc("lightColor"   ), 0.5,0.45,0.4 );
    glUniform3f(sh1->getUloc("diffuseColor" ), 1.0,1.0,1.0  );
    glUniform3f(sh1->getUloc("ambientColor" ), 0.2,0.25,0.3 );
    glUniform3f(sh1->getUloc("specularColor"), 2.0,2.0,2.0  );
    glUniform3f(sh1->getUloc("lightPos"     ), 10.0,-10.0,-10.0 );
    //glUniform3f(sh1->getUloc("lightPos"     ), 10.0,10.0,10.0 );





    glmesh = new GLMesh();
    glmesh->init_hardface( Solids::Cube );

    //msh_normals = mshbuild.normals2GLmesh(0.1);

    Camera& cam = screens[0]->cam;
    cam.zmin = 1.0;
    cam.zmax = 1000.0;
    cam.zoom = 5.00f;
    cam.aspect = screens[0]->HEIGHT/(float)screens[0]->WIDTH;

    // fix look-at point
    //screen->camLookAt = new Vec3f();
    //(*screen->camLookAt)={0.0,0.0,0.0};

    glCamMesh = cam2mesh( cam );


    float sc = 1000.0;
    float z = sc*1.20710678119; // octagon    a/2 + a*sqrt(2)
    float x = sc;
    float y = sc*0.5;
    float Bilboard_verts[] = {
        x*-0.5f,y*-0.5f,z*1.0f,   x*+0.5f,y*-0.5f,z*1.0f,   x*-0.5f,y*+0.5f,z*1.0f,
        x*+0.5f,y*+0.5f,z*1.0f,   x*+0.5f,y*-0.5f,z*1.0f,   x*-0.5f,y*+0.5f,z*1.0f
    };

    float Bilboard_UVs[] = {
        0.0f,0.0f,   1.0f,0.0f,   0.0f,1.0f,
        1.0f,1.0f,   1.0f,0.0f,   0.0f,1.0f
    };

    glScreenMesh = new GLMesh();
    glScreenMesh->init( 6, 0,  NULL, Bilboard_verts, NULL, NULL, Bilboard_UVs );

    frameBuff1.init( 2048, 1024 );

    // ===== Prepare texture by rendering
    //glBindFramebuffer(GL_FRAMEBUFFER, frameBuff1.buff );
    frameBuff1.bind();

    glClearColor( 0.6f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glEnable(GL_DEPTH_TEST);

    sh1->use();
    cam.lookAt( Vec3dZ, 20.0 );
    setCamera( *sh1, cam );
    sh1->setModelPoseT( Vec3dZero, Mat3dIdentity );
    drawMeshArray( {20,200}, {3.0,3.0} , *glmesh, *sh1,  100.0,  true );


};

//void AeroCraftGUI::update(){
//    AppSDL2OGL3::update();
//}

void AeroCraftGUI::draw( Camera& cam ){

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glViewport(0,0,screen->WIDTH,screen->HEIGHT);
    glScissor (0,0,screen->WIDTH,screen->HEIGHT);

    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);

    sh1->use();
    //cam.lookAt( Vec3dZ, 20.0 );
    setCamera( *sh1, cam );
    //setCameraOrtho(*sh1,cam);
    //sh1->setModelPose( myCraft->pos, myCraft->rotMat );
    sh1->setModelPoseT( Vec3dZero, Mat3dIdentity );

    GLuint ucolor = sh1->getUloc("baseColor");
    glUniform4f( ucolor, 0.0, 0.0, 0.0, 1.0 );
    glmesh->draw();
    drawMeshArray( {20,200}, {3.0,3.0} , *glmesh, *sh1,  100.0,  false );

    glUniform4f( ucolor, 1.0, 0.0, 0.0, 1.0 );
    sh1->setModelPoseT( Vec3dZero, Mat3dIdentity );
    glCamMesh->draw();

    //glmesh->draw(GL_TRIANGLES);
    //printf( "%f %f %f \n", cam.pos.x, cam.pos.y, cam.pos.z );
    shTx->use();

    glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, frameBuff1.texRGB );

    setCamera(*shTx, cam);
    shTx->setModelPoseT( Vec3dZero, Mat3dIdentity );
    glScreenMesh->draw();

    //DEBUG_draw(cam,myCraft->pos,myCraft->rotMat);
    //DEBUG_draw(cam,Vec3dZero,Mat3dIdentity);

};

void AeroCraftGUI::eventHandling( const SDL_Event& event  ){
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

AeroCraftGUI * app;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	SDL_ShowCursor(SDL_DISABLE);
    SDL_DisplayMode dm;
    SDL_GetDesktopDisplayMode(0, &dm);
    app = new AeroCraftGUI( dm.w-150, dm.h-100 );
    app->loop( 1000000 );
    app->quit();
    return 0;
}



