
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
#include "GLobjects.h"
#include "GLInstances.h"
#include "IO_utils.h"
#include "Shader.h"


#include "DrawOGL3.h"

#include "Solids.h"
#include "CMesh.h"


#include "CMesh.h"

#include "TerrainOGL3.h"
#include "HorizontOGL3.h"



#include "testUtils.h"

// =============== Global variables

Shader *shBranches,*shLeafs,*shTex,*shTexView;

//GLInstances instances;
//GLBillboards bilboards;
//GLuint texture_1;

GLMesh   *mshBranches,*mshLeafs,*msh1;

GLuint       texTest;
GLMesh       *glSprite,*glAtlas,*glAtlas_;
GLBillboards bilboards;
FrameBuffer  frameBuff1;


int frameCount = 0;
double lastTime = 0.0;

GLuint vao;     // vertex array object
int delay = 10;
int VSync = 0;
SDL_Window * window     = NULL;
SDL_GLContext   context = NULL;
int WIDTH  = 800;
int HEIGHT = 800;
float ASPECT_RATIO = HEIGHT/WIDTH;

int mouseX, mouseY;
Quat4f qCamera;
Mat3f  mouseMat;
Vec3f  camPos = (Vec3f){ 0.0f, 0.0f, -10.0f };

bool bTransparent  = false;

long lastCPUtick = 0;
double ticks_per_second=0;

// =============== Functions

int setup(){

    shBranches=new Shader();
    //shBranches->init( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/color3D.glslf"   );
    shBranches->init( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/normal2color.glslf"   );
    shBranches->getDefaultUniformLocation();

    shTexView=new Shader();
    shTexView->init( "common_resources/shaders/texture3D.glslv",   "common_resources/shaders/texture.glslf"   );
    shTexView->getDefaultUniformLocation();

    //ticks_per_second = calibrate_timer(100);
    //lastTime = glfwGetTime();
    lastTime = 0.0;
    qCamera.setOne();


    //glquad =new GLMesh();
    //glquad->init( 6, 0, NULL, DEFAULT_Bilboard_verts, NULL, NULL, DEFAULT_Bilboard_UVs );
    //glquad->init( 6, 0, NULL, DEFAULT_Bilboard_verts, DEFAULT_Bilboard_verts, NULL, NULL );

    //texTest = makeTestTextureRGBA( 256, 256);

    //frameBuff1.init( 2048, 64 );
    //frameBuff1.init( 2048, 128 );
    frameBuff1.init( 4096, 256 );

    // ===== Prepare texture by rendering

    //glBindFramebuffer(GL_FRAMEBUFFER, frameBuff1.buff );
    frameBuff1.bind();

    //glClearColor(0.0, 0.0, 0.8, 1.0);
    glClearColor(0.9, 0.9, 0.9, 0.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );

    Mat4f camMat;
    //camera( 1/16.0, camMat);
    //camera( 1.0, camMat);
    //renderPhases( camMat);
    //renderPhases2D( camMat, 16 );
    //renderPhasesOct( camMat, 8 );

    // ==== END   : RENDER TO TEXTURE

    return 0;
};

void draw( ){
    glBindFramebuffer(GL_FRAMEBUFFER, 0); glViewport(0,0,WIDTH,HEIGHT);
    glClearColor(0.8, 0.8, 0.8, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );
	// Simulate all particles

	Mat4f camMat;
	// camera( 1.0, camMat);
	//renderPhases(camMat);

    Mat3f modelMat; modelMat.setOne(); modelMat.mul(0.5f);
    //Vec3f modelPos;

    //glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, texTest );
    glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, frameBuff1.texRGB );

    // ======= Bare Texture Atlas
    Shader * sh;


    //sh = shTex;
    sh = shTexView;
    sh->use();
    sh->set_modelMat( (GLfloat*)&modelMat );
    sh->set_camPos  ( (GLfloat*)&camPos   );
    sh->set_camMat  ( (GLfloat*)&camMat   );
    glUniform1i(sh->getUloc("texture_1"), 0);

    sh->set_modelPos( (const GLfloat[]){-4.0f,-16.0f,5.0f} );
    glAtlas->draw(GL_TRIANGLES); glAtlas->drawPoints(10.0);

    sh->set_modelPos( (const GLfloat[]){-4.0f-8.0f,-24.0f,5.0f} );
    glAtlas_->draw(GL_TRIANGLES); glAtlas_->drawPoints(10.0);

    // ======= Animated Texture
    //sh = shTex;


    sh = shTex;
    sh->use();
    sh->set_modelMat( (GLfloat*)&modelMat );
    sh->set_camPos  ( (GLfloat*)&camPos   );
    sh->set_camMat  ( (GLfloat*)&camMat   );

    glUniform1i(sh->getUloc("texture_1"), 0);
    //glUniform1f(sh->getUloc("nPhases"), 32.0f);

    //glUniform2f(sh->getUloc("uv0s"), 0.0   ,0.0);


    /*
    //float the = +1.5;
    float the = -0.2;
    float phi = frameCount*0.01;
    Vec2f cst; cst.fromAngle(the);
    Vec2f uv = selectPhaseOct( (Vec3f){cst.a*cos(phi),cst.b,cst.a*sin(phi)}, 8 );
    glUniform2f(sh->getUloc("uv0"), uv.x, uv.y );
    glUniform2f(sh->getUloc("du" ), 1/16.0,0.0 );
    glUniform2f(sh->getUloc("dv" ), 0.0   ,1/16.0);

    sh->set_modelPos( (const GLfloat[]){0.0f,0.0,10.0f} );
    glSprite->draw(GL_TRIANGLES);
    glSprite->drawPoints(10.0);


    sh = shBranches;
    sh->use();
    sh->set_camPos  ( (GLfloat*)&camPos );
    sh->set_camMat  ( (GLfloat*)&camMat );
    sh->set_modelMat( (GLfloat*)&modelMat );

    //camMat.m
    Vec3f mpos = (Vec3f){uv.x*16-4+0.5,uv.y*16-16+0.5,10.0};
    sh->set_modelPos( (GLfloat*)&mpos );
    //mshBranches->draw(GL_TRIANGLES);


    */

    sh = shBranches;
    sh->use();
    sh->set_camPos  ( (GLfloat*)&camPos );
    sh->set_camMat  ( (GLfloat*)&camMat );
    sh->set_modelMat( (GLfloat*)&modelMat );

    msh1->draw(GL_TRIANGLES);



}

void init();
void quit();
void die( char const *msg );
void inputHanding();

int main(int argc, char *argv[]){
    init();
	setup();
    for ( frameCount=1; frameCount<1000000; frameCount++)    {

        draw(); SDL_GL_SwapWindow(window);
 		//if( !STOP ) draw();
		inputHanding();
        SDL_Delay(delay);
    }
    quit();
    return 0;
}

// FUNCTION ======	inputHanding
void inputHanding(){

    //float posstep = 0.1f; if(RayTerrain){ posstep = 2.0f; }
    float step          = 0.1f;
    float keyRotSpeed   = 0.01f;

    const Uint8 *keys = SDL_GetKeyboardState(NULL);

    if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.roll2  (  (float)keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.roll2  ( -(float)keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.dpitch2(  (float)keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.dpitch2( -(float)keyRotSpeed ); }

    if( keys[ SDL_SCANCODE_Q  ] ){ camPos.add_mul( mouseMat.c, +step ); }
	if( keys[ SDL_SCANCODE_E  ] ){ camPos.add_mul( mouseMat.c, -step ); }
    if( keys[ SDL_SCANCODE_W  ] ){ camPos.add_mul( mouseMat.b, +step ); }
	if( keys[ SDL_SCANCODE_S  ] ){ camPos.add_mul( mouseMat.b, -step ); }
	if( keys[ SDL_SCANCODE_A  ] ){ camPos.add_mul( mouseMat.a, -step ); }
	if( keys[ SDL_SCANCODE_D  ] ){ camPos.add_mul( mouseMat.a, +step ); }

	//if( keys[SDL_SCANCODE_W]||keys[SDL_SCANCODE_S]||keys[SDL_SCANCODE_A]||keys[ SDL_SCANCODE_D  ]  )printf( "camPos (%g,%g,%g)\n", camPos.x, camPos.y, camPos.z );

	SDL_Event event;
	while(SDL_PollEvent(&event)){
		if( event.type == SDL_KEYDOWN ){
            switch( event.key.keysym.sym ){
                case SDLK_ESCAPE: quit(); break;
                case SDLK_t: bTransparent   =!bTransparent;    break;
                //case SDLK_KP_PLUS:  terrain_size[0] *=1.1; terrain_size[2] *=1.1; break;
                //case SDLK_KP_MINUS: terrain_size[0] /=1.1; terrain_size[2] /=1.1; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            printf( "" );
		}
		if( event.type == SDL_QUIT){ quit();  };
	}

	int dmx,dmy;
	SDL_GetMouseState( &mouseX, &mouseY );
    Uint32 buttons = SDL_GetRelativeMouseState( &dmx, &dmy);
    //printf( " %i %i \n", mx,my );
    float mouseRotSpeed = 0.002;
    if ( buttons & SDL_BUTTON(SDL_BUTTON_RIGHT)) {
        Quat4f q; q.fromTrackball( 0, 0, -dmx*mouseRotSpeed, dmy*mouseRotSpeed ); qCamera.qmul_T( q );
        //qCamera.dyaw2(-dmx*mouseRotSpeed); qCamera.dpitch2(-dmy*mouseRotSpeed);
        //qCamera.dpitch2(-dmy*mouseRotSpeed); qCamera.dyaw2(-dmx*mouseRotSpeed);
        //qCamera.normalize();

        //pitch +=  dmy*mouseRotSpeed;
        //yaw   +=  dmx*mouseRotSpeed;
        qCamera.toMatrix(mouseMat);
        printf("mouseMat:\n");
        mouseMat.print();
    }
}

void init(){
    if (SDL_Init(SDL_INIT_VIDEO) < 0) die( "Unable to initialize SDL" );
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3); // Opengl 3.2
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3); // Opengl 3.2
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    //SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 32);
    window = SDL_CreateWindow("Tutorial2", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WIDTH, HEIGHT, SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN);
    if ( !window ) die("Unable to create window");
    context = SDL_GL_CreateContext( window );
    //SDL_GL_SetSwapInterval(1); // VSync On
    SDL_GL_SetSwapInterval(VSync);

    glewExperimental = true; // Needed for core profile
	if (glewInit() != GLEW_OK) {
		fprintf(stderr, "Failed to initialize GLEW\n");
		getchar();
		quit();
		//return -1;
	}

	// vertex array object
	glGenVertexArrays(1, &vao);  				// Allocate and assign a Vertex Array Object to our handle
	glBindVertexArray(vao); 					// Bind our Vertex Array Object as the current used object
}

void quit(){
	glDeleteVertexArrays(1, &vao);
    if( context != NULL ) SDL_GL_DeleteContext( context );
    if( window  != NULL ) SDL_DestroyWindow   ( window  );
    SDL_Quit();
	exit(0);
};

void die( char const *msg ){
    printf("%s: %s\n", msg, SDL_GetError());
    SDL_Quit();
    exit(1);
}
