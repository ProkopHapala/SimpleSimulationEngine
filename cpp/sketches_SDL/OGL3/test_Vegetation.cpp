
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
//#include "GLInstances.h"
#include "IO_utils.h"
#include "Shader.h"

#include "Solids.h"
#include "CMesh.h"

#include "testUtils.h"

// =============== Global variables

Shader *shBranches,*shLeafs;

//GLInstances instances;
//GLBillboards bilboards;
//GLuint texture_1;

GLMesh   *mshBranches,*mshLeafs;


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
Vec3f  camPos = (Vec3f){ 0.0f, 0.0f, 0.0f };

bool bTransparent  = false;

long lastCPUtick = 0;
double ticks_per_second=0;

// =============== Functions

void tree_step( int level, Vec3f pos, Vec3f dir, std::vector<Vec3f>& branches, std::vector<Vec3f>& leafs ){
    static const float drnd = 0.6;
    //dir.x *= randf(1.0-drnd,1.0);
    //dir.y *= randf(1.0-drnd,1.0);
    //dir.z *= randf(1.0-drnd,1.0);
    float l = dir.norm();
    dir.add( randf(-drnd,drnd)*l, randf(-drnd,drnd)*l, randf(-drnd,drnd)*l );
    dir.mul( randf(0.5,0.9) );
    Vec3f pos_ = pos + dir;
    branches.push_back(pos );
    branches.push_back(pos_);
    if( level==0 ){
        leafs.push_back(pos_);
    }else{
        level--;
        tree_step( level, pos_, dir,  branches, leafs );
        tree_step( level, pos_, dir,  branches, leafs );
    }
}

int setup(){

    shBranches=new Shader();
    shBranches->init( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/color3D.glslf"   );
    shBranches->getDefaultUniformLocation();

    shLeafs=new Shader();
    shLeafs->init( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/pointSprite.glslf"   );
    shLeafs->getDefaultUniformLocation();

    std::vector<Vec3f> branches;
    std::vector<Vec3f> leafs;
    tree_step( 5, (Vec3f){0.0f,0.0f,0.0f}, (Vec3f){0.0f,1.0f,0.0f}, branches, leafs );

    //for(int i=0; i<branches.size(); i++){ printf("%i (%g,%g,%g)\n", i, branches[i].x, branches[i].y, branches[i].z ); }

    mshBranches = new GLMesh();
    mshBranches->init( branches.size(), 0, NULL, &branches[0],  NULL, NULL, NULL );

    mshLeafs = new GLMesh();
    mshLeafs->init( leafs.size(), 0, NULL, &leafs[0],  NULL, NULL, NULL );

    //ticks_per_second = calibrate_timer(100);
    //lastTime = glfwGetTime();
    lastTime = 0.0;
    qCamera.setOne();
    return 0;
};

void draw( ){
    glClearColor(0.8, 0.8, 0.8, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );
	// Simulate all particles

    Mat4f camMat,mRot,mPersp;
    qCamera.toMatrix(mouseMat);
    mRot.setOne(); mRot.setRot(mouseMat);
    float fov = 3.0;
    mPersp.setPerspective( fov, fov*ASPECT_RATIO, 1.0, 1000.0 );
    camMat.set_mmul_TN( mRot, mPersp );
    //Mat3f objRot; objRot.setOne();

    Mat3f modelMat; modelMat.setOne();
    Vec3f modelPos; modelPos.set(0.0f,0.0f,10.0f);

    Shader * sh;

    sh = shBranches;
    sh->use();
    sh->set_camPos  ( (GLfloat*)&camPos );
    sh->set_camMat  ( (GLfloat*)&camMat );
    sh->set_modelPos( (GLfloat*)&modelPos );
    sh->set_modelMat( (GLfloat*)&modelMat );

    glLineWidth(3.0); mshBranches->draw(GL_LINES);
    //mshBranch->drawPoints(10.0);

    sh = shLeafs;
    sh->use();
    sh->set_camPos  ( (GLfloat*)&camPos );
    sh->set_camMat  ( (GLfloat*)&camMat );
    sh->set_modelPos( (GLfloat*)&modelPos );
    sh->set_modelMat( (GLfloat*)&modelMat );
    //glEnable( GL_PROGRAM_POINT_SIZE ); // somehow does not work ... perhaps look inside the shader
    glPointSize(30.0);
    mshLeafs->draw( GL_POINTS );

    //glUniformMatrix3fv( sh->getUloc( "camRot" ), 1, GL_FALSE, (GLfloat*)&mouseMat );
    //glUniform4fv( sh->getUloc( "keyColor" ), 1, (const GLfloat[]){1.0f,0.0f,1.0f,10.0f} );

    //uploadArrayBuffer( instances.pose_Up, instances.nInstances*3*sizeof(GLfloat), instance_Up );
    //bilboards.draw( GL_TRIANGLES );

    //glPointSize(10.0); bilboards.draw( GL_POINTS );

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

    if( keys[ SDL_SCANCODE_W  ] ){ camPos.add_mul( mouseMat.c, +step ); }
	if( keys[ SDL_SCANCODE_S  ] ){ camPos.add_mul( mouseMat.c, -step ); }
	if( keys[ SDL_SCANCODE_A  ] ){ camPos.add_mul( mouseMat.a, -step ); }
	if( keys[ SDL_SCANCODE_D  ] ){ camPos.add_mul( mouseMat.a, +step ); }

	if( keys[SDL_SCANCODE_W]||keys[SDL_SCANCODE_S]||keys[SDL_SCANCODE_A]||keys[ SDL_SCANCODE_D  ]  )printf( "camPos (%g,%g,%g)\n", camPos.x, camPos.y, camPos.z );

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
