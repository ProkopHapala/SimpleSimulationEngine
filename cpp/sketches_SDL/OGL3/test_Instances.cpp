
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


#include "GLInstances.h"
#include "Shader.h"

// =============== Global variables

//const int MaxParticles = 100000;
const int MaxParticles = 16;

GLuint VertexArrayID;
GLfloat* instance_pos = new GLfloat[MaxParticles * 4];
GLubyte* instance_color         = new GLubyte[MaxParticles * 4];

//GLuint programID;
Shader* shaderParticle;
//GLuint CameraRight_worldspace_ID;
//GLuint CameraUp_worldspace_ID;
//GLuint ViewProjMatrixID;

GLInstances instances;
//GLuint billboard_vertex_buffer;
//GLuint particles_color_buffer;
//GLuint particles_position_buffer;

int nVerts = 0;
int ParticlesCount = 0;
int frameCount = 0;
double lastTime = 0.0;

GLuint vao;     // vertex array object
int delay = 1;
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

// =============== Functions

int setup(){

	shaderParticle=new Shader();
	shaderParticle->init( "common_resources/shaders/Particle.glslv",   "common_resources/shaders/Particle.glslf"   );

    shaderParticle->use();
    shaderParticle->getDefaultUniformLocation();

	instance_pos = new GLfloat[MaxParticles * 4];
	instance_color         = new GLubyte[MaxParticles * 4];

	nVerts = 50;
    float dx = 1.0f/nVerts;
	Vec3f g_vertex_buffer_data[nVerts*3];
	for( int i=0; i<nVerts; i++){
        float x = i*dx;
        g_vertex_buffer_data[i].x = x*sin(x*10.0);
        g_vertex_buffer_data[i].y = x*cos(x*10.0);
        g_vertex_buffer_data[i].z = x;
	}

    float time = frameCount * 0.005;
    ParticlesCount = MaxParticles;
	for( int i=0; i<ParticlesCount; i++ ){
        instance_pos[4*i+0] = sin(time+i);
        instance_pos[4*i+1] = cos(time+i);
        instance_pos[4*i+2] = -5.0;
        instance_pos[4*i+3] =  0.5f;

        // actually direction
        instance_color[4*i+0] = rand() % 0xFF;
        instance_color[4*i+1] = rand() % 0xFF;
        instance_color[4*i+2] = rand() % 0xFF;
        instance_color[4*i+3] = 0xFF;
	}

	instances.init( MaxParticles, nVerts, g_vertex_buffer_data, instance_pos, instance_color );

    //lastTime = glfwGetTime();
    lastTime = 0.0;
    qCamera.setOne();
    return 0;
};

void physics(){
    float time = frameCount * 0.005;
    ParticlesCount = MaxParticles;
	for( int i=0; i<ParticlesCount-2; i++ ){
        instance_pos[4*i+0] = sin(time+i);
        instance_pos[4*i+1] = cos(time+i);
        instance_pos[4*i+2] =  5.0;
        instance_pos[4*i+3] =  0.5f;
	}
}

void draw( ){

    glClearColor(0.8, 0.8, 0.8, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );
	// Simulate all particles

	glDisable(GL_ALPHA);

	//physics();

	//glEnable(GL_BLEND);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glDisable(GL_BLEND);

    Mat4f camMat,mRot,mPersp;
    qCamera.toMatrix(mouseMat);
    mRot.setOne(); mRot.setRot(mouseMat);
    float fov = 3.0;
    mPersp.setPerspective( fov, fov*ASPECT_RATIO, 1.0, 1000.0 );
    camMat.set_mmul_TN( mRot, mPersp );
    //Mat3f objRot; objRot.setOne();

    //printf("======\n");
    //printf( "qCamera (%g,%g,%g,%g)\n", qCamera.x, qCamera.y, qCamera.z, qCamera.w );
    //printf( "camPos (%g,%g,%g)\n", camPos.x, camPos.y, camPos.z );
    //mouseMat.print();
    //camMat.print();

    shaderParticle->use();
    shaderParticle->set_camPos( (GLfloat*)&camPos );
    shaderParticle->set_camMat( (GLfloat*)&camMat );

    // Update the buffers that OpenGL uses for rendering.
	// There are much more sophisticated means to stream data from the CPU to the GPU,
	// but this is outside the scope of this tutorial.
	// http://www.opengl.org/wiki/Buffer_Object_Streaming
    //instances.upload_pos   ( ParticlesCount, instance_pos   );
    //instances.upload_colors( ParticlesCount, instance_color );
    instances.draw();

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
    float step          = 0.005f;
    float keyRotSpeed   = 0.002f;

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
