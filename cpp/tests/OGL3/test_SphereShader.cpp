#include <stdlib.h>
#include <stdio.h>

//#define GL3_PROTOTYPES 1
#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>
#include <SDL2/SDL.h>

#include "Vec2.h"
#include "Vec3.h"

#include "GLObject.h"
#include "Shader.h"

Shader   * shader1;
GLObject * object1;
GLObject * object2;

int mouseX, mouseY;

SDL_Window * window     = NULL;
SDL_GLContext   context = NULL;

GLuint vao;     // vertex array object

/*
GLfloat vertexes[4][2] = {
	{  -1.0f,  -1.0f  },
	{  -1.0f,   1.0f  },
	{   1.0f,  -1.0f  },
	{   1.0f,   1.0f  } };
*/


GLfloat vertexes[4][2] = {
	{  -0.2f,  -0.2f  },
	{  -0.2f,   0.2f  },
	{   0.2f,  -0.2f  },
	{   0.2f,   0.2f  } };

GLfloat  afineMat[4] = {
  1.0, 0.0,
  0.0, 1.0
};

GLfloat  origin[2] = {  0.0, 0.0 };


int WIDTH  = 800;
int HEIGHT = 800;

GLfloat resolution[2];
GLfloat sphere[4];
GLfloat light_dir[4];

// ============= FUNCTIONS

int frameCount = 0;
bool STOP = false;

void quit();
void die ( char const *msg );
void inputHanding ();
void init();
void draw();
void loop( int niters );

void setup(){

	object1 = new GLObject( );
	object1->nVert   = 4;
	object1->vertDim = 2;
	object1->vertexes = &vertexes[0][0];
	object1->init();

	shader1=new Shader();
	//shader1->init( "shaders/plain_vert.c", "shaders/sphere_frag.c" );
	shader1->init( "shaders/afine2D_vert.c", "shaders/sphere_frag.c" );

    sphere[0] = 0.0;
    sphere[1] = 0.0;
    sphere[2] = 0.0;
    sphere[3] = 0.5;

    Vec3f light_dir_; light_dir_.set( 1, 1, 2 ); light_dir_.normalize();
    light_dir[0]=light_dir_.x;
    light_dir[1]=light_dir_.x;
    light_dir[2]=light_dir_.x;

    resolution[0] = (float)WIDTH;
    resolution[1] = (float)HEIGHT;

	glUseProgram(shader1->shaderprogram);

    GLuint uloc;
    uloc = glGetUniformLocation( shader1->shaderprogram, "resolution" );	glUniform2fv(uloc, 1, resolution  );
    uloc = glGetUniformLocation( shader1->shaderprogram, "sphere"     );	glUniform4fv(uloc, 1, sphere      );
    uloc = glGetUniformLocation( shader1->shaderprogram, "light_dir"  );	glUniform3fv(uloc, 1, light_dir   );

    uloc = glGetUniformLocation( shader1->shaderprogram, "afineMat"   );	glUniformMatrix2fv(uloc, 1, GL_FALSE, afineMat  );
    uloc = glGetUniformLocation( shader1->shaderprogram, "origin"     );	glUniform2fv      (uloc, 1, origin   );

}

void draw(){
    glClearColor(0.0, 0.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );

    glEnable( GL_DEPTH_TEST );
    glDepthFunc( GL_LESS );

    GLuint uloc;

    // ---- constant sphere in center
    sphere[0] = 0.0;
    sphere[1] = 0.0;
    origin[0] = sphere[0];
    origin[1] = sphere[0];
    uloc = glGetUniformLocation( shader1->shaderprogram, "origin"     );	glUniform2fv      (uloc, 1, origin   );

    uloc = glGetUniformLocation( shader1->shaderprogram, "sphere"     );	glUniform4fv(uloc, 1, sphere      );
    glEnableVertexAttribArray(0); object1->draw();

    // ---- shifted sphere in mouse position
    sphere[0] = -2*(mouseX*2-WIDTH)/resolution[0];
    sphere[1] =  2*(mouseY*2-HEIGHT)/resolution[1];
    origin[0] = -sphere[0]*0.4;
    origin[1] = -sphere[1]*0.4;
    uloc = glGetUniformLocation( shader1->shaderprogram, "origin"     );	glUniform2fv      (uloc, 1, origin   );

    //printf( "  %i %i %f %f \n", mouseX, mouseY, sphere[0], sphere[1] );
    uloc = glGetUniformLocation( shader1->shaderprogram, "sphere"     );	glUniform4fv(uloc, 1, sphere      );
    glEnableVertexAttribArray(0); object1->draw();

    SDL_GL_SwapWindow(window);

}

// FUNCTION ======	inputHanding
void inputHanding(){
	SDL_Event event;
	while(SDL_PollEvent(&event)){
		if( event.type == SDL_KEYDOWN ){
			if(event.key.keysym.sym == SDLK_ESCAPE ) { quit(); }
			if(event.key.keysym.sym == SDLK_SPACE  ) { STOP=!STOP; }

		}
		if( event.type == SDL_QUIT){ quit();  };
	}
	SDL_GetMouseState( &mouseX, &mouseY );
}

void loop( int nframes ){
    for ( int iframe=1; iframe<nframes; iframe++)    {
 		if( !STOP ) draw();
		inputHanding();
		frameCount++;
        SDL_Delay(10);
    }
}

int main(int argc, char *argv[]){
    init();
	setup();
	loop( 100000 );
    quit();
    return 0;
}

void init(){
    if (SDL_Init(SDL_INIT_VIDEO) < 0) die( "Unable to initialize SDL" );
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3); // Opengl 3.2
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 2); // Opengl 3.2
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    window = SDL_CreateWindow("Tutorial2", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WIDTH, HEIGHT, SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN);
    if ( !window ) die("Unable to create window");
    context = SDL_GL_CreateContext( window );
    SDL_GL_SetSwapInterval(1);

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
