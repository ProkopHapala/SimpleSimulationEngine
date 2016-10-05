#include <stdlib.h>
#include <stdio.h>

//#define GL3_PROTOTYPES 1
#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>
#include <SDL2/SDL.h>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"

//#include "IO_utils.h"
#include "Shader.h"
#include "GLObject.h"

Shader   * shader1;
GLObject * object1;
GLObject * object2;

GLfloat vertexes[4][2] = {
	{  -0.8f,  -0.8f  },
	{  -0.8f,   0.8f  },
	{   0.8f,  -0.8f  },
	{   0.8f,   0.8f  } };

SDL_Window * window     = NULL;
SDL_GLContext   context = NULL;

GLuint vao;     // vertex array object

const int WIDTH  = 800;
const int HEIGHT = 800;

GLfloat resolution[2];
//GLfloat light_dir [3] = { -0.35707142142f, -0.35707142142f, 0.7f };
GLfloat light_dir [3] = { -0.35707142142f, 0.7f, 0.35707142142f };

const int natoms_max = 100;
int natoms = natoms_max;

GLfloat atoms[natoms_max*4];
GLfloat coefs[natoms_max*4];

bool redraw = true;

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
	shader1->init( "shaders/plain_vert.c", "shaders/orbitalRayMarch.c" );

	glUseProgram(shader1->shaderprogram);

	natoms=10;
	for( int i=0; i<natoms; i++ ){
		int ii = i<<2;
		atoms[ii   ] = randf( -1.0,1.0  );
		atoms[ii+1 ] = randf( -1.0,1.0  );
		atoms[ii+2 ] = randf(  1.0,5.0 );
		atoms[ii+3 ] = 2.0;
		printf( "atom %i %i (%3.3f,%3.3f,%3.3f) %3.3f \n", i, ii, atoms[ii], atoms[ii]+1, atoms[ii]+2, atoms[ii+3]  ); 
	}

}

void draw(){
	if( redraw  ){
		glClearColor(0.0, 0.0, 1.0, 1.0);
		glClear(GL_COLOR_BUFFER_BIT);

		resolution[0] = (float)WIDTH;
    	resolution[1] = (float)HEIGHT;

		GLuint uloc;
    	uloc = glGetUniformLocation( shader1->shaderprogram, "resolution" );	glUniform2fv(uloc, 1, resolution  );
		uloc = glGetUniformLocation( shader1->shaderprogram, "light_dir" );	    glUniform3fv(uloc, 1, light_dir   );

		uloc = glGetUniformLocation( shader1->shaderprogram, "natoms"  );		glUniform1iv(uloc, 1, &natoms    );
		uloc = glGetUniformLocation( shader1->shaderprogram, "atoms"    );		glUniform4fv(uloc, natoms, atoms );
		uloc = glGetUniformLocation( shader1->shaderprogram, "coefs"    );		glUniform4fv(uloc, natoms, coefs );

        glEnableVertexAttribArray(0);
		object1->draw();

		SDL_GL_SwapWindow(window);

		redraw = false;
	}
}

// FUNCTION ======	inputHanding
void inputHanding(){
	SDL_Event event;
	while(SDL_PollEvent(&event)){
		if( event.type == SDL_KEYDOWN ){
			if(event.key.keysym.sym == SDLK_ESCAPE ) { quit(); }
			if(event.key.keysym.sym == SDLK_SPACE  ) { STOP=!STOP; }
			redraw = true;
		}
		if( event.type == SDL_QUIT){ quit();  };
	}
}

void loop( int nframes ){
    for ( int iframe=1; iframe<nframes; iframe++)    {
 		if( !STOP ) draw();
		inputHanding();
		frameCount++;
        SDL_Delay(100);
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
