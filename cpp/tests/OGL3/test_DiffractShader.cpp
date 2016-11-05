#include <stdlib.h>
#include <stdio.h>

//#define GL3_PROTOTYPES 1
#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>
#include <SDL2/SDL.h>


const float  INV_RAND_MAX = 1.0f/RAND_MAX;
inline float randf(){ return INV_RAND_MAX*rand(); }

//#include "IO_utils.h"
#include "Shader.h"
#include "GLObject.h"


Shader   * shader1;
GLObject * object1;
GLObject * object2;

SDL_Window * window     = NULL;
SDL_GLContext   context = NULL;

GLuint vao;     // vertex array object

GLfloat vertexes[4][2] = {
	{  -1.0f,  -1.0f  },
	{  -1.0f,   1.0f  },
	{   1.0f,  -1.0f  },
	{   1.0f,   1.0f  } };


const int window_size = 800;

float aperture    = 1.5;
float wave_length = 5e-7;
float distance    = 100e+6;

int nx = 3;
int ny = 3;
float phase = 0;

float Iscale    = 0.1;
float scale     = 0.2;

const int nsource_max = 100;
int nsource = nsource_max;

float spacing = 5.0;
float random_pos_scale = 0.1;

bool redraw = true;


GLfloat sources[nsource_max*4];

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
    object1->draw_mode = GL_TRIANGLE_STRIP;
	object1->nVert   = 4;
    object1->buffs[0].setup(0,2,GL_FALSE,&vertexes[0][0],'v'); // vertexes
	object1->init();

	shader1=new Shader();
	shader1->init( "shaders/plain_vert.c", "shaders/diffract_frag.c" );

	glUseProgram(shader1->shaderprogram);
}

void draw(){
	if( redraw  ){
		glClearColor(0.0, 0.0, 1.0, 1.0);
		glClear(GL_COLOR_BUFFER_BIT);

		// https://en.wikipedia.org/wiki/Diffraction
		// sin Theta = 1.22 * wavelength / aperture
		//  spot_size / distance = wavelength / aperture
		//  spot_size * aperture = wavelength * distance
		//  spot_size [m] * aperture [m] = wavelength [ micron ] * distance [ 1000 km ]

		//float aperture2  = aperture*aperture;
		float freq       = 6.28318530718/wave_length;
		float spot_size  = wave_length * distance / aperture;
		//float rho_scale2 = 1/(2*spot_size*spot_size) ;

		//printf( " aperture2: %f freq: %f spot_size: %f rho_scale2: %f \n", aperture2, freq, spot_size, rho_scale2 );
		printf( "distance: %f wave_length: %f aperture: %f freq: %f spot_size: %f \n", distance, wave_length, aperture, freq, spot_size );

		float dist2  = distance * distance;
		//float freq   = 6.28318530718/wave_length;
		nsource=nx*ny;
		int i = 0;
		for( int iy=0; iy<ny; iy++ ){
			for( int ix=0; ix<nx; ix++ ){
				int ii = i<<2;
				float x      = ( ix-nx/2 + (randf() - 0.5)*random_pos_scale )*spacing;
				float y      = ( iy-ny/2 + (randf() - 0.5)*random_pos_scale )*spacing;
				float r2     = x*x + y*y;
				float a      = r2/dist2;
				float dr     = distance*a*( 0.5 + a*( -0.125 + a*( 0.0625 + a* -0.0390625 ) ) );
				float phi    = freq * dr;
				float wave_x = cos( phi );
				float wave_y = sin( phi );
				sources[ii   ] = x;
				sources[ii+1 ] = y;
				sources[ii+2 ] = wave_x;
				sources[ii+3 ] = -wave_y;
				i++;
			}
		}

		GLuint uloc;
		uloc = glGetUniformLocation( shader1->shaderprogram, "window_size"     );	glUniform1iv(uloc, 1, &window_size  );

		uloc = glGetUniformLocation( shader1->shaderprogram, "nsource"     );		glUniform1iv(uloc, 1, &nsource     );
		uloc = glGetUniformLocation( shader1->shaderprogram, "wave_length" );	    glUniform1fv(uloc, 1, &wave_length );
		uloc = glGetUniformLocation( shader1->shaderprogram, "distance"    );	    glUniform1fv(uloc, 1, &distance    );
		uloc = glGetUniformLocation( shader1->shaderprogram, "aperture"    );	    glUniform1fv(uloc, 1, &aperture    );

		uloc = glGetUniformLocation( shader1->shaderprogram, "Iscale"     );		glUniform1fv(uloc, 1, &Iscale   );
		uloc = glGetUniformLocation( shader1->shaderprogram, "scale"      );		glUniform1fv(uloc, 1, &scale    );
		uloc = glGetUniformLocation( shader1->shaderprogram, "phase"      );		glUniform1fv(uloc, 1, &phase   );
		uloc = glGetUniformLocation( shader1->shaderprogram, "sources"    );		glUniform4fv(uloc, nsource, sources );

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

			if(event.key.keysym.sym == SDLK_LEFTBRACKET  ) { nx=(nx+1)%10; ny=(ny+1)%10;   redraw = true;  }
			if(event.key.keysym.sym == SDLK_RIGHTBRACKET ) { nx=(nx-1)%10; ny=(ny-1)%10;   redraw = true; }

			if(event.key.keysym.sym == SDLK_PAGEDOWN     ) { distance/=1.1; redraw = true; }
			if(event.key.keysym.sym == SDLK_PAGEUP       ) { distance*=1.1; redraw = true; }

			if(event.key.keysym.sym == SDLK_HOME         ) { wave_length/=1.1; redraw = true; }
			if(event.key.keysym.sym == SDLK_END          ) { wave_length*=1.1; redraw = true; }

			if(event.key.keysym.sym == SDLK_INSERT       ) { aperture/=1.1; redraw = true; }
			if(event.key.keysym.sym == SDLK_DELETE       ) { aperture*=1.1; redraw = true; }

			if(event.key.keysym.sym == SDLK_KP_PLUS      ) { spacing/=1.1; redraw = true; }
			if(event.key.keysym.sym == SDLK_KP_MINUS     ) { spacing*=1.1; redraw = true; }

			if(event.key.keysym.sym == SDLK_KP_DIVIDE    ) { random_pos_scale/=1.1; redraw = true; }
			if(event.key.keysym.sym == SDLK_KP_MULTIPLY  ) { random_pos_scale*=1.1; redraw = true; }

			if(event.key.keysym.sym == SDLK_UP           ) { Iscale/=1.1; redraw = true; }
			if(event.key.keysym.sym == SDLK_DOWN         ) { Iscale*=1.1; redraw = true; }

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
    window = SDL_CreateWindow("Tutorial2", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, window_size, window_size, SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN);
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
