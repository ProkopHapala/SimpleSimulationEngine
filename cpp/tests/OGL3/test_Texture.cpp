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


GLfloat vertexes[4][2] = {
	{  -0.9f,  -0.9f  },
	{  -0.9f,   0.9f  },
	{   0.9f,  -0.9f  },
	{   0.9f,   0.9f  } };

GLfloat  origin[2] = {  0.0, 0.0 };

const int imgW = 256;
const int imgH = 256;

unsigned int imgData [imgW*imgH];

int WIDTH  = 800;
int HEIGHT = 800;

GLfloat resolution[2];
GLfloat sphere[4];
GLfloat light_dir[4];

GLuint textureID;

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

    // ------------- object

	object1 = new GLObject( );
	object1->nVert   = 4;
	object1->vertDim = 2;
	object1->vertexes = &vertexes[0][0];
	object1->init();

    // ------------- shader

	shader1=new Shader();
	//shader1->init( "shaders/plain_vert.c", "shaders/sphere_frag.c" );
	//shader1->init( "shaders/afine2D_vert.c", "shaders/sphere_frag.c" );
	shader1->init( "shaders/plain_vert.c", "shaders/texture_frag.c" );

    resolution[0] = (float)WIDTH;
    resolution[1] = (float)HEIGHT;

	glUseProgram(shader1->shaderprogram);

    GLuint uloc;
    uloc = glGetUniformLocation( shader1->shaderprogram, "resolution" );	glUniform2fv(uloc, 1, resolution  );

    // ------------- texture

    for( int iy=0; iy<imgH; iy++ ){
        for( int ix=0; ix<imgW; ix++ ){
            float r = sin(ix*0.1) + 1.0f;
            float g = sin(iy*0.1) + 1.0f;
            float b = sin((ix+iy)*0.1) + 1.0f;
            imgData[ iy*imgW + ix ] =  ((int)(127*r) <<16) | ((int)(127*g)<<8) | ((int)(127*b));
            //imgData[ iy*imgW + ix ]   = (( 0xFF & (ix^iy) )<<16 ) | ( 0xFF&(ix^(-iy) )<<8 ); // | ( 0xFF&(ix*iy) );
            //imgData[ iy*imgW + ix ]   = (( 0xFF & (ix^iy) )<<8 ) | ( 0xFF&(ix^(-iy)) ); // | ( 0xFF&(ix*iy) );
            //imgData[ iy*imgW + ix ]   = ix^iy;
        }
    }
    glGenTextures(1, &textureID);    // Create one OpenGL texture
    glBindTexture(GL_TEXTURE_2D, textureID); // "Bind" the newly created texture : all future texture functions will modify this texture
    glTexImage2D(GL_TEXTURE_2D, 0,GL_RGBA, imgW, imgH, 0, GL_RGBA, GL_UNSIGNED_BYTE, imgData);   // Give the image to OpenGL
    //glTexImage2D(GL_TEXTURE_2D, 0, GL_ARGB, imgW, imgH, 0, GL_ARGB, GL_UNSIGNED_BYTE, imgData);   // Give the image to OpenGL
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    //glActiveTexture(GL_TEXTURE0 );
    //glBindTexture(GL_TEXTURE_2D, textureID );
    //glBindSampler(0, uloc);

}

void draw(){
    glClearColor(0.0, 0.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );

    GLuint uloc;

    uloc = glGetUniformLocation( shader1->shaderprogram, "texture1");     glUniform1i(uloc, 0);
    //glActiveTexture(GL_TEXTURE0 );
    //glBindTexture(GL_TEXTURE_2D, textureID );
    //glBindSampler(0, uloc);


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
