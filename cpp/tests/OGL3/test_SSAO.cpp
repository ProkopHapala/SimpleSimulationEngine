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

Shader   * shader_pre, * shader_post;
GLObject * object1;
GLObject * object2;

int mouseX, mouseY;

SDL_Window * window     = NULL;
SDL_GLContext   context = NULL;

GLuint vao;     // vertex array object

bool depth_instead_color = true;

GLfloat vertexes[4][2] = {
	{  -0.9f,  -0.9f  },
	{  -0.9f,   0.9f  },
	{   0.9f,  -0.9f  },
	{   0.9f,   0.9f  } };

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

GLuint textureID;
GLuint FramebufferName = 0;
GLuint depthrenderbuffer;

// ============= FUNCTIONS

int frameCount = 0;
bool STOP = false;

void quit();
void die ( char const *msg );
void inputHanding ();
void init();
void draw();
void loop( int niters );


bool checkFramebufferStatus(){
    // check FBO status
    GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
    switch(status)
    {
    case GL_FRAMEBUFFER_COMPLETE:
        printf( "Framebuffer complete.\n" );
        return true;

    case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT:
        printf( "[ERROR] Framebuffer incomplete: Attachment is NOT complete.\n" );
        return false;

    case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT:
        printf( "[ERROR] Framebuffer incomplete: No image is attached to FBO.\n" );
        return false;
/*
    case GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS:
        std::cout << "[ERROR] Framebuffer incomplete: Attached images have different dimensions." << std::endl;
        return false;

    case GL_FRAMEBUFFER_INCOMPLETE_FORMATS:
        std::cout << "[ERROR] Framebuffer incomplete: Color attached images have different internal formats." << std::endl;
        return false;
*/
    case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER:
        printf( "[ERROR] Framebuffer incomplete: Draw buffer.\n" );
        return false;

    case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER:
        printf( "[ERROR] Framebuffer incomplete: Read buffer.\n" );
        return false;

    case GL_FRAMEBUFFER_UNSUPPORTED:
        printf( "[ERROR] Framebuffer incomplete: Unsupported by FBO implementation.\n" );
        return false;

    default:
        printf( "[ERROR] Framebuffer incomplete: Unknown error.\n" );
        return false;
    }
}

void setup(){

    // ------------- object

	object1 = new GLObject( );
	object1->nVert   = 4;
	object1->vertDim = 2;
	object1->vertexes = &vertexes[0][0];
	object1->init();

    // ------------- shader for rendering geometry

    shader_pre=new Shader();
	shader_pre->init( "shaders/afine2D_vert.c", "shaders/sphere_frag.c" );

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

	glUseProgram(shader_pre->shaderprogram);

    GLuint uloc;
    uloc = glGetUniformLocation( shader_pre->shaderprogram, "resolution" );	glUniform2fv(uloc, 1, resolution  );
    uloc = glGetUniformLocation( shader_pre->shaderprogram, "sphere"     );	glUniform4fv(uloc, 1, sphere      );
    uloc = glGetUniformLocation( shader_pre->shaderprogram, "light_dir"  );	glUniform3fv(uloc, 1, light_dir   );

    uloc = glGetUniformLocation( shader_pre->shaderprogram, "afineMat"   );	glUniformMatrix2fv(uloc, 1, GL_FALSE, afineMat  );
    uloc = glGetUniformLocation( shader_pre->shaderprogram, "origin"     );	glUniform2fv      (uloc, 1, origin   );

    // ------------- shader for blitting from texture

    shader_post=new Shader();
	shader_post->init( "shaders/plain_vert.c", "shaders/texture_frag.c" );

    resolution[0] = (float)WIDTH;
    resolution[1] = (float)HEIGHT;

	glUseProgram(shader_post->shaderprogram);

    uloc = glGetUniformLocation( shader_post->shaderprogram, "resolution" );	glUniform2fv(uloc, 1, resolution  );

    // ------------- texture

    glGenTextures(1, &textureID);    // Create one OpenGL texture
    glBindTexture(GL_TEXTURE_2D, textureID); // "Bind" the newly created texture : all future texture functions will modify this texture
    //glTexImage2D(GL_TEXTURE_2D, 0,GL_RGBA, imgW, imgH, 0, GL_RGBA, GL_UNSIGNED_BYTE, imgData);   // Give the image to OpenGL
    //if( depth_instead_color ){ glTexImage2D(GL_TEXTURE_2D, 0,GL_DEPTH_COMPONENT24, WIDTH, HEIGHT, 0, GL_DEPTH_COMPONENT24, GL_FLOAT,         0); }
    //else                     { glTexImage2D(GL_TEXTURE_2D, 0,GL_RGB,               WIDTH, HEIGHT, 0, GL_RGB,               GL_UNSIGNED_BYTE, 0); }
    if( depth_instead_color ){ glTexImage2D(GL_TEXTURE_2D, 0,GL_DEPTH_COMPONENT, WIDTH, HEIGHT, 0, GL_DEPTH_COMPONENT, GL_FLOAT,         0); }
    else                     { glTexImage2D(GL_TEXTURE_2D, 0,GL_RGB,             WIDTH, HEIGHT, 0, GL_RGB,               GL_UNSIGNED_BYTE, 0); }
    //glTexImage2D(GL_TEXTURE_2D, 0, GL_ARGB, imgW, imgH, 0, GL_ARGB, GL_UNSIGNED_BYTE, imgData);   // Give the image to OpenGL
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    //glActiveTexture(GL_TEXTURE0 );
    //glBindTexture(GL_TEXTURE_2D, textureID );
    //glBindSampler(0, uloc);

    // ------------- frameBuffer

    glGenFramebuffers(1, &FramebufferName);
    glBindFramebuffer(GL_FRAMEBUFFER, FramebufferName);
    // The depth buffer

    glGenRenderbuffers(1, &depthrenderbuffer);
    glBindRenderbuffer(GL_RENDERBUFFER, depthrenderbuffer);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, WIDTH, HEIGHT );
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthrenderbuffer);

    // Set "renderedTexture" as our colour attachement #0

    GLenum DrawBuffers[1];
    if( depth_instead_color ){
        glFramebufferTexture(GL_FRAMEBUFFER,  GL_DEPTH_ATTACHMENT, textureID, 0);
        DrawBuffers[0] = {GL_DEPTH_ATTACHMENT};
    }else{
        glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, textureID, 0);
        DrawBuffers[0] = {GL_COLOR_ATTACHMENT0};
    }
    glDrawBuffers(1, DrawBuffers); // "1" is the size of DrawBuffers

    if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE){
        printf(" problem in FBO ! \n ");
        checkFramebufferStatus();
    }


}

void draw(){
    glClearColor(0.0, 0.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );

    GLuint uloc;

    // ------- render to texture

    //glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glBindFramebuffer(GL_FRAMEBUFFER, FramebufferName);
    glClearColor(0.0, 0.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );
    glEnable( GL_DEPTH_TEST );
    glDepthFunc( GL_LESS );

    glUseProgram(shader_pre->shaderprogram);
    origin[0] = 0;    origin[1] = 0;
    uloc = glGetUniformLocation( shader_pre->shaderprogram, "origin"     );	glUniform2fv(uloc, 0, origin   );
    for( int i; i<10; i++ ){
        //sphere[0] = randf()-0.5f;    sphere[1] = randf()-0.5f;    origin[0] = randf()-0.5f;    origin[1] = sphere[0];
        float d = 0.2;
        sphere[0] = d*i-1.0;    sphere[1] = d*i-1.0;     sphere[2] = d*i-1.0;
        //sphere[0] = 0.0;    sphere[1] = 0.0;    origin[0] = sphere[0];    origin[1] = sphere[0];
        uloc = glGetUniformLocation( shader_pre->shaderprogram, "sphere"     );	glUniform4fv(uloc, 1, sphere   );
        glEnableVertexAttribArray(0); object1->draw();
    }

    // -------- from texture to Screen

    glUseProgram(shader_post->shaderprogram);
    uloc = glGetUniformLocation( shader_post->shaderprogram, "texture1");
    glUniform1i(uloc, 0);

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    //glViewport(0,0,WIDTH,HEIGHT);
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
