
// read this tutorial
// http://www.opengl-tutorial.org/beginners-tutorials/tutorial-4-a-colored-cube/

// http://gamedev.stackexchange.com/questions/93055/getting-the-real-fragment-depth-in-glsl

// TODO see later vbo indexing to avoid duplicating triangles
// http://www.opengl-tutorial.org/intermediate-tutorials/tutorial-9-vbo-indexing/

#include <stdlib.h>
#include <stdio.h>

//#define GL3_PROTOTYPES 1
#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>
#include <SDL2/SDL.h>

#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

#include "Mesh.h"

#include "GLObject.h"
#include "Shader.h"


Shader   * shader1;
GLObject * object1;


Mesh mesh;



GLuint vao;     // vertex array object



/*
const int nVert = 3;
GLfloat vertexes[nVert*3] = {
  1.0f,  0.0f,  0.0f,
  0.0f,  1.0f,  0.0f,
  0.0f,  0.0f,  1.0f,
};
GLfloat colors[nVert*3] = {
  1.0f,  0.0f,  0.0f,
  0.0f,  1.0f,  0.0f,
  0.0f,  0.0f,  1.0f,
};

GLfloat UVs[nVert*2] = {
  0.0f,  0.0f,
  1.0f,  0.0f,
  0.0f,  1.0f
};
*/


//http://www.opengl-tutorial.org/beginners-tutorials/tutorial-4-a-colored-cube/
// Our vertices. Three consecutive floats give a 3D vertex; Three consecutive vertices give a triangle.
// A cube has 6 faces with 2 triangles each, so this makes 6*2=12 triangles, and 12*3 vertices

const int nVert = 36;
GLfloat vertexes[nVert*3] = {
    -1.0f,-1.0f,-1.0f, // triangle 1 : begin
    -1.0f,-1.0f, 1.0f,
    -1.0f, 1.0f, 1.0f, // triangle 1 : end
    1.0f, 1.0f,-1.0f, // triangle 2 : begin
    -1.0f,-1.0f,-1.0f,
    -1.0f, 1.0f,-1.0f, // triangle 2 : end
    1.0f,-1.0f, 1.0f,
    -1.0f,-1.0f,-1.0f,
    1.0f,-1.0f,-1.0f,
    1.0f, 1.0f,-1.0f,
    1.0f,-1.0f,-1.0f,
    -1.0f,-1.0f,-1.0f,
    -1.0f,-1.0f,-1.0f,
    -1.0f, 1.0f, 1.0f,
    -1.0f, 1.0f,-1.0f,
    1.0f,-1.0f, 1.0f,
    -1.0f,-1.0f, 1.0f,
    -1.0f,-1.0f,-1.0f,
    -1.0f, 1.0f, 1.0f,
    -1.0f,-1.0f, 1.0f,
    1.0f,-1.0f, 1.0f,
    1.0f, 1.0f, 1.0f,
    1.0f,-1.0f,-1.0f,
    1.0f, 1.0f,-1.0f,
    1.0f,-1.0f,-1.0f,
    1.0f, 1.0f, 1.0f,
    1.0f,-1.0f, 1.0f,
    1.0f, 1.0f, 1.0f,
    1.0f, 1.0f,-1.0f,
    -1.0f, 1.0f,-1.0f,
    1.0f, 1.0f, 1.0f,
    -1.0f, 1.0f,-1.0f,
    -1.0f, 1.0f, 1.0f,
    1.0f, 1.0f, 1.0f,
    -1.0f, 1.0f, 1.0f,
    1.0f,-1.0f, 1.0f
};
// One color for each vertex. They were generated randomly.
GLfloat colors[nVert*3] = {
    0.583f,  0.771f,  0.014f,
    0.609f,  0.115f,  0.436f,
    0.327f,  0.483f,  0.844f,
    0.822f,  0.569f,  0.201f,
    0.435f,  0.602f,  0.223f,
    0.310f,  0.747f,  0.185f,
    0.597f,  0.770f,  0.761f,
    0.559f,  0.436f,  0.730f,
    0.359f,  0.583f,  0.152f,
    0.483f,  0.596f,  0.789f,
    0.559f,  0.861f,  0.639f,
    0.195f,  0.548f,  0.859f,
    0.014f,  0.184f,  0.576f,
    0.771f,  0.328f,  0.970f,
    0.406f,  0.615f,  0.116f,
    0.676f,  0.977f,  0.133f,
    0.971f,  0.572f,  0.833f,
    0.140f,  0.616f,  0.489f,
    0.997f,  0.513f,  0.064f,
    0.945f,  0.719f,  0.592f,
    0.543f,  0.021f,  0.978f,
    0.279f,  0.317f,  0.505f,
    0.167f,  0.620f,  0.077f,
    0.347f,  0.857f,  0.137f,
    0.055f,  0.953f,  0.042f,
    0.714f,  0.505f,  0.345f,
    0.783f,  0.290f,  0.734f,
    0.722f,  0.645f,  0.174f,
    0.302f,  0.455f,  0.848f,
    0.225f,  0.587f,  0.040f,
    0.517f,  0.713f,  0.338f,
    0.053f,  0.959f,  0.120f,
    0.393f,  0.621f,  0.362f,
    0.673f,  0.211f,  0.457f,
    0.820f,  0.883f,  0.371f,
    0.982f,  0.099f,  0.879f
};


GLfloat camPos[3] = { 0.0f,  0.0f,  -10.0f };
GLfloat camMat[9] = {
  1.0f,  0.0f,  0.0f,
  0.0f,  1.0f,  0.0f,
  0.0f,  0.0f,  1.0f
};
GLfloat modelPos[3] = { 0.0f,  0.0f,  -30.0f };
GLfloat modelMat[9] = {
  1.0f,  0.0f,  0.0f,
  0.0f,  1.0f,  0.0f,
  0.0f,  0.0f,  1.0f
};
GLfloat light_dir[3] = { 1.0f,  1.0f,  1.0f };



int WIDTH  = 800;
int HEIGHT = 800;


int mouseX, mouseY;
SDL_Window * window     = NULL;
SDL_GLContext   context = NULL;
Quat4f qCamera;



//   http://www.songho.ca/opengl/gl_projectionmatrix.html

void getPerspectiveMatrix( float xmin, float xmax, float ymin, float ymax, float zmin, float zmax, float * mat ){
    float invdx = xmax-xmin;
    float invdy = ymax-ymin;
    float invdz = zmax-zmin;
    mat[0 ]  = 2*zmin*invdx; mat[1 ] = 0;            mat[2 ] =  (xmax+xmin)*invdx;  mat[3 ] = 0;
    mat[4 ]  = 0;            mat[5 ] = 2*zmin*invdy; mat[6 ] =  (ymin+ymax)*invdy;  mat[7 ] = 0;
    mat[8 ]  = 0;            mat[9 ] = 0;            mat[10] = -(zmin+zmax)*invdz;  mat[11] = -2*zmax*zmin*invdz;
    mat[12]  = 0;            mat[13] = 0;            mat[14] = -1;                  mat[15] = 0;
}

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
	object1->nVert    = nVert;
	object1->vertDim  = 3;
	object1->clrDim   = 3;
	object1->vertexes = &vertexes[0];
	object1->colors   = &vertexes[0];
	//object1->normals  = &normals [0];
	//object1->vertexes = &colors[0];
	//object1->colors   = &colors  [0];
	object1->init();

	shader1=new Shader();
	//shader1->init( "shaders/plain_vert.c", "shaders/sphere_frag.c" );
	shader1->init( "shaders/minimal3D_vert.c", "shaders/minimal3D_frag.c" );
	//shader1->init( "shaders/default_vert.c", "shaders/defaut_frag.c" );
	glUseProgram(shader1->shaderprogram);

	mesh.fromFileOBJ("common_resources/turret.obj");

	qCamera.setOne();
}

void draw(){
    glClearColor(0.0, 0.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );

    glEnable( GL_DEPTH_TEST );
    glDepthFunc( GL_LESS );

    getPerspectiveMatrix( -WIDTH, WIDTH, -HEIGHT, HEIGHT, 1.0, 10.0, camMat );

    //Quat4f qCamera_; convert(qCamera,qCamera_);
    Mat3f mouseMat; qCamera.toMatrix(mouseMat);
    printf( " (%3.3f,%3.3f,%3.3f,%3.3f) \n", qCamera.x,qCamera.y,qCamera.z,qCamera.w );
    //printf( " (%3.3f,%3.3f,%3.3f,%3.3f) \n", qCamera.x,qCamera.y,qCamera.z,qCamera.w );
    printf( "mouseMat (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) \n", mouseMat.ax, mouseMat.ay, mouseMat.az,  mouseMat.bx, mouseMat.by, mouseMat.bz,   mouseMat.cx, mouseMat.cy, mouseMat.cz );

    GLuint uloc;
    //uloc = glGetUniformLocation( shader1->shaderprogram, "MVC"   ); glUniformMatrix4fv(uloc, 1, GL_FALSE, MVC );
    uloc = glGetUniformLocation( shader1->shaderprogram, "modelPos" ); glUniform3fv      (uloc, 1, modelPos );
    //uloc = glGetUniformLocation( shader1->shaderprogram, "modelMat" ); glUniformMatrix3fv(uloc, 1, GL_FALSE, modelMat );
    uloc = glGetUniformLocation( shader1->shaderprogram, "modelMat" ); glUniformMatrix3fv(uloc, 1, GL_FALSE, (float*)&mouseMat );
    uloc = glGetUniformLocation( shader1->shaderprogram, "camMat"   ); glUniformMatrix4fv(uloc, 1, GL_FALSE, camMat   );

    uloc = glGetUniformLocation( shader1->shaderprogram, "light_dir"); glUniform3fv(uloc, 1, light_dir  );
    object1->draw();

    SDL_GL_SwapWindow(window);

}


// =============================================================
// ========== BORING DEFAULT RUTINES
// =============================================================

// FUNCTION ======	inputHanding
void inputHanding(){
    float posstep = 0.1;
	SDL_Event event;
	while(SDL_PollEvent(&event)){
		if( event.type == SDL_KEYDOWN ){
            switch( event.key.keysym.sym ){
                case SDLK_ESCAPE: quit(); break;
                case SDLK_w: modelPos[1] +=posstep; break;
                case SDLK_s: modelPos[1] -=posstep; break;
                case SDLK_a: modelPos[0] +=posstep; break;
                case SDLK_d: modelPos[0] -=posstep; break;
                case SDLK_q: modelPos[2] +=posstep*10; break;
                case SDLK_e: modelPos[2] -=posstep*10; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
		}
		if( event.type == SDL_QUIT){ quit();  };
	}

	int dmx,dmy;
	SDL_GetMouseState( &mouseX, &mouseY );
    Uint32 buttons = SDL_GetRelativeMouseState( &dmx, &dmy);
    //printf( " %i %i \n", mx,my );
    float mouseRotSpeed = 0.01;
    if ( buttons & SDL_BUTTON(SDL_BUTTON_RIGHT)) {
        Quat4f q; q.fromTrackball( 0, 0, -dmx*mouseRotSpeed, dmy*mouseRotSpeed );
        printf( " %i %i  (%3.3f,%3.3f,%3.3f,%3.3f) \n", dmx,dmy, q.x,q.y,q.z,q.w );
        qCamera.qmul_T( q );
        //qCamera.normalize();
    }
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
