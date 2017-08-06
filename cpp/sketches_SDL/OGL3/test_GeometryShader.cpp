
// read this tutorial
// https://open.gl/geometry

#include <stdlib.h>
#include <stdio.h>

#include <GL/glew.h>
//#define GL3_PROTOTYPES 1
//#define GL_GLEXT_PROTOTYPES
//#include <GL/gl.h>
#include <SDL2/SDL.h>

#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

#include "Mesh.h"

#include "GLfunctions.h"
#include "GLobjects.h"

#include "GLObject.h"
#include "Shader.h"

int WIDTH  = 800;
int HEIGHT = 800;
SDL_Window   * window  = NULL;
SDL_GLContext  context = NULL;

GLuint vao;
Shader   *shader1,*shConst;
//GLObject * object1;
GLMesh   *glmesh, *gledges;

Mesh mesh;

int mouseX, mouseY;
Quat4f qCamera;
Mat3f  mouseMat;
Vec3f  camPos   = (Vec3f){ 0.0f,0.0f, -10.0f };
Vec3f  modelPos = (Vec3f){ 0.0f,0.0f, -30.0f };


GLuint shaderProgram,vertexShader,fragmentShader,geometryShader;
GLuint vbo;
int npoints =20;

// ============= FUNCTIONS

int frameCount = 0;
bool STOP = false;

void quit();
void die ( char const *msg );
void inputHanding ();
void init();
void draw();
void loop( int niters );

int render_type = 1;

const char geometryShaderSrc[]= GLSL(150 core,
layout(points) in;
layout(line_strip, max_vertices = 33) out;

vec2 monopoleForce(vec2 dp){
    //float r = length(dp);
    return dp/dot(dp,dp);
}

vec2 getForce( vec2 p ){
    vec2 f = vec2(1.0,0.0);
    f += monopoleForce( p-vec2( -0.5,0.0 ) )*+0.1;
    f += monopoleForce( p-vec2(  0.5,0.0 ) )*-0.1;
    //p*=5.0;
    //vec2 f = vec2( sin(p.x)*cos(p.y), -cos(p.x)*sin(p.y) );
    return f;
}

const int n = 33;
const vec2 dq = vec2( cos(6.28318530718/n), sin( 6.28318530718/n ) );
void main(){
    //gl_Position = gl_in[0].gl_Position + vec4(-0.1, 0.0, 0.0, 0.0); EmitVertex();
    //gl_Position = gl_in[0].gl_Position + vec4(0.1, 0.0, 0.0, 0.0);  EmitVertex();
    //vec2 p  = gl_in[0].gl_Position.xy;
    /*
    vec2 q  = vec2( 0.02, 0.0 );
    for(int i=0; i<(n+1); i++){
        gl_Position = gl_in[0].gl_Position + vec4( q, 0.0, 0.0 );
        EmitVertex();
        q = vec2( q.x*dq.x - q.y*dq.y, q.x*dq.y + q.y*dq.x );
    }
    */
    const float dt = 0.1;
    vec2 p      = gl_in[0].gl_Position.xy;
    gl_Position = gl_in[0].gl_Position;
    EmitVertex();
    for(int i=0; i<n; i++){
        p += getForce( p ) * dt;
        gl_Position.xy = p;
        EmitVertex();
    }
    EndPrimitive();
}
);

const char vertexShaderSrc[]= GLSL(150 core,
    in vec2 pos;
    void main(){
        gl_Position = vec4(pos, 0.0, 1.0);
    }
);

const char fragmentShaderSrc[]= GLSL(150 core,
    out vec4 outColor;
    void main(){
        outColor = vec4(1.0, 0.0, 0.0, 1.0);
    }
);

GLuint createShader(GLenum type, const GLchar* src) {
    GLuint shader = glCreateShader(type);
    glShaderSource(shader, 1, &src, nullptr);
    glCompileShader(shader);
    GLint isCompiled;
    glGetShaderiv( shader, GL_COMPILE_STATUS, &isCompiled );
    if( isCompiled == false)    {
        int maxLength;
        glGetShaderiv( shader, GL_INFO_LOG_LENGTH, &maxLength );
        char * errLog = (char *)malloc(maxLength); 								// The maxLength includes the NULL character
        glGetShaderInfoLog( shader, maxLength, &maxLength, errLog );
        printf( " Error in compilation of shader : \n"  );
        printf( " %s \n", errLog );
    }
    return shader;
}

void setup(){

GLuint vertexShader   = createShader(GL_VERTEX_SHADER, vertexShaderSrc);
GLuint fragmentShader = createShader(GL_FRAGMENT_SHADER, fragmentShaderSrc);
GLuint shaderProgram  = glCreateProgram();
glAttachShader(shaderProgram, vertexShader);
glAttachShader(shaderProgram, fragmentShader);

geometryShader = createShader(GL_GEOMETRY_SHADER, geometryShaderSrc);
glAttachShader(shaderProgram, geometryShader);

glLinkProgram (shaderProgram);
glUseProgram  (shaderProgram);

glGenBuffers(1, &vbo);

float *points = new float[2*npoints];
for( int i=0; i<npoints; i++ ){
    int i2 = i<<1;
    points[i2+0]=-0.9f;
    points[i2+1]=-0.5f + i*1.0/npoints;
}

glBindBuffer(GL_ARRAY_BUFFER, vbo);
glBufferData(GL_ARRAY_BUFFER, npoints*2*sizeof(float), points, GL_STATIC_DRAW);

// Create VAO
//GLuint vao;
//glGenVertexArrays(1, &vao);
//glBindVertexArray(vao);

// Specify layout of point data
GLint posAttrib = glGetAttribLocation(shaderProgram, "pos");
glEnableVertexAttribArray(posAttrib);
glVertexAttribPointer(posAttrib, 2, GL_FLOAT, GL_FALSE, 0, 0);

}

void draw(){
    glClearColor(0.0, 0.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );

    glEnable( GL_DEPTH_TEST );
    glDepthFunc( GL_LESS );

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    glDrawArrays(GL_POINTS, 0, npoints);

    SDL_GL_SwapWindow(window);

}


// =============================================================
// ========== BORING DEFAULT RUTINES
// =============================================================

// FUNCTION ======	inputHanding
void inputHanding(){
    float posstep = 0.1;
	SDL_Event event;

	const Uint8 *keys = SDL_GetKeyboardState(NULL);
    if( keys[ SDL_SCANCODE_W  ] ){ modelPos.y += posstep; }
	if( keys[ SDL_SCANCODE_S  ] ){ modelPos.y -= posstep; }
	if( keys[ SDL_SCANCODE_A  ] ){ modelPos.x -= posstep; }
	if( keys[ SDL_SCANCODE_D  ] ){ modelPos.x += posstep; }
    if( keys[ SDL_SCANCODE_Q  ] ){ modelPos.z += posstep; }
	if( keys[ SDL_SCANCODE_E  ] ){ modelPos.z -= posstep; }

	while(SDL_PollEvent(&event)){
		if( event.type == SDL_KEYDOWN ){
            switch( event.key.keysym.sym ){
                case SDLK_ESCAPE: quit(); break;
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
        //printf( " %i %i  (%3.3f,%3.3f,%3.3f,%3.3f) \n", dmx,dmy, q.x,q.y,q.z,q.w );
        //qCamera.qmul_T( q );
        qCamera.qmul( q );
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
