#include <stdlib.h>
#include <stdio.h>

#include <GL/glew.h>
//#define GL3_PROTOTYPES 1
//#define GL_GLEXT_PROTOTYPES
//#include <GL/gl.h>
#include <SDL2/SDL.h>

#include "Vec2.h"
#include "Vec3.h"

#include "GLObject.h"
#include "Shader.h"
#include "GLfunctions.h"
#include "GLobjects.h"

#include "Mesh.h"

// ============= GLOBAL VARIABLES

const int WIDTH  = 800;
const int HEIGHT = 800;
SDL_Window * window     = NULL;
SDL_GLContext   context = NULL;

GLuint    vao;
Shader   *shader_pre,*shader_post;
GLObject *object1;

GLMesh      *glmesh;
FrameBuffer  frameBuff1;

int mouseX, mouseY;
Quat4f qCamera;
Mat3f  mouseMat;
Vec3f  camPos   = (Vec3f){ 0.0f,0.0f, -10.0f };
Vec3f  modelPos = (Vec3f){ 0.0f,0.0f, -30.0f };

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

    Mesh mesh;
    mesh.fromFileOBJ( "common_resources/turret.obj" );
    mesh.polygonsToTriangles(false);
    mesh.tris2normals(true);
    mesh.findEdges( );

    glmesh = new GLMesh();
    glmesh->init_d( mesh.points.size(), mesh.triangles.size()*3, ((int*)&mesh.triangles[0]), (double*)&(mesh.points [0]), (double*)&(mesh.normals[0]), NULL, NULL );

    shader_pre=new Shader();
	shader_pre->init( "common_resources/shaders/shade3D.glslv",   "common_resources/shaders/shade3D.glslf"   );
	shader_pre->getDefaultUniformLocation();

    shader_pre->use();
    glUniform3fv( shader_pre->getUloc("light_pos"    ), 1, (const float[]){ 1.0f,  1.0f, -1.0f } );
    glUniform3fv( shader_pre->getUloc("lightColor"   ), 1, (const float[]){ 1.0f,  0.9f,  0.8f } );
    glUniform3fv( shader_pre->getUloc("diffuseColor" ), 1, (const float[]){ 1.0f,  1.0f,  1.0f } );
    glUniform3fv( shader_pre->getUloc("ambientColor" ), 1, (const float[]){ 0.2f,  0.2f,  0.3f } );
    glUniform3fv( shader_pre->getUloc("specularColor"), 1, (const float[]){ 0.0f,  0.0f,  0.0f } );

    qCamera.setOne();

    // ------------- object

    static const GLfloat vertexes[4][2] = {
        {  -0.9f,  -0.9f  },
        {  -0.9f,   0.9f  },
        {   0.9f,  -0.9f  },
        {   0.9f,   0.9f  } };

	object1 = new GLObject( );
	object1->draw_mode = GL_TRIANGLE_STRIP;
	object1->nVert   = 4;
    object1->buffs[0].setup(0,2,GL_FALSE,&vertexes[0][0],'v'); // vertexes
	object1->init();

    // ------------- shader for blitting from texture

    shader_post=new Shader();
	shader_post->init( "shaders/plain_vert.c", "shaders/SSAO_frag.c" );

	shader_post->use();
    glUniform2fv( shader_post->getUloc( "resolution" ), 1, (const float[]){WIDTH,HEIGHT} );

    // ---- prepare FrameBuffer

    //newTexture2D   ( texRGB, WIDTH, HEIGHT, NULL, GL_RGB,             GL_UNSIGNED_BYTE );
    //newTexture2D   ( texZ  , WIDTH, HEIGHT, NULL, GL_DEPTH_COMPONENT, GL_FLOAT         );
    //frameBuff1.init( texRGB, texZ, WIDTH, HEIGHT );

    frameBuff1.init( WIDTH, HEIGHT );

    printf( "========== SETUP DONE ! \n");

}

void draw(){
    glClearColor(0.0, 0.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );

    // ------- render to texture

    glBindFramebuffer(GL_FRAMEBUFFER, frameBuff1.buff );

    glClearColor(0.0, 0.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );
    glEnable( GL_DEPTH_TEST );
    glDepthFunc( GL_LESS );

    qCamera.toMatrix(mouseMat);
    Mat4f camMat;   camMat.setPerspective( 20.0, 20.0, -2.0, -1000.0 );

    shader_pre->use();
    shader_pre->set_modelPos( (GLfloat*)&modelPos );
    shader_pre->set_modelMat( (GLfloat*)&mouseMat );
    //shader1->set_camPos  ( (GLfloat*)&camPos );
    shader_pre->set_camMat  ( (GLfloat*)&camMat );

    glmesh->draw();

    // -------- from texture to Screen

    shader_post->use();
    glUniform1i(shader_post->getUloc("texZ"  ), 0);
    glUniform1i(shader_post->getUloc("texRGB"), 1);

    glActiveTexture(GL_TEXTURE0);    glBindTexture  (GL_TEXTURE_2D, frameBuff1.texZ   );
    glActiveTexture(GL_TEXTURE1);    glBindTexture  (GL_TEXTURE_2D, frameBuff1.texRGB );

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    //glViewport(0,0,WIDTH,HEIGHT);
    glEnableVertexAttribArray(0);
    object1->draw();

    SDL_GL_SwapWindow(window);
}

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
