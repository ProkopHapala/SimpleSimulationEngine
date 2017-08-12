
// read this tutorial
// http://www.opengl-tutorial.org/beginners-tutorials/tutorial-4-a-colored-cube/

// http://gamedev.stackexchange.com/questions/93055/getting-the-real-fragment-depth-in-glsl

// TODO see later vbo indexing to avoid duplicating triangles
// http://www.opengl-tutorial.org/intermediate-tutorials/tutorial-9-vbo-indexing/

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

#include "IO_utils.h"

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

void setup(){

    char* names[] = { "shade3D.frag", "color3D.vert", "color3D.frag" };

    //char* src = fileGetSection( "common_resources/shaders/Basic_small.glslf", "//>>const3D.vert", "//<<" );
    //printf("src : \n %s\n", src);

    int nkeys = 3;
    char ** srcs = fileGetSections( "common_resources/shaders/Basic.glslf", nkeys, names, "//>>" );
    for(int ikey=0; ikey<nkeys; ikey++){ printf("##### shader: %s\n%s\n", names[ikey], srcs[ikey] ); }
    //exit(0);

    mesh.fromFileOBJ( "common_resources/turret.obj" );
    mesh.polygonsToTriangles(false);
    mesh.tris2normals(true);
    mesh.findEdges( );

    /*
    object1 = new GLObject( );
    object1->nVert   = mesh.points.size();
    object1->setIndexes( mesh.triangles.size()*3, (int*)&mesh.triangles[0] );
    object1->buffs[0].setup(0,3,GL_FALSE,(double*)&(mesh.points [0]), object1->nVert, 'v'); // vertexes
    object1->buffs[1].setup(1,3,GL_FALSE,(double*)&(mesh.normals[0]), object1->nVert, 'n'); // normals
    object1->init();
    */

    glmesh = new GLMesh();
    glmesh->init_d( mesh.points.size(), mesh.triangles.size()*3, ((int*)&mesh.triangles[0]), (double*)&(mesh.points [0]), (double*)&(mesh.normals[0]), NULL, NULL );
    //glmesh->init_d( mesh.points.size(), 0, NULL, (double*)&(mesh.points [0]), (double*)&(mesh.normals[0]), NULL, NULL );

    gledges = new GLMesh();
    *gledges = *glmesh;
    gledges->nInds = mesh.edges.size()*2;
    Vec2i * evts = mesh.exportEdgeVs();
    newArrayBuffer( gledges->inds, gledges->nInds*3*sizeof(GLuint), (int*)evts, GL_STATIC_DRAW );
    delete [] evts;

	shConst=new Shader();
	shConst->init_default();
	//shConst->init( "common_resources/shaders/const3D.glslv",   "common_resources/shaders/const3D.glslf"   );
	//shConst->init( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/color3D.glslf"   );
	//shConst->init( "common_resources/shaders/shade3D.glslv",   "common_resources/shaders/shade3D.glslf"   );
	shConst->getDefaultUniformLocation();

	shader1=new Shader();
	char*  shader_names[] = { "shade3D.vert", "shade3D.frag" };
	char** shader_srcs    = fileGetSections( "common_resources/shaders/Basic.glslf", 2, shader_names, "//>>" );
	for(int i=0; i<2; i++){ printf("##### shader: %s\n%s\n", shader_names[i], shader_srcs[i] ); }
	shader1->init_str(shader_srcs[0],shader_srcs[1],NULL);
	//saveStr("shade3D.vert.glslv", shader_srcs[0]);
	//saveStr("shade3D.frag.glslv", shader_srcs[1]);
    //shader1->init( "shade3D.vert.glslv",   "shade3D.frag.glslv"   );
	//exit(0);
	//shader1->init( "common_resources/shaders/shade3D.glslv",   "common_resources/shaders/shade3D.glslf"   );
	shader1->getDefaultUniformLocation();

	shader1->use();
    glUniform3fv( shader1->getUloc("light_pos"    ), 1, (const float[]){ 1.0f,  1.0f, -1.0f } );
    glUniform3fv( shader1->getUloc("lightColor"   ), 1, (const float[]){ 1.0f,  0.9f,  0.8f } );
    glUniform3fv( shader1->getUloc("diffuseColor" ), 1, (const float[]){ 1.0f,  1.0f,  1.0f } );
    glUniform3fv( shader1->getUloc("ambientColor" ), 1, (const float[]){ 0.2f,  0.2f,  0.3f } );
    glUniform3fv( shader1->getUloc("specularColor"), 1, (const float[]){ 0.0f,  0.0f,  0.0f } );

    qCamera.setOne();

}

void draw(){
    glClearColor(0.0, 0.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );

    glEnable( GL_DEPTH_TEST );
    glDepthFunc( GL_LESS );

    qCamera.toMatrix(mouseMat);
    Mat4f camMat;   camMat.setPerspective( 20.0, 20.0, -2.0, -1000.0 );

    shader1->use();
    shader1->set_modelPos( (GLfloat*)&modelPos );
    shader1->set_modelMat( (GLfloat*)&mouseMat );
    //shader1->set_camPos  ( (GLfloat*)&camPos );
    shader1->set_camMat  ( (GLfloat*)&camMat );

    // shading
    //uloc = glGetUniformLocation( shader1->shaderprogram, "light_dir"); glUniform3fv(uloc, 1, light_dir  );
    //object1->draw();

    glmesh->draw();
    //glmesh->draw( GL_TRIANGLES );
    //glPointSize(10.0); glmesh->draw( GL_POINTS );
    //glLineWidth(30.0); glmesh->draw( GL_LINE_STRIP );

    shConst->use();
    shConst->set_modelPos( (GLfloat*)&modelPos );
    shConst->set_modelMat( (GLfloat*)&mouseMat );
    shConst->set_camMat  ( (GLfloat*)&camMat );
    glUniform4fv( shConst->getUloc("baseColor"), 1, (const float[]){1.0f,0.0f,0.0f,1.0f} ); gledges->drawPoints(10.0);
    glUniform4fv( shConst->getUloc("baseColor"), 1, (const float[]){0.0f,1.0f,0.0f,1.0f} ); glLineWidth(3.0); gledges->draw(GL_LINES);

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
