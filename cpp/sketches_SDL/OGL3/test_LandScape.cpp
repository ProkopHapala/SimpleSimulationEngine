
// http://outerra.blogspot.nl/2012/11/maximizing-depth-buffer-range-and.html

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
#include "Mat4.h"

#include "Solids.h"
#include "Noise.h"
#include "Mesh.h"

#include "GL3Utils.h"
#include "GLObject.h"
#include "Shader.h"

#include "testUtils.h"

//============ Globals

Shader   * shader1;
GLObject * object1,*obj_terrain;

Mesh mesh;

GLuint vao;     // vertex array object
GLuint textureID;
GLuint uloc;

GLfloat modelPos[3] = { 0.0f,  0.0f,  -5.0f };
GLfloat modelMat[9] = {
  1.0f,  0.0f,  0.0f,
  0.0f,  1.0f,  0.0f,
  0.0f,  0.0f,  1.0f
};

Vec3f cam_pos = (Vec3f){ 0.0f,  0.0f,  -10.0f };
//GLfloat cam_pos      [3] = { 0.0f,  0.0f,  -10.0f };
GLfloat light_pos    [3] = { 1.0f,  1.0f,   -1.0f };
GLfloat lightColor   [3] = { 1.0f,  0.9f,   0.8f  };
//GLfloat lightColor   [3] = { 1.0f,  1.0f,   1.0f  };
GLfloat diffuseColor [3] = { 1.0f,  1.0f,   1.0f  };
//GLfloat diffuseColor [3] = { 0.0f,  0.0f,   0.0f  };
GLfloat ambientColor [3] = { 0.2f,  0.2f,   0.3f  };
GLfloat specularColor[3] = { 1.0f,  1.0f,   1.0f  };
//GLfloat specularColor[3] = { 0.0f,  0.0f,   0.0f  };

float resolution  [2] = {800,800};
float terrain_0   [2] = {128.0,128.0};
float terrain_size[2] = {256.0,256.0};

int ninstancs = 10;
//GLfloat * instance_points;

Vec3f * instance_points;

int WIDTH  = 800;
int HEIGHT = 800;

int mouseX, mouseY;
SDL_Window * window     = NULL;
SDL_GLContext   context = NULL;
Quat4f qCamera;

int frameCount = 0;
bool STOP = false;

void quit();
void die ( char const *msg );
void inputHanding ();
void init();
void draw();
void loop( int niters );

int  render_type  = 1;
bool terrain_mode = true;
bool RayTerrain   = false;

// speed test
int delay = 1; int VSync = 0;
//int delay = 10; int VSync = 1;
//float camMat[16];

Mat4f camMat;

Mat3f mouseMat;

int   Ter_nquads = 100+1;
float Ter_tg     = 0.4;
float Ter_z0     = 1.0;
float Ter_dx     = 1.0/(Ter_nquads-1);
float Ter_dz     = 2.0*Ter_dx;
float Ter_fsc    = 1+Ter_dz*Ter_tg;

// ===============================================
// ======================= Functions
// ===============================================

void setup(){

    if ( render_type == 0      ){
        // --- vertex const color
        shader1=new Shader();
        shader1->init( "shaders/basicColor3D_vert.c", "shaders/basicColor3D_frag.c" );
    }else if ( render_type == 1 ){
        // --- shading
        shader1=new Shader();
        shader1->init( "shaders/basicShading3D_vert.c", "shaders/basicShading3D_frag.c" );
        glUseProgram(shader1->shaderprogram);
    };
	//mesh.fromFileOBJ("common_resources/turret.obj");

    int nVert = countVerts( Solids::Icosahedron_nfaces, Solids::Icosahedron_ngons );
    GLfloat * verts   = new GLfloat[nVert*3];
    GLfloat * normals = new GLfloat[nVert*3];
    hardFace( Solids::Icosahedron_nfaces, Solids::Icosahedron_ngons, Solids::Icosahedron_faces, Solids::Icosahedron_verts, verts, normals );

    object1 = new GLObject( );
    object1->nVert    = nVert;
    object1->buffs[0].setup(0,3,GL_FALSE,verts,  'v'); // vertexes
    object1->buffs[1].setup(1,3,GL_FALSE,normals,'n'); // normals
    object1->init();

        // shading
    if ( render_type == 1 ){
        uloc = glGetUniformLocation( shader1->shaderprogram, "cam_pos"       ); glUniform3fv      (uloc, 1, (GLfloat*)&cam_pos      );
        uloc = glGetUniformLocation( shader1->shaderprogram, "light_pos"     ); glUniform3fv      (uloc, 1, light_pos     );
        uloc = glGetUniformLocation( shader1->shaderprogram, "lightColor"    ); glUniform3fv      (uloc, 1, lightColor    );
        uloc = glGetUniformLocation( shader1->shaderprogram, "diffuseColor"  ); glUniform3fv      (uloc, 1, diffuseColor  );
        uloc = glGetUniformLocation( shader1->shaderprogram, "ambientColor"  ); glUniform3fv      (uloc, 1, ambientColor  );
        uloc = glGetUniformLocation( shader1->shaderprogram, "specularColor" ); glUniform3fv      (uloc, 1, specularColor );
    };

    ninstancs = 100; // 30 ms/frame

    instance_points = new Vec3f[ninstancs];
    for (int i=0; i<ninstancs; i++){ instance_points[i] = (Vec3f){randf(-15.0,15.0),randf(-15.0,15.0),randf(-20.0,-500.0)};};

    Vec3f * strip = new Vec3f[Ter_nquads*2];
    for(int i=0; i<Ter_nquads; i++){
        int i2 = i<<1;
        float x  = (i-0.5*Ter_nquads)*Ter_dx;
        strip[i].set(-1.0, randf(), x);
        strip[i].set( 1.0, randf(), x);
    }
    obj_terrain = new GLObject();
    obj_terrain->draw_mode = GL_TRIANGLE_STRIP;
    obj_terrain->nVert     = Ter_nquads*2;
    obj_terrain->buffs[0].setup(0,3,GL_FALSE,strip,'v');
    obj_terrain->init();

	qCamera.setOne();

	delay = 1;
}

//=========================
void draw(){

    long time_start = getCPUticks();

    glClearColor(0.0, 0.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );

    glEnable( GL_DEPTH_TEST );
    glDepthFunc( GL_LESS );

    long time_0 = getCPUticks();

    qCamera.toMatrix(mouseMat);


    Mat4f m1; m1.setOne(); m1.set(mouseMat);
    //Mat4f m2; m2.getPerspectiveMatrix( -WIDTH, WIDTH, -HEIGHT, HEIGHT, 1.0, 20.0 );
    Mat4f m2; m2.getPerspectiveMatrix( -WIDTH*0.0001, WIDTH*0.0001, -HEIGHT*0.0001, HEIGHT*0.0001, 5.0, 1000.0 );

    camMat.set_mmul_TN( m1, m2 );

    // ============= Terrain

    long time_1 = getCPUticks();

    // ============= Objects
    glUseProgram(shader1->shaderprogram);

    uloc = glGetUniformLocation( shader1->shaderprogram, "camMat"   ); glUniformMatrix4fv(uloc, 1, GL_FALSE, (float*)&camMat   );
    uloc = glGetUniformLocation( shader1->shaderprogram, "modelMat" ); glUniformMatrix3fv(uloc, 1, GL_FALSE, (float*)&mouseMat );
    uloc = glGetUniformLocation( shader1->shaderprogram, "modelPos" ); // glUniform3fv      (uloc, 1, modelPos );

    //object1->draw();
    object1->preDraw();


    for(int i=0; i<ninstancs; i++){
        //glUniform3fv( uloc, 1, instance_points+i*3 );
        Vec3f p = instance_points[i] - cam_pos;
        glUniform3fv( uloc, 1, (GLfloat*)&p );
        //glDrawArrays( object1->draw_mode, 0, object1->nVert);
        object1->draw_instance();
    }
    object1->afterDraw();

    obj_terrain->preDraw();
    glUniform3fv (uloc, 1, (GLfloat*)modelPos );
    obj_terrain->draw_instance();
    obj_terrain->afterDraw();

    long time_2 = getCPUticks();

    SDL_GL_SwapWindow(window);

    long time_3 = getCPUticks();


    double Ttot     = (time_3-time_start)*1e-6;
    double Tterrain = (time_1-time_0)*1e-6;
    double Tobject  = (time_2-time_1)*1e-6;
    double Tswap    = (time_3-time_2)*1e-6;

    //printf("camPos (%g,%g,%g) \n", cam_pos[0], cam_pos[1], cam_pos[2]);
    printf("camPos (%g,%g,%g) \n", cam_pos.x, cam_pos.y, cam_pos.z );
    //printf( "Ttot %3.2f terrain %3.2f objects %3.2f swap %3.2f \n", Ttot, Tterrain, Tobject, Tswap );

}

// =============================================================
// ========== BORING DEFAULT RUTINES
// =============================================================

// FUNCTION ======	inputHanding
void inputHanding(){

    float posstep = 0.1; if(RayTerrain){ posstep = 2.0; }
    double step = 1.0;

	SDL_Event event;
	while(SDL_PollEvent(&event)){
		if( event.type == SDL_KEYDOWN ){
            switch( event.key.keysym.sym ){
                case SDLK_ESCAPE: quit(); break;

                case SDLK_w: cam_pos.add_mul( mouseMat.c,  step ); break;
                case SDLK_s: cam_pos.add_mul( mouseMat.c, -step ); break;
                case SDLK_a: cam_pos.add_mul( mouseMat.a,  step ); break;
                case SDLK_d: cam_pos.add_mul( mouseMat.a, -step ); break;

                case SDLK_KP_PLUS:  terrain_size[0] *=1.1; terrain_size[2] *=1.1; break;
                case SDLK_KP_MINUS: terrain_size[0] /=1.1; terrain_size[2] /=1.1; break;
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
        SDL_Delay(delay);
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
    //SDL_GL_SetSwapInterval(1); // VSync On
    SDL_GL_SetSwapInterval(VSync);

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
