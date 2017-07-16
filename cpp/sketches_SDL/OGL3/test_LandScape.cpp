
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
float ASPECT_RATIO = HEIGHT/WIDTH;

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

Vec3f camPos = (Vec3f){ 0.0f, 0.0f, 0.0f };
Mat4f camMat,mRot,mPersp;
Mat3f mouseMat;

float pitch=0,yaw=0;

// ===============================================
// ======================= Functions
// ===============================================

Vec3d terrainFunc( Vec2d p ){ return (Vec3d){p.x*10.0,sin(p.x)*sin(p.y*0.5)*10.0,p.y*10.0}; };

void setup(){

    shader1=new Shader();

    //shader1->init( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/color3D.glslf"   );
    shader1->init( "common_resources/shaders/const3D.glslv",   "common_resources/shaders/const3D.glslf"   );
    //shader1->init( "common_resources/shaders/const3D.glslv",   "common_resources/shaders/pointSprite.glslf"   );

    /*
    shader1->init( "common_resources/shaders/shade3D.glslv",   "common_resources/shaders/shade3D.glslf"   );
    GLuint uloc;
    uloc = glGetUniformLocation( shader1->shaderprogram, "cam_pos"       ); glUniform3fv      (uloc, 1, (GLfloat*)&camPos );
    uloc = glGetUniformLocation( shader1->shaderprogram, "light_pos"     ); glUniform3fv      (uloc, 1, light_pos         );
    uloc = glGetUniformLocation( shader1->shaderprogram, "lightColor"    ); glUniform3fv      (uloc, 1, lightColor        );
    uloc = glGetUniformLocation( shader1->shaderprogram, "diffuseColor"  ); glUniform3fv      (uloc, 1, diffuseColor      );
    uloc = glGetUniformLocation( shader1->shaderprogram, "ambientColor"  ); glUniform3fv      (uloc, 1, ambientColor      );
    uloc = glGetUniformLocation( shader1->shaderprogram, "specularColor" ); glUniform3fv      (uloc, 1, specularColor     );
    */

    /*
    int nVert = countVerts( Solids::Icosahedron_nfaces, Solids::Icosahedron_ngons );
    GLfloat * verts   = new GLfloat[nVert*3];
    GLfloat * normals = new GLfloat[nVert*3];
    hardFace( Solids::Icosahedron_nfaces, Solids::Icosahedron_ngons, Solids::Icosahedron_faces, Solids::Icosahedron_verts, verts, normals );

    object1 = new GLObject( );
    object1->nVert    = nVert;
    object1->buffs[0].setup(0,3,GL_FALSE,verts,  'v'); // vertexes
    object1->buffs[1].setup(1,3,GL_FALSE,normals,'n'); // normals
    object1->init();
    */

    /*
    object1 = new GLObject( );
    object1->setup( countVerts( Solids::Icosahedron_nfaces, Solids::Icosahedron_ngons ) );
    hardFace( Solids::Icosahedron_nfaces, Solids::Icosahedron_ngons, Solids::Icosahedron_faces, Solids::Icosahedron_verts, object1->buffs[0].cbuff, object1->buffs[1].cbuff );
    object1->init();
    */
    /*
    object1 = new GLObject( );
    object1->setup( countVerts( Solids::Cube_nfaces, Solids::Cube_ngons ) );
    hardFace( Solids::Cube_nfaces, Solids::Cube_ngons, Solids::Cube_faces, Solids::Cube_verts, object1->buffs[0].cbuff, object1->buffs[1].cbuff );
    object1->init();
    */

    object1 = new GLObject( );
    object1->setup( countVerts( Solids::Octahedron_nfaces, Solids::Octahedron_ngons ) );
    hardFace( Solids::Octahedron_nfaces, Solids::Octahedron_ngons, Solids::Octahedron_faces, Solids::Octahedron_verts, object1->buffs[0].cbuff, object1->buffs[1].cbuff );
    object1->init();

    obj_terrain = qaudPatchHard( 100, (Vec2d){-50.0,-50.0}, (Vec2d){1.0,0.0}, (Vec2d){0.0,1.0}, terrainFunc );

    printf( "obj_terrain nVert %i", obj_terrain->nVert);
    //exit(0);
    /*
    for(int i=0; i<obj_terrain->nVert; i++){
        Vec3f& p  = ((Vec3f*)obj_terrain->buffs[0].cbuff)[i];
        Vec3f& nv = ((Vec3f*)obj_terrain->buffs[1].cbuff)[i];
        printf( " %i (%g,%g,%g) (%g,%g,%g) \n", i, p.x, p.y, p.z,   nv.x, nv.y, nv.z );
    }
    */
    //exit(0);

    ninstancs = 100; // 30 ms/frame
    instance_points = new Vec3f[ninstancs];
    //for (int i=0; i<ninstancs; i++){ instance_points[i] = (Vec3f){randf(-15.0,15.0),randf(-15.0,15.0),randf(5.0,100.0)};};
    for (int i=0; i<ninstancs; i++){ instance_points[i] = (Vec3f){0.0,0.0,5.0*i};};


    shader1->getDefaultUniformLocation();
	qCamera.setOne();
	delay = 1;
}

//=========================
void draw(){

    long time_start = getCPUticks();

    glClearColor(0.8, 0.8, 0.8, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );

    glEnable( GL_DEPTH_TEST );
    glDepthFunc   ( GL_LESS );

    qCamera.toMatrix(mouseMat);

    //mouseMat.fromEuler( yaw, pitch, 0.0 );
    //mouseMat.fromEuler( 0.0, pitch, yaw );

    //printf("====\n"); mouseMat.print();

    float fov = 3.0;
    mRot.setOne(); mRot.set(mouseMat);
    mPersp.getPerspectiveMatrix( fov, fov*ASPECT_RATIO, -1.0, -1000.0 );
    //camMat.set_mmul_TN( m1, mPersp );
    //camMat.set_mmul( m1, mPersp );
    mRot = mRot.transposed( );
    camMat.set_mmul( mRot, mPersp );
    //camMat.set_mmul( mPersp, mRot );
    Mat3f objRot; objRot.setOne();

    // ============= Objects
    glUseProgram(shader1->shaderprogram);

    shader1->set_camMat  ( (float*)&camMat );
    shader1->set_modelMat( (float*)&objRot );

    Vec3f  p;
    Quat4f c;

    //object1->draw();
    object1->preDraw();
    for(int i=0; i<ninstancs; i++){
        p = instance_points[i] - camPos;
        shader1->set_modelPos ( (GLfloat*)&p );
        c.set( i*0.01f, 0.5f,1-i*0.01f, 1.0f  );
        //c.set( 1.0f, 0.0f, 1.0f, 1.0f );
        shader1->set_baseColor( (GLfloat*)&c );
        object1->draw_instance();
    }
    object1->afterDraw();

    p = (Vec3f){0.0,0.0,0.0} - camPos;

    shader1->set_modelPos( (GLfloat*)&p );
    obj_terrain->draw_mode = GL_POINTS; glPointSize( 15.0 ); obj_terrain->draw_default();
    obj_terrain->draw_mode = GL_LINES;  glLineWidth( 5.0 ); obj_terrain->draw_default();

    /*
    // https://www.khronos.org/opengl/wiki/Primitive#Point_primitives
    obj_terrain->preDraw();
    //obj_terrain->draw_mode = GL_TRIANGLES;
    obj_terrain->draw_instance();
    //p.add(0.0,1.0,0.0);  glUniform3fv( uloc_pos, 1, (GLfloat*)&p );
    //obj_terrain->draw_mode = GL_POINTS;    obj_terrain->draw_instance();
    //obj_terrain->draw_mode = GL_LINES;     obj_terrain->draw_instance();
    obj_terrain->afterDraw();
    */

    SDL_GL_SwapWindow(window);
}

// =============================================================
// ========== BORING DEFAULT RUTINES
// =============================================================

// FUNCTION ======	inputHanding
void inputHanding(){

    float posstep = 0.1f; if(RayTerrain){ posstep = 2.0f; }
    float step          = 0.25f;
    float keyRotSpeed   = 0.002f;

    const Uint8 *keys = SDL_GetKeyboardState(NULL);

    //if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.dyaw  (  (float)keyRotSpeed ); }
	//if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.dyaw  ( -(float)keyRotSpeed ); }
    if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.roll2  (  (float)keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.roll2  ( -(float)keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.dpitch2(  (float)keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.dpitch2( -(float)keyRotSpeed ); }

    if( keys[ SDL_SCANCODE_W  ] ){ camPos.add_mul( mouseMat.c, +step ); }
	if( keys[ SDL_SCANCODE_S  ] ){ camPos.add_mul( mouseMat.c, -step ); }
	if( keys[ SDL_SCANCODE_A  ] ){ camPos.add_mul( mouseMat.a, -step ); }
	if( keys[ SDL_SCANCODE_D  ] ){ camPos.add_mul( mouseMat.a, +step ); }

	SDL_Event event;
	while(SDL_PollEvent(&event)){
		if( event.type == SDL_KEYDOWN ){
            switch( event.key.keysym.sym ){
                case SDLK_ESCAPE: quit(); break;
                case SDLK_KP_PLUS:  terrain_size[0] *=1.1; terrain_size[2] *=1.1; break;
                case SDLK_KP_MINUS: terrain_size[0] /=1.1; terrain_size[2] /=1.1; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            printf( "" );
		}
		if( event.type == SDL_QUIT){ quit();  };
	}

	int dmx,dmy;
	SDL_GetMouseState( &mouseX, &mouseY );
    Uint32 buttons = SDL_GetRelativeMouseState( &dmx, &dmy);
    //printf( " %i %i \n", mx,my );
    float mouseRotSpeed = 0.01;
    if ( buttons & SDL_BUTTON(SDL_BUTTON_RIGHT)) {
        Quat4f q; q.fromTrackball( 0, 0, -dmx*mouseRotSpeed, dmy*mouseRotSpeed ); qCamera.qmul_T( q );
        //qCamera.dyaw2(-dmx*mouseRotSpeed); qCamera.dpitch2(-dmy*mouseRotSpeed);
        //qCamera.dpitch2(-dmy*mouseRotSpeed); qCamera.dyaw2(-dmx*mouseRotSpeed);
        //qCamera.normalize();

        //pitch +=  dmy*mouseRotSpeed;
        //yaw   +=  dmx*mouseRotSpeed;

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
	loop( 10000000 );
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
