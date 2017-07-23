
// Copied tutorial from
//  http://www.opengl-tutorial.org/intermediate-tutorials/billboards-particles/particles-instancing/
//  https://github.com/opengl-tutorials/ogl/tree/master/tutorial18_billboards_and_particles

#include <stdlib.h>
#include <stdio.h>

#include <vector>
#include <algorithm>

#include <GL/glew.h>
#include <SDL2/SDL.h>

#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "Mat4.h"

#include "GL3Utils.h"
#include "GLObject.h"
#include "GLfunctions.h"
#include "GLInstances.h"
#include "IO_utils.h"
#include "Shader.h"

#include "Solids.h"
#include "CMesh.h"

#include "testUtils.h"

// =============== Global variables

//const int MaxParticles = 100000;
const int MaxParticles = 80000;
//const int MaxParticles = 16;
//const int MaxParticles = 256;
//const int MaxParticles = 256*256;

//GLuint VertexArrayID;
GLfloat* instance_pos;
GLfloat* instance_dir;
GLfloat* instance_Up;
GLfloat* instance_sc;

//GLuint programID;
Shader *shGeom,*shFraged,*shFragedDepth;
//GLuint CameraRight_worldspace_ID;
//GLuint CameraUp_worldspace_ID;
//GLuint ViewProjMatrixID;

GLInstances instances;
//GLuint billboard_vertex_buffer;
//GLuint particles_color_buffer;
//GLuint particles_position_buffer;

int nVerts = 0;
int ParticlesCount = 0;
int frameCount = 0;
double lastTime = 0.0;

GLuint vao;     // vertex array object
int delay = 10;
int VSync = 0;
SDL_Window * window     = NULL;
SDL_GLContext   context = NULL;
int WIDTH  = 800;
int HEIGHT = 800;
float ASPECT_RATIO = HEIGHT/WIDTH;

int mouseX, mouseY;
Quat4f qCamera;
Mat3f  mouseMat;
Vec3f  camPos = (Vec3f){ 0.0f, 0.0f, 0.0f };

bool viewGeom = false;
bool customDepth = false;

long lastCPUtick = 0;
double ticks_per_second=0;

// =============== Functions

const char str_glslf_sin[]= GLSL(330,
in        vec3 fpos_world;
out       vec4 gl_FragColor;
void main(){ gl_FragColor = vec4( sin( fpos_world*30.0 ), 1.0 ); }
);

char* strconcat( char* str1, char * str2 ){
    int n1 = strlen(str1);
    int n2 = strlen(str2);
    //printf("%s\n", str1);
    //printf("%s\n", str2);
    char* str12 = new char[n1+n2+1];
    strcpy(str12,str1);
    strcat(str12,str2);
    return str12;
}

void replaceStr_inpl(char* str, char* mask, char* mod){
    char * found  = strstr(str,mask);
    if(found) strncpy(str,mod,strlen(mask));
}

char* replaceStr(char* str, char* mask, char* mod){
    char * str_  = new char[strlen(str)];
    strcpy(str_,str);
    char * found = strstr(str_,mask);
    if(found) strncpy(found,mod,strlen(mask));
    return str_;
}

double calibrate_timer(int delay){
    long t1 = getCPUticks();
    SDL_Delay(delay);
    long t2 = getCPUticks();
    return (t2-t1)/(delay*0.001d);
}

int setup(){

    shFraged=new Shader();
    //shFraged->init( "common_resources/shaders/Instance3D.glslv",   "common_resources/shaders/Sphere3D.glslf" );
    shFraged->init( "common_resources/shaders/Instance3D.glslv",   "common_resources/shaders/Sphere3D.glslf" );
    shFraged->getDefaultUniformLocation();

    char* str_glslv_Instance3D     = filetobuf( "common_resources/shaders/Instance3D.glslv"  );
    char* str_glslf_Sphere3D       = filetobuf( "common_resources/shaders/Sphere3D.glslf"    );
    char* str_glslf_Sphere3D_depth = replaceStr( str_glslf_Sphere3D, "#define CUSTOM_DEPTH_0", "#define CUSTOM_DEPTH_1");
    printf("str_glslf_Sphere3D_depth:>>%s<<\n", str_glslf_Sphere3D_depth );

    shGeom=new Shader();
	shGeom-> init_str ( str_glslv_Instance3D, str_glslf_sin );
    shGeom->getDefaultUniformLocation();

    shFragedDepth=new Shader();
    shFragedDepth->init_str( str_glslv_Instance3D, str_glslf_Sphere3D_depth );
    shFragedDepth->getDefaultUniformLocation();

    delete [] str_glslv_Instance3D; delete [] str_glslf_Sphere3D; delete [] str_glslf_Sphere3D_depth;

	instance_pos    = new GLfloat[MaxParticles*4];
	instance_dir    = new GLfloat[MaxParticles*4];
	instance_Up     = new GLfloat[MaxParticles*4];
	instance_sc     = new GLfloat[MaxParticles*4];

	float span = 20.0f;
    float time = frameCount * 0.005;
    ParticlesCount = MaxParticles;

    int nside = int( pow( ParticlesCount, 1.0/3.0) );

	for( int i=0; i<ParticlesCount; i++ ){
        Vec3f pos, dir, up, sc;
        // pos.set( randf(-span,span), randf(-span,span), randf(-span,span) );
        // sc.set( randf(0.5,1.5),randf(0.5,1.5),randf(0.5,1.5) );
        pos.set( i%nside , (i/nside )%nside , (i/(nside*nside))%nside ); pos.mul(2.5);
        //sc.set( 1.0,1.0,1.0 );
        float sz = randf(0.5,1.5); sc.set( sz,sz,sz );

        dir.set( randf(-1.0f,1.0f), randf(-1.0f,1.0f), randf(-1.0f,1.0f) ); dir.normalize();
        up .set( randf(-1.0f,1.0f), randf(-1.0f,1.0f), randf(-1.0f,1.0f) ); up.makeOrthoU( dir ); up.normalize();
        //printf( "%i %g  %g %g \n", i, up.dot(dir),    dir.norm2(), up.norm2()  );
        *((Vec3f*)(instance_pos+(3*i))) = pos;
        *((Vec3f*)(instance_dir+(3*i))) = dir;
        *((Vec3f*)(instance_Up +(3*i))) = up;
        *((Vec3f*)(instance_sc +(3*i))) = sc;
	}

    //CMesh mesh = Solids::Octahedron;
    CMesh mesh = Solids::Icosahedron;
    nVerts = countVerts( mesh.nfaces, mesh.ngons );
    Vec3f * model_vpos = new Vec3f[nVerts];
	Vec3f * model_vnor = new Vec3f[nVerts];
	hardFace( mesh.nfaces, mesh.ngons, mesh.faces, mesh.verts, (GLfloat*)model_vpos, (GLfloat*)model_vnor );

	instances.init( MaxParticles, nVerts, model_vpos, model_vnor, instance_pos, instance_dir, instance_Up, instance_sc );

	delete [] model_vpos;
	delete [] model_vnor;

    ticks_per_second = calibrate_timer(100);

    //lastTime = glfwGetTime();
    lastTime = 0.0;
    qCamera.setOne();
    return 0;
};

void physics(){
    float time = frameCount * 0.005;
    ParticlesCount = MaxParticles;

    float dphi = 0.001;
    double ca  = cos(dphi);
    float  sa  = sin(dphi);

	for( int i=0; i<ParticlesCount-2; i++ ){
        //Vec3f pos, dir, up;
        //pos.set( randf(-span,span), randf(-span,span), randf(-span,span) );
        //dir.set( randf(-1.0f,1.0f), randf(-1.0f,1.0f), randf(-1.0f,1.0f) ); dir.normalize();

        ((Vec3f*)(instance_Up+(3*i)))->rotate_csa( ca, sa, *(Vec3f*)(instance_dir+(3*i)) );
        //up.makeOrthoU( dir ); up.normalize();

	}
}

void draw( ){

    if((frameCount%100)==0){
        long t2     = getCPUticks();
        double lag  = (t2-lastCPUtick)/100.0d;
        printf( "%f [Mtick/frame] %f fps\n", lag*1e-6, ticks_per_second/lag );
        lastCPUtick = t2;
    };

    glClearColor(0.8, 0.8, 0.8, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );
	// Simulate all particles

	physics();

	//glEnable(GL_BLEND);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glDisable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);

    Mat4f camMat,mRot,mPersp;
    qCamera.toMatrix(mouseMat);
    mRot.setOne(); mRot.setRot(mouseMat);
    float fov = 3.0;
    mPersp.setPerspective( fov, fov*ASPECT_RATIO, 1.0, 1000.0 );
    camMat.set_mmul_TN( mRot, mPersp );
    //Mat3f objRot; objRot.setOne();

    //printf("======\n");
    //printf( "qCamera (%g,%g,%g,%g)\n", qCamera.x, qCamera.y, qCamera.z, qCamera.w );
    //printf( "camPos (%g,%g,%g)\n", camPos.x, camPos.y, camPos.z );
    //mouseMat.print();
    //camMat.print();

    Shader       * sh;
    if      ( viewGeom    ){ sh = shGeom; }
    else if ( customDepth ){ sh = shFragedDepth; }else{ sh = shFraged; };
    sh->use();
    sh->set_camPos( (GLfloat*)&camPos );
    sh->set_camMat( (GLfloat*)&camMat );

    uploadArrayBuffer( instances.pose_Up, instances.nInstances*3*sizeof(GLfloat), instance_Up );
    instances.draw( GL_TRIANGLES );

}

void init();
void quit();
void die( char const *msg );
void inputHanding();

int main(int argc, char *argv[]){
    init();
	setup();
    for ( frameCount=1; frameCount<1000000; frameCount++)    {

        draw(); SDL_GL_SwapWindow(window);
 		//if( !STOP ) draw();
		inputHanding();
        SDL_Delay(delay);
    }
    quit();
    return 0;
}

// FUNCTION ======	inputHanding
void inputHanding(){

    //float posstep = 0.1f; if(RayTerrain){ posstep = 2.0f; }
    float step          = 0.5f;
    float keyRotSpeed   = 0.01f;

    const Uint8 *keys = SDL_GetKeyboardState(NULL);

    if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.roll2  (  (float)keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.roll2  ( -(float)keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.dpitch2(  (float)keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.dpitch2( -(float)keyRotSpeed ); }

    if( keys[ SDL_SCANCODE_W  ] ){ camPos.add_mul( mouseMat.c, +step ); }
	if( keys[ SDL_SCANCODE_S  ] ){ camPos.add_mul( mouseMat.c, -step ); }
	if( keys[ SDL_SCANCODE_A  ] ){ camPos.add_mul( mouseMat.a, -step ); }
	if( keys[ SDL_SCANCODE_D  ] ){ camPos.add_mul( mouseMat.a, +step ); }

	if( keys[SDL_SCANCODE_W]||keys[SDL_SCANCODE_S]||keys[SDL_SCANCODE_A]||keys[ SDL_SCANCODE_D  ]  )printf( "camPos (%g,%g,%g)\n", camPos.x, camPos.y, camPos.z );

	SDL_Event event;
	while(SDL_PollEvent(&event)){
		if( event.type == SDL_KEYDOWN ){
            switch( event.key.keysym.sym ){
                case SDLK_ESCAPE: quit(); break;
                case SDLK_g: viewGeom   =!viewGeom;    break;
                case SDLK_f: customDepth=!customDepth; break;
                //case SDLK_KP_PLUS:  terrain_size[0] *=1.1; terrain_size[2] *=1.1; break;
                //case SDLK_KP_MINUS: terrain_size[0] /=1.1; terrain_size[2] /=1.1; break;
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
    float mouseRotSpeed = 0.002;
    if ( buttons & SDL_BUTTON(SDL_BUTTON_RIGHT)) {
        Quat4f q; q.fromTrackball( 0, 0, -dmx*mouseRotSpeed, dmy*mouseRotSpeed ); qCamera.qmul_T( q );
        //qCamera.dyaw2(-dmx*mouseRotSpeed); qCamera.dpitch2(-dmy*mouseRotSpeed);
        //qCamera.dpitch2(-dmy*mouseRotSpeed); qCamera.dyaw2(-dmx*mouseRotSpeed);
        //qCamera.normalize();

        //pitch +=  dmy*mouseRotSpeed;
        //yaw   +=  dmx*mouseRotSpeed;
        qCamera.toMatrix(mouseMat);
        printf("mouseMat:\n");
        mouseMat.print();
    }
}

void init(){
    if (SDL_Init(SDL_INIT_VIDEO) < 0) die( "Unable to initialize SDL" );
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3); // Opengl 3.2
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3); // Opengl 3.2
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    //SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 32);
    window = SDL_CreateWindow("Tutorial2", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WIDTH, HEIGHT, SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN);
    if ( !window ) die("Unable to create window");
    context = SDL_GL_CreateContext( window );
    //SDL_GL_SetSwapInterval(1); // VSync On
    SDL_GL_SetSwapInterval(VSync);

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
