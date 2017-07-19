
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

Shader   *shader1, *shaderParticle;
GLObject *object1,*obj_terrain,*obj_flicker;

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


const int MaxParticles  = 16;
GLfloat g_particule_position_size_data [MaxParticles * 4];
GLubyte g_particule_color_data         [MaxParticles * 4];

GLuint particles_color_buffer;
GLuint particles_position_buffer;
GLuint billboard_vertex_buffer;

GLuint CameraRight_worldspace_ID, CameraUp_worldspace_ID, ViewProjMatrixID;
//GLuint TextureID                  = glGetUniformLocation(programID, "myTextureSampler");

// ===============================================
// ======================= Functions
// ===============================================

/*
Vec3d terrainFunc( Vec2d p ){ return (Vec3d){p.x*10.0,sin(p.x)*sin(p.y*0.5)*10.0,p.y*10.0}; };
*/

double arr_func( int n, const double * xs ){ };
struct S { int a, b, c, d, e; };

GLObject * makeOgl_flat( const CMesh& mesh ){
    GLObject * ogl = new GLObject();
    ogl->setup( countVerts( mesh.nfaces, mesh.ngons ) );
    hardFace( mesh.nfaces, mesh.ngons, mesh.faces, mesh.verts, ogl->buffs[0].cbuff, ogl->buffs[1].cbuff );
    ogl->init();
    return ogl;
}

void setup(){

    arr_func( 3, (const double[]){1.0,2.0,3.0} );
    //struct S s = { .c = 3, .d=4.0 }; // works only in C99 not in C++11

    shader1=new Shader();

    shader1->init( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/color3D.glslf"   );
    //shader1->init( "common_resources/shaders/color3D_depth.glslv",   "common_resources/shaders/color3D_depth.glslf"   );
    //shader1->init( "common_resources/shaders/Particle.glslv",   "common_resources/shaders/Particle.glslf"   );
    //shader1->init( "common_resources/shaders/const3D.glslv",   "common_resources/shaders/const3D.glslf"   );
    //shader1->init( "common_resources/shaders/const3D.glslv",   "common_resources/shaders/pointSprite.glslf"   );
    //shader1->init( "common_resources/shaders/pos3D.glslv",   "common_resources/shaders/pos3D.glslf"   );
    shader1->getDefaultUniformLocation();

    shaderParticle=new Shader();
    shaderParticle->init( "common_resources/shaders/Particle.glslv",   "common_resources/shaders/Particle.glslf"   );
    //shaderParticle->getDefaultUniformLocation();

    //object1 = makeOgl_flat( Solids::Tetrahedron );
    object1 = makeOgl_flat( Solids::Octahedron );
    //object1 = makeOgl_flat( Solids::Icosahedron );

    ninstancs = 100; // 30 ms/frame
    instance_points = new Vec3f[ninstancs];
    //for (int i=0; i<ninstancs; i++){ instance_points[i] = (Vec3f){randf(-15.0,15.0),randf(-15.0,15.0),randf(5.0,100.0)};};
    for (int i=0; i<ninstancs; i++){ instance_points[i] = (Vec3f){0.0,0.0,5.0*i};};



	qCamera.setOne();
	delay = 1;


	// ============ Particles

		// Vertex shader
	GLuint CameraRight_worldspace_ID  = glGetUniformLocation(shaderParticle->shaderprogram, "CameraRight_worldspace");
	GLuint CameraUp_worldspace_ID     = glGetUniformLocation(shaderParticle->shaderprogram, "CameraUp_worldspace");
	GLuint ViewProjMatrixID           = glGetUniformLocation(shaderParticle->shaderprogram, "VP");
	//GLuint TextureID                  = glGetUniformLocation(programID, "myTextureSampler");

    // The VBO containing the 4 vertices of the particles.
	// Thanks to instancing, they will be shared by all particles.
	static const GLfloat g_vertex_buffer_data[] = {
		 -0.5f, -0.5f, 0.0f,
		  0.5f, -0.5f, 0.0f,
		 -0.5f,  0.5f, 0.0f,
		  0.5f,  0.5f, 0.0f,
	};
	glGenBuffers(1, &billboard_vertex_buffer);
	glBindBuffer(GL_ARRAY_BUFFER, billboard_vertex_buffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data), g_vertex_buffer_data, GL_STATIC_DRAW);


    // The VBO containing the positions and sizes of the particles
	glGenBuffers(1, &particles_position_buffer);
	glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
	// Initialize with empty (NULL) buffer : it will be updated later, each frame.
	glBufferData(GL_ARRAY_BUFFER, MaxParticles * 4 * sizeof(GLfloat), NULL, GL_STREAM_DRAW);

	// The VBO containing the colors of the particles
	glGenBuffers(1, &particles_color_buffer);
	glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
	// Initialize with empty (NULL) buffer : it will be updated later, each frame.
	glBufferData(GL_ARRAY_BUFFER, MaxParticles * 4 * sizeof(GLubyte), NULL, GL_STREAM_DRAW);

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

    mRot.setOne(); mRot.setRot(mouseMat);


    mPersp.setPerspective( fov, fov*ASPECT_RATIO, -1.0, -1000.0 );
    //camMat.set_mmul_TN( m1, mPersp );
    //camMat.set_mmul( m1, mPersp );
    //mRot = mRot.transposed( );
    camMat.set_mmul_TN( mRot, mPersp );
    //camMat.set_mmul( mPersp, mRot );
    Mat3f objRot; objRot.setOne();

    // ============= Objects

    shader1->use();
    shader1->set_camPos  ( (float*)&camPos );
    shader1->set_camMat  ( (float*)&camMat );
    shader1->set_modelMat( (float*)&objRot );


    Vec3f  p;
    Quat4f c;


    //object1->draw();
    object1->preDraw();
    for(int i=0; i<ninstancs; i++){
        p = instance_points[i]; p.y+=5.0;
        shader1->set_modelPos ( (GLfloat*)&p );
        c.set( i*0.01f, 0.5f,1-i*0.01f, 1.0f  );
        //c.set( 1.0f, 0.0f, 1.0f, 1.0f );
        shader1->set_baseColor( (GLfloat*)&c );
        object1->draw_instance();
    }
    object1->afterDraw();


    // =========== Particles

    if( true ){ // particle on

    int ParticlesCount = MaxParticles;

    for(int ip=0; ip<ParticlesCount; ip++){

        float t = frameCount * 0.00154f;
        g_particule_position_size_data[4*ParticlesCount+0] = (float)(20.0*sin(ip+t));
        g_particule_position_size_data[4*ParticlesCount+1] = (float)(20.0*sin(ip+t*2.0));
        g_particule_position_size_data[4*ParticlesCount+2] = (float)(20.0*sin(ip+t*3.0));

        g_particule_position_size_data[4*ParticlesCount+3] = (float)(1.0*sin(ip+t + 1215.0));

        g_particule_color_data[4*ParticlesCount+0] = (float)(sin(ip+t*15.0));
        g_particule_color_data[4*ParticlesCount+1] = (float)(sin(ip+t*13.0));
        g_particule_color_data[4*ParticlesCount+2] = (float)(sin(ip+t*15.0+0.5));
        g_particule_color_data[4*ParticlesCount+3] = (float)(sin(ip+t*10.0));

    }


    glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
    glBufferData(GL_ARRAY_BUFFER, MaxParticles * 4 * sizeof(GLfloat), NULL, GL_STREAM_DRAW); // Buffer orphaning, a common way to improve streaming perf. See above link for details.
    glBufferSubData(GL_ARRAY_BUFFER, 0, ParticlesCount * sizeof(GLfloat) * 4, g_particule_position_size_data);

    glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
    glBufferData(GL_ARRAY_BUFFER,  MaxParticles * 4 * sizeof(GLubyte), NULL, GL_STREAM_DRAW); // Buffer orphaning, a common way to improve streaming perf. See above link for details.
    glBufferSubData(GL_ARRAY_BUFFER, 0, ParticlesCount * sizeof(GLubyte) * 4, g_particule_color_data);

    shaderParticle->use();

    //glEnable(GL_BLEND);
    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Same as the billboards tutorial
    //glUniform3f(CameraRight_worldspace_ID, ViewMatrix[0][0], ViewMatrix[1][0], ViewMatrix[2][0] );
    //glUniform3f(CameraUp_worldspace_ID   , ViewMatrix[0][1], ViewMatrix[1][1], ViewMatrix[2][1] );
    //glUniformMatrix4fv(ViewProjMatrixID, 1, GL_FALSE, &ViewProjectionMatrix[0][0]);

    glUniform3f(CameraRight_worldspace_ID, mouseMat.a.x, mouseMat.b.x, mouseMat.c.x );
    glUniform3f(CameraUp_worldspace_ID   , mouseMat.a.y, mouseMat.b.y, mouseMat.c.y );
    glUniformMatrix4fv(ViewProjMatrixID, 1, GL_FALSE, (GLfloat*)&mPersp );

    // 1rst attribute buffer : vertices
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, billboard_vertex_buffer);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0 );

    // 2nd attribute buffer : positions of particles' centers
    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
    glVertexAttribPointer(1,4,GL_FLOAT,GL_FALSE,0,(void*)0);

    // 3rd attribute buffer : particles' colors
    glEnableVertexAttribArray(2);
    glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
    glVertexAttribPointer(2, 4,GL_UNSIGNED_BYTE,GL_TRUE,0,(void*)0);


    glVertexAttribDivisor(0, 0); // particles vertices : always reuse the same 4 vertices -> 0
    glVertexAttribDivisor(1, 1); // positions : one per quad (its center)                 -> 1
    glVertexAttribDivisor(2, 1); // color : one per quad                                  -> 1
    glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, ParticlesCount);

    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);

    } // particle on

    SDL_GL_SwapWindow(window);
}

// =============================================================
// ========== BORING DEFAULT RUTINES
// =============================================================

// FUNCTION ======	inputHanding
void inputHanding(){

    float posstep = 0.1f; if(RayTerrain){ posstep = 2.0f; }
    float step          = 0.1f;
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

    //if( keys[ SDL_SCANCODE_W  ] ){ camPos.z +=step ; }
	//if( keys[ SDL_SCANCODE_S  ] ){ camPos.z -=step ; }
	//if( keys[ SDL_SCANCODE_A  ] ){ camPos.x -=step ; }
	//if( keys[ SDL_SCANCODE_D  ] ){ camPos.x +=step ; }

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
    float mouseRotSpeed = 0.002;
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
    //SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 32);
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
