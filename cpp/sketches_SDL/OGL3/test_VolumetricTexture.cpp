

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
#include "GLobjects.h"
#include "GLInstances.h"
#include "IO_utils.h"
#include "Shader.h"

#include "Solids.h"
#include "CMesh.h"

#include "testUtils.h"

// =============== Global variables



// http://www.real-time-volume-graphics.org/?page_id=28
// http://moddb.wikia.com/wiki/OpenGL:Tutorials:3D_Textures

Shader *sh1;
GLuint  tx3D_1;

//GLInstances instances;
//GLBillboards bilboards;
//GLuint texture_1;

//GLMesh   *mshBranches,*mshLeafs;
//GLMesh  *msh1;

//GLuint       texTest;

//GLuint txHi=0,txLo=0,txDist=0;

GLMesh       *glSprite;
//GLBillboards bilboards;
//FrameBuffer  frameBuff1;
//FrameBuffer  frameBuff2;
//FrameBuffer  frameBuff3;


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
Vec3f  camPos = (Vec3f){ 0.0f, 0.0f, -10.0f };

bool bTransparent  = false;

long lastCPUtick = 0;
double ticks_per_second=0;

// =============== Functions

GLMesh* makeQuad3D( Vec2f p0, Vec2f p1, Vec2f u0, Vec2f u1 ){
    GLfloat verts[]  = { p0.x,p0.y,0.0, p1.x,p1.y,0.0, p0.x,p1.y,0.0,  p0.x,p0.y,0.0, p1.x,p1.y,0.0, p1.x,p0.y,0.0 };
    //Vec2f vUVs   [] = { u0.x,u0.y,     u1.x,u1.y,     u0.x,u1.y,       u0.x,u0.y,     u1.x,u1.y,     u1.x,u0.y     };
    GLfloat colors[] = { u0.x,u0.y,0.0, u1.x,u1.y,0.0, u0.x,u1.y,0.0,  u0.x,u0.y,0.0, u1.x,u1.y,0.0, u1.x,u0.y,0.0 };
    GLMesh* glquad = new GLMesh();
    //glquad->init( 6, 0, NULL, verts, NULL, NULL, vUVs );
    glquad->init( 6, 0, NULL, verts, NULL, colors, NULL );
    return glquad;
}

uint8_t* makeImage3D( int nx, int ny, int nz ){
    int i = 0;
    uint8_t * buff = new uint8_t[nx*ny*nz*4];
    for(int iz=0; iz<nz; iz++){
        for(int iy=0; iy<ny; iy++){
            for(int ix=0; ix<nx; ix++){
                int ioff=i*4;
                //buff[ioff+0]=(int)255/(1+(ix-nx/2)*(ix-nx/2)*0.1) + 32*randf(-1,1);
                //buff[ioff+1]=(int)255/(1+(iy-ny/2)*(iy-ny/2)*0.1) + 32*randf(-1,1);
                //buff[ioff+2]=(int)255/(1+(iz-nz/2)*(iz-nz/2)*0.1) + 32*randf(-1,1);

                uint8_t b = ix^iy^iz;
                //uint8_t b = ix*iz + iy*iz;
                buff[ioff+0]=(int)255/(1+(ix-nx/2)*(ix-nx/2)*0.1) + b;
                buff[ioff+1]=(int)255/(1+(iy-ny/2)*(iy-ny/2)*0.1) + b;
                buff[ioff+2]=(int)255/(1+(iz-nz/2)*(iz-nz/2)*0.1) + b;
                buff[ioff+3]=255;
                i++;
            }
        }
    }
    return buff;
}

void camera(float aspect, Mat4f& camMat){
    Mat4f mRot,mProj;
    qCamera.toMatrix(mouseMat);
    mRot.setOne(); mRot.setRot(mouseMat);
    // float fov = 3.0; mPersp.setPerspective( fov, fov*ASPECT_RATIO, 1.0, 1000.0 );
    mProj.setOrthographic( 8.0, 8.0*aspect, -1.0, -1000.0);
    //mProj.setOne();
    camMat.set_mmul_TN( mRot, mProj );
}

int setup(){

    // from here: http://moddb.wikia.com/wiki/OpenGL:Tutorials:3D_Textures
    glEnable(GL_TEXTURE_3D);
    //int NX=16,NY=16,NZ=16;
    int NX=64,NY=64,NZ=16;
    uint8_t* buff = makeImage3D( NX, NY, NZ );

    glGenTextures(1, &tx3D_1);
    glBindTexture(GL_TEXTURE_3D, tx3D_1);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    //glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    //glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_REPEAT);
    glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA8, NX, NY, NZ, 0, GL_RGBA, GL_UNSIGNED_BYTE, buff);
    //glTexImage3D(GL_TEXTURE_3D,0,GL_INTENSITY,NX,NY,NZ,0,GL_LUMINANCE,GL_UNSIGNED_BYTE,buff);

    //glBindTexture(GL_TEXTURE_3D, tx3D_1);
    delete [] buff;


    /*
    // from here:https://www.opengl.org/discussion_boards/showthread.php/175065-Loading-Volume-Data-to-3D-texture
    //assuming that the data at hand is a 256x256x256 unsigned byte data
    int XDIM=256, YDIM=256, ZDIM=256;
    const int size = XDIM*YDIM*ZDIM;
    //bool LoadVolumeFromFile(const char* fileName) {
    //   FILE *pFile = fopen(fileName,"rb");
    //   if(NULL == pFile) {
    //   return false;
    //   }
    //   GLubyte* pVolume=new GLubyte[size];
    //   fread(pVolume,sizeof(GLubyte),size,pFile);
    //   fclose(pFile);

       //load data into a 3D texture
       glGenTextures(1, &amp;textureID);
       glBindTexture(GL_TEXTURE_3D, textureID);

       // set the texture parameters
       glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP);
       glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP);
       glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP);
       glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
       glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
       glTexImage3D(GL_TEXTURE_3D,0,GL_INTENSITY,NX,NY,NZ,0,GL_LUMINANCE,GL_UNSIGNED_BYTE,pVolume);
       delete [] pVolume;
       return true;
    //}
    */

    sh1=new Shader();
    //sh1->init( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/color3D.glslf"   );
    sh1->init( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/cut3DTexture.glslf"   );
    sh1->getDefaultUniformLocation();

    //ticks_per_second = calibrate_timer(100);
    //lastTime = glfwGetTime();
    lastTime = 0.0;
    qCamera.setOne();

    //glAtlas  = makeQuad3D( {0.0f,0.0f}, {32.0f,2.0f}, {0.0f,0.0f}, {1.f,1.0f} );
    //glquad = makeQuad3D( {0.0f,0.0f}, {0.5f,1.0f}, {0.0f,0.0f}, {0.125f,1.0f} );
    glSprite = makeQuad3D( {0.0f,0.0f}, {1.0f,1.0f}, {0.0f,0.0f}, {1.0f,1.0f} );

    printf( "SETUP DONE \n" );

    return 0;
};

void draw( ){
    glBindFramebuffer(GL_FRAMEBUFFER, 0); glViewport(0,0,WIDTH,HEIGHT);
    glClearColor(0.8, 0.8, 0.8, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );
	// Simulate all particles

	Mat4f camMat;  camera( 1.0, camMat);
	//renderPhases(camMat);

    Mat3f modelMat; modelMat.setOne(); modelMat.mul(4.0f);
    //Vec3f modelPos;

    //printf( "DEBUG 1 \n" );

    glEnable    (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable   ( GL_ALPHA_TEST );

    //printf( "DEBUG 2 \n" );

    Shader * sh;
    sh = sh1;
    sh->use();
    sh->set_modelMat( (GLfloat*)&modelMat );
    sh->set_camPos  ( (GLfloat*)&camPos   );
    sh->set_camMat  ( (GLfloat*)&camMat   );

    glActiveTexture(GL_TEXTURE0);
    //glBindTexture(GL_TEXTURE_2D, frameBuff1.texRGB );
    glBindTexture(GL_TEXTURE_3D, tx3D_1);
    glUniform1i(sh->getUloc("texture_1"), 0);
    glUniform3f(sh->getUloc("txOffset"), 0.0,0.0,frameCount*0.01);
    //glUniform1i(sh->getUloc("baseColor"), 0);
    //printf( "DEBUG 3 \n" );


    //msh1->draw(GL_TRIANGLES);

    sh->set_camPos  ( (const GLfloat[]){4.5,0.0,-10.0}   );
    //glBindTexture(GL_TEXTURE_2D, frameBuff1.texRGB );
    glSprite->draw(GL_TRIANGLES);


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
    float step          = 0.1f;
    float keyRotSpeed   = 0.01f;

    const Uint8 *keys = SDL_GetKeyboardState(NULL);

    if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.roll2  (  (float)keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.roll2  ( -(float)keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.dpitch2(  (float)keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.dpitch2( -(float)keyRotSpeed ); }

    if( keys[ SDL_SCANCODE_Q  ] ){ camPos.add_mul( mouseMat.c, +step ); }
	if( keys[ SDL_SCANCODE_E  ] ){ camPos.add_mul( mouseMat.c, -step ); }
    if( keys[ SDL_SCANCODE_W  ] ){ camPos.add_mul( mouseMat.b, +step ); }
	if( keys[ SDL_SCANCODE_S  ] ){ camPos.add_mul( mouseMat.b, -step ); }
	if( keys[ SDL_SCANCODE_A  ] ){ camPos.add_mul( mouseMat.a, -step ); }
	if( keys[ SDL_SCANCODE_D  ] ){ camPos.add_mul( mouseMat.a, +step ); }

	if( keys[SDL_SCANCODE_W]||keys[SDL_SCANCODE_S]||keys[SDL_SCANCODE_A]||keys[ SDL_SCANCODE_D  ]  )printf( "camPos (%g,%g,%g)\n", camPos.x, camPos.y, camPos.z );

	SDL_Event event;
	while(SDL_PollEvent(&event)){
		if( event.type == SDL_KEYDOWN ){
            switch( event.key.keysym.sym ){
                case SDLK_ESCAPE: quit(); break;
                case SDLK_t: bTransparent   =!bTransparent;    break;
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
