
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
#include "GLobjects.h"
#include "GLInstances.h"
#include "IO_utils.h"
#include "Shader.h"

#include "Solids.h"
#include "CMesh.h"

#include "testUtils.h"

// =============== Global variables

// http://wiki.polycount.com/wiki/Transparency_map
// http://www.valvesoftware.com/publications/2007/SIGGRAPH2007_AlphaTestedMagnification.pdf
// /home/prokop/Dropbox/MyDevSW/Python/ImageProcessing/DistanceField.py

Shader *shModel,*shTex,*shGenSDFTex,*shViewSDFTex;

//GLInstances instances;
//GLBillboards bilboards;
//GLuint texture_1;

//GLMesh   *mshBranches,*mshLeafs;
GLMesh  *msh1;

//GLuint       texTest;

GLuint txHi=0,txLo=0,txDist=0;

GLMesh       *glSprite,*glAtlas;
GLBillboards bilboards;
FrameBuffer  frameBuff1;
FrameBuffer  frameBuff2;
FrameBuffer  frameBuff3;


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

/*
GLMesh* makeQuad3D( Vec2f p0, Vec2f p1, Vec2f u0, Vec2f u1 ){
    GLfloat verts[] = { p0.x,p0.y,0.0, p1.x,p1.y,0.0, p0.x,p1.y,0.0,   p0.x,p0.y,0.0, p1.x,p1.y,0.0, p1.x,p0.y,0.0 };
    Vec2f vUVs   [] = { u0.x,u0.y,     u1.x,u1.y,     u0.x,u1.y,       u0.x,u0.y,     u1.x,u1.y,     u1.x,u0.y     };
    GLMesh* glquad = new GLMesh();
    glquad->init( 6, 0, NULL, verts, NULL, NULL, vUVs );
    return glquad;
}
*/

float* makeImage( int n ){
    float * out = new float[n*n];
    float dx = 1.0/n;
    float dy = 1.0/n;

    const int np = 3;
    float ps[] ={
       0.5, 0.5,   0.4 , 0.4,  1.0,
       0.2,0.2,  0.1 , 0.1,  1.0,
       0.75,0.75,  0.1 , 0.5,  1.0,
    };

    for(int iy=0; iy<n; iy++){
        for(int ix=0; ix<n; ix++){
            int i = iy*n+ix;
            bool b = true;
            for( int ip=0; ip<np; ip++ ){
                int iip  = ip*5;
                float x  = (ix-n*ps[iip+0])*dx/ps[iip+2];
                float y  = (iy-n*ps[iip+1])*dy/ps[iip+3];
                bool bi  = (x*x + y*y) < 1.0;
                if( ps[iip+4]>0 ){ b&=!bi; }else{ b|=bi; }
            }
            if(b){ out[i] = 1.0; }else{  out[i] = 0.0; }
            /*
            if(r2 > 0.1){ out[i] = 1.0;
            }else{
                out[i] = 0.0;
            }
            */
        }
    }
    return out;
}

float* img2distanceField( int n, int m, float* in ){
    int l = n/m;
    float* out = new float[l*l];
    float dnorm2 = 2.0/(m*m);
    float vmax   = sqrt(0.5);
    for(int iy=0; iy<l; iy++){
        for(int ix=0; ix<l; ix++){
            double dOn=vmax,dOff=vmax;
            for(int jy=0; jy<m; jy++){
                for(int jx=0; jx<m; jx++){
                    float d = (sq(jx-m*0.5f) + sq(jy-m*0.5f))*dnorm2;
                    int iin = (iy*m+jy)*n + (ix*m+jx);
                    if( in[iin] > 0.5 ){ dOff=fmin(dOff,d); }else{ dOn=fmin(dOn,d); }
                } // jx
            } // jy
            out[ iy*l+ix ] = (  ( sqrt(dOn) - sqrt(dOff) )*0.5 + 0.5 );
        } // ix
    } // iy
    return out;
}

GLuint img2texture(  int W, int H, float* in){
    //uint32_t * c_img1 = new uint32_t[H*W];
    uint8_t * cimg = new uint8_t[H*W];
    for( int iy=0; iy<H; iy++ ){
        for( int ix=0; ix<W; ix++ ){
            //c_img1[ iy*W + ix ] = (a<<24) | (b<<16) | (g<<8) | (r);
            int i = iy*W + ix;
            cimg[i] = (int)(255.0*in[i]);
        }
    }
    GLuint texID;
    //newTexture2D( texID, W, H, cimg, GL_RGBA, GL_UNSIGNED_BYTE );
    newTexture2D( texID, W, H, cimg, GL_RED, GL_UNSIGNED_BYTE );
    delete [] cimg;
    return texID;
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

    //ticks_per_second = calibrate_timer(100);
    //lastTime = glfwGetTime();
    lastTime = 0.0;
    qCamera.setOne();

    // ==== BEGIN : RENDER TO TEXTURE
    shModel=new Shader();
    shModel->init( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/color3D.glslf"   );
    shModel->getDefaultUniformLocation();

    shTex=new Shader();
    shTex->init( "common_resources/shaders/texture3D.glslv",   "common_resources/shaders/texture.glslf"   );
    shTex->getDefaultUniformLocation();

    shViewSDFTex=new Shader();
    shViewSDFTex->init( "common_resources/shaders/texture3D.glslv",   "common_resources/shaders/texSDF.glslf"   );
    shViewSDFTex->getDefaultUniformLocation();

    shGenSDFTex=new Shader();
    shGenSDFTex->init( "common_resources/shaders/texture3D.glslv",   "common_resources/shaders/genTexSDF.glslf"   );
    shGenSDFTex->getDefaultUniformLocation();

    //glAtlas  = makeQuad3D( {0.0f,0.0f}, {32.0f,2.0f}, {0.0f,0.0f}, {1.f,1.0f} );
    //glquad = makeQuad3D( {0.0f,0.0f}, {0.5f,1.0f}, {0.0f,0.0f}, {0.125f,1.0f} );
    glSprite = makeQuad3D( {0.0f,0.0f}, {1.0f,1.0f}, {0.0f,0.0f}, {1.0f,1.0f} );

    //glquad =new GLMesh();
    //glquad->init( 6, 0, NULL, DEFAULT_Bilboard_verts, NULL, NULL, DEFAULT_Bilboard_UVs );
    //glquad->init( 6, 0, NULL, DEFAULT_Bilboard_verts, DEFAULT_Bilboard_verts, NULL, NULL );

    //texTest = makeTestTextureRGBA( 256, 256);


    int nHi = 1024;

    float* fHi = makeImage( nHi );
    txHi = img2texture(  nHi, nHi, fHi);

    int nLo = 16;
    //int nLo = 32;
    float* fLo = makeImage( nLo );
    txLo = img2texture(  nLo, nLo, fLo);

    int m   = nHi/nLo;
    float* fDist = img2distanceField( nHi, m, fHi );
    txDist  = img2texture(  nLo, nLo, fDist);

    delete [] fHi; delete [] fLo; delete [] fDist;



    //frameBuff1.init( 2048, 64 );
    nHi = 1024;
    frameBuff1.init( nHi, nHi );
    // ===== Prepare texture by rendering
    //glBindFramebuffer(GL_FRAMEBUFFER, frameBuff1.buff );
    frameBuff1.bind();

    //glClearColor(0.0, 0.0, 0.8, 1.0);
    glClearColor(0.9, 0.9, 0.9, 0.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );
    Mat4f camMat; camera( 1/16.0, camMat);
    //renderPhases( camMat);

    std::vector<Vec3f> verts;
    std::vector<Vec3f> normals;
    //pushTris_Treestep( 5, (Vec3f){0.0f,0.0f,0.0f}, (Vec3f){0.0f,1.0f,0.0f}, (Vec3f){1.0f,0.0f,0.0f}, verts, &normals );
    //for(int i=0; i<branches.size(); i++){ printf("%i (%g,%g,%g)\n", i, branches[i].x, branches[i].y, branches[i].z ); }

    float sz=1.0;
    for(int i=0; i<10; i++){
        verts.push_back({randf(-sz,sz),randf(-sz,sz),0.0}); normals.push_back( {0.0,0.0,1.0} );
        verts.push_back({randf(-sz,sz),randf(-sz,sz),0.0}); normals.push_back( {0.0,1.0,0.0} );
        verts.push_back({randf(-sz,sz),randf(-sz,sz),0.0}); normals.push_back( {1.0,0.0,0.0} );
    }
    msh1 = new GLMesh();
    msh1->init( verts.size(), 0, NULL, &verts[0],  &normals[0], NULL, NULL );
    msh1->draw(GL_TRIANGLES);

    // downsample
    nLo = 32;
    frameBuff2.init( nLo, nLo );
    frameBuff2.bind();

    //glActiveTexture(GL_TEXTURE0);
    //glBindTexture(GL_TEXTURE_2D, frameBuff1.texRGB );
    //glUniform1i(sh->getUloc("texture_1"), 0);
    //glUniform1i(sh->getUloc("baseColor"), 0);
    //printf( "DEBUG 3 \n" );

    //msh1->draw(GL_TRIANGLES);

    //Mat4f camMat;
    camera( 1.0, camMat);
    Mat3f modelMat; modelMat.setOne(); modelMat.mul(16.0f);
    glEnable    (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable   (GL_ALPHA_TEST );
    Shader * sh;
    sh = shGenSDFTex;
    sh->use();
    sh->set_modelMat( (GLfloat*)&modelMat );
    sh->set_camPos  ( (const float[]){+8.0,+8.0,-10.0} );
    sh->set_camMat  ( (GLfloat*)&camMat   );

    int msub = 16;
    sh->setUniformi("msub",  msub            );
    //sh->setUniformf("lpix",  1.4*(1.0f/nHi)/(nHi/((float)(nLo*msub)))  );
    //sh->setUniformf("lpix",  1.4*(1.0f/nHi)*(((float)(nLo*msub))/nHi)  );
    //sh->setUniformf("lpix",  1.4*(((float)(nLo*msub))/(nHi*nHi))  );
    sh->setUniformf("lpix",  1.4/(((float)(nLo*msub)))  );

    glActiveTexture(GL_TEXTURE0);
    glUniform1i(sh->getUloc("texture_1"), 0);
    glBindTexture(GL_TEXTURE_2D, frameBuff1.texRGB );   glSprite->draw(GL_TRIANGLES);

    frameBuff3.init( 32, 32 );
    frameBuff3.bind();

    sh = shTex;
    sh->use();
    sh->set_modelMat( (GLfloat*)&modelMat );
    sh->set_camPos  ( (const float[]){+8.0,+8.0,-10.0} );
    sh->set_camMat  ( (GLfloat*)&camMat   );

    glBindTexture(GL_TEXTURE_2D, frameBuff1.texRGB );   glSprite->draw(GL_TRIANGLES);

    // ==== END   : RENDER TO TEXTURE

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
    sh = shViewSDFTex;
    sh->use();
    sh->set_modelMat( (GLfloat*)&modelMat );
    sh->set_camPos  ( (GLfloat*)&camPos   );
    sh->set_camMat  ( (GLfloat*)&camMat   );

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, frameBuff1.texRGB );
    glUniform1i(sh->getUloc("texture_1"), 0);
    //glUniform1i(sh->getUloc("baseColor"), 0);
    //printf( "DEBUG 3 \n" );


    //msh1->draw(GL_TRIANGLES);

    sh->set_camPos  ( (const GLfloat[]){4.5,0.0,-10.0}   );
    glBindTexture(GL_TEXTURE_2D, frameBuff1.texRGB );   glSprite->draw(GL_TRIANGLES);

    sh->set_camPos  ( (const GLfloat[]){0.0,0.0,-10.0}   );
    glBindTexture(GL_TEXTURE_2D, frameBuff2.texRGB );   glSprite->draw(GL_TRIANGLES);

    sh->set_camPos  ( (const GLfloat[]){4.5,4.5,-10.0}   );
    glBindTexture(GL_TEXTURE_2D, frameBuff3.texRGB );   glSprite->draw(GL_TRIANGLES);

    /*
    camPos = {-3.0,0.0,-10.0};
    sh->set_camPos  ( (GLfloat*)&camPos   );
    glBindTexture(GL_TEXTURE_2D, txHi );   glSprite->draw(GL_TRIANGLES);

    camPos.x += 4.1;
    sh->set_camPos  ( (GLfloat*)&camPos   );
    glBindTexture(GL_TEXTURE_2D, txLo );   glSprite->draw(GL_TRIANGLES);

    camPos.x += 4.1;
    sh->set_camPos  ( (GLfloat*)&camPos   );
    glBindTexture(GL_TEXTURE_2D, txDist );   glSprite->draw(GL_TRIANGLES);
    */

    //glSprite->drawPoints(10.0);
    //printf( "DEBUG 4 \n" );



    /*
    glUniform1i(sh->getUloc("texture_1"), 0);
    glUniform1f(sh->getUloc("nPhases"), 32.0f);

    glUniform1f(sh->getUloc("angle0"), 0 );
    sh->set_modelPos( (const GLfloat[]){-4.0f,0.5f,10.0f} );
    glAtlas->draw(GL_TRIANGLES);
    glAtlas->drawPoints(10.0);

    glUniform1f(sh->getUloc("angle0"), frameCount*0.005 );
    sh->set_modelPos( (const GLfloat[]){0.0f,2.5f,10.0f} );
    glSprite->draw(GL_TRIANGLES);
    glSprite->drawPoints(10.0);
    */

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
