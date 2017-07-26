
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

Shader *shBranches,*shLeafs,*shTex;

//GLInstances instances;
//GLBillboards bilboards;
//GLuint texture_1;

GLMesh   *mshBranches,*mshLeafs;

GLuint       texTest;
GLMesh       *glquad;
GLBillboards bilboards;
FrameBuffer  frameBuff1;


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

GLuint makeTestTextureRGBA( int W, int H ){
    double dx = 1.0d/W;
    double dy = 1.0d/H;
    uint32_t * c_img1 = new uint32_t[H*W];
    for( int iy=0; iy<H; iy++ ){
        for( int ix=0; ix<W; ix++ ){
            uint8_t r,g,b,a;
            r=ix; g=ix^iy; b=iy; a=255;
            c_img1[ iy*W + ix ] = (a<<24) | (b<<16) | (g<<8) | (r);
        }
    }
    GLuint texID;
    newTexture2D( texID, W, H, c_img1, GL_RGBA, GL_UNSIGNED_BYTE );
    return texID;
}

int pushCylinderTris( int n, float r1, float r2, Vec3f base, Vec3f tip, Vec3f up, std::vector<Vec3f>& verts, std::vector<Vec3f>* normals ){
	int nvert=0;

	Vec3f dir,left;
	dir.set_sub( tip, base );
	dir.normalize();
    left = dir.getOrtho(up);

    float alfa = 2*M_PI/n;
    Vec2f rot,drot;
    rot .set(1.0f,0.0f);
    drot.set( cos( alfa ), sin( alfa ) );

	Vec3f q; q.set(dir); q.add_mul( up, -(r1-r2) );
	float pnab =  dir.dot( q )/q.norm();
	float pnc  =  sqrt( 1 - pnab*pnab );

	Vec3f op,opn;
	for(int i=0; i<=n; i++ ){
		Vec3f p,pn;
		p .set( rot.x*up.x + rot.y*left.x, rot.x*up.y + rot.y*left.y, rot.x*up.z + rot.y*left.z );
		pn.set( pnab*p.x   + pnc*dir.x   , pnab*p.y   + pnc*dir.y   , pnab*p.z   + pnc*dir.z    );
		if( i>0 ){
            verts  .push_back(base+op*r1); verts  .push_back(tip+p*r2); verts  .push_back(tip+op*r2);
            verts  .push_back(base+op*r1); verts  .push_back(tip+p*r2); verts  .push_back(base+p*r1);
            if(normals){
                normals->push_back(opn); normals->push_back(pn); normals->push_back(opn);
                normals->push_back(opn); normals->push_back(pn); normals->push_back(pn);
            }
		}
		op=p; opn=pn;
        rot.mul_cmplx( drot );
	}
	return nvert;
};

void pushTris_Treestep( int level, Vec3f pos, Vec3f dir, Vec3f up, std::vector<Vec3f>& verts, std::vector<Vec3f>* normals ){
    float l = dir.norm();
    static const float drnd = 0.15;
    dir.add( randf(-drnd,drnd)*l, randf(-drnd,drnd)*l, randf(-drnd,drnd)*l );
    float f = randf(0.5,0.9);
    dir.mul( f );
    up=cross(dir,up); up.normalize();
    Vec3f pos_ = pos + dir;
    pushCylinderTris( 6, 0.1*l, 0.1*l*f, pos, pos_, (Vec3f){0.0f,1.0f,0.0f},verts,normals);
    //pushCylinderTris(verts,normals);
    if( level>0 ){
        level--;
        pushTris_Treestep( level, pos_, dir+up*0.3*l, up, verts, normals );
        pushTris_Treestep( level, pos_, dir-up*0.3*l, up, verts, normals );
    }
}

void tree_step( int level, Vec3f pos, Vec3f dir, std::vector<Vec3f>& branches, std::vector<Vec3f>& leafs ){
    static const float drnd = 0.6;
    //dir.x *= randf(1.0-drnd,1.0);
    //dir.y *= randf(1.0-drnd,1.0);
    //dir.z *= randf(1.0-drnd,1.0);
    float l = dir.norm();
    dir.add( randf(-drnd,drnd)*l, randf(-drnd,drnd)*l, randf(-drnd,drnd)*l );
    dir.mul( randf(0.5,0.9) );
    Vec3f pos_ = pos + dir;
    branches.push_back(pos );
    branches.push_back(pos_);
    if( level==0 ){
        leafs.push_back(pos_);
    }else{
        level--;
        tree_step( level, pos_, dir,  branches, leafs );
        tree_step( level, pos_, dir,  branches, leafs );
    }
}

void camera(Mat4f& camMat){
    Mat4f mRot,mPersp;
    qCamera.toMatrix(mouseMat);
    mRot.setOne(); mRot.setRot(mouseMat);
    float fov = 3.0;
    mPersp.setPerspective( fov, fov*ASPECT_RATIO, 1.0, 1000.0 );
    camMat.set_mmul_TN( mRot, mPersp );
}

void renderPhases(){

    Mat4f camMat;  camera(camMat);

    int nPhases = 8;
    Mat3f modelMat; modelMat.setOne();
    Vec3f modelPos; modelPos.set(-7.0f,0.0f,10.0f);
    for(int iph=0; iph<nPhases; iph++){

        modelMat.rotate( 2*M_PI/16, (Vec3f){0.0f,1.0f,0.0f});
        modelPos.add( 1.5f, 0.0f, 0.0f );

        Shader * sh;
        sh = shBranches;
        sh->use();
        sh->set_camPos  ( (GLfloat*)&camPos );
        sh->set_camMat  ( (GLfloat*)&camMat );
        sh->set_modelPos( (GLfloat*)&modelPos );
        sh->set_modelMat( (GLfloat*)&modelMat );

        glEnable(GL_DEPTH_TEST);
        mshBranches->draw(GL_TRIANGLES);
    }
}



int setup(){

    shBranches=new Shader();
    shBranches->init( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/color3D.glslf"   );
    shBranches->getDefaultUniformLocation();

    shLeafs=new Shader();
    shLeafs->init( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/pointSprite.glslf"   );
    shLeafs->getDefaultUniformLocation();

    std::vector<Vec3f> branches;
    std::vector<Vec3f> leafs;
    tree_step( 5, (Vec3f){0.0f,0.0f,0.0f}, (Vec3f){0.0f,1.0f,0.0f}, branches, leafs );

    std::vector<Vec3f> verts;
    std::vector<Vec3f> normals;
    pushTris_Treestep( 5, (Vec3f){0.0f,0.0f,0.0f}, (Vec3f){0.0f,1.0f,0.0f}, (Vec3f){1.0f,0.0f,0.0f}, verts, &normals );

    //for(int i=0; i<branches.size(); i++){ printf("%i (%g,%g,%g)\n", i, branches[i].x, branches[i].y, branches[i].z ); }

    mshBranches = new GLMesh();
    //mshBranches->init( branches.size(), 0, NULL, &branches[0],  NULL, NULL, NULL );
    mshBranches->init( verts.size(), 0, NULL, &verts[0],  &normals[0], NULL, NULL );

    mshLeafs = new GLMesh();
    mshLeafs->init( leafs.size(), 0, NULL, &leafs[0],  NULL, NULL, NULL );

    //ticks_per_second = calibrate_timer(100);
    //lastTime = glfwGetTime();
    lastTime = 0.0;
    qCamera.setOne();

    // ==== BEGIN : RENDER TO TEXTURE

    shTex=new Shader();
    //shTex->init( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/color3D.glslf"   );
    //shTex->init( "common_resources/shaders/texture3D.glslv",   "common_resources/shaders/texture.glslf"   );
    //shTex->init( "common_resources/shaders/texture3D_anim.glslv",   "common_resources/shaders/texture.glslf"   );
    shTex->init( "common_resources/shaders/texture3D_anim.glslv",   "common_resources/shaders/texture_anim.glslf"   );
    shTex->getDefaultUniformLocation();

    glquad =new GLMesh();
    glquad->init( 6, 0, NULL, DEFAULT_Bilboard_verts, NULL, NULL, DEFAULT_Bilboard_UVs );
    //glquad->init( 6, 0, NULL, DEFAULT_Bilboard_verts, DEFAULT_Bilboard_verts, NULL, NULL );

    texTest = makeTestTextureRGBA( 256, 256);

    frameBuff1.init( 800, 800 );

    // ===== Prepare texture by rendering

    glBindFramebuffer(GL_FRAMEBUFFER, frameBuff1.buff );

    //glClearColor(0.0, 0.0, 0.8, 1.0);
    glClearColor(0.9, 0.9, 0.9, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );

    renderPhases();

    // ==== END   : RENDER TO TEXTURE

    return 0;
};

void draw( ){
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glClearColor(0.8, 0.8, 0.8, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );
	// Simulate all particles

	renderPhases();

    Mat3f modelMat; modelMat.setOne(); modelMat.mul(5.0f);
    Vec3f modelPos; modelPos.set(0.0f,0.0f,10.0f);
    Mat4f camMat; camera(camMat);

    //glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, texTest );
    glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, frameBuff1.texRGB );

    Shader * sh;
    sh = shTex;
    sh->use();
    sh->set_modelPos( (GLfloat*)&modelPos );
    sh->set_modelMat( (GLfloat*)&modelMat );
    sh->set_camPos  ( (GLfloat*)&camPos   );
    sh->set_camMat  ( (GLfloat*)&camMat   );

    glUniform1f(sh->getUloc("angle0"), frameCount*0.01 );
    glUniform1i(sh->getUloc("texture_1"), 0);

    glquad->draw(GL_TRIANGLES);
    glquad->drawPoints(10.0);

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
