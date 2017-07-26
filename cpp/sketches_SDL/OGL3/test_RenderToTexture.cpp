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
#include "GLInstances.h"

#include "Mesh.h"

// ============= GLOBAL VARIABLES

const int WIDTH  = 800;
const int HEIGHT = 800;
SDL_Window * window     = NULL;
SDL_GLContext   context = NULL;
int frameCount = 0;
bool STOP = false;

GLuint    vao;
Shader   *shader_pre,*shader_post,*shTex;
GLObject *object1;

GLuint       texTest;
GLMesh       *glmesh,*glquad;
GLBillboards bilboards;
FrameBuffer  frameBuff1;

int mouseX, mouseY;
Quat4f qCamera;
Mat3f  mouseMat;
Vec3f  camPos   = (Vec3f){ 0.0f,0.0f,  10.0f };
Vec3f  modelPos = (Vec3f){ 0.0f,0.0f,  30.0f };

bool preRender = true;

// ============= FUNCTIONS

void quit();
void die ( char const *msg );
void inputHanding ();
void init();
void draw();
void loop( int niters );

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

void setup(){

    Mesh mesh;
    mesh.fromFileOBJ( "common_resources/turret.obj" );
    mesh.polygonsToTriangles(false);
    mesh.tris2normals(true);
    mesh.findEdges( );

    glmesh = new GLMesh();
    glmesh->init_d( mesh.points.size(), mesh.triangles.size()*3, ((int*)&mesh.triangles[0]), (double*)&(mesh.points [0]), (double*)&(mesh.normals[0]), NULL, NULL );

    shader_pre=new Shader();
	//shader_pre->init( "common_resources/shaders/shade3D.glslv",   "common_resources/shaders/shade3D.glslf"   );
	shader_pre->init( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/color3D.glslf"   );
	shader_pre->getDefaultUniformLocation();

    shader_pre->use();
    glUniform3fv( shader_pre->getUloc("light_pos"    ), 1, (const float[]){ 1.0f,  1.0f, -1.0f } );
    glUniform3fv( shader_pre->getUloc("lightColor"   ), 1, (const float[]){ 1.0f,  0.9f,  0.8f } );
    glUniform3fv( shader_pre->getUloc("diffuseColor" ), 1, (const float[]){ 1.0f,  1.0f,  1.0f } );
    glUniform3fv( shader_pre->getUloc("ambientColor" ), 1, (const float[]){ 0.2f,  0.2f,  0.3f } );
    glUniform3fv( shader_pre->getUloc("specularColor"), 1, (const float[]){ 0.0f,  0.0f,  0.0f } );


    shTex=new Shader();
    //shTex->init( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/color3D.glslf"   );
    shTex->init( "common_resources/shaders/texture3D.glslv",   "common_resources/shaders/texture.glslf"   );
    shTex->getDefaultUniformLocation();


    glquad =new GLMesh();
    glquad->init( 6, 0, NULL, DEFAULT_Bilboard_verts, NULL, NULL, DEFAULT_Bilboard_UVs );
    //glquad->init( 6, 0, NULL, DEFAULT_Bilboard_verts, DEFAULT_Bilboard_verts, NULL, NULL );

    texTest = makeTestTextureRGBA( 256, 256);


    qCamera.setOne();

    printf("DEBUG 1 \n");

    // ------------- object

    /*
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
	*/

	/*
	int MaxParticles = 1600;
    Vec3f * instance_pos = new Vec3f[MaxParticles];
	Vec2f * instance_sc  = new Vec2f[MaxParticles];
	float span     = 20.0f;
    float time     = frameCount * 0.005;
    int ParticlesCount = MaxParticles;
    int nside = int( pow( ParticlesCount, 1.0/3.0) );
	for( int i=0; i<ParticlesCount; i++ ){
        Vec3f pos;
        Vec2f sc;
        // pos.set( randf(-span,span), randf(-span,span), randf(-span,span) );
        // sc.set( randf(0.5,1.5),randf(0.5,1.5),randf(0.5,1.5) );
        pos.set( i%nside , (i/nside)%nside , (i/(nside*nside)) ); pos.mul(1.2);
        //sc.set( 1.0,1.0,1.0 );
        //float sz = randf(0.5,1.5); sc.set( sz,sz );
        sc.set( randf(0.5,1.5),randf(0.5,1.5) );

        instance_pos[i] = pos;
        instance_sc[i]  = sc;
	}
	bilboards.init( MaxParticles, 6, DEFAULT_Bilboard_UVs, instance_pos, instance_sc );
    delete [] instance_pos;
    delete [] instance_sc;
    */

    printf("DEBUG 2 \n");


    // ------------- shader for blitting from texture
    /*
    shader_post=new Shader();
    //shader_post->init( "common_resources/shaders/const3D.glslv", "common_resources/shaders/hardSprite.glslf" );
    shader_post->init( "common_resources/shaders/Bilboard3D.glslv", "common_resources/shaders/hardSprite.glslf" );
    //shader_post->init( "common_resources/shaders/Bilboard3D.glslv", "common_resources/shaders/texture.glslf" );
    //shader_post->init( "common_resources/shaders/Bilboard3D.glslv", "common_resources/shaders/const3D.glslf" );
    shader_post->getDefaultUniformLocation();

    //shader_post=new Shader();
	//shader_post->init( "shaders/plain_vert.c", "shaders/SSAO_frag.c" );
	//shader_post->use();
    //glUniform2fv( shader_post->getUloc( "resolution" ), 1, (const float[]){WIDTH,HEIGHT} );

    printf("DEBUG 3 \n");
    */

    // ---- prepare FrameBuffer
    frameBuff1.init( WIDTH, HEIGHT );

    printf("DEBUG 4 \n");


    printf( "========== SETUP DONE ! \n");
}

void draw(){
    glClearColor(0.8, 0.8, 0.8, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );

    // ------- render to texture

    Mat4f camMat;

    if( preRender ){
        qCamera.toMatrix(mouseMat);
        camMat;   camMat.setPerspective( 20.0, 20.0, 2.0, 1000.0 );

        glBindFramebuffer(GL_FRAMEBUFFER, frameBuff1.buff );
        //glBindFramebuffer(GL_FRAMEBUFFER, 0);

        glEnable( GL_DEPTH_TEST ); //glDepthFunc( GL_LESS );

        shader_pre->use();
        shader_pre->set_modelPos( (GLfloat*)&modelPos );
        shader_pre->set_modelMat( (GLfloat*)&mouseMat );
        //shader1->set_camPos  ( (GLfloat*)&camPos );
        shader_pre->set_camMat  ( (GLfloat*)&camMat );

        glmesh->draw();
        //glmesh->drawPoints(10.0);

        //mouseMat.print();
        //camMat.print();

        glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, frameBuff1.texRGB );
    } else {
        glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, texTest );
    };

    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    //camMat.setPerspective( 20.0, 20.0, 1.0, 1000.0 );
    //camMat.setPerspective( 30.0, 30.0, 1.0, 1000.0 );
    camMat.setPerspective( 40.0, 40.0, 2.0, 1000.0 );
    //camMat.setOne();
    qCamera.toMatrix(mouseMat);

    Shader *sh = shTex;
    sh->use();
    sh->set_modelPos( (GLfloat*)&modelPos );
    sh->set_modelMat( (GLfloat*)&mouseMat );
    //sh->set_camPos  ( (GLfloat*)&camPos   );
    sh->set_camMat  ( (GLfloat*)&camMat   );
    glUniform1i(sh->getUloc("texture_1"), 0);

    glquad->draw(GL_TRIANGLES);
    glquad->drawPoints(10.0);

    // -------- from texture to Screen

    /*
    shader_post->use();
    glUniform1i(shader_post->getUloc("texZ"  ), 0);
    glUniform1i(shader_post->getUloc("texRGB"), 1);

    glActiveTexture(GL_TEXTURE0);    glBindTexture  (GL_TEXTURE_2D, frameBuff1.texZ   );
    glActiveTexture(GL_TEXTURE1);    glBindTexture  (GL_TEXTURE_2D, frameBuff1.texRGB );

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    //glViewport(0,0,WIDTH,HEIGHT);
    glEnableVertexAttribArray(0);
    object1->draw();
    */


    /*

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    //glActiveTexture(GL_TEXTURE0);  glBindTexture  (GL_TEXTURE_2D, frameBuff1.texRGB );

    glClearColor(0.0, 1.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );

    shader_pre->use();
    glUniform3fv( shader_pre->getUloc("light_pos"    ), 1, (const float[]){ 1.0f,  1.0f, -1.0f } );
    glUniform3fv( shader_pre->getUloc("lightColor"   ), 1, (const float[]){ 1.0f,  0.9f,  0.8f } );
    glUniform3fv( shader_pre->getUloc("diffuseColor" ), 1, (const float[]){ 1.0f,  1.0f,  1.0f } );
    glUniform3fv( shader_pre->getUloc("ambientColor" ), 1, (const float[]){ 0.2f,  0.2f,  0.3f } );
    glUniform3fv( shader_pre->getUloc("specularColor"), 1, (const float[]){ 0.0f,  0.0f,  0.0f } );
    glmesh->draw();

    // TODO:
    // copied from test_Sprites.cpp
    // why the fuck does it not work ?

    Mat4f mRot,mPersp;
    qCamera.toMatrix(mouseMat);
    mRot.setOne(); mRot.setRot(mouseMat);
    float fov = 3.0;
    mPersp.setPerspective( fov, fov, 1.0, 1000.0 );
    camMat.set_mmul_TN( mRot, mPersp );
    //Mat3f objRot; objRot.setOne();

    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);

    //Shader * sh = shader_post;
    Shader * sh = shader_pre;
    sh->use();
    sh->set_camPos( (GLfloat*)&camPos );
    sh->set_camMat( (GLfloat*)&camMat );

    //glUniformMatrix3fv( sh->getUloc( "camRot"   ), 1, GL_FALSE, (GLfloat*)&mouseMat );
    //glUniform4fv      ( sh->getUloc( "keyColor" ), 1, (const GLfloat[]){0.0f,0.0f,1.0f,10.0f} );

    //glmesh->draw();

    //uploadArrayBuffer( instances.pose_Up, instances.nInstances*3*sizeof(GLfloat), instance_Up );
    bilboards.draw( GL_TRIANGLES );
    */

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
                case SDLK_SPACE: preRender=!preRender; break;
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
