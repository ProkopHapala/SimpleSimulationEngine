
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
//const int MaxParticles = 1600;
const int MaxParticles = 16;
//const int MaxParticles = 256;
//const int MaxParticles = 256*256;

GLfloat *instance_pos,*instance_sc;

Shader *shSprite,*shSpriteBlend;

//GLInstances instances;
GLBillboards bilboards;
GLuint texture_1;

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

bool bTransparent  = false;

long lastCPUtick = 0;
double ticks_per_second=0;

// =============== Functions

const char str_glslf_sin[]= GLSL(330,
    in        vec3 fpos_world;
    out       vec4 gl_FragColor;
    void main(){ gl_FragColor = vec4( sin( fpos_world*30.0 ), 1.0 ); }
    //void main(){ gl_FragColor = vec4( 0.0,0.0,0.0, 1.0 ); }
);

double calibrate_timer(int delay){
    long t1 = getCPUticks();
    SDL_Delay(delay);
    long t2 = getCPUticks();
    return (t2-t1)/(delay*0.001d);
}

int juliaPoint( double x, double y, int maxIters ){
    //const int maxIters = 64;
    const double cX=-0.73,cY=0.27015;
    double x_;
    for(int i=0; i<maxIters; i++){
        x_ = x*x-y*y  + cX;
        y  = 2.0d*x*y + cY;
        x = x_;
        if((x*x+y*y)>4.0d) return i;
    }
    return maxIters;
};

int setup(){

    shSprite=new Shader();
    shSprite->init( "common_resources/shaders/Bilboard3D.glslv", "common_resources/shaders/hardSprite.glslf" );
    shSprite->getDefaultUniformLocation();

    shSpriteBlend=new Shader();
    shSpriteBlend->init( "common_resources/shaders/Bilboard3D.glslv", "common_resources/shaders/texture.glslf" );
    shSpriteBlend->getDefaultUniformLocation();

    /*
    char* str_glslv_Bilboard3D = filetobuf( "common_resources/shaders/Bilboard3D.glslv"  ); if(str_glslv_Bilboard3D==NULL){ printf("fail to load shader!\n"); exit(1); };
    shSprite=new Shader();
    shSprite->init_str( str_glslv_Bilboard3D, str_glslf_sin );
    shSprite->getDefaultUniformLocation();
    delete [] str_glslv_Bilboard3D;
    */

	instance_pos    = new GLfloat[MaxParticles*3];
	instance_sc     = new GLfloat[MaxParticles*2];

	float span = 20.0f;
    float time = frameCount * 0.005;
    ParticlesCount = MaxParticles;

    int nside = int( pow( ParticlesCount, 1.0/3.0) );

	for( int i=0; i<ParticlesCount; i++ ){
        Vec3f pos;
        Vec2f sc;
        // pos.set( randf(-span,span), randf(-span,span), randf(-span,span) );
        // sc.set( randf(0.5,1.5),randf(0.5,1.5),randf(0.5,1.5) );
        pos.set( i%nside , (i/nside)%nside , (i/(nside*nside)) ); pos.mul(1.2);
        //sc.set( 1.0,1.0,1.0 );
        //float sz = randf(0.5,1.5); sc.set( sz,sz );

        //sc.set( randf(0.5,1.5),randf(0.5,1.5) );

        sc.set( 1.0, 1.0 );

        *((Vec3f*)(instance_pos+(3*i))) = pos;
        *((Vec2f*)(instance_sc +(2*i))) = sc;
	}

	bilboards.init( MaxParticles, 6, DEFAULT_Bilboard_UVs, instance_pos, instance_sc );

    int imgH  = 512; int imgW  = 512;
    double dx = 1.0d/(imgW);
    double dy = 1.0d/(imgH);
    uint32_t * c_img1 = new uint32_t[imgH*imgW];
    for( int iy=0; iy<imgH; iy++ ){
        for( int ix=0; ix<imgW; ix++ ){
            uint8_t r,g,b,a;
            //r=0; g=ix^iy; b=255; a=ix^iy;
            double x = (ix*2-imgW)*dx;
            double y = (iy*2-imgH)*dy;

            int nMaxIter = 64;

            /*
            int iters = (nMaxIter-juliaPoint(x*2,y,nMaxIter))*4; if(iters>255)iters=255;
            r=iters; g=0; b=iters*2.0; a=255-iters;
            */

            uint8_t  mask = 31;
            uint8_t  cut  = 4;
            //a=(((iy&mask)<cut)|((ix&mask)<cut))*255;
            a=((((ix+iy)&mask)<cut)|(((ix-iy+mask)&mask)<cut) |  ((ix&mask)<cut)|((iy&mask)<cut))*255;
            r=ix/2;  g=iy/2; b=(ix+iy)/4;
            c_img1[ iy*imgW + ix ] = (a<<24) | (b<<16) | (g<<8) | (r);
        }
    }
    newTexture2D( texture_1, imgW, imgH, c_img1, GL_RGBA, GL_UNSIGNED_BYTE );

    ticks_per_second = calibrate_timer(100);
    //lastTime = glfwGetTime();
    lastTime = 0.0;
    qCamera.setOne();
    return 0;
};

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

    Mat4f camMat,mRot,mPersp;
    qCamera.toMatrix(mouseMat);
    mRot.setOne(); mRot.setRot(mouseMat);
    float fov = 3.0;
    mPersp.setPerspective( fov, fov*ASPECT_RATIO, 1.0, 1000.0 );
    camMat.set_mmul_TN( mRot, mPersp );
    //Mat3f objRot; objRot.setOne();

    Shader * sh=shSprite;
    if( bTransparent ){
    	glDisable(GL_BLEND);
        glEnable(GL_DEPTH_TEST);
    }else{
        glDisable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        sh=shSpriteBlend;
    }

    sh->use();
    sh->set_camPos( (GLfloat*)&camPos );
    sh->set_camMat( (GLfloat*)&camMat );

    glUniformMatrix3fv( sh->getUloc( "camRot" ), 1, GL_FALSE, (GLfloat*)&mouseMat );
    glUniform4fv( sh->getUloc( "keyColor" ), 1, (const GLfloat[]){1.0f,0.0f,1.0f,10.0f} );

    //uploadArrayBuffer( instances.pose_Up, instances.nInstances*3*sizeof(GLfloat), instance_Up );
    bilboards.draw( GL_TRIANGLES );

    //glPointSize(10.0); bilboards.draw( GL_POINTS );

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


    SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 1);
    //SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 2);
    //SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 8);
    SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 16);
    glEnable(GL_MULTISAMPLE);


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
