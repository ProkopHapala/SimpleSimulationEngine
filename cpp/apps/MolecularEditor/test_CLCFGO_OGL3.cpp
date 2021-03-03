#include <stdlib.h>
#include <stdio.h>

#include <GL/glew.h>
//#define GL3_PROTOTYPES 1
//#define GL_GLEXT_PROTOTYPES
//#include <GL/gl.h>
#include <SDL2/SDL.h>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"

//#include "IO_utils.h"
#include "Shader.h"
#include "GLObject.h"

long timeStart;
Vec3d DEBUG_dQdp;
int DEBUG_iter     = 0;
int DEBUG_log_iter = 0;

#include "Grid.h"
#include "GaussianBasis.h"
#include "CLCFGO.h"
#include "CLCFGO_tests.h"

// ===== Globals

bool bRun = false;
double dt;
CLCFGO eff;

Shader   * shader1, *shClear;
GLObject * object1;
GLObject * object2;

const int WIDTH  = 800;
const int HEIGHT = 800;
Vec2f resolution{WIDTH,HEIGHT};
Vec3f light_dir {  0.1, 0.2, 0.5 };
Vec3f lookAt    { 0.0f, 0.f,  0.f};
Vec3f camPos    { 0.0f, 0.f, -5.f};
float zoom = 5.0;
float iso  = 0.001;

Quat4f qCamera = Quat4fIdentity;

GLfloat vertexes[4][2] = {
	{  -0.99f,  -0.99f  },
	{  -0.99f,   0.99f  },
	{   0.99f,  -0.99f  },
	{   0.99f,   0.99f  } };

SDL_Window * window     = NULL;
SDL_GLContext   context = NULL;

GLuint vao;     // vertex array object





const int natoms_max = 100;
int natoms = natoms_max;

GLfloat atoms[natoms_max*4];
GLfloat coefs[natoms_max*4];

bool redraw = true;

// ============= FUNCTIONS

int frameCount = 0;
bool STOP = false;

void quit();
void die ( char const *msg );
void inputHanding ();
void init();
void draw();
void loop( int niters );

void setup(){


    // ==== Init eFF CLC-FGO

    //eff.loadFromFile( "data/H2.fgo", true);
    eff.loadFromFile( "inputs/H2.fgo", true);
    dt = 0.001;
    //exit(0);

    eff.turnAllSwitches(false);
    eff.bNormalize     = 1;
    eff.bEvalAE        = 1;
    //eff.bEvalAECoulomb = 1;
    //eff.bEvalCoulomb   = 1;
    eff.bEvalPauli     = 1;
    eff.bEvalKinetic   = 1;

    eff.bOptEPos = 1;
    eff.bOptSize = 1;

    eff.iPauliModel = 0;


    {
        // MO 0
        float sz=0.5;
        float l = 0.25;
        eff.ecoef[0]= 1.0; eff.epos[0].set(-l, 1.,0.);  eff.esize[0] = sz;
        eff.ecoef[1]= 1.0; eff.epos[1].set( l, 1.,0.);  eff.esize[1] = sz;
        // MO 1
        eff.ecoef[2]= 1.0; eff.epos[2].set(-l,-1.,0.);  eff.esize[2] = sz;
        eff.ecoef[3]=-1.0; eff.epos[3].set( l,-1.,0.);  eff.esize[3] = sz;
    }

    light_dir.normalize();


    object1 = new GLObject( );
    object1->draw_mode = GL_TRIANGLE_STRIP;
	object1->nVert   = 4;
    object1->buffs[0].setup(0,2,GL_FALSE,&vertexes[0][0],'v'); // vertexes
	object1->init();

    shClear=new Shader();
	shClear->init( "common_resources/shaders/plain_vert.glslv", "common_resources/shaders/clear.glslf" );

	shader1=new Shader();
	shader1->init( "common_resources/shaders/plain_vert.glslv", "common_resources/shaders/CLCFGO.glslf" );
	//shader1->init( "common_resources/shaders/plain_vert.glslv", "common_resources/shaders/rayMarchZero.glslf" );

}

void draw(){
    //glDisable( GL_DEPTH_TEST );
    //glColorMask (GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
    //glDepthMask (GL_TRUE);
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );

    //double E = eff.eval();
    //float F2 = eff.moveGD(dt);

    //bool bAlphaBlend = true;
    bool bAlphaBlend = false;
    if(bAlphaBlend){
        glDisable( GL_DEPTH_TEST );
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    }else{
        //glDisable( GL_DEPTH_TEST );
        glEnable( GL_DEPTH_TEST );
        glDepthFunc( GL_LESS );
    }
    Mat3f mouseMat;
    qCamera.toMatrix(mouseMat);
    camPos = lookAt  + mouseMat.c * -5.;

    GLuint uloc;

    //shClear->use();
    //uloc = glGetUniformLocation( shader1->shaderprogram, "iso" );	        glUniform1fv(uloc, 1, &iso     );
    //uloc = glGetUniformLocation( shader1->shaderprogram, "zoom" );	        glUniform1fv(uloc, 1, &zoom     );
    //object1->draw();

    shader1->use();
    uloc = glGetUniformLocation( shader1->shaderprogram, "iso" );	        glUniform1fv(uloc, 1, &iso     );
    uloc = glGetUniformLocation( shader1->shaderprogram, "zoom" );	        glUniform1fv(uloc, 1, &zoom     );
    uloc = glGetUniformLocation( shader1->shaderprogram, "lookAt" );	    glUniform3fv(uloc, 1, (float*)&lookAt      );
    uloc = glGetUniformLocation( shader1->shaderprogram, "camPos" );	    glUniform3fv(uloc, 1, (float*)&camPos      );
    uloc = glGetUniformLocation( shader1->shaderprogram, "resolution" );	glUniform2fv(uloc, 1, (float*)&resolution  );
    uloc = glGetUniformLocation( shader1->shaderprogram, "light_dir" );	    glUniform3fv(uloc, 1, (float*)&light_dir   );
    uloc = glGetUniformLocation( shader1->shaderprogram, "natoms"  );       glUniform1iv(uloc, 1, &eff.perOrb );
    GLuint uloc_atoms  = glGetUniformLocation( shader1->shaderprogram, "atoms"   );
    GLuint uloc_coefs  = glGetUniformLocation( shader1->shaderprogram, "coefs"   );
    Quat4f atoms[eff.perOrb];
    Quat4f coefs[eff.perOrb];
    for(int i=0; i<eff.nOrb; i++){
    //for(int i=0; i<1; i++){
        int io0 = eff.getOrbOffset(i);
        //for(int j=0; j<eff.perOrb; j++){
        for(int j=0; j<2; j++){

            atoms[j].f = (Vec3f)eff.epos [io0+j];
            atoms[j].e = eff.esize[io0+j];
            coefs[j].f = Vec3fZero;
            coefs[j].e = eff.ecoef[io0+j];

            //atoms[j].set(j*0.25,0,0,0.1);
            //coefs[j].set(0,0,0,1.0 - 2*(j%2) );
        }
        glUniform4fv(uloc_atoms, natoms, (float*)atoms );
        glUniform4fv(uloc_coefs, natoms, (float*)coefs );
        glEnableVertexAttribArray(0);
        object1->draw();
    }

    SDL_GL_SwapWindow(window);
    //redraw = false;

}

/*
// FUNCTION ======	inputHanding
void inputHanding(){
	SDL_Event event;
	while(SDL_PollEvent(&event)){
		if( event.type == SDL_KEYDOWN ){
			if(event.key.keysym.sym == SDLK_ESCAPE ) { quit(); }
			if(event.key.keysym.sym == SDLK_SPACE  ) { STOP=!STOP; }
			redraw = true;
		}
		if( event.type == SDL_QUIT){ quit();  };
	}
}
*/

// FUNCTION ======	inputHanding
void inputHanding(){
    float posstep = 0.1;
	SDL_Event event;

	const Uint8 *keys = SDL_GetKeyboardState(NULL);
	//Vec3d& v = lookAt;
	Vec3f& v = camPos;
    if( keys[ SDL_SCANCODE_W  ] ){ v.y += posstep; }
	if( keys[ SDL_SCANCODE_S  ] ){ v.y -= posstep; }
	if( keys[ SDL_SCANCODE_A  ] ){ v.x -= posstep; }
	if( keys[ SDL_SCANCODE_D  ] ){ v.x += posstep; }
    if( keys[ SDL_SCANCODE_Q  ] ){ v.z += posstep; }
	if( keys[ SDL_SCANCODE_E  ] ){ v.z -= posstep; }

    //if( keys[ SDL_SCANCODE_KP_PLUS   ] ){ resolution.mul(  1.1); }
	//if( keys[ SDL_SCANCODE_KP_MINUS  ] ){ resolution.mul(1/1.1); }

    if( keys[ SDL_SCANCODE_KP_PLUS   ] ){ zoom /= 1.1; }
	if( keys[ SDL_SCANCODE_KP_MINUS  ] ){ zoom *= 1.1; }

	while(SDL_PollEvent(&event)){
		if( event.type == SDL_KEYDOWN ){
            switch( event.key.keysym.sym ){
                case SDLK_ESCAPE: quit(); break;
            }
		}
		if( event.type == SDL_QUIT){ quit();  };
	}

	int mouseX,mouseY;
    float mouseRotSpeed = 0.01;
	int dmx,dmy;
	SDL_GetMouseState( &mouseX, &mouseY );
    Uint32 buttons = SDL_GetRelativeMouseState( &dmx, &dmy);
    //printf( " %i %i \n", mx,my );
    if ( buttons & SDL_BUTTON(SDL_BUTTON_RIGHT)) {
        Quat4f q; q.fromTrackball( 0, 0, dmx*mouseRotSpeed, dmy*mouseRotSpeed );
        //printf( " %i %i  (%3.3f,%3.3f,%3.3f,%3.3f) \n", dmx,dmy, q.x,q.y,q.z,q.w );
        //qCamera.qmul_T( q );
        qCamera.qmul( q );
        //qCamera.normalize();
    }
    //printf( "input handling camPos(%g,%g,%g) lookAt(%g,%g,%g) \n", camPos.x,camPos.y,camPos.z,   lookAt.x,lookAt.y,lookAt.z );

}

void loop( int nframes ){
    for ( int iframe=1; iframe<nframes; iframe++)    {
 		if( !STOP ) draw();
		inputHanding();
		frameCount++;
        //SDL_Delay(100);
        SDL_Delay(20);
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
    //SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 0); // Helps with glClear()
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
