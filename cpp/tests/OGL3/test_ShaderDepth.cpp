
// read this tutorial
// http://www.opengl-tutorial.org/beginners-tutorials/tutorial-4-a-colored-cube/

// http://gamedev.stackexchange.com/questions/93055/getting-the-real-fragment-depth-in-glsl

// TODO see later vbo indexing to avoid duplicating triangles
// http://www.opengl-tutorial.org/intermediate-tutorials/tutorial-9-vbo-indexing/

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

#include "Solids.h"
#include "GL3Utils.h"
#include "Mesh.h"

#include "GLObject.h"
#include "Shader.h"


Shader   * shader1;
GLObject * object1;

Mesh mesh;

GLuint vao;     // vertex array object

GLfloat modelPos[3] = { 0.0f,  0.0f,  -30.0f };
GLfloat modelMat[9] = {
  1.0f,  0.0f,  0.0f,
  0.0f,  1.0f,  0.0f,
  0.0f,  0.0f,  1.0f
};

GLfloat camRot[9] = {
  1.0f,  0.0f,  0.0f,
  0.0f,  1.0f,  0.0f,
  0.0f,  0.0f,  1.0f
};


GLfloat cam_pos      [3] = { 0.0f,  0.0f,  -10.0f };
GLfloat light_pos    [3] = { 1.0f,  1.0f,   -1.0f };

GLfloat lightColor   [3] = { 1.0f,  0.9f,   0.8f  };
//GLfloat lightColor   [3] = { 1.0f,  1.0f,   1.0f  };
GLfloat diffuseColor [3] = { 1.0f,  1.0f,   1.0f  };
//GLfloat diffuseColor [3] = { 0.0f,  0.0f,   0.0f  };
GLfloat ambientColor [3] = { 0.2f,  0.2f,   0.3f  };
GLfloat specularColor[3] = { 1.0f,  1.0f,   1.0f  };
//GLfloat specularColor[3] = { 0.0f,  0.0f,   0.0f  };


int WIDTH  = 800;
int HEIGHT = 800;


int mouseX, mouseY;
SDL_Window * window     = NULL;
SDL_GLContext   context = NULL;
Quat4f qCamera;

//   http://www.songho.ca/opengl/gl_projectionmatrix.html

void getPerspectiveMatrix( float xmin, float xmax, float ymin, float ymax, float zmin, float zmax, float * mat ){
    float invdx = xmax-xmin;
    float invdy = ymax-ymin;
    float invdz = zmax-zmin;
    mat[0 ]  = 2*zmin*invdx; mat[1 ] = 0;            mat[2 ] =  (xmax+xmin)*invdx;  mat[3 ] = 0;
    mat[4 ]  = 0;            mat[5 ] = 2*zmin*invdy; mat[6 ] =  (ymin+ymax)*invdy;  mat[7 ] = 0;
    mat[8 ]  = 0;            mat[9 ] = 0;            mat[10] = -(zmin+zmax)*invdz;  mat[11] = -2*zmax*zmin*invdz;
    mat[12]  = 0;            mat[13] = 0;            mat[14] = -1;                  mat[15] = 0;
}

// ============= FUNCTIONS

int frameCount = 0;
bool STOP = false;

void quit();
void die ( char const *msg );
void inputHanding ();
void init();
void draw();
void loop( int niters );


int render_type = 1;

void setup(){


    if ( render_type == 0      ){
        // --- vertex const color
        shader1=new Shader();
        shader1->init( "shaders/basicColor3D_vert.c", "shaders/basicColor3D_frag.c" );
        glUseProgram(shader1->shaderprogram);

    }else if ( render_type == 1 ){
        // --- shading
        shader1=new Shader();
        shader1->init( "shaders/basicShading3D_vert.c", "shaders/basicShading3D_frag.c" );
        glUseProgram(shader1->shaderprogram);
    };
	//mesh.fromFileOBJ("common_resources/turret.obj");


	/*
    object1 = new GLObject( );
    object1->nVert    = nVert;
    object1->buffs[0].setup(0,3,GL_FALSE,vertexes,'v'); // vertexes
    object1->buffs[1].setup(1,3,GL_FALSE,vertexes,'n'); // normals
    object1->init();
    */

    /*
    object1 = new GLObject( );
    object1->nVert    = nVert;
    object1->buffs[0].setup(0,3,GL_FALSE,vertexes,'v'); // vertexes
    object1->buffs[1].setup(1,3,GL_FALSE,vertexes,'n'); // normals
    object1->init();
    */

    /*
    object1 = new GLObject( );
    object1->nVert    = Solids::Octahedron_nverts;
    object1->setIndexes( Solids::Octahedron_nfaces*3, Solids::Octahedron_faces );
    printf( "nVert %i nInd %i \n", object1->nVert, object1->nInd);
    object1->buffs[0].setup(0,3,GL_FALSE,(double*)Solids::Octahedron_verts, object1->nVert, 'v'); // vertexes
    object1->buffs[1].setup(1,3,GL_FALSE,(double*)Solids::Octahedron_verts, object1->nVert, 'n'); // normals
    object1->init();
    */

    /*
    object1 = new GLObject( );
    object1->nVert    = Solids::Tetrahedron_nverts;
    object1->setIndexes( Solids::Tetrahedron_nfaces*3, Solids::Tetrahedron_faces );
    printf( "nVert %i nInd %i \n", object1->nVert, object1->nInd);
    object1->buffs[0].setup(0,3,GL_FALSE,(double*)Solids::Tetrahedron_verts, object1->nVert, 'v'); // vertexes
    object1->buffs[1].setup(1,3,GL_FALSE,(double*)Solids::Tetrahedron_verts, object1->nVert, 'n'); // normals
    object1->init();
    */

    /*
    object1 = new GLObject( );
    object1->nVert    = Solids::Icosahedron_nverts;
    object1->setIndexes( Solids::Icosahedron_nfaces*3, Solids::Icosahedron_faces );
    printf( "nVert %i nInd %i \n", object1->nVert, object1->nInd);
    object1->buffs[0].setup(0,3,GL_FALSE,(double*)Solids::Icosahedron_verts, object1->nVert, 'v'); // vertexes
    object1->buffs[1].setup(1,3,GL_FALSE,(double*)Solids::Icosahedron_verts, object1->nVert, 'n'); // normals
    object1->init();
    */

    /*
    int nVert = 3*countTris( Solids::Tetrahedron_nfaces, Solids::Tetrahedron_ngons );
    GLfloat * verts   = new GLfloat[nVert];
    GLfloat * normals = new GLfloat[nVert];
    hardFace( Solids::Tetrahedron_nfaces, Solids::Tetrahedron_ngons, Solids::Tetrahedron_faces, Solids::Tetrahedron_verts, verts, normals );
    */

    int nVert = countVerts( Solids::Icosahedron_nfaces, Solids::Icosahedron_ngons );
    GLfloat * verts   = new GLfloat[nVert*3];
    GLfloat * normals = new GLfloat[nVert*3];
    hardFace( Solids::Icosahedron_nfaces, Solids::Icosahedron_ngons, Solids::Icosahedron_faces, Solids::Icosahedron_verts, verts, normals );

    object1 = new GLObject( );
    object1->nVert    = nVert;
    object1->buffs[0].setup(0,3,GL_FALSE,verts,'v'); // vertexes
    object1->buffs[1].setup(1,3,GL_FALSE,normals,'n'); // normals
    object1->init();

	qCamera.setOne();
}

void draw(){
    glClearColor(0.0, 0.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );

    glEnable( GL_DEPTH_TEST );
    glDepthFunc( GL_LESS );

    float camMat[16];
    getPerspectiveMatrix( -WIDTH, WIDTH, -HEIGHT, HEIGHT, 1.0, 10.0, camMat );

    //Quat4f qCamera_; convert(qCamera,qCamera_);
    Mat3f mouseMat; qCamera.toMatrix(mouseMat);
    //printf( " (%3.3f,%3.3f,%3.3f,%3.3f) \n", qCamera.x,qCamera.y,qCamera.z,qCamera.w );
    //printf( " (%3.3f,%3.3f,%3.3f,%3.3f) \n", qCamera.x,qCamera.y,qCamera.z,qCamera.w );
    //printf( "mouseMat (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) \n", mouseMat.ax, mouseMat.ay, mouseMat.az,  mouseMat.bx, mouseMat.by, mouseMat.bz,   mouseMat.cx, mouseMat.cy, mouseMat.cz );

    GLuint uloc;
    uloc = glGetUniformLocation( shader1->shaderprogram, "modelPos" ); glUniform3fv      (uloc, 1, modelPos );
    //uloc = glGetUniformLocation( shader1->shaderprogram, "modelMat" ); glUniformMatrix3fv(uloc, 1, GL_FALSE, modelMat );
    uloc = glGetUniformLocation( shader1->shaderprogram, "modelMat" ); glUniformMatrix3fv(uloc, 1, GL_FALSE, (float*)&mouseMat );
    uloc = glGetUniformLocation( shader1->shaderprogram, "camMat"   ); glUniformMatrix4fv(uloc, 1, GL_FALSE, camMat   );

    // shading
    if ( render_type == 1 ){
        uloc = glGetUniformLocation( shader1->shaderprogram, "cam_pos"       ); glUniform3fv      (uloc, 1, cam_pos      );
        uloc = glGetUniformLocation( shader1->shaderprogram, "light_pos"     ); glUniform3fv      (uloc, 1, light_pos     );
        uloc = glGetUniformLocation( shader1->shaderprogram, "lightColor"    ); glUniform3fv      (uloc, 1, lightColor    );
        uloc = glGetUniformLocation( shader1->shaderprogram, "diffuseColor"  ); glUniform3fv      (uloc, 1, diffuseColor  );
        uloc = glGetUniformLocation( shader1->shaderprogram, "ambientColor"  ); glUniform3fv      (uloc, 1, ambientColor  );
        uloc = glGetUniformLocation( shader1->shaderprogram, "specularColor" ); glUniform3fv      (uloc, 1, specularColor );
    };

    //uloc = glGetUniformLocation( shader1->shaderprogram, "light_dir"); glUniform3fv(uloc, 1, light_dir  );

    object1->draw();

    SDL_GL_SwapWindow(window);

}


// =============================================================
// ========== BORING DEFAULT RUTINES
// =============================================================

// FUNCTION ======	inputHanding
void inputHanding(){
    float posstep = 0.1;
	SDL_Event event;
	while(SDL_PollEvent(&event)){
		if( event.type == SDL_KEYDOWN ){
            switch( event.key.keysym.sym ){
                case SDLK_ESCAPE: quit(); break;
                case SDLK_w: modelPos[1] +=posstep; break;
                case SDLK_s: modelPos[1] -=posstep; break;
                case SDLK_a: modelPos[0] +=posstep; break;
                case SDLK_d: modelPos[0] -=posstep; break;
                case SDLK_q: modelPos[2] +=posstep*10; break;
                case SDLK_e: modelPos[2] -=posstep*10; break;
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
