
// read this tutorial
// http://prideout.net/blog/?p=61      The Little Grasshopper : Tron, Volumetric Lines, and Meshless Tubes

#include <stdlib.h>
#include <stdio.h>

#include <GL/glew.h>
//#define GL3_PROTOTYPES 1
//#define GL_GLEXT_PROTOTYPES
//#include <GL/gl.h>
#include <SDL2/SDL.h>

#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

//#include "Mesh.h"
#include "Solids.h"

#include "GLfunctions.h"
#include "GLobjects.h"

#include "GLObject.h"
#include "Shader.h"
#include "IO_utils.h"

int WIDTH  = 800;
int HEIGHT = 800;
SDL_Window   * window  = NULL;
SDL_GLContext  context = NULL;

GLuint vao;
Shader   *shader1,*shConst;

GLuint p_sh;
//GLObject * object1;
//GLMesh   *glmesh, *gledges;
//Mesh mesh;

GLuint curve_vbo;

int mouseX, mouseY;
Quat4f qCamera;
Mat3f  mouseMat;
Vec3f  camPos   = (Vec3f){ 0.0f,0.0f, -10.0f };
Vec3f  modelPos = (Vec3f){ 0.0f,0.0f, 0.0f };

//GLuint shaderProgram,vertexShader,fragmentShader,geometryShader;
//GLuint vbo;

int NPoints = 40;

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

/*
void RenderCurve(Curve curve){
    glBindBuffer(GL_ARRAY_BUFFER, curve.Vbo);
    glEnableVertexAttribArray(PositionSlot  );
    glEnableVertexAttribArray(NormalSlot    );
    glEnableVertexAttribArray(PathCoordSlot );

    GLsizei stride = sizeof(float) * 7;
    const GLvoid* normalOffset    = (GLvoid*) (sizeof(float) * 3);
    const GLvoid* pathCoordOffset = (GLvoid*) (sizeof(float) * 6);

    glVertexAttribPointer(PositionSlot,  3, GL_FLOAT, GL_FALSE, stride, 0              );
    glVertexAttribPointer(NormalSlot,    3, GL_FLOAT, GL_FALSE, stride, normalOffset   );
    glVertexAttribPointer(PathCoordSlot, 1, GL_FLOAT, GL_FALSE, stride, pathCoordOffset);

    glDrawArrays(GL_LINE_STRIP_ADJACENCY_EXT, 0, curve.Count);

    glDisableVertexAttribArray(PositionSlot);
    glDisableVertexAttribArray(NormalSlot);
    glDisableVertexAttribArray(PathCoordSlot);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}
*/

const int PositionSlot  = 0;
const int NormalSlot    = 1;
const int PathCoordSlot = 2;


void checkShaderError( char * errLog, const char * comment ){
    if( errLog != NULL ){
        //printf( " Error in Vertex Shader: \n" );
        printf( "%s %s \n", comment, errLog );
        free( errLog );
        //return -1;
        exit(-1);
    }
}

GLuint LoadProgram(const char* vsName, const char* gsName, const char* fsName, GLenum prim){
    GLint compileSuccess, linkSuccess;
    GLchar compilerSpew[256];

    GLuint programHandle = glCreateProgram();

    printf("DEBUG 0 \n");

    //const char* vsSource = glswGetShader(vsKey);
    const char* vsSource = filetobuf( vsName );
    if(!vsSource) printf(  "Can't find vertex shader: %s\n", vsName);
    GLuint vsHandle = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vsHandle, 1, &vsSource, 0);
    glCompileShader(vsHandle);
    glGetShaderiv(vsHandle, GL_COMPILE_STATUS, &compileSuccess);
    glGetShaderInfoLog(vsHandle, sizeof(compilerSpew), 0, compilerSpew);
    if(!compileSuccess) printf( "Can't compile %s:\n%s", vsName, compilerSpew );
    glAttachShader(programHandle, vsHandle);

    printf("DEBUG 1 \n");

    //const char* fsSource = glswGetShader(fsKey);
    const char* fsSource = filetobuf( fsName );
    if(!fsSource) printf(  "Can't find fragment shader: %s\n", fsName);
    GLuint fsHandle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fsHandle, 1, &fsSource, 0);
    glCompileShader(fsHandle);
    glGetShaderiv(fsHandle, GL_COMPILE_STATUS, &compileSuccess);
    glGetShaderInfoLog(fsHandle, sizeof(compilerSpew), 0, compilerSpew);
    if(!compileSuccess) printf( "Can't compile %s:\n%s", fsName, compilerSpew);
    glAttachShader(programHandle, fsHandle);

    printf("DEBUG 2 \n");

    if (gsName) {
        //const char* gsSource = glswGetShader(gsKey);
        const char* gsSource = filetobuf( gsName );
        if(!gsSource) printf( "Can't find geometry shader: %s\n", gsName);
        GLuint gsHandle = glCreateShader(GL_GEOMETRY_SHADER);
        glShaderSource(gsHandle, 1, &gsSource, 0);
        glCompileShader(gsHandle);
        glGetShaderiv(gsHandle, GL_COMPILE_STATUS, &compileSuccess);
        glGetShaderInfoLog(gsHandle, sizeof(compilerSpew), 0, compilerSpew);
        if(!compileSuccess) printf("Can't compile %s:\n%s", gsName, compilerSpew);
        glAttachShader(programHandle, gsHandle);

        if (prim == GL_TRIANGLE_STRIP) {
            glProgramParameteriEXT(programHandle, GL_GEOMETRY_INPUT_TYPE_EXT, GL_LINES_ADJACENCY_EXT);
            glProgramParameteriEXT(programHandle, GL_GEOMETRY_OUTPUT_TYPE_EXT, GL_TRIANGLE_STRIP);
            glProgramParameteriEXT(programHandle, GL_GEOMETRY_VERTICES_OUT_EXT, 24);
        } else {
            glProgramParameteriEXT(programHandle, GL_GEOMETRY_INPUT_TYPE_EXT, GL_LINES);
            glProgramParameteriEXT(programHandle, GL_GEOMETRY_OUTPUT_TYPE_EXT, GL_TRIANGLE_STRIP);
            glProgramParameteriEXT(programHandle, GL_GEOMETRY_VERTICES_OUT_EXT, 4);
        }
    }

    printf("DEBUG 3 \n");

    glBindAttribLocation(programHandle, PositionSlot,  "Position"  );
    glBindAttribLocation(programHandle, NormalSlot,    "Normal"    );
    glBindAttribLocation(programHandle, PathCoordSlot, "PathCoord" );

    printf("DEBUG 4 \n");

    glLinkProgram(programHandle);
    glGetProgramiv(programHandle, GL_LINK_STATUS, &linkSuccess);
    glGetProgramInfoLog(programHandle, sizeof(compilerSpew), 0, compilerSpew);
    if(!linkSuccess) printf( "error linking shader %s",compilerSpew);

    printf("DEBUG 5 \n");

    glUseProgram(programHandle);

    printf("DEBUG 6 \n");
    return programHandle;
}



float * makeLines( const CMesh& cmesh){
    NPoints = cmesh.nedge * 2;
    float * points = new float[NPoints*7];
    float * pi = points;
    for (size_t i=0; i<NPoints; i++ ) {
        int ip  = ((int*)cmesh.edges)[i];
        Vec3d p = cmesh.verts[ip];
        *pi++ = p.x; *pi++ = p.y; *pi++ = p.z;
        *pi++ = 0;   *pi++ = 0;   *pi++ = 1;
        *pi++ = (float)i/float(NPoints);
    }
    return points;
}


void setup(){

    /*
    shader1 = new Shader();
    shader1->init( "shaders/Tubes.glslv", "shaders/Tubes.glslf", "shaders/Tubes_GS.glslv" );
    shader1->use();
    shader1->setUniformVec3f("DiffuseMaterial", (Vec3f){1.0f, 0.5f, 0.125f}   );
    shader1->setUniformVec3f("AmbientMaterial", (Vec3f){0.125f, 0.125f, 0.0f} );
    shader1->setUniformVec3f("SpecularMaterial",(Vec3f){0.5f, 0.5f, 0.5f}     );
    shader1->setUniformf    ("Shininess", 50.0f);
    */

    //p_sh =  LoadProgram(  "shaders/Tubes.glslv",  "shaders/Tubes_GS.glslv", "shaders/Tubes.glslf", GL_TRIANGLE_STRIP );
    p_sh =  LoadProgram(  "shaders/Lines2Tubes.glslv",  "shaders/Lines2Tubes_GS.glslv", "shaders/Lines2Tubes.glslf", GL_TRIANGLE_STRIP );
    //p_sh =  LoadProgram(  "common_resources/shaders/color3D.glslv",  NULL, "common_resources/shaders/color3D.glslf", GL_TRIANGLE_STRIP );

    //glUniform1iv(, 1, (GLint*  )&i           ); };
    glUniform3fv(glGetUniformLocation(p_sh,"DiffuseMaterial" ), 1, (const float[]){1.0f,   0.5f,   0.125f} );
    glUniform3fv(glGetUniformLocation(p_sh,"AmbientMaterial" ), 1, (const float[]){0.125f, 0.125f, 0.0f}   );
    glUniform3fv(glGetUniformLocation(p_sh,"SpecularMaterial"), 1, (const float[]){0.5f,   0.5f,   0.5f}   );
    glUniform1f (glGetUniformLocation(p_sh,"Shininess"       ), 50.0f );

    //uniform vec3 modelPos;
    //uniform mat3 modelMat;
    //uniform vec3 camPos;
    //uniform mat4 camMat;
    //uniform vec3 light_dir;

    printf("DEBUG 7 \n");

    //glmesh = new GLMesh();
    //glmesh->init_d( mesh.points.size(), mesh.triangles.size()*3, ((int*)&mesh.triangles[0]), (double*)&(mesh.points [0]), (double*)&(mesh.normals[0]), NULL, NULL );

    //float * points = makeLines( Solids::Tetrahedron );
    float * points = makeLines( Solids::Icosahedron );

    printf("DEBUG 8 \n");

    glGenBuffers(1, &curve_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, curve_vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*7*NPoints, points, GL_STATIC_DRAW);
    delete [] points;

    printf("DEBUG 9 \n");

    qCamera.setOne();
}

void draw(){
    glClearColor(0.5f, 0.5f, 0.5f, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );

    glEnable( GL_DEPTH_TEST );
    glDisable(GL_BLEND);

    glUseProgram(p_sh);

    Vec3f light_dir = (Vec3f){1.0,1.0,1.0}; light_dir.normalize();

    // debug - simple shader



    Mat4f camMat;   camMat.setOne();    camMat.setPerspective( 10.0, 10.0, 2.0, 1000.0 );
    Mat3f modelMat; modelMat.setOne(); qCamera.toMatrix( modelMat);
    glUniformMatrix4fv(glGetUniformLocation(p_sh,"camMat"   ), 1, GL_FALSE, (GLfloat*)&camMat   );
    glUniformMatrix3fv(glGetUniformLocation(p_sh,"modelMat" ), 1, GL_FALSE, (GLfloat*)&modelMat );
    glUniform3fv      (glGetUniformLocation(p_sh,"modelPos" ), 1, (GLfloat*)&modelPos );
    glUniform3fv      (glGetUniformLocation(p_sh,"camPos"   ), 1, (GLfloat*)&camPos   );
    glUniform3fv      (glGetUniformLocation(p_sh,"light_dir"), 1, (GLfloat*)&light_dir );





    glBindBuffer(GL_ARRAY_BUFFER, curve_vbo);
    glEnableVertexAttribArray(PositionSlot);
    glEnableVertexAttribArray(NormalSlot);
    glEnableVertexAttribArray(PathCoordSlot);

    GLsizei stride = sizeof(float) * 7;
    const GLvoid* normalOffset    = (GLvoid*) (sizeof(float) * 3);
    const GLvoid* pathCoordOffset = (GLvoid*) (sizeof(float) * 6);

    glVertexAttribPointer(PositionSlot,  3, GL_FLOAT, GL_FALSE, stride, 0);
    glVertexAttribPointer(NormalSlot,    3, GL_FLOAT, GL_FALSE, stride, normalOffset);
    glVertexAttribPointer(PathCoordSlot, 1, GL_FLOAT, GL_FALSE, stride, pathCoordOffset);

    //glDrawArrays(GL_LINE_STRIP_ADJACENCY_EXT, 0, NPoints );
    //glPointSize(100.0); glDrawArrays(GL_POINTS, 0, NPoints );
    glDrawArrays(GL_LINES, 0, NPoints );

    glDisableVertexAttribArray(PositionSlot);
    glDisableVertexAttribArray(NormalSlot);
    glDisableVertexAttribArray(PathCoordSlot);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    SDL_GL_SwapWindow(window);

}


// =============================================================
// ========== BORING DEFAULT RUTINES
// =============================================================

// FUNCTION ======	inputHanding
void inputHanding(){
    float posstep = 0.1;
	SDL_Event event;

	const Uint8 *keys = SDL_GetKeyboardState(NULL);
	/*
    if( keys[ SDL_SCANCODE_W  ] ){ modelPos.y += posstep; }
	if( keys[ SDL_SCANCODE_S  ] ){ modelPos.y -= posstep; }
	if( keys[ SDL_SCANCODE_A  ] ){ modelPos.x -= posstep; }
	if( keys[ SDL_SCANCODE_D  ] ){ modelPos.x += posstep; }
    if( keys[ SDL_SCANCODE_Q  ] ){ modelPos.z += posstep; }
	if( keys[ SDL_SCANCODE_E  ] ){ modelPos.z -= posstep; }
	*/

    if( keys[ SDL_SCANCODE_W  ] ){ camPos.y += posstep; }
	if( keys[ SDL_SCANCODE_S  ] ){ camPos.y -= posstep; }
	if( keys[ SDL_SCANCODE_A  ] ){ camPos.x -= posstep; }
	if( keys[ SDL_SCANCODE_D  ] ){ camPos.x += posstep; }
    if( keys[ SDL_SCANCODE_Q  ] ){ camPos.z += posstep; }
	if( keys[ SDL_SCANCODE_E  ] ){ camPos.z -= posstep; }


	while(SDL_PollEvent(&event)){
		if( event.type == SDL_KEYDOWN ){
            switch( event.key.keysym.sym ){
                case SDLK_ESCAPE: quit(); break;
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
