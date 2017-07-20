
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

#include "Shader.h"

// =============== Types

// CPU representation of a particle
struct Particle{
	//glm::vec3 pos, speed;
	Vec3f pos, speed;
	unsigned char r,g,b,a; // Color
	float size, angle, weight;
	float life; // Remaining life of the particle. if <0 : dead and unused.
	float cameradistance; // *Squared* distance to the camera. if dead : -1.0f

	bool operator<(const Particle& that) const {
		// Sort in reverse order : far particles drawn first.
		return this->cameradistance > that.cameradistance;
	}
};

// =============== Global variables

//const int MaxParticles = 100000;
const int MaxParticles = 16;
Particle ParticlesContainer[MaxParticles];
int LastUsedParticle = 0;

GLuint VertexArrayID;
GLfloat* g_particule_position_size_data = new GLfloat[MaxParticles * 4];
GLubyte* g_particule_color_data         = new GLubyte[MaxParticles * 4];

//GLuint programID;
Shader* shaderParticle;

GLuint CameraRight_worldspace_ID;
GLuint CameraUp_worldspace_ID;
GLuint ViewProjMatrixID;

GLuint billboard_vertex_buffer;
GLuint particles_color_buffer;
GLuint particles_position_buffer;

int nVerts = 0;
int ParticlesCount = 0;
int frameCount = 0;
double lastTime = 0.0;

GLuint vao;     // vertex array object
int delay = 1;
int VSync = 0;
SDL_Window * window     = NULL;
SDL_GLContext   context = NULL;
int WIDTH  = 800;
int HEIGHT = 800;
float ASPECT_RATIO = HEIGHT/WIDTH;

int mouseX, mouseY;
Quat4f qCamera;

// =============== Functions

// Finds a Particle in ParticlesContainer which isn't used yet.
// (i.e. life < 0);
int FindUnusedParticle(){
	for(int i=LastUsedParticle; i<MaxParticles; i++){
		if (ParticlesContainer[i].life < 0){
			LastUsedParticle = i;
			return i;
		}
	}
	for(int i=0; i<LastUsedParticle; i++){
		if (ParticlesContainer[i].life < 0){
			LastUsedParticle = i;
			return i;
		}
	}
	return 0; // All particles are taken, override the first one
}

void SortParticles(){
	std::sort(&ParticlesContainer[0], &ParticlesContainer[MaxParticles]);
}

int setup(){

	// Create and compile our GLSL program from the shaders
	//GLuint programID = LoadShaders( "Particle.vertexshader", "Particle.fragmentshader" );
	//programID = LoadShaders( "Particle.glslv",  "Particle.glslf" );

	shaderParticle=new Shader();
	//common_resources/shaders/
	shaderParticle->init( "common_resources/shaders/Particle.glslv",   "common_resources/shaders/Particle.glslf"   );
    //shaderParticle->init( "Particle.glslv",   "Particle.glslf"   );

	CameraRight_worldspace_ID  = glGetUniformLocation(shaderParticle->shaderprogram, "CameraRight_worldspace");
	CameraUp_worldspace_ID     = glGetUniformLocation(shaderParticle->shaderprogram, "CameraUp_worldspace");
	ViewProjMatrixID           = glGetUniformLocation(shaderParticle->shaderprogram, "VP");

	// fragment shader
	//GLuint TextureID  = glGetUniformLocation(programID, "myTextureSampler");

	g_particule_position_size_data = new GLfloat[MaxParticles * 4];
	g_particule_color_data         = new GLubyte[MaxParticles * 4];

	for(int i=0; i<MaxParticles; i++){
		ParticlesContainer[i].life = -1.0f;
		ParticlesContainer[i].cameradistance = -1.0f;
	}

	/*
	nVerts = 8;
	//static const
	GLfloat g_vertex_buffer_data[nVerts*3] = {
        -0.5f, -0.5f, 0.0f,
         0.5f, -0.5f, 0.0f,
        -0.5f,  0.5f, 0.0f,
         0.5f,  0.5f, 0.0f,

        -1.5f, -1.5f, 0.0f,
         1.5f, -1.5f, 0.0f,
        -1.5f,  1.5f, 0.0f,
         1.5f,  1.5f, 0.0f,
	};
	*/

	nVerts = 50;
    float dx = 1.0f/nVerts;
	Vec3f g_vertex_buffer_data[nVerts*3];
	for( int i=0; i<nVerts; i++){
        float x = i*dx;
        g_vertex_buffer_data[i].x = x*sin(x*10.0);
        g_vertex_buffer_data[i].y = x*cos(x*10.0);
        g_vertex_buffer_data[i].z = x;
	}


	glGenBuffers(1, &billboard_vertex_buffer);
	glBindBuffer(GL_ARRAY_BUFFER, billboard_vertex_buffer);
	glBufferData(GL_ARRAY_BUFFER, nVerts*3* sizeof(GLfloat), (GLfloat*)g_vertex_buffer_data, GL_STATIC_DRAW);

	// The VBO containing the positions and sizes of the particles
	glGenBuffers(1, &particles_position_buffer);
	glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
	// Initialize with empty (NULL) buffer : it will be updated later, each frame.
	glBufferData(GL_ARRAY_BUFFER, MaxParticles * 4 * sizeof(GLfloat), NULL, GL_STREAM_DRAW);

	// The VBO containing the colors of the particle
	glGenBuffers(1, &particles_color_buffer);
	glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
	// Initialize with empty (NULL) buffer : it will be updated later, each frame.
	glBufferData(GL_ARRAY_BUFFER, MaxParticles * 4 * sizeof(GLubyte), NULL, GL_STREAM_DRAW);

    //lastTime = glfwGetTime();
    lastTime = 0.0;
    return 0;
};



void physics(){

		//double currentTime = glfwGetTime();
		double currentTime = frameCount * 0.005;
		double delta = currentTime - lastTime;
		lastTime = currentTime;

		printf( "time: %g delta: %g \n", currentTime, delta );
		int newparticles = (int)(delta*10000.0);
		if (newparticles > (int)(0.016f*10000.0))
			newparticles = (int)(0.016f*10000.0);

		for(int i=0; i<newparticles; i++){
			int particleIndex = FindUnusedParticle();
			ParticlesContainer[particleIndex].life = 5.0f; // This particle will live 5 seconds.
			//ParticlesContainer[particleIndex].pos = glm::vec3(0,0,-20.0f);
			ParticlesContainer[particleIndex].pos = (Vec3f){0.0f,0.0f,-20.0f};
			float spread = 1.5f;
			Vec3f maindir = (Vec3f){0.0f, 10.0f, 0.0f};
			// Very bad way to generate a random direction;
			// See for instance http://stackoverflow.com/questions/5408276/python-uniform-spherical-distribution instead,
			// combined with some user-controlled parameters (main direction, spread, etc)
			//glm::vec3 randomdir = glm::vec3(
			Vec3f randomdir = (Vec3f){
				(rand()%2000 - 1000.0f)/1000.0f,
				(rand()%2000 - 1000.0f)/1000.0f,
				(rand()%2000 - 1000.0f)/1000.0f
			};

			ParticlesContainer[particleIndex].speed = maindir + randomdir*spread;
			// Very bad way to generate a random color
			ParticlesContainer[particleIndex].r = rand() % 256;
			ParticlesContainer[particleIndex].g = rand() % 256;
			ParticlesContainer[particleIndex].b = rand() % 256;
			ParticlesContainer[particleIndex].a = (rand() % 256) / 3;
			ParticlesContainer[particleIndex].size = (rand()%1000)/2000.0f + 0.1f;
		}

    ParticlesCount = 0;
	for(int i=0; i<MaxParticles; i++){
		Particle& p = ParticlesContainer[i]; // shortcut
		if(p.life > 0.0f){
			// Decrease life
			p.life -= delta;
			if (p.life > 0.0f){
				p.speed.add( (Vec3f){0.0f,-9.81f, 0.0f} * (float)delta * 0.5f );
				p.pos  .add( p.speed * (float)delta );
				//p.cameradistance = glm::length2( p.pos );

				// Fill the GPU buffer
				g_particule_position_size_data[4*ParticlesCount+0] = p.pos.x;
				g_particule_position_size_data[4*ParticlesCount+1] = p.pos.y;
				g_particule_position_size_data[4*ParticlesCount+2] = p.pos.z;

				g_particule_position_size_data[4*ParticlesCount+3] = p.size;

				g_particule_color_data[4*ParticlesCount+0] = p.r;
				g_particule_color_data[4*ParticlesCount+1] = p.g;
				g_particule_color_data[4*ParticlesCount+2] = p.b;
				g_particule_color_data[4*ParticlesCount+3] = p.a;

			}else{
				// Particles that just died will be put at the end of the buffer in SortParticles();
				p.cameradistance = -1.0f;
			}

			ParticlesCount++;

		}
	}
}


void draw( ){

    glClearColor(0.8, 0.8, 0.8, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );
	// Simulate all particles

	glDisable(GL_ALPHA);

	physics();
	SortParticles();

	//printf("%d ",ParticlesCount);

	// Update the buffers that OpenGL uses for rendering.
	// There are much more sophisticated means to stream data from the CPU to the GPU,
	// but this is outside the scope of this tutorial.
	// http://www.opengl.org/wiki/Buffer_Object_Streaming


	glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
	glBufferData(GL_ARRAY_BUFFER, MaxParticles * 4 * sizeof(GLfloat), NULL, GL_STREAM_DRAW); // Buffer orphaning, a common way to improve streaming perf. See above link for details.
	glBufferSubData(GL_ARRAY_BUFFER, 0, ParticlesCount * sizeof(GLfloat) * 4, g_particule_position_size_data);

	glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
	glBufferData(GL_ARRAY_BUFFER, MaxParticles * 4 * sizeof(GLubyte), NULL, GL_STREAM_DRAW); // Buffer orphaning, a common way to improve streaming perf. See above link for details.
	glBufferSubData(GL_ARRAY_BUFFER, 0, ParticlesCount * sizeof(GLubyte) * 4, g_particule_color_data);

	//glEnable(GL_BLEND);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glDisable(GL_BLEND);

	// Use our shader
	//glUseProgram(programID);
	glUseProgram(shaderParticle->shaderprogram);

	// Bind our texture in Texture Unit 0
	//glActiveTexture(GL_TEXTURE0);
	//glBindTexture(GL_TEXTURE_2D, Texture);
	// Set our "myTextureSampler" sampler to user Texture Unit 0
	//glUniform1i(TextureID, 0);

	// Same as the billboards tutorial
	//glUniform3f(CameraRight_worldspace_ID, ViewMatrix[0][0], ViewMatrix[1][0], ViewMatrix[2][0]);
	//glUniform3f(CameraUp_worldspace_ID   , ViewMatrix[0][1], ViewMatrix[1][1], ViewMatrix[2][1]);
	glUniform3f(CameraRight_worldspace_ID, 1.0f, 0.0f, 0.0f);
	glUniform3f(CameraUp_worldspace_ID   ,  0.0f, 1.0f, 0.0f );
	//glUniformMatrix4fv(ViewProjMatrixID, 1, GL_FALSE, &ViewProjectionMatrix[0][0]);

	static GLfloat VPM[16] = {
	 1.81066, 0.000000, 0.00159, 0.00159255,
     0.00000, 2.414210, 0.00000, 0.00000000,
     0.00288, 0.000000,-1.00200,-0.99999900,
    -0.01440, 0.000000, 4.80980, 4.99999000
    };
    glUniformMatrix4fv(ViewProjMatrixID, 1, GL_FALSE, VPM );


	// 1rst attribute buffer : vertices
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, billboard_vertex_buffer);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0,(void*)0 );

	// 2nd attribute buffer : positions of particles' centers
	glEnableVertexAttribArray(1);
	glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
	glVertexAttribPointer(1,4,GL_FLOAT,GL_FALSE,0,(void*)0);

	// 3rd attribute buffer : particles' colors
	glEnableVertexAttribArray(2);
	glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
	glVertexAttribPointer(2,4,GL_UNSIGNED_BYTE,GL_TRUE,0,(void*)0);

	// These functions are specific to glDrawArrays*Instanced*.
	// The first parameter is the attribute buffer we're talking about.
	// The second parameter is the "rate at which generic vertex attributes advance when rendering multiple instances"
	// http://www.opengl.org/sdk/docs/man/xhtml/glVertexAttribDivisor.xml
	glVertexAttribDivisor(0, 0); // particles vertices : always reuse the same 4 vertices -> 0
	glVertexAttribDivisor(1, 1); // positions : one per quad (its center)                 -> 1
	glVertexAttribDivisor(2, 1); // color : one per quad                                  -> 1
	//glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, nVerts, ParticlesCount);
	glDrawArraysInstanced(GL_TRIANGLE_FAN, 0, nVerts, ParticlesCount);
	//glLineWidth(5); glDrawArraysInstanced(GL_LINE_STRIP, 0, nVerts, ParticlesCount);
    //glDrawArraysInstanced(GL_TRIANGLES, 0, nVerts, ParticlesCount);

	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
	glDisableVertexAttribArray(2);

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
    float keyRotSpeed   = 0.002f;

    const Uint8 *keys = SDL_GetKeyboardState(NULL);

    //if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.roll2  (  (float)keyRotSpeed ); }
	//if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.roll2  ( -(float)keyRotSpeed ); }
	//if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.dpitch2(  (float)keyRotSpeed ); }
	//if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.dpitch2( -(float)keyRotSpeed ); }

    //if( keys[ SDL_SCANCODE_W  ] ){ camPos.add_mul( mouseMat.c, +step ); }
	//if( keys[ SDL_SCANCODE_S  ] ){ camPos.add_mul( mouseMat.c, -step ); }
	//if( keys[ SDL_SCANCODE_A  ] ){ camPos.add_mul( mouseMat.a, -step ); }
	//if( keys[ SDL_SCANCODE_D  ] ){ camPos.add_mul( mouseMat.a, +step ); }

	SDL_Event event;
	while(SDL_PollEvent(&event)){
		if( event.type == SDL_KEYDOWN ){
            switch( event.key.keysym.sym ){
                case SDLK_ESCAPE: quit(); break;
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
