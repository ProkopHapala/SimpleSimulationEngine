
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "Vec2.h"

#include "testUtils.h"
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "AppSDL2OGL.h"

#include "OCL.h"


#define R2SAFE  1.0e-2
#define F2MAX   10.0

//int n = 1024;
//constexpr int n = 32;
//constexpr int n = 256;
float time_step = 0.05;
constexpr int n = 1024;
float damp      = 0.5;
Vec2f pos   [n];
Vec2f vel   [n];
Vec2f force [n];
Vec2f force_[n];

inline void acum_force( const Vec2f& p1, const Vec2f& p2, Vec2f& f ){
    Vec2f dp; dp.set_sub( p2, p1 );
    float ir2 = 1/( dp.norm2() + R2SAFE );
    //float ir  = sqrt(ir2);
    float ir6 = ir2*ir2*ir2;
    float fr  = ir6 - ir6*ir6;
    f.add_mul( dp, fr );
}

void add_ineraction_forces(){
    for(int i=0; i<n; i++){
        Vec2f f; f.set(0.0);
        for(int j=0; j<n; j++){
            acum_force( pos[i], pos[j], f );
        }
        force[i] = f;
    }
}

void clean_force( ){ for(int i=0; i<n; i++){ force[i].set(0.0); } }

void add_confine_force( const Vec2f& center, float strength ){
    for(int i=0; i<n; i++){ force[i].add_mul( pos[i], strength ); }
}

void add_external_force( const Vec2f& center, float strength, float width ){
    float w2 = width*width;
    for(int i=0; i<n; i++){
        Vec2f dp; dp.set_sub( pos[i], center );
        float ir2 = 1/(dp.norm2() + w2);
        float ir  = sqrt(ir2);
        force[i].add_mul( dp, strength*ir2*ir );
    }
}

float move_leap_frog( float dt ){
    float cdamp = 1 - damp*dt;
    double f2max = 0;
    for(int i=0; i<n; i++){
        float f2 = force[i].norm2();
        if(f2>f2max) f2max = f2;
    }
    if(f2max > F2MAX ) dt *= sqrt( sqrt(F2MAX/f2max) );
    for(int i=0; i<n; i++){
        vel[i].mul(cdamp);
        vel[i].add_mul( force[i], dt );
        pos[i].add_mul( vel[i]  , dt );
    }
    return f2max;
}

// ==================== Declaration

class TestApp_clNBody2D : public AppSDL2OGL {
	public:
	int per_frame = 10;

	int err;
	OCLsystem cl;
	OCLtask *task1,*task2;

	//bool use_GPU = false;
	int method = 2;

	//virtual void loop( int nframes );
	virtual void eventHandling( const SDL_Event& event );
	virtual void draw   ();
	//virtual void drawHUD();

	TestApp_clNBody2D( int& id, int WIDTH_, int HEIGHT_ );

};

// ==================== Implementation

TestApp_clNBody2D::TestApp_clNBody2D( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    float span = 5.0;
    for(int i=0; i<n; i++){
        pos[i].set(randf(-span,span),randf(-span,span));
        vel[i].set(0.0);
    }

    // --- OpenCL
    cl.init();
    cl.newBuffer( "pos",   n*2, sizeof(float), (float*)pos,    CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR );
    cl.newBuffer( "vel",   n*2, sizeof(float), (float*)vel,    CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR );
    cl.newBuffer( "force", n*2, sizeof(float), (float*)force_, CL_MEM_READ_WRITE );
    //cl.buffers[1].read_on_finish = true;
    err = cl.buildProgram( "cl/nbody.cl" );                                       OCL_checkError(err, "cl.buildProgram");

    task1 = new OCLtask( &cl, cl.newKernel("NBody_force"), 1, n, 0 );
    task1->args = { INTarg(n), BUFFarg(0), BUFFarg(2) };
    task1->print_arg_list();
    //exit(0);

    task2 = new OCLtask( &cl, cl.newKernel("NBody_step"), 1, n, 0 );
    task2->args = { INTarg(n), BUFFarg(0), BUFFarg(1), FLOATarg(time_step), FLOATarg(damp), INTarg(per_frame) };

    printf("DEBUG = 5\n");

    //char ret[1024];
    //size_t nbytes=0;
    //err = clGetKernelArgInfo ( cl.kernels[0], 0, CL_KERNEL_ARG_TYPE_NAME, 1024 , ret , &nbytes );  printf("CL_KERNEL_ARG_TYPE_NAME  >>%s<<\n", ret );      OCL_checkError(err, "clGetKernelArgInfo");
    //err = clGetKernelArgInfo ( cl.kernels[0], 0, CL_KERNEL_ARG_NAME     , 1024 , ret , &nbytes );  printf("CL_KERNEL_ARG_NAME       >>%s<<\n", ret );      OCL_checkError(err, "clGetKernelArgInfo");

    // --- Test OpenCL force
    //task1->enque();
    task1->enque();
    //cl.finish();
    clFinish(cl.commands); cl.download(2);

    clean_force( );
    add_ineraction_forces();
    for(int i=0;  i<n; i++ ){  printf( "i (%g,%g) (%g,%g) \n", i,  force_[i].x, force_[i].y,   force[i].x, force[i].y ); }

    cl.buffers[2].p_cpu = (float*)force;    // bind Force output array

    //exit(0);

}

void TestApp_clNBody2D::draw(){
    long t1,t2,t3;
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    //delay = 100;
    double f2err;
    t1=getCPUticks();
    switch(method){
        case 3:
            task2->enque();
            clFinish(cl.commands);
            cl.download(0);
            break;
        case 2:
            for(int itr=0; itr<per_frame; itr++){
                cl.upload(0);
                task1->enque();
                clFinish(cl.commands);
                cl.download(2);
                if( LMB ) add_external_force( {mouse_begin_x, mouse_begin_y}, -40.0, 1.0 );
                f2err = move_leap_frog( time_step );
            }
            break;
        case 1:
            for(int itr=0; itr<per_frame; itr++){
                add_ineraction_forces( );
                if( LMB ) add_external_force( {mouse_begin_x, mouse_begin_y}, -40.0, 1.0 );
                f2err = move_leap_frog( time_step );
            }
            break;
    }

	double fticks = getCPUticks() - t1;
    printf( "%i METHOD %i Ferr= %g T= %g [Mtick] %g [t/op]\n", frameCount*per_frame, method, sqrt(f2err), fticks*1e-6, fticks/(n*n) );

	glBegin(GL_POINTS);
    for(int i=0; i<n; i++){  glVertex3f( pos[i].x, pos[i].y, 0.0 ); }
    glEnd();


};


void TestApp_clNBody2D::eventHandling( const SDL_Event& event ){
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_f: use_GPU=!use_GPU; break;
                //case SDLK_g: fullGPU=!fullGPU; break;
                case SDLK_KP_1:  method=1; break;
                case SDLK_KP_2:  method=2; break;
                case SDLK_KP_3:  method=3; break;
            } break;
    };
    AppSDL2OGL::eventHandling( event );
};


// ===================== MAIN

TestApp_clNBody2D * app;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	app = new TestApp_clNBody2D( junk , 800, 600 );
	app->loop( 1000000 );
	return 0;
}
















