
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "fastmath.h"
#include "Vec3.h"

#include "testUtils.h"
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw2D.h"
#include "AppSDL2OGL.h"

#include "spaceFrame2D.h"
#include "OCL.h"

#define R2SAFE  1.0e-2
#define F2MAX   10.0

void drawSquareFrame(int nx, int ny, float * pos){
    for(int iy=0;iy<ny;iy++){
        glBegin(GL_LINE_STRIP);
        for(int ix=0;ix<nx;ix++){
            int i4  = (iy*nx + ix)<<2;
            glVertex3f( pos[i4], pos[i4+1], pos[i4+2] );
            //printf( "(%i,%i) (%g,%g,%g)\n", ix,iy, pos[i4], pos[i4+1], pos[i4+2] );
        }
        glEnd();
    }
    for(int ix=0;ix<nx;ix++){
        glBegin(GL_LINE_STRIP);
        for(int iy=0;iy<ny;iy++){
            int i4  = (iy*nx + ix)<<2;
            glVertex3f( pos[i4], pos[i4+1], pos[i4+2] );
            //printf( "(%i,%i) (%g,%g,%g)\n", ix,iy, pos[i4], pos[i4+1], pos[i4+2] );
        }
        glEnd();
    }
}

void drawVecs(int n, float * pos, float * vs, float sc){
    glBegin(GL_LINES);
    for(int i=0;i<n;i++){
        int i4  = i<<2;
        glVertex3f(pos[i4]          ,pos[i4+1]            ,pos[i4+2]);
        glVertex3f(pos[i4]+vs[i4]*sc,pos[i4+1]+vs[i4+1]*sc,pos[i4+2]+vs[i4+2]*sc);
    }
    glEnd();
}


class TestApp_clConvolve2D : public AppSDL2OGL {
	public:
	int per_frame = 64;

	int err;
	OCLsystem cl;
	OCLtask *task1,*task2;
	OCLtask *taskMove, *taskConstr;

	//bool use_GPU = false;
	int method = 1;

	GLuint texture1 = 0;

	int i_picked = -1;

	float dt   = 0.1;
	float damp = 0.999;
	float K_constr = 10.0;

	//virtual void loop( int nframes );
	virtual void eventHandling( const SDL_Event& event );
	virtual void draw   ();
	//virtual void drawHUD();

	void switchMethod( int i );

	TestApp_clConvolve2D( int& id, int WIDTH_, int HEIGHT_ );

};

// ==================== Implementation

TestApp_clConvolve2D::TestApp_clConvolve2D( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    genPos0(nx, ny, pos, l0*1.1, 0.1 );

    // --- OpenCL
    cl.init();
    cl.newBuffer( "pos",     ntot, sizeof(float)*4, pos  , CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR );
    cl.newBuffer( "vel",     ntot, sizeof(float)*4, vel  , CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR );
    cl.newBuffer( "force",   ntot, sizeof(float)*4, force, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR );

    cl.newBuffer( "iConstrains",   nConstrMax, sizeof(int),     iConstrains, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR );
    cl.newBuffer( "constrains" ,   nConstrMax, sizeof(float)*4,  constrains, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR );
    err = cl.buildProgram( "cl/cloth.cl" );     OCL_checkError(err, "cl.buildProgram");

    task1 = new OCLtask( &cl, cl.newKernel("cloth_force"), 2, 0, 0 );
    task1->global[0] = nx; task1->global[1] = ny;
    task1->local [0] = 16; task1->local [1] = 16;
    task1->args = { BUFFarg(0), BUFFarg(2) };
    task1->print_arg_list();

    task2 = new OCLtask( &cl, cl.newKernel("cloth_dynamics"), 2, 0, 0 );
    task2->global[0] = nx; task2->global[1] = ny;
    task2->local [0] = 16; task2->local [1] = 16;
    task2->args = { BUFFarg(0), BUFFarg(1), FLOATarg(dt), FLOATarg(damp) };
    task2->print_arg_list();

    taskMove = new OCLtask( &cl, cl.newKernel("move_leapfrog"), 1, ntot-nx, 32 );
    taskMove->args = { BUFFarg(0), BUFFarg(1), BUFFarg(2), FLOATarg(dt), FLOATarg(damp) };
    taskMove->print_arg_list();

    taskConstr = new OCLtask( &cl, cl.newKernel("harmonic_constr"), 1, 1, 0 );
    taskConstr->args = { BUFFarg(3), BUFFarg(4), BUFFarg(0),  BUFFarg(2) };
    taskConstr->print_arg_list();
    //exit(0);
    //switchMethod( 1 );

    set_array(4*ntot, force , 0.0f);
    set_array(4*ntot, force_, 0.0f);

    printf( " --- Test force OpenCL N2 \n" );
    cl.upload(0);
    task1->enque();
    clFinish(cl.commands);
    cl.download(2);
    evalForce(nx, ny, pos, force_ );
    //for(int i=0; i<ntot; i++ ){ int i4=i<<2; printf( "%i (%g,%g,%g) (%g,%g,%g) \n", i, force[i4], force[i4+1], force[i4+2], force_[i4], force_[i4+1], force_[i4+2] ); }
    checkDiff( ntot, force, force_ );

    //exit(0);

    camX0 = 32.0f;
    camY0 = 32.0f;
    zoom  = 48.0f;

}

void TestApp_clConvolve2D::draw(){
    long t1,t2,t3;
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    //delay = 100;
    double f2err;
    t1=getCPUticks();
    switch(method){
        case 1:
            for(int itr=0; itr<per_frame; itr++){
                evalForce(nx, ny, pos, force );
                if( i_picked >= 0 ){
                    int i4 = (i_picked<<2);
                    addHarmonicForce( *((Vec3f*)(pos+i4)), {mouse_begin_x,mouse_begin_y,0.0}, *((Vec3f*)(force+i4)), K_constr );
                }
                move_leapfrog( ntot-nx, pos, vel, force, dt, damp );
            }
            break;
        case 2:
            for(int itr=0; itr<per_frame; itr++){
                cl.upload(0);
                task1->enque();
                clFinish(cl.commands);
                cl.download(2);
                if( i_picked >= 0 ){
                    int i4 = (i_picked<<2);
                    addHarmonicForce( *((Vec3f*)(pos+i4)), {mouse_begin_x,mouse_begin_y,0.0}, *((Vec3f*)(force+i4)), K_constr );
                }
                move_leapfrog( ntot-nx, pos, vel, force, dt, damp );
            }
            break;
        case 3:
            cl.upload(0);
            cl.upload(1);
            //cl.upload(1);
            taskConstr->global[0] = 1;
            iConstrains       [0] = i_picked;
            constrains        [0] = mouse_begin_x;
            constrains        [1] = mouse_begin_y;
            constrains        [2] = 0.0;
            constrains        [3] = K_constr;
            cl.upload(3);
            cl.upload(4);
            for(int itr=0; itr<per_frame; itr++){
                task1     ->enque();
                if(i_picked >= 0)taskConstr->enque();
                taskMove  ->enque();
            }
            clFinish(cl.commands);
            cl.download(0);
            cl.download(1);
            break;
        case 4:
            cl.upload(0);
            cl.upload(1);
            for(int itr=0; itr<per_frame; itr++){
                task2   ->enque();
            }
            clFinish(cl.commands);
            cl.download(0);
            cl.download(1);
            break;
    }
    t1 = getCPUticks() - t1;

    printf( " %i method %i Time %g [Mticks] %g [ticks/pix]\n", frameCount, method,  t1*1e-6, ((double)t1)/(per_frame*ntot) );

    glColor3f(0.0f,0.0f,0.0f); drawSquareFrame( nx, ny, pos);
    //glColor3f(1.0f,0.0f,0.0f); drawVecs(ntot, pos, force, 10.0);
    //exit(0);

};


void TestApp_clConvolve2D::switchMethod( int i ){
    method=i;
    //genXOR2D(nx, ny, buff);
    //genZero(nx, ny, buff);
    //cl.upload(0);
}

void TestApp_clConvolve2D::eventHandling( const SDL_Event& event ){
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_KP_1:  switchMethod( 1 ); break;
                case SDLK_KP_2:  switchMethod( 2 ); break;
                case SDLK_KP_3:  switchMethod( 3 ); break;
                case SDLK_KP_4:  switchMethod( 4 ); break;
            } break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    i_picked = findNearest( {mouse_begin_x,mouse_begin_y,0.0}, ntot, pos );
                    break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    i_picked = -1;
                    break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
};


// ===================== MAIN

TestApp_clConvolve2D * app;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	app = new TestApp_clConvolve2D( junk , 800, 600 );
	app->loop( 1000000 );
	return 0;
}
















