
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "Vec2.h"

#include "testUtils.h"
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw2D.h"
#include "AppSDL2OGL.h"

#include "NBody2D.h"

void drawAtomsByCell( Vec2f* pos, int* c2p ){
    //glBegin(GL_POINTS);
    for(int iy=0; iy<ny; iy++){ for(int ix=0; ix<nx; ix++){ int i = nx*iy+ix;
    //for(int i=0; i<ncell; i++){
            int iStart = c2p[i  ];
            int iEnd   = c2p[i+1];
            if(iStart<iEnd){
            //printf( "(%i,%i) %i (%i,%i) \n", ix, iy, i, iStart, iEnd );
            Draw::setRGB( rand_hash2( i + 15464) );
            for(int i=iStart; i<iEnd; i++){
                //glVertex3f( pos[i].x, pos[i].y, 0.0 );
                Draw2D::drawCircle(pos[i],0.5,8,false);
            }
            }
    //}
    } }
    //glEnd();
    //exit(0);
}

// ==================== Declaration

class TestApp_clNBody2DTiled : public AppSDL2OGL {
	public:
	int per_frame = 1;

	int err;
	OCLsystem cl;
	OCLtask *task1,*task2;

	//bool use_GPU = false;
	int method = 1;

	//virtual void loop( int nframes );
	virtual void eventHandling( const SDL_Event& event );
	virtual void draw   ();
	//virtual void drawHUD();

	TestApp_clNBody2DTiled( int& id, int WIDTH_, int HEIGHT_ );

};

// ==================== Implementation

TestApp_clNBody2DTiled::TestApp_clNBody2DTiled( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    /*
    float span = 5.0;
    for(int i=0; i<n; i++){
        pos[i].set(randf(-span,span),randf(-span,span));
        vel[i].set(0.0);
    }
    */

    initParticles( 1.5, 0.0 );

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

void TestApp_clNBody2DTiled::draw(){
    long t1,t2,t3;
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable(GL_DEPTH_TEST);



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
                if( LMB ) add_external_force( {mouse_begin_x, mouse_begin_y}, -4.0, 1.0 );
                f2err = move_leap_frog( time_step );
            }
            break;
        case 1:
            for(int itr=0; itr<per_frame; itr++){
                atomsToCells( );
                //add_ineraction_forces( );

                clean_force( );
                add_ineraction_forces_cells( pos, cell2pos );

                if( LMB ) add_external_force( {mouse_begin_x, mouse_begin_y}, -4.0, 1.0 );
                f2err = move_leap_frog( time_step );
            }
            break;
    }


	double fticks = getCPUticks() - t1;
    printf( "%i METHOD %i Ferr= %g T= %g [Mtick] %g [t/op]\n", frameCount*per_frame, method, sqrt(f2err), fticks*1e-6, fticks/(n*n) );



    /*
    printf("frame %i \n", frameCount);
    glColor3f(0.48,0.48,0.48);
    Draw2D::drawGrid( grid_orig[0], grid_orig[1], grid_orig[0]+cell_size*nx+0.1, grid_orig[0]+cell_size*ny+0.1, cell_size, cell_size );

    atomsToCells( );
    glColor3f(0.5,0.45,0.45);
    add_ineraction_forces_cells( pos_, cell2pos );
    drawAtomsByCell( pos_, cell2pos );
    */

    //exit(0);

    drawAtomsByCell( pos_, cell2pos );

    /*
    glColor3f(0.9,0.9,0.9);
	glBegin(GL_POINTS);
    for(int i=0; i<n; i++){
        //printf( "%i %i\n", i, atom2cell[i] );
        Draw::setRGB( rand_hash2( atom2cell[i]+15454 ) );
        glVertex3f( pos[i].x, pos[i].y, 0.0 );
    }
    glEnd();
    */


    //exit(0);


};


void TestApp_clNBody2DTiled::eventHandling( const SDL_Event& event ){
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

TestApp_clNBody2DTiled * app;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	app = new TestApp_clNBody2DTiled( junk , 800, 600 );
	app->loop( 1000000 );
	return 0;
}
















