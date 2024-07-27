
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

#include "OCL.h"

void drawAtomsByCell(int n, Vec2f* pos, int* c2p ){
    //glBegin(GL_POINTS);
    for(int iy=0; iy<ny; iy++){ for(int ix=0; ix<nx; ix++){ int i = nx*iy+ix;
    //for(int i=0; i<ncell; i++){
            int iStart = c2p[i  ];
            int iEnd   = c2p[i+1];
            if(iStart<iEnd){
            //printf( "(%i,%i) %i (%i,%i) \n", ix, iy, i, iStart, iEnd );
            if(c2p) Draw::setRGB( rand_hash2( i + 15464) );
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
	int per_frame = 10;

	int err;
	OCLsystem cl;
	OCLtask *taskN2,*taskTiled,*taskSorted;

	//bool use_GPU = false;
	int method = 3;
	int nloc   = 16;
	int niters = 0;
    double t_interact = 0;
    double t_arrange  = 0;
    double t_move     = 0;

	//virtual void loop( int nframes );
	virtual void eventHandling( const SDL_Event& event );
	virtual void draw();
	virtual void quit();
	//virtual void drawHUD();

    void switchMethod( int i );

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

    //Vec2i v2i_test[10];
    //for(int i=0; i<10; i++){ v2i_test[i].set(i,10-i); }
    //for(int i=0; i<10; i++){ printf( "%i   : %i %i \n", i, v2i_test[i].x, v2i_test[i].y ); }
    //exit(0);

    initParticles( 1.5, 0.25 );

    // --- OpenCL
    printf(" ::: Init OpenCL .... \n");
    cl.initOCL();

    // http://stackoverflow.com/questions/5237181/is-there-a-limit-to-opencl-local-memory
    //cl_ulong local_mem_size;
    //clGetDeviceInfo(cl.device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &local_mem_size, 0); printf("local_mem_size %li \n",local_mem_size); exit(0);

    cl.newBuffer( "pos",   n*2, sizeof(float), pos,    CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR );
    cl.newBuffer( "vel",   n*2, sizeof(float), vel,    CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR );
    cl.newBuffer( "force", n*2, sizeof(float), force_, CL_MEM_READ_WRITE );

    cl.newBuffer( "bounds", ncell*2, sizeof(int), cellBounds,     CL_MEM_READ_ONLY  );
    cl.newBuffer( "order",  ncell,   sizeof(int), sortedCells,    CL_MEM_READ_ONLY  );
    //cl.newBuffer( "debug",      ncell,   sizeof(float), (float*)DEBUG_int_buff, CL_MEM_WRITE_ONLY );
    //cl.newBuffer( "neighCells", 9,       sizeof(float), (float*)neighCells,     CL_MEM_READ_ONLY  );

    //cl.buffers[1].read_on_finish = true;
    err = cl.buildProgram( "cl/nbody_SR.cl" );                                       OCL_checkError(err, "cl.buildProgram");
    taskN2 = new OCLtask( &cl, cl.newKernel("NBody_force"), 1, n, 64 );
    taskN2->args = { INTarg(n), BUFFarg(0), BUFFarg(2) };
    taskN2->print_arg_list();
    //exit(0);

    nloc = 32;
    taskTiled = new OCLtask( &cl, cl.newKernel("force_Tiled"), 1, ncell*nloc, nloc );
    taskTiled->args = { INTarg(nx), BUFFarg(0), BUFFarg(2), BUFFarg(3) };

    taskSorted = new OCLtask( &cl, cl.newKernel("force_Tiled_Sorted"), 1, ncell*nloc, nloc );
    taskSorted->args = { INTarg(nx), BUFFarg(0), BUFFarg(2), BUFFarg(3), BUFFarg(4) };

    printf(" ::: Init OpenCL .... DONE \n");

    //char ret[1024];
    //size_t nbytes=0;
    //err = clGetKernelArgInfo ( cl.kernels[0], 0, CL_KERNEL_ARG_TYPE_NAME, 1024 , ret , &nbytes );  printf("CL_KERNEL_ARG_TYPE_NAME  >>%s<<\n", ret );      OCL_checkError(err, "clGetKernelArgInfo");
    //err = clGetKernelArgInfo ( cl.kernels[0], 0, CL_KERNEL_ARG_NAME     , 1024 , ret , &nbytes );  printf("CL_KERNEL_ARG_NAME       >>%s<<\n", ret );      OCL_checkError(err, "clGetKernelArgInfo");

    printf( " --- Test Force tiles \n" );
    atomsToCells();
    sortIntersBins();
    //for(int i=0; i<ncell; i++){ printf("%i : %i %i %i \n", i, cellBounds[i].x, cellBounds[i].y, DEBUG_int_buff[i] ); }; exit(0);
    clean_array( n, force  ); add_ineraction_forces      ( n, pos, force );
    clean_array( n, force_ ); add_ineraction_forces_cells( n, pos, force_, cell2pos );
    printf("check CPU-Tiled vs CPU-naieve\n"); checkDiff ( n, force_, force );
    clean_array( n, force_ ); add_forces_Inters          ( n, pos, force_ );
    printf("check CPU-Inter vs CPU-naieve\n"); checkDiff ( n, force_, force );



    printf( " --- Test force OpenCL N2 \n" );
    cl.upload(0); taskN2->enque(); // cl.finish();
    clFinish(cl.commands); cl.download(2);
    //for(int i=0;  i<n; i++ ){  printf( "i (%g,%g) (%g,%g) \n", i,  force_[i].x, force_[i].y,   force[i].x, force[i].y ); }
    checkDiff( n, force_, force );

    printf( " --- Test force OpenCL Tiled \n" );
    atomsToCells( );
    prepareCellBounds();
    clean_array( n, force_ );
    cl.upload(0); cl.upload(2); cl.upload(3); //cl.upload(5);
    taskTiled->enque();
    clFinish(cl.commands);
    cl.download(2); //cl.download(4);
    //for(int i=0; i<ncell; i++){ printf("%i %i %i \n", i, DEBUG_int_buff[i], DEBUG_int_buff_[i] ); }
    //for(int i=0;  i<n; i++ ){  printf( "%i (%g,%g) (%g,%g) \n", i,  force_[i].x, force_[i].y,   force[i].x, force[i].y ); }
    checkDiff( n, force_, force );
    //exit(0);

    printf( " --- Test force OpenCL Sorted \n" );
    atomsToCells( );
    sortCells();
    printf("nFilledCells %i\n", nFilledCells );
    taskSorted->global.x = nloc*nFilledCells;
    taskSorted->local .x = nloc;
    //taskTiled->global.x = nloc*nFilledCells;
    //taskTiled->local .x = nloc;
    clean_array( n, force_ );
    cl.upload(0); cl.upload(2); cl.upload(3); cl.upload(4); //cl.upload(5);
    taskSorted->enque();
    //taskTiled->enque();
    clFinish(cl.commands);
    cl.download(2); //cl.download(4);
    //for(int i=0; i<ncell; i++){ printf("%i %i %i \n", i, DEBUG_int_buff[i], DEBUG_int_buff_[i] ); }
    //for(int i=0;  i<n; i++ ){  printf( "%i (%g,%g) (%g,%g) \n", i,  force_[i].x, force_[i].y,   force[i].x, force[i].y ); }
    checkDiff( n, force_, force );
    //exit(0);

    cl.buffers[2].p_cpu = (float*)force;    // bind Force output array

    //exit(0);

}

void TestApp_clNBody2DTiled::draw(){
    long t1,t2,t3;
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable(GL_DEPTH_TEST);

    float touchStrength =  -0.5;
    float touchRadius   =   3.0;

    //delay = 100;
    double f2err;
    t1=getCPUticks();
    switch(method){
        case 1: // ---- N^2 implementation on CPU
            for(int itr=0; itr<per_frame; itr++){
                niters++;
                //atomsToCells( );
                long t_1 =  getCPUticks();
                    add_ineraction_forces( n, pos, force );
                long t_2 = getCPUticks();
                    if( LMB ) add_external_force( n, force, {mouse_begin_x, mouse_begin_y}, touchStrength, touchRadius );
                    f2err = move_leap_frog( n, pos, vel, force, time_step );
                long t_3 = getCPUticks();
                t_interact+=(t_2-t_1); t_move+=(t_3-t_2);
            }
            break;
        case 2: // ---- Tilled implementation on CPU
            for(int itr=0; itr<per_frame; itr++){
                niters++;
                long t_1 =  getCPUticks();
                    atomsToCells();
                    clean_array(n,force);
                long t_2 = getCPUticks();
                    add_ineraction_forces_cells( n, pos, force, cell2pos );
                    if( LMB ) add_external_force( n, force, {mouse_begin_x, mouse_begin_y}, touchStrength, touchRadius );
                long t_3 = getCPUticks();
                    f2err = move_leap_frog( n, pos, vel, force, time_step );
                long t_4 = getCPUticks();
                t_arrange+=(t_2-t_1); t_interact+=(t_3-t_2); t_move+=(t_4-t_3);
            }
            break;
        case 3: // ---- Inters implementation on CPU
            for(int itr=0; itr<per_frame; itr++){
                niters++;
                long t_1 =  getCPUticks();
                    atomsToCells();
                    sortIntersBins();
                    clean_array(n,force);
                long t_2 = getCPUticks();
                    add_forces_Inters( n, pos, force );
                    if( LMB ) add_external_force( n, force, {mouse_begin_x, mouse_begin_y}, touchStrength, touchRadius );
                long t_3 = getCPUticks();
                    f2err = move_leap_frog( n, pos, vel, force, time_step );
                long t_4 = getCPUticks();
                t_arrange+=(t_2-t_1); t_interact+=(t_3-t_2); t_move+=(t_4-t_3);
            }
            break;
        case 4:// ---- N^2 implementation on GPU
            for(int itr=0; itr<per_frame; itr++){
                niters++;
                long t_1 =  getCPUticks();
                    cl.upload(0);
                    taskN2->enque();
                    clFinish(cl.commands);
                    cl.download(2);
                long t_2 = getCPUticks();
                    if( LMB ) add_external_force( n, force, {mouse_begin_x, mouse_begin_y}, touchStrength, touchRadius );
                    f2err = move_leap_frog( n, pos, vel, force, time_step );
                long t_3 = getCPUticks();
                t_interact+=(t_2-t_1); t_move+=(t_3-t_2);
            }
            break;
        case 5:// ---- Tiled implementation on GPU
            for(int itr=0; itr<per_frame; itr++){
                niters++;
                long t_1 =  getCPUticks();
                    atomsToCells();
                    prepareCellBounds();
                long t_2 =  getCPUticks();
                    cl.upload(0); cl.upload(3);// cl.upload(5);
                    taskTiled->enque();
                    clFinish(cl.commands);
                    cl.download(2);// cl.download(4);
                long t_3 = getCPUticks();
                    if( LMB ) add_external_force( n, force, {mouse_begin_x, mouse_begin_y}, touchStrength, touchRadius );
                    f2err = move_leap_frog( n, pos, vel, force, time_step );
                long t_4 = getCPUticks();
                t_arrange+=(t_2-t_1); t_interact+=(t_3-t_2); t_move+=(t_4-t_3);
            }
            break;
        case 6:// ---- Tiled implementation on GPU
            for(int itr=0; itr<per_frame; itr++){
                niters++;
                long t_1 =  getCPUticks();
                    atomsToCells();
                    prepareCellBounds();
                    sortCells();
                long t_2 =  getCPUticks();
                    cl.upload(0); cl.upload(3); cl.upload(4);
                    taskSorted->global.x = nloc*nFilledCells;
                    taskSorted->local .y = nloc;
                    taskSorted->enque();
                    clFinish(cl.commands);
                    cl.download(2);// cl.download(4);
                long t_3 = getCPUticks();
                    if( LMB ) add_external_force( n, force, {mouse_begin_x, mouse_begin_y}, touchStrength, touchRadius );
                    f2err = move_leap_frog( n, pos, vel, force, time_step );
                long t_4 = getCPUticks();
                t_arrange+=(t_2-t_1); t_interact+=(t_3-t_2); t_move+=(t_4-t_3);
            }
            break;
    }

    glColor3f(0.5f,0.0f,0.5f); Draw2D::drawCircle( {mouse_begin_x, mouse_begin_y}, touchRadius, 16, false );

	double fticks = getCPUticks() - t1;
    //printf( "%i METHOD %i Ferr= %g T= %g [Mtick] %g [t/op]\n", frameCount*per_frame, method, sqrt(f2err), fticks*1e-6, fticks/(n*n) );
    printf( "%i METHOD %i T %g (%g,%g,%g) [Mticks]\n", frameCount, method, (t_arrange+t_interact+t_move)*1e-6/niters,  t_arrange*1e-6/niters, t_interact*1e-6/niters, t_move*1e-6/niters );



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

    Vec2f a = {1.0,1.0};

    glColor3f(0.48,0.48,0.48); Draw2D::drawGrid( grid_orig[0], grid_orig[1], grid_orig[0]+cell_size*nx+0.1, grid_orig[0]+cell_size*ny+0.1, cell_size, cell_size );
    glColor3f(0.15,0.15,0.15); Draw2D::drawRectangle( {-xspan,-yspan}, {xspan,yspan}, false );
    //drawAtomsByCell( pos_, cell2pos );
    glColor3f(0.15,0.15,0.15); for(int i=0; i<n; i++){  Draw2D::drawCircle(pos[i],0.5,8,false); };

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


void TestApp_clNBody2DTiled::switchMethod( int i ){
    method=i;
    niters     = 0;
    t_interact = 0;
    t_arrange  = 0;
    t_move     = 0;
}

void TestApp_clNBody2DTiled::quit( ){
    printf("    n  nx  ny   ncell nloc niters CellSz  t_arrange t_interact t_move\n");
    printf("%5i %3i %3i %6i %6i   %2i  %2.3f  %5.5f   %5.5f   %5.5f\n", n, nx, ny, ncell, nloc, niters, cell_size, t_arrange*1e-6/niters, t_interact*1e-6/niters, t_move*1e-6/niters );
    exit(0);
}

void TestApp_clNBody2DTiled::eventHandling( const SDL_Event& event ){
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_f: use_GPU=!use_GPU; break;
                //case SDLK_g: fullGPU=!fullGPU; break;
                case SDLK_KP_1:  switchMethod( 1 ); break;
                case SDLK_KP_2:  switchMethod( 2 ); break;
                case SDLK_KP_3:  switchMethod( 3 ); break;
                case SDLK_KP_4:  switchMethod( 4 ); break;
                case SDLK_KP_5:  switchMethod( 5 ); break;
                case SDLK_KP_6:  switchMethod( 6 ); break;
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
















