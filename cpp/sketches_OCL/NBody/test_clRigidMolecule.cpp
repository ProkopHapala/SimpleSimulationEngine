
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "testUtils.h"
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"

#include "RigidMolecule.h"
#include "OCL.h"

// ==================== Declaration

class TestApp_clRigidMolecule : public AppSDL2OGL_3D {
	public:
	int per_frame = 10;

	int err;
	OCLsystem cl;
	OCLtask *task1,*task2;

	//bool use_GPU = false;
	int method = 2;

	int atomView = 0;

	//virtual void loop( int nframes );
	virtual void eventHandling( const SDL_Event& event );
	virtual void draw   ();
	//virtual void drawHUD();

	TestApp_clRigidMolecule( int& id, int WIDTH_, int HEIGHT_ );

};

// ==================== Implementation

TestApp_clRigidMolecule::TestApp_clRigidMolecule( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    atomView=glGenLists(1);
	glNewList( atomView, GL_COMPILE );
		//glDisable ( GL_LIGHTING );
		//glColor3f ( 1.0f, 1.0f, 1.0f );
		Draw3D::drawSphere_oct( 3, 0.95f, {0.0d,0.0d,0.0d} );
		//Draw3D::drawSphere_oct(4,1.0,{0.0,0.0,0.0});
	glEndList();

	initParticles( nMols, pos, 6.0, 0.0 );

    /*
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

    */

}

void TestApp_clRigidMolecule::draw(){
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
    }

	double fticks = getCPUticks() - t1;
    //printf( "%i METHOD %i Ferr= %g T= %g [Mtick] %g [t/op]\n", frameCount*per_frame, method, sqrt(f2err), fticks*1e-6, fticks/(n*n) );

    glDisable(GL_LIGHTING);
    transformAllAtoms( nMols, nAtoms, pos, molAtoms, atomsT );
    Mat3f mrot; mrot.setOne();
    glColor3f(1.0f,0.0f,0.0f);
    for(int i=0; i<nMols; i++){
        //printf( " %i (%g,%g,%g)  \n", i, atomsT[i].x, atomsT[i].y, atomsT[i].z );
        //Draw3D::drawShape(atomsT[i],mrot,atomView);
        Draw3D::drawPointCross(*((Vec3f*)(pos+(i<<3))),0.1);
    }
    glEnable(GL_LIGHTING);
    //glColor3f(0.0f,1.0f,0.0f);
    glColor3f(0.8f,0.8f,0.8f);
    for(int i=0; i<nMols*nAtoms; i++){
        //printf( " %i (%g,%g,%g)  \n", i, atomsT[i].x, atomsT[i].y, atomsT[i].z );
        Draw3D::drawShape(atomsT[i],mrot,atomView);
        //Draw3D::drawPointCross(atomsT[i],0.1);
    }
    //exit(0);

};


void TestApp_clRigidMolecule::eventHandling( const SDL_Event& event ){
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

TestApp_clRigidMolecule * app;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	app = new TestApp_clRigidMolecule( junk , 800, 600 );
	app->loop( 1000000 );
	return 0;
}
















