
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "testUtils.h"
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"

#include "RigidMolecule.h"
#include "OCL.h"

// ==================== Declaration

void checkDiff( int n,  float * a_,  Vec3f * a ){
    for(int i=0; i<n; i++){
        int i4 = i<<2;
        printf( "%i (%g,%g,%g) (%g,%g,%g) \n",i,  a_[i4],a_[i4+1],a_[i4+2], a[i].x,a[i].y,a[i].z );
        float err2 = sq(a_[i4]-a[i].x) + sq(a_[i4+1]-a[i].y)  + sq(a_[i4+2]-a[i].z);
        if(err2 > 1e-8){
            printf(" incorrect force on molecule %i err2 = %g \n", i, err2 );
            exit(0);
        }
    }
}

class TestApp_clRigidMolecule : public AppSDL2OGL_3D {
	public:
	int per_frame = 5;

	int err;
	OCLsystem cl;
	OCLtask *task1,*task2;

	float dt   = 0.01;
	float damp = 0.998;

	//bool use_GPU = false;
	int method = 1;

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


	// test quaternion
	/*
	Quat4f q;
	Mat3f  m;
	Vec3f  v,vT,vT_;
	q.fromUniformS3({randf(),randf(),randf()});
	q.toMatrix(m);
	v.set(randf(),randf(),randf()); printVec(v  );
	printf( "q.norm2() %g \n", q.norm2() );

	m.dot_to_T    (v,  vT  );   printVec(vT );
	//m.dot_to    (v,  vT  );
	//m.dot_to     (vT, vT_ );
	q.transformVec  (v,  vT_); printVec(vT );
	q.untransformVec(vT, vT_); printVec(vT_);
    exit(0);
    */

	initParticles( nMols, pos, 4.0, 0.0 );
	setArray     ( nMols*8, vel,   0.0f );
	setArray     ( nMols*8, force, 0.0f );

	for(int i=0; i<nAtoms; i++){
        int i4 = i<<2;
        molAtoms_[i4  ]=molAtoms[i].x;
        molAtoms_[i4+1]=molAtoms[i].y;
        molAtoms_[i4+2]=molAtoms[i].z;
        molAtoms_[i4+3]=0.0f;
        printf("%i (%g,%g,%g)\n", i, molAtoms[i].x,molAtoms[i].y,molAtoms[i].z);
	}

    // --- OpenCL
    cl.init();
    cl.newBuffer( "molAtoms",   nAtoms*4,       sizeof(float), (float*)molAtoms_, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR );
    cl.newBuffer( "pos",        nMols*8,        sizeof(float), (float*)pos,       CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR );
    cl.newBuffer( "force",      nMols*8,        sizeof(float), (float*)force,     CL_MEM_READ_WRITE );
    cl.newBuffer( "atomsT",     nAtoms*nMols*4, sizeof(float), (float*)atomsT_,   CL_MEM_READ_WRITE );
    //cl.buffers[1].read_on_finish = true;
    err = cl.buildProgram( "cl/RigidMolecule.cl" );                                       OCL_checkError(err, "cl.buildProgram");

    /*
    //  check transform
    task1 = new OCLtask( &cl, cl.newKernel("RigidMolTransform"), 1, nMols*nAtoms, nAtoms );
    task1->args = { INTarg(nAtoms), BUFFarg(0), BUFFarg(1), BUFFarg(3) };
    task1->print_arg_list();
    printf( "==== transformAllAtoms CPU \n" );
    transformAllAtoms( nMols, nAtoms, pos, molAtoms, atomsT );
    printf( "==== transformAllAtoms GPU \n" );
    cl.upload(0); cl.upload(1);
    task1->enque();
    clFinish(cl.commands);
    cl.download(3);
    printf( "==== check results \n" );
    checkDiff( nMols*nAtoms,  atomsT_, atomsT );
    exit(0);
    */

    /*
    //  check atom-wise forces
    task1 = new OCLtask( &cl, cl.newKernel("RigidMolAtomForce"), 1, nMols*nAtoms, nAtoms );
    task1->args = { INTarg(nAtoms), BUFFarg(0), BUFFarg(1), BUFFarg(3) };
    task1->print_arg_list();
    printf( "==== transformAllAtoms CPU \n" );
    transformAllAtoms ( nMols, nAtoms, pos, molAtoms, atomsT );
    setArray          ( nMols*8, force_, 0.0f );
    RBodyForce        ( nMols, nAtoms, pos, force_, atomsT, molAtoms );
    printf( "==== transformAllAtoms GPU \n" );
    cl.upload(0); cl.upload(1);
    task1->enque();
    clFinish(cl.commands);
    cl.download(3);
    printf( "==== check results \n" );
    checkDiff( nMols*nAtoms,  atomsT_, forceAtomT );
    exit(0);
    */


    task1 = new OCLtask( &cl, cl.newKernel("RigidMolForce"), 1, nMols*nAtoms, nAtoms );
    task1->args = { INTarg(nAtoms), BUFFarg(0), BUFFarg(1), BUFFarg(2) };
    task1->print_arg_list();

    setArray         ( nMols*8, force , 0.0f );
    setArray         ( nMols*8, force_, 0.0f );
    transformAllAtoms( nMols, nAtoms, pos, molAtoms, atomsT );
    RBodyForce       ( nMols, nAtoms, pos, force_, atomsT, molAtoms );

    cl.upload(0); cl.upload(1);
    task1->enque();
    clFinish(cl.commands);
    cl.download(2);

    for(int i=0; i<nMols; i++){
        int i8 = i<<3;
        printf( "%i (%g,%g,%g,%g, %g,%g,%g,%g) (%g,%g,%g,%g, %g,%g,%g,%g) \n",i,
            force [i8],force [i8+1],force [i8+2],force [i8+3],force [i8+4],force [i8+5],force [i8+6],force [i8+7],
            force_[i8],force_[i8+1],force_[i8+2],force_[i8+3],force_[i8+4],force_[i8+5],force_[i8+6],force_[i8+7]  );
        float err2 = sq(force [i8  ]-force_[i8  ]) + sq(force [i8+1]-force_[i8+1]) + sq(force [i8+2]-force_[i8+2]) + sq(force [i8+3]-force_[i8+3])
                   + sq(force [i8+4]-force_[i8+4]) + sq(force [i8+5]-force_[i8+5]) + sq(force [i8+6]-force_[i8+6]) + sq(force [i8+7]-force_[i8+7]);
        if(err2 > 1e-8){
            printf(" incorrect force on molecule %i err2 = %g \n", i, err2 );
            exit(0);
        }
    }


    //exit(0);


}

void TestApp_clRigidMolecule::draw(){
    long t1,t2,t3;
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    //delay = 100;
    double f2err;
    t1=getCPUticks();
    switch(method){
        case 1:
            for(int itr=0;itr<per_frame; itr++){
                transformAllAtoms( nMols, nAtoms, pos, molAtoms, atomsT );
                setArray      ( nMols*8, force, 0.0f );
                RBodyForce    ( nMols, nAtoms, pos, force, atomsT, molAtoms );
                move_leap_frog( nMols, pos, vel, force, dt, damp );
                //task2->enque();
                //clFinish(cl.commands);
                //cl.download(0);
            }
            break;
        case 2:
            for(int itr=0;itr<per_frame; itr++){
                cl.upload(0); cl.upload(1);
                task1->enque();
                clFinish(cl.commands);
                cl.download(2);
                move_leap_frog( nMols, pos, vel, force, dt, damp );
            }
            break;
    }
    double fticks = getCPUticks() - t1;
    printf( "%i METHOD %i Ferr= %g T= %g [Mtick] %g [t/op]\n", frameCount*per_frame, method, sqrt(f2err), fticks*1e-6, fticks/(nMols*nMols*nAtoms*nAtoms*per_frame) );

    transformAllAtoms( nMols, nAtoms, pos, molAtoms, atomsT );  // just for visualization

    glDisable(GL_LIGHTING);
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
    int i=0;
    for(int imol=0; imol<nMols; imol++){
        Draw::color_of_hash(i+154);
        for(int iatom=0; iatom<nAtoms; iatom++){
            //printf( " %i (%g,%g,%g)  \n", i, atomsT[i].x, atomsT[i].y, atomsT[i].z );
            Draw3D::drawShape(atomsT[i],mrot,atomView);
            //Draw3D::drawPointCross(atomsT[i],0.1);
           // Draw3D::drawVecInPos(forceAtomT[i]*10.0f,atomsT[i]);
           i++;
    }}
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
















