
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

#include "interpolate2D.h"
#include "OCL.h"

#define R2SAFE  1.0e-2
#define F2MAX   10.0

class TestApp_clInterpolate2D : public AppSDL2OGL {
	public:
	int per_frame = 50;
	float dt=0.5,damp=0.99;

	int err;
	OCLsystem cl;
	OCLtask *task1,*task2,*task3;

	//bool use_GPU = false;
	int  method = 1;
	char job    = 'r';

	GLuint texture1 = 0;

	//virtual void loop( int nframes );
	virtual void eventHandling( const SDL_Event& event );
	virtual void draw   ();
	//virtual void drawHUD();

	void switchMethod( int i );
	void switchJob( char ch );

	void job_InterpolateLine();
	void job_InterpolateDeriv();
	void job_relax();

	TestApp_clInterpolate2D( int& id, int WIDTH_, int HEIGHT_ );

};

// ==================== Implementation

TestApp_clInterpolate2D::TestApp_clInterpolate2D( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {
    int err = CL_SUCCESS;
    //genXOR2D(nx, ny, hf);
    genRCOS(nx, ny, hf, 0.1);
    arrayDerivs2D( nx, ny, hf, Dhf, 1.0f, 1.0f );
    texture1 = Draw::makeTexture<Draw::float2RGBA>( nx, ny, hf );
    pointsOnLine( {50.0f,20.0f}, {100.0f/nPoints,200.0f/nPoints}, nPoints, points );

    // --- OpenCL
    cl.init();

    // copy input data
    //ret = clEnqueueWriteImage(
    //  instancePtr->queue, instancePtr->image_x, CL_FALSE, image_x_origin, image_x_region,
    //  0,                               // row pitch
    //  0,                               // slice pitch
    //  x, 0, NULL, &wrEv[numWrEv++]);
    //checkResult(ret);

    // cl_image_format imageFormat = {CL_RGBA, CL_FLOAT};  // this will probably make float4 image (?)

    cl.newBuffer(    "points",   nPoints, 2*sizeof(float),   points, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR );
    cl.newBuffer(      "vals",   nPoints,   sizeof(float),     vals, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR );

    cl.newBufferImage2D(  "hf", nx, ny,   sizeof(float),      hf, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR , {CL_A,  CL_FLOAT} );
    cl.newBufferImage2D( "Dhf", nx, ny,   sizeof(float),     Dhf, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR , {CL_RG, CL_FLOAT} );
    cl.newBuffer(      "Dvals", nPoints,2*sizeof(float),   Dvals, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR );

    err = cl.buildProgram( "cl/interpolate2D.cl" );     OCL_checkError(err, "cl.buildProgram");

    task1 = new OCLtask( &cl, cl.newKernel("getValInPoints"), 1, nPoints, 32 );
    task1->args = { BUFFarg(2), BUFFarg(0), BUFFarg(1) };
    task1->print_arg_list();

    task2 = new OCLtask( &cl, cl.newKernel("getF2InPoints"), 1, nPoints, 32 );
    task2->args = { BUFFarg(3), BUFFarg(0), BUFFarg(4) };
    task2->print_arg_list();

    task3 = new OCLtask( &cl, cl.newKernel("relaxPoints"), 1, nPoints, 32 );
    task3->args = { BUFFarg(3), BUFFarg(0), BUFFarg(4), INTarg(per_frame), FLOATarg(dt), FLOATarg(damp) };
    task3->print_arg_list();

    /*
    //cl.upload(1);
    task1->enque();
    clFinish(cl.commands);
    cl.download(1);
    lerp(nx, ny, hf, nPoints, points, vals_ );
    //for(int i=0; i<nPoints; i++ ){ printf( "%i : %g  |  %g \n", i, vals[i], vals_[i] ); }
    //checkDiff( nPoints, vals, vals_ );
    printf("==== CHECK PASSED ====");
    //exit(0);
    */
}

void TestApp_clInterpolate2D::job_InterpolateLine(){
    long t1=getCPUticks();
    switch(method){
        case 1:
            lerp(nx, ny, hf, nPoints, points, vals );
            break;
        case 2:
            task1->enque();
            clFinish(cl.commands);
            cl.download(1);
            break;
    }
    t1 = getCPUticks() - t1;
    printf( " %i method %i Time %g [Mticks] %g [ticks/op]\n", frameCount, method,  t1*1e-6, ((double)t1)/(per_frame*nPoints) );

    glColor3f(1.0f,0.0f,1.0f);
    glBegin(GL_LINE_STRIP);
    for(int i=0; i<nPoints; i++){ glVertex3f( 10.0f*(points[i].x/nx-0.5f), vals[i]-5.0f, 10.0f ); }
    glEnd();
    glBegin(GL_LINE_STRIP);
    for(int i=0; i<nPoints; i++){ glVertex3f( vals[i]-5.0f, 10.0f*(points[i].y/ny-0.5f), 10.0f ); }
    glEnd();
    glColor3f(0.0f,1.0f,0.0f);
    glBegin(GL_LINE_STRIP);
    for(int i=0; i<nPoints; i++){ glVertex3f( 10.0f*(points[i].x/nx-0.5f), 10.0f*(points[i].y/ny-0.5f), 10.0f ); }
    glEnd();
}

void TestApp_clInterpolate2D::job_InterpolateDeriv(){
    long t1=getCPUticks();
    switch(method){
        case 1:
            lerpD( nx, ny, (Vec2f*)Dhf, nPoints, points, Dvals );
            break;
        case 2:
            task2->enque();
            clFinish(cl.commands);
            cl.download(4);
            break;
    }
    t1 = getCPUticks() - t1;
    printf( " %i method %i Time %g [Mticks] %g [ticks/op]\n", frameCount, method,  t1*1e-6, ((double)t1)/(per_frame*nPoints) );

    float fscale = -10.0;
    glColor3f(1.0f,0.0f,1.0f);
    glBegin(GL_LINES);
    for(int i=0; i<nPoints; i++){
        Vec2f p; p.set( 10.0f*(points[i].x/nx-0.5f), 10.0f*(points[i].y/ny-0.5f) );
        //printf( " %i (%g,%g) (%g,%g)\n", i, p.x, p.y, Dvals[i].x, Dvals[i].y );
        glVertex3f( p.x,                   p.y,                   10.0f );
        //glVertex3f( p.x+5, p.y+5, 10.0f );
        glVertex3f( p.x+Dvals[i].x*fscale, p.y+Dvals[i].y*fscale, 10.0f );
    }
    glEnd();

    //exit(0);
}


void TestApp_clInterpolate2D::job_relax(){
    long t1=getCPUticks();
    //pointsOnLine( {50.0f,20.0f}, {0.79f,0.35f}, nPoints, points );
    //pointsOnLine( {50.0f,20.0f}, {100.0f/nPoints,200.0f/nPoints}, nPoints, points );
    genPointsHash( nPoints, points, {20.0f,20.0f}, {230.0f,230.0f}, 156468 );

    switch(method){
        case 1:
            relaxPoints( nx, ny, (Vec2f*)Dhf, nPoints, points, Dvals, per_frame, dt, damp );
            //printf("relaxPoint %i %f %f \n", per_frame, damp, dt );
            break;
        case 2:
            cl.upload(0);
            task3->enque();
            clFinish(cl.commands);
            cl.download(0);
            //exit(0);
            break;
    }
    t1 = getCPUticks() - t1;
    printf( " %i method %i Time %g [Mticks] %g [ticks/op]\n", frameCount, method,  t1*1e-6, ((double)t1)/(per_frame*nPoints) );

    glColor3f(1.0f,0.0f,1.0f);
    glBegin(GL_POINTS);
    for(int i=0; i<nPoints; i++){
        glVertex3f( 10.0f*(points[i].x/nx-0.5f), 10.0f*(points[i].y/ny-0.5f), 10.0f );
    }
    glEnd();

    //exit(0);
}



void TestApp_clInterpolate2D::draw(){
    long t1,t2,t3;
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable(GL_LIGHTING);


    glDeleteTextures( 1, &texture1 );
    Draw::makeTexture<Draw::float2RGBA>( nx, ny, hf );
    if(texture1) Draw2D::renderImage(texture1,{-5.0,-5.0,5.0,5.0});

    switch(job){
        case 'l': job_InterpolateLine(); break;
        case 'd': job_InterpolateDeriv(); break;
        case 'r': job_relax(); break;
    }


};

void TestApp_clInterpolate2D::switchMethod( int i ){
    method=i;
    setArray( nPoints  ,          vals, 0 );
    setArray( nPoints*2, (float*)Dvals, 0 );
    //genXOR2D(nx, ny, hf);
    //genZero(nx, ny, buff);
    //cl.upload(0);
}

void TestApp_clInterpolate2D::switchJob( char ch ){
    job=ch;
    //setArray( nPoints, vals, 0 );
    //genXOR2D(nx, ny, hf);
    //genZero(nx, ny, buff);
    //cl.upload(0);
}

void TestApp_clInterpolate2D::eventHandling( const SDL_Event& event ){
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_KP_1:  switchMethod( 1 ); break;
                case SDLK_KP_2:  switchMethod( 2 ); break;
                //case SDLK_KP_3:  switchMethod( 3 ); break;
                //case SDLK_KP_4:  switchMethod( 4 ); break;
                case SDLK_l:  switchJob( 'l' ); break;
                case SDLK_d:  switchJob( 'd' ); break;
                case SDLK_r:  switchJob( 'r' ); break;
            } break;
    };
    AppSDL2OGL::eventHandling( event );
};


// ===================== MAIN

TestApp_clInterpolate2D * app;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	app = new TestApp_clInterpolate2D( junk , 800, 600 );
	app->loop( 1000000 );
	return 0;
}
















