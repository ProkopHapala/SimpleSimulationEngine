
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
	int per_frame = 1;

	int err;
	OCLsystem cl;
	OCLtask *task1,*task2,*task3;

	//bool use_GPU = false;
	int method = 1;

	GLuint texture1 = 0;

	//virtual void loop( int nframes );
	virtual void eventHandling( const SDL_Event& event );
	virtual void draw   ();
	//virtual void drawHUD();

	void switchMethod( int i );

	TestApp_clInterpolate2D( int& id, int WIDTH_, int HEIGHT_ );

};

// ==================== Implementation

TestApp_clInterpolate2D::TestApp_clInterpolate2D( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {
    int err = CL_SUCCESS;
    //genXOR2D(nx, ny, hf);
    genRCOS(nx, ny, hf, 0.1);
    texture1 = Draw::makeTexture<Draw::float2RGBA>( nx, ny, hf );
    pointsOnLine( {2.0f,2.0f}, {1.0f,1.0f}, nPoints, points );

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
    printf("DEBUG 1 \n");
    cl_image_format imageFormat = {CL_A, CL_FLOAT};
    printf("DEBUG 2 \n");
    cl.newBufferImage2D( "hf", nx, ny,   sizeof(float),       hf, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR , imageFormat );
    printf("DEBUG 3 \n");
    cl.newBuffer(    "points",   nPoints, 2*sizeof(float),   points, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR );
    printf("DEBUG 4 \n");
    cl.newBuffer(      "vals",   nPoints,   sizeof(float),     vals, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR );
    printf("DEBUG 5 \n");
    err = cl.buildProgram( "cl/interpolate2D.cl" );     OCL_checkError(err, "cl.buildProgram");
    printf("DEBUG 6 \n");
    task1 = new OCLtask( &cl, cl.newKernel("getValInPoints"), 1, nPoints, 1 );
    printf("DEBUG 7 \n");
    task1->args = { BUFFarg(0), BUFFarg(1), BUFFarg(2) };
    printf("DEBUG 8 \n");
    task1->print_arg_list();
    printf("DEBUG 9 \n");
    //exit(0);
    //switchMethod( 4 );

    //cl.upload(1);
    task1->enque();
    clFinish(cl.commands);
    cl.download(2);

    lerp(nx, ny, hf, nPoints, points, vals_ );
    for(int i=0; i<nPoints; i++ ){ printf( "%i : %g  |  %g \n", i, vals[i], vals_[i] ); }
    //checkDiff( nPoints, vals, vals_ );
    printf("==== CHECK PASSED ====");
    //exit(0);
}

void TestApp_clInterpolate2D::draw(){
    long t1,t2,t3;
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable(GL_LIGHTING);

    //delay = 100;
    double f2err;
    t1=getCPUticks();
    switch(method){
        case 1:
            lerp(nx, ny, hf, nPoints, points, vals );
            break;
        case 2:
            task1->enque();
            clFinish(cl.commands);
            cl.download(2);
            break;
    }
    t1 = getCPUticks() - t1;

    printf( " %i method %i Time %g [Mticks] %g [ticks/op]\n", frameCount, method,  t1*1e-6, ((double)t1)/(per_frame*nPoints) );

    glDeleteTextures( 1, &texture1 );
    Draw::makeTexture<Draw::float2RGBA>( nx, ny, hf );
    if(texture1) Draw2D::renderImage(texture1,{-5.0,-5.0,5.0,5.0});

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
};

void TestApp_clInterpolate2D::switchMethod( int i ){
    method=i;
    setArray( nPoints, vals, 0 );
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
















