
#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include <SDL2/SDL_opengl.h>

#include "drawMath.h"

#include "AeroCraft.h" // THE HEADER

void AeroCraft::render(){
	glPushMatrix();
	float glmat[16];
	toGLMat( pos, rotMat, glmat);
	glMultMatrixf( glmat );
	glDisable (GL_LIGHTING);
	glColor3f( 0.3f,0.3f,0.3f );
	drawLine( {0,0,0},wingLeft.lpos);
	drawLine( {0,0,0},wingRight.lpos);
	drawLine( {0,0,0},elevator.lpos);
	wingLeft.render();
	wingRight.render();
	rudder.render();
	elevator.render();
	glPopMatrix();
};


