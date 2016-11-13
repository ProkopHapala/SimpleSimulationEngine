
#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include <SDL2/SDL_opengl.h>

#include "Draw3D.h"

#include "AeroCraft.h" // THE HEADER

void AeroCraft::render(){
	glPushMatrix();
        Mat3d camMat;
        camMat.setT(rotMat);
        float glmat[16];
        Draw3D::toGLMat( pos, camMat, glmat);
        glMultMatrixf( glmat );
        glDisable (GL_LIGHTING);
        glColor3f( 0.3f,0.3f,0.3f );
        Draw3D::drawSphere_oct(1,0.2,{0.0,0.0,0.0});
        Draw3D::drawLine( {0,0,0},wingLeft .lpos);
        Draw3D::drawLine( {0,0,0},wingRight.lpos);
        Draw3D::drawLine( {0,0,0},elevator .lpos);
        wingLeft .render();
        wingRight.render();
        rudder   .render();
        elevator .render();
	glPopMatrix();
};


