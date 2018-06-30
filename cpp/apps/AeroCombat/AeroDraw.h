
#ifndef AeroDraw_h
#define AeroDraw_h

#include "AeroCraft.h"
#include "AeroTest.h"

#include "Draw.h"
#include "Draw2D.h"
#include "Draw3D.h"

void renderAeroCraft ( const AeroCraft& craft, bool bPos ){
	glPushMatrix();
        Mat3d camMat;
        camMat.set(craft.rotMat);
        float glmat[16];

        Vec3d pos;
        if   (bPos){ pos = craft.pos;            }
        else       { pos = (Vec3d){0.0,0.0,0.0}; }

        Draw3D::toGLMat( pos, camMat, glmat);
        //printf( "%g %g %g\n", pos.x, pos.y, pos.z);
        //printf( "%g %g %g\n", camMat.ax, camMat.ay, camMat.az);
        //printf( "%g %g %g\n", camMat.bx, camMat.by, camMat.bz);
        //printf( "%g %g %g\n", camMat.cx, camMat.cy, camMat.cz);
        glMultMatrixf( glmat );
        //glDisable (GL_LIGHTING);
        //glColor3f( 0.3f,0.3f,0.3f );

        Draw3D::drawSphere_oct(1,0.2,{0.0,0.0,0.0});

        //Draw3D::drawLine( {0,0,0},wingLeft .lpos);
        //Draw3D::drawLine( {0,0,0},wingRight.lpos);
        //Draw3D::drawLine( {0,0,0},elevator .lpos);
        //wingLeft .render();
        //wingRight.render();
        //rudder   .render();
        //elevator .render();


		for( int i=0; i<craft.nPanels; i++ ){
            Draw3D::drawLine( {0,0,0}, craft.panels[i].lpos);
            //printf( "%g %g %g\n", panels[i].lpos.x, panels[i].lpos.y, panels[i].lpos.z);
            //panels[i].render( );
            //Draw3D::drawKite( craft.panels[i].lpos, craft.panels[i].lrot, sqrt(craft.panels[i].area) );
            double sz = sqrt(craft.panels[i].area);
            Draw3D:: drawPanel( craft.panels[i].lpos, craft.panels[i].lrot, {sz,sz*0.25} );
		}

	glPopMatrix();
};

void renderAeroCraftDebug ( const AeroCraft& craft, bool bPos ){
	glPushMatrix();
        Mat3d camMat;
        camMat.setT(craft.rotMat);
        float glmat[16];

        Vec3d pos;
        if   (bPos){ pos = craft.pos;            }
        else       { pos = (Vec3d){0.0,0.0,0.0}; }

        Draw3D::toGLMat( pos, camMat, glmat);
        //printf( "%g %g %g\n", pos.x, pos.y, pos.z);
        //printf( "%g %g %g\n", camMat.ax, camMat.ay, camMat.az);
        //printf( "%g %g %g\n", camMat.bx, camMat.by, camMat.bz);
        //printf( "%g %g %g\n", camMat.cx, camMat.cy, camMat.cz);
        glMultMatrixf( glmat );
        //glDisable (GL_LIGHTING);
        //glColor3f( 0.3f,0.3f,0.3f );

        Draw3D::drawSphere_oct(1,0.2,{0.0,0.0,0.0});

        //Draw3D::drawLine( {0,0,0},wingLeft .lpos);
        //Draw3D::drawLine( {0,0,0},wingRight.lpos);
        //Draw3D::drawLine( {0,0,0},elevator .lpos);
        //wingLeft .render();
        //wingRight.render();
        //rudder   .render();
        //elevator .render();


		for( int i=0; i<craft.nPanels; i++ ){
            Draw3D::drawLine( {0,0,0}, craft.panels[i].lpos);
            //printf( "%g %g %g\n", panels[i].lpos.x, panels[i].lpos.y, panels[i].lpos.z);
            //panels[i].render( );
            //Draw3D::drawKite( craft.panels[i].lpos, craft.panels[i].lrot, sqrt(craft.panels[i].area) );
            double sz = sqrt(craft.panels[i].area);
            Draw3D:: drawPanel( craft.panels[i].lpos, craft.panels[i].lrot, {sz,sz*0.25} );
		}

	glPopMatrix();
};

#endif
