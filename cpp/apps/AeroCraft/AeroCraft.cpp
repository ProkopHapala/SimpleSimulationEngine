
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
        //printf( "%g %g %g\n", pos.x, pos.y, pos.z);
        //printf( "%g %g %g\n", camMat.ax, camMat.ay, camMat.az);
        //printf( "%g %g %g\n", camMat.bx, camMat.by, camMat.bz);
        //printf( "%g %g %g\n", camMat.cx, camMat.cy, camMat.cz);
        glMultMatrixf( glmat );
        glDisable (GL_LIGHTING);
        glColor3f( 0.3f,0.3f,0.3f );

        Draw3D::drawSphere_oct(1,0.2,{0.0,0.0,0.0});

        /*
        Draw3D::drawLine( {0,0,0},wingLeft .lpos);
        Draw3D::drawLine( {0,0,0},wingRight.lpos);
        Draw3D::drawLine( {0,0,0},elevator .lpos);
        wingLeft .render();
        wingRight.render();
        rudder   .render();
        elevator .render();
        */

		for( int i=0; i<nPanels; i++ ){
            Draw3D::drawLine( {0,0,0}, panels[i].lpos);
            //printf( "%g %g %g\n", panels[i].lpos.x, panels[i].lpos.y, panels[i].lpos.z);
            panels[i].render( );
		}

	glPopMatrix();
};


	int AeroCraft::fromFile( const char * fname ){
        char buf  [1024];
        FILE * pFile;
        pFile = fopen (fname,"r");
        printf(" loading molTypes from: >>%s<<\n", fname );

        fscanf (pFile, " %lf %lf %lf %lf\n", &mass, &Ibody.xx, &Ibody.yy, &Ibody.zz );
        printf( " %lf %lf %lf %lf\n", mass, Ibody.xx, Ibody.yy, Ibody.zz );

        fscanf ( pFile, "%i\n", &nPanels);
        panels = new AeroSurface[nPanels];
        for(int i=0; i<nPanels; i++){
            fgets( buf, 256, pFile);
            //AeroSurface * a = new AeroSurface();
            panels[i].fromString( buf );
            panels[i].craft = this;
        }

        Ibody.invert_to( invIbody );
        qrot.setOne();
        qrot.toMatrix(rotMat);
        L.set(0,0,0);
        setMass( mass );
        vel.set(0,0,0);
        pos.set(0,500,0);
        clean_temp();

        printf("AeroCraft loaded\n");
        fclose(pFile);
        return nPanels;
    }


