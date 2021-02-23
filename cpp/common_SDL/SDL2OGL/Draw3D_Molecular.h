
#ifndef  Draw3D_Molecular_h
#define  Draw3D_Molecular_h

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include <math.h>
#include <cstdlib>
#include <stdio.h>
#include <stdint.h>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

#include "Draw.h"
//#include "Mesh.h"
#include "Draw3D.h"

//#include <SDL2/SDL_opengl.h>

namespace Draw3D{

void plotSurfPlane( Vec3d normal, double c0, Vec2d d, Vec2i n ){
    Vec3d da,db;
    normal.getSomeOrtho( da,db );
    da.mul( d.a/da.norm() );
    db.mul( d.b/db.norm() );
    //glColor3f(1.0f,0.0f,0.0f); Draw3D::drawVecInPos(normal, {0.0,0.0,0.0} );
    //glColor3f(0.0f,1.0f,0.0f); Draw3D::drawVecInPos(da*10, {0.0,0.0,0.0} );
    //glColor3f(0.0f,0.0f,1.0f); Draw3D::drawVecInPos(db*10, {0.0,0.0,0.0} );
    Draw3D::drawRectGridLines( n*2, (da*-n.a)+(db*-n.b) + normal*c0, da, db );
}

void makeSphereOgl( int& ogl, int nsub, float sz ){
    ogl = Draw::list(ogl);
    //glNewList( ogl, GL_COMPILE );
        //glEnable( GL_LIGHTING );
        //glColor3f( 0.8f, 0.8f, 0.8f );
        //Draw3D::drawSphere_oct(3, 0.5, {0.0,0.0,0.0} );
        Draw3D::drawSphere_oct( nsub, sz, {0.0,0.0,0.0} );
    glEndList();
}


void atomsREQ( int n, Vec3d* ps, Vec3d* REQs, int ogl_sph, float qsc=1, float Rsc=1, float Rsub=0 ){
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);
    for(int i=0; i<n; i++){
        float q = (float)REQs[i].z*qsc;
        glColor3f(1-fmax(0,-q),1-fmax(q,-q),1-fmax(0,+q));
        //Draw3D::drawShape( ogl_sph, ps[i], Mat3dIdentity*(REQs[i].x*Rsc) );
        Draw3D::drawShape( ogl_sph, ps[i], Mat3dIdentity*((REQs[i].x-Rsub)*Rsc) );
    }
}

void bondLabels( int n, const Vec2i* b2a, const Vec3d* apos, int fontTex, float sz=0.02 ){
    char str[256];
    for(int i=0; i<n; i++){
        /*
        Vec2i ib = b2a[i];
        //glColor3f(0.0f,0.0f,0.0f);
        if(i==ibpicked) glColor3f(1.0f,0.0f,0.0f); ;
        Draw3D::drawLine(world.apos[ib.x],world.apos[ib.y]);
        sprintf(str,"%i\0",i);
        */
        //Draw3D::drawText(str, (world.apos[ib.x]+world.apos[ib.y])*0.5, fontTex, 0.02, 0,0);
        Vec2i ib = b2a[i];
        sprintf(str,"%i\0",i);
        Draw3D::drawText(str, (apos[ib.x]+apos[ib.y])*0.5, fontTex, sz, 0);
    }
}

}; // namespace Draw3D

#endif

