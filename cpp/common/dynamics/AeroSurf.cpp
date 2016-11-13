
#include <cstdio>
#include <cstdlib>

#include <SDL2/SDL_opengl.h>

//#include <Draw3D.h>

#include "AeroSurf.h" // THE HEADER

void AeroSurface::render( ){
	//drawLine( const Vec3d& p1, const Vec3d& p2 );
	double sz = 3*sqrt( C.z );
	glEnable (GL_LIGHTING);
	glColor3f( 0.5f,0.5f,0.5f );
	glBegin  (GL_QUADS);
		glNormal3d( lrot.b.x, lrot.b.y, lrot.b.z );
		glVertex3d( lpos.x-sz*lrot.a.x, lpos.y-sz*lrot.a.y, lpos.z-sz*lrot.a.z );
		glVertex3d( lpos.x-sz*lrot.c.x, lpos.y-sz*lrot.c.y, lpos.z-sz*lrot.c.z );
		glVertex3d( lpos.x+sz*lrot.a.x, lpos.y+sz*lrot.a.y, lpos.z+sz*lrot.a.z );
		glVertex3d( lpos.x+sz*lrot.c.x, lpos.y+sz*lrot.c.y, lpos.z+sz*lrot.c.z );
	glEnd();
};

void AeroSurface::fromString( const char * str ){
    sscanf ( str, " %lf %lf %lf    %lf %lf %lf    %lf %lf %lf    %lf %lf %lf",
                &lpos.x, &lpos.y, &lpos.z,
                &lrot.bx, &lrot.by, &lrot.bz,
                &lrot.cx, &lrot.cy, &lrot.cz,
                &C.a, &C.b, &C.c
            );
    //printf ( " %lf %lf %lf    %lf %lf %lf    %lf %lf %lf    %lf %lf %lf\n", lpos.x, lpos.y, lpos.z, lrot.bx, lrot.by, lrot.bz, lrot.cx, lrot.cy, lrot.cz, C.a, C.b, C.c );

    printf ( " %lf %lf %lf   %lf %lf %lf\n", lpos.x, lpos.y, lpos.z, C.a, C.b, C.c );

    lrot.b.normalize();
    double bc = lrot.b.dot(lrot.c);
    lrot.c.add_mul( lrot.b, -bc );
    lrot.a.set_cross( lrot.b, lrot.c );
    lrot.a.normalize();
};







