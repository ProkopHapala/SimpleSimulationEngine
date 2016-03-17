
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








