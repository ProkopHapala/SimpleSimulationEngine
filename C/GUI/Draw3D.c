#ifndef  Draw3D_c
#define  Draw3D_c

#include "VectorTypes.h"

#include <SDL2/SDL_opengl.h>

inline void Draw3f_toGLMat( const float3 pos, const float3x3 rot, float* glMat ){
    //printf("pos (%3.3f,%3.3f,%3.3f)\n", pos.x,pos.y,pos.z);
    glMat[0 ] = rot.ax;   glMat[1 ] = rot.ay;   glMat[2 ] = rot.az;   glMat[3 ]  = 0;
    glMat[4 ] = rot.bx;   glMat[5 ] = rot.by;   glMat[6 ] = rot.bz;   glMat[7 ]  = 0;
    glMat[8 ] = rot.cx;   glMat[9 ] = rot.cy;   glMat[10] = rot.cz;   glMat[11]  = 0;
    glMat[12] = pos. x;   glMat[13] = pos. y;   glMat[14] = pos. z;   glMat[15]  = 1;
}

inline void Draw3f_toGLMatCam( const float3 pos, const float3x3 rot, float* glMat ){
	glMat[0 ] = rot.ax;   glMat[1 ] = rot.bx;   glMat[2 ] = -rot.cx;   glMat[3 ]  = 0;
	glMat[4 ] = rot.ay;   glMat[5 ] = rot.by;   glMat[6 ] = -rot.cy;   glMat[7 ]  = 0;
	glMat[8 ] = rot.az;   glMat[9 ] = rot.bz;   glMat[10] = -rot.cz;   glMat[11]  = 0;
	glMat[12] = -pos. x;  glMat[13] = -pos. y;  glMat[14] = -pos. z;   glMat[15]  = 1;
};


void Draw3f_point( float3 vec ){
	//glDisable (GL_LIGHTING);
	glBegin   (GL_POINTS);
		glVertex3d( vec.x, vec.y, vec.z );
	glEnd();
};

void Draw3f_pointCross( float3 vec, float sz ){
	//glDisable (GL_LIGHTING);
	glBegin   (GL_LINES);
		glVertex3f( vec.x-sz, vec.y, vec.z ); glVertex3f( vec.x+sz, vec.y, vec.z );
		glVertex3f( vec.x, vec.y-sz, vec.z ); glVertex3f( vec.x, vec.y+sz, vec.z );
		glVertex3f( vec.x, vec.y, vec.z-sz ); glVertex3f( vec.x, vec.y, vec.z+sz );
	glEnd();
};

void Draw3f_drawVec( float3 vec ){
	//glDisable (GL_LIGHTING);
	glBegin   (GL_LINES);
		glVertex3d( 0, 0, 0 ); glVertex3d( vec.x, vec.y, vec.z );
	glEnd();
};

void Draw3f_drawVecInPos( float3 v, float3 pos ){
	//glDisable (GL_LIGHTING);
	glBegin   (GL_LINES);
		glVertex3d( pos.x, pos.y, pos.z ); glVertex3d( pos.x+v.x, pos.y+v.y, pos.z+v.z );
	glEnd();
};

void Draw3f_line( float3 p1, float3 p2 ){
	//glDisable (GL_LIGHTING);
	glBegin   (GL_LINES);
		glVertex3d( p1.x, p1.y, p1.z ); glVertex3d( p2.x, p2.y, p2.z );
	glEnd();
};


void Draw3f_arrow( float3 p1, float3 p2, float sz ){
	//glDisable (GL_LIGHTING);
	float3 up,lf,p;
    float3 fw = f3_sub(p2,p1); 
    f3_normalize(&fw);
    f3_someOrtho(fw,&up,&lf);
    f3_mulf_p(&fw,sz); 
    f3_mulf_p(&lf,sz); 
    f3_mulf_p(&up,sz);
	glBegin   (GL_LINES);
		glVertex3f( p1.x, p1.y, p1.z ); glVertex3f( p2.x, p2.y, p2.z );
		p = f3_add(f3_sub(p2,fw),up); glVertex3f( p.x, p.y, p.z ); glVertex3f( p2.x, p2.y, p2.z );
        p = f3_sub(f3_sub(p2,fw),up); glVertex3f( p.x, p.y, p.z ); glVertex3f( p2.x, p2.y, p2.z );
        p = f3_add(f3_sub(p2,fw),lf); glVertex3f( p.x, p.y, p.z ); glVertex3f( p2.x, p2.y, p2.z );
        p = f3_sub(f3_sub(p2,fw),lf); glVertex3f( p.x, p.y, p.z ); glVertex3f( p2.x, p2.y, p2.z );
	glEnd();
};


/*
void drawPolyLine( int n, Vec3d * ps, bool closed ){   // closed=false
    //printf("%i %i\n", n, closed );
    if(closed){ glBegin(GL_LINE_LOOP); }else{ glBegin(GL_LINE_STRIP); }
    for(int i=0; i<n; i++){
        //printf("%i (%3.3f,%3.3f,%3.3f)\n", i, ps[i].x, ps[i].y, ps[i].z );
        glVertex3d( ps[i].x, ps[i].y, ps[i].z );
    };
    glEnd();
};
*/


void Draw3f_box( float x0, float x1, float y0, float y1, float z0, float z1, float r, float g, float b ){
	glBegin(GL_QUADS);
		glColor3f( r, g, b );
		glNormal3f(0,0,-1); glVertex3f( x0, y0, z0 ); glVertex3f( x1, y0, z0 ); glVertex3f( x1, y1, z0 ); glVertex3f( x0, y1, z0 );
		glNormal3f(0,-1,0); glVertex3f( x0, y0, z0 ); glVertex3f( x1, y0, z0 ); glVertex3f( x1, y0, z1 ); glVertex3f( x0, y0, z1 );
		glNormal3f(-1,0,0); glVertex3f( x0, y0, z0 ); glVertex3f( x0, y1, z0 ); glVertex3f( x0, y1, z1 ); glVertex3f( x0, y0, z1 );
		glNormal3f(0,0,+1); glVertex3f( x1, y1, z1 ); glVertex3f( x0, y1, z1 ); glVertex3f( x0, y0, z1 ); glVertex3f( x1, y0, z1 );
		glNormal3f(0,+1,1); glVertex3f( x1, y1, z1 ); glVertex3f( x0, y1, z1 ); glVertex3f( x0, y1, z0 ); glVertex3f( x1, y1, z0 );
		glNormal3f(+1,0,0); glVertex3f( x1, y1, z1 ); glVertex3f( x1, y0, z1 ); glVertex3f( x1, y0, z0 ); glVertex3f( x1, y1, z0 );
	glEnd();
}

void Draw3f_axis( float sc ){
	//glDisable (GL_LIGHTING);
	glBegin   (GL_LINES);
		glColor3f( 1, 0, 0 ); glVertex3f( 0, 0, 0 ); glVertex3f( 1*sc, 0, 0 );
		glColor3f( 0, 1, 0 ); glVertex3f( 0, 0, 0 ); glVertex3f( 0, 1*sc, 0 );
		glColor3f( 0, 0, 1 ); glVertex3f( 0, 0, 0 ); glVertex3f( 0, 0, 1*sc );
	glEnd();
}

#endif