
#include <SDL2/SDL_opengl.h>

#include "drawMath2D.h"  // THE HEADER

namespace Draw2D{

void drawPoint( const Vec2f& vec ){
	glDisable (GL_LIGHTING);
	glBegin   (GL_POINTS);
		glVertex3f( vec.x, vec.y, 0 );
	glEnd();
};

void drawPointCross( const Vec2f& vec, float d ){
	glDisable (GL_LIGHTING);
	glBegin   (GL_LINES);
		glVertex3f( vec.x-d, vec.y,   0   );  glVertex3f( vec.x+d, vec.y,   0 );
		glVertex3f( vec.x,   vec.y-d, 0   );  glVertex3f( vec.x,   vec.y+d, 0  );
	glEnd();
};

void drawVec( const Vec2f& vec ){
	glDisable (GL_LIGHTING);
	glBegin   (GL_LINES);
		glVertex3f( 0, 0, 0 ); glVertex3f( vec.x, vec.y, 0 );
	glEnd();
};

void drawVecInPos( const Vec2f& v, const Vec2f& pos ){
	glDisable (GL_LIGHTING);
	glBegin   (GL_LINES);
		glVertex3f( pos.x, pos.y, 0 ); glVertex3f( pos.x+v.x, pos.y+v.y, 0 );
	glEnd();
};

void drawLine( const Vec2f& p1, const Vec2f& p2 ){
	glDisable (GL_LIGHTING);
	glBegin   (GL_LINES);
		glVertex3f( p1.x, p1.y, 0 ); glVertex3f( p2.x, p2.y, 0 );
	glEnd();
};

void drawTriangle( const Vec2f& p1, const Vec2f& p2, const Vec2f& p3 ){
	glBegin   (GL_TRIANGLES);
		glNormal3f( 0, 0, 1 );
		glVertex3f( p1.x, p1.y, 0 );
		glVertex3f( p2.x, p2.y, 0 );
		glVertex3f( p3.x, p3.y, 0 );
	glEnd();
};


void drawRectangle( float p1x, float p1y, float p2x, float p2y, bool filled ){
	if( filled){ glBegin(GL_QUADS); }else{ glBegin(GL_LINE_LOOP); };
	glBegin   (GL_QUADS);
		glVertex3f( p1x, p1y, 0 );
		glVertex3f( p1x, p2y, 0 );
		glVertex3f( p2x, p2y, 0 );
		glVertex3f( p2x, p1y, 0 );
	glEnd();
};

void drawRectangle( const Vec2f& p1, const Vec2f& p2, bool filled ){
	drawRectangle( p1.x, p1.y, p2.x, p2.y, filled );
};

void drawPoint_d( const Vec2d& vec ){
	Vec2f vec_;    convert( vec, vec_ ); drawPoint     ( vec_ );
};

void drawPointCross_d ( const Vec2d& vec, double d         ){
	Vec2f vec_;    convert( vec, vec_ ); drawPointCross( vec_, (float)d);
};

void drawVec_d( const Vec2d& vec ){
	Vec2f vec_;
	convert( vec, vec_ );
	drawVec( vec_ );
};

void drawVecInPos_d( const Vec2d& v,   const Vec2d& pos ){
	Vec2f v_,pos_; convert( v, v_ );     convert( pos, pos_ ); drawVecInPos( v_, pos_);
};

void drawLine_d( const Vec2d& p1,  const Vec2d& p2  ){
	Vec2f p1_,p2_; convert( p1, p1_ );   convert( p2, p2_ );   drawLine( p1_, p2_);
};

void drawRectangle_d( const Vec2d& p1,  const Vec2d& p2, bool filled ){
	Vec2f p1_,p2_; convert( p1, p1_ );   convert( p2, p2_ );   drawRectangle( p1_, p2_, filled );
};

void drawTriangle_d( const Vec2d& p1,  const Vec2d& p2, const Vec2d& p3 ){
	Vec2f p1_,p2_,p3_;  convert( p1, p1_ );  convert( p2, p2_ ); convert( p3, p3_ ); drawTriangle( p1_, p2_, p3_ );
};

void drawCircle( const Vec2f& center, float radius, int n ){
	glBegin   (GL_LINES);
	float dphi =  6.28318530718f / n;
	Vec2f drot; drot.fromAngle( dphi );
	Vec2f v;    v.set( radius, 0.0f );
	for ( int i=0; i<=n; i++ ){
		glVertex3f( center.x + v.x, center.y + v.y, 0 );
		v.mul_cmplx( drot );
	}
	glEnd();
};

void drawCircle_d( const Vec2d& center_, float radius, int n ){
	Vec2f center; convert( center_, center );
	drawCircle( center, radius, n );
};

void drawPoints( int npoints, Vec2d * points ){
	glBegin   (GL_POINTS);
	for( int i=0; i<npoints; i++ ){
		Vec2f p; convert( points[i], p );  glVertex3f( (float)p.x, (float)p.y, 0.0f );
	}
	glEnd();
};

void drawLines( int nlinks, int * links, Vec2d * points ){
	int n2 = nlinks<<1;
	glBegin   (GL_LINES);
	for( int i=0; i<n2; i+=2 ){
		Vec2f p1, p2;
		convert( points[links[i]],   p1 );  glVertex3f( (float)p1.x, (float)p1.y, 0.0f );
		convert( points[links[i+1]], p2 );  glVertex3f( (float)p2.x, (float)p2.y, 0.0f );
		//drawLine_d( points[links[i]], points[links[i+1]] );
		//printf ( " %i %i %i %f %f \n", i, links[i], links[i+1], points[links[i]].x, points[links[i+1]].x );
	}
	glEnd();
};

void drawConvexPolygon( int n, Vec2d * points, bool filled ){
	if( filled ){
		glBegin   ( GL_TRIANGLE_FAN );
	}else{
		glBegin   ( GL_LINE_LOOP   );
	}
	for( int i=0; i<n; i++ ){
		glVertex3f( (float)points[i].x, (float)points[i].y, 0 );
	}
	glEnd();
};

void drawPolarFunc( double x0, double y0, double fscale, int n, double phi0, double * data ){
		double dphi = M_PI_2/n;
		glBegin(GL_LINE_STRIP);
		for( int i=-1; i<n; i++ ){
			int ii = i;	if( i<0 ) ii=n-1;
			double cd  = data[ii];
			double phi = dphi*i + phi0;
			double x   = cos( phi );
			double y   = sin( phi );
			glVertex3f( (float)( x0 + fscale*x ), (float)(y0 + fscale*y), 0 );
		}
		glEnd();
};



void drawShape( const Vec2d& pos, const Vec2d& rot, int shape ){
	glPushMatrix();
	//glTranslatef( pos.x, pos.y , 0 );
	//glRotatef( phi*(180/M_PI), 0, 0, 1 );
	float glMat[16];
	toGLMat( pos, rot, glMat );
	glMultMatrixf( glMat );
	glCallList( shape );
	glPopMatrix();
};


}; // namespace Draw2D
