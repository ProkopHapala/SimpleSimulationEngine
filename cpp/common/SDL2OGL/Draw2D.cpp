
#include <SDL2/SDL_opengl.h>

#include "Draw.h"
#include "Draw2D.h"  // THE HEADER

//namespace Draw2D{

float Draw2D::z_layer = 0.0f; // should be initialized like this http://stackoverflow.com/questions/19136059/namespace-global-variable-losing-value-c

void Draw2D::drawPoint( const Vec2f& vec ){
	glDisable (GL_LIGHTING);
	glBegin   (GL_POINTS);
		glVertex3f( vec.x, vec.y, z_layer );
	glEnd();
};

void Draw2D::drawPointCross( const Vec2f& vec, float d ){
	glDisable (GL_LIGHTING);
	glBegin   (GL_LINES);
		glVertex3f( vec.x-d, vec.y,   z_layer   );  glVertex3f( vec.x+d, vec.y,   z_layer );
		glVertex3f( vec.x,   vec.y-d, z_layer   );  glVertex3f( vec.x,   vec.y+d, z_layer  );
	glEnd();
};

void Draw2D::drawVec( const Vec2f& vec ){
	glDisable (GL_LIGHTING);
	glBegin   (GL_LINES);
		glVertex3f( 0, 0, z_layer ); glVertex3f( vec.x, vec.y, z_layer );
	glEnd();
};

void Draw2D::drawVecInPos( const Vec2f& v, const Vec2f& pos ){
	glDisable (GL_LIGHTING);
	glBegin   (GL_LINES);
		glVertex3f( pos.x, pos.y, z_layer ); glVertex3f( pos.x+v.x, pos.y+v.y, z_layer );
	glEnd();
};

void Draw2D::drawBody2d( const Vec2f& rot, const Vec2f& pos, float l1, float l2 ){
    //printf( "(%3.3f,%3.3f) (%3.3f,%3.3f)\n",  rot.x, rot.y, pos.x, pos.y );
	glDisable (GL_LIGHTING);
	glBegin   (GL_LINES);
		glVertex3f( pos.x, pos.y, z_layer ); glVertex3f( pos.x+rot.x*l1, pos.y+rot.y*l1, z_layer );
		glVertex3f( pos.x, pos.y, z_layer ); glVertex3f( pos.x+rot.y*l2, pos.y-rot.x*l2, z_layer );
	glEnd();
};

void Draw2D::drawLine( const Vec2f& p1, const Vec2f& p2 ){
	glDisable (GL_LIGHTING);
	glBegin   (GL_LINES);
		glVertex3f( p1.x, p1.y, z_layer ); glVertex3f( p2.x, p2.y, z_layer );
	glEnd();
};

void Draw2D::drawTriangle( const Vec2f& p1, const Vec2f& p2, const Vec2f& p3 ){
	glBegin   (GL_TRIANGLES);
		glNormal3f( 0, 0, 1 );
		glVertex3f( p1.x, p1.y, z_layer );
		glVertex3f( p2.x, p2.y, z_layer );
		glVertex3f( p3.x, p3.y, z_layer );
	glEnd();
};


void Draw2D::drawRectangle( float p1x, float p1y, float p2x, float p2y, bool filled ){
	if( filled){ glBegin(GL_QUADS); }else{ glBegin(GL_LINE_LOOP); };
	glBegin   (GL_QUADS);
		glVertex3f( p1x, p1y, z_layer );
		glVertex3f( p1x, p2y, z_layer );
		glVertex3f( p2x, p2y, z_layer );
		glVertex3f( p2x, p1y, z_layer );
	glEnd();
};

void Draw2D::drawRectangle( const Vec2f& p1, const Vec2f& p2, bool filled ){
	drawRectangle( p1.x, p1.y, p2.x, p2.y, filled );
};

void Draw2D::drawPoint_d( const Vec2d& vec ){
	Vec2f vec_;    convert( vec, vec_ ); drawPoint     ( vec_ );
};

void Draw2D::drawPointCross_d ( const Vec2d& vec, double d         ){
	Vec2f vec_;    convert( vec, vec_ ); drawPointCross( vec_, (float)d);
};

void Draw2D::drawVec_d( const Vec2d& vec ){
	Vec2f vec_;
	convert( vec, vec_ );
	drawVec( vec_ );
};

void Draw2D::drawVecInPos_d( const Vec2d& v,   const Vec2d& pos ){
	Vec2f v_,pos_; convert( v, v_ );     convert( pos, pos_ ); drawVecInPos( v_, pos_);
};

void Draw2D::drawBody2d_d( const Vec2d& rot,   const Vec2d& pos, float l1, float l2 ){
	Vec2f rot_,pos_; convert( rot, rot_ );     convert( pos, pos_ ); drawBody2d( rot_, pos_, l1, l2 );
};

void Draw2D::drawLine_d( const Vec2d& p1,  const Vec2d& p2  ){
	Vec2f p1_,p2_; convert( p1, p1_ );   convert( p2, p2_ );   drawLine( p1_, p2_);
};

void Draw2D::drawRectangle_d( const Vec2d& p1,  const Vec2d& p2, bool filled ){
	Vec2f p1_,p2_; convert( p1, p1_ );   convert( p2, p2_ );   drawRectangle( p1_, p2_, filled );
};

void Draw2D::drawTriangle_d( const Vec2d& p1,  const Vec2d& p2, const Vec2d& p3 ){
	Vec2f p1_,p2_,p3_;  convert( p1, p1_ );  convert( p2, p2_ ); convert( p3, p3_ ); drawTriangle( p1_, p2_, p3_ );
};

void Draw2D::drawCircle( const Vec2f& center, float radius, int n, bool filled ){
    //printf( " z_layer %3.3f \n", z_layer );
	if( filled){ glBegin(GL_TRIANGLE_FAN); }else{ glBegin(GL_LINE_LOOP); };
	float dphi =  6.28318530718f / n;
	Vec2f drot; drot.fromAngle( dphi );
	Vec2f v;    v.set( radius, 0.0f );
	for ( int i=0; i<n; i++ ){
		glVertex3f( center.x + v.x, center.y + v.y, z_layer );
		v.mul_cmplx( drot );
	}
	glEnd();
};

void Draw2D::drawCircle_d( const Vec2d& center_, float radius, int n, bool filled ){
	Vec2f center; convert( center_, center );
	drawCircle( center, radius, n , filled );
};




void Draw2D::drawPoints( int npoints, Vec2d * points ){
	glBegin   (GL_POINTS);
	for( int i=0; i<npoints; i++ ){
		Vec2f p; convert( points[i], p );
		glVertex3f( p.x, p.y, z_layer );
	}
	glEnd();
};

void Draw2D::drawPoints( int npoints, Vec2d * points, float sc  ){
	glBegin   (GL_LINES);
	for( int i=0; i<npoints; i++ ){
		Vec2f p; convert( points[i], p );
		glVertex3f( p.x-sc, p.y, z_layer );
		glVertex3f( p.x+sc, p.y, z_layer );
        glVertex3f( p.x, p.y-sc, z_layer );
		glVertex3f( p.x, p.y+sc, z_layer );
	}
	glEnd();
};

void Draw2D::drawLines( int n, Vec2d * points ){
	glBegin   (GL_LINE_STRIP);
	for( int i=0; i<n; i++ ){
		Vec2f p1;
		convert( points[i-1], p1 );
		glVertex3f( p1.x, p1.y, z_layer );
		//drawLine_d( points[links[i]], points[links[i+1]] );
		//printf ( " %i %i %i %f %f \n", i, links[i], links[i+1], points[links[i]].x, points[links[i+1]].x );
	}
	glEnd();
};

void Draw2D::drawLines( int nlinks, int * links, Vec2d * points ){
	int n2 = nlinks<<1;
	glBegin   (GL_LINES);
	for( int i=0; i<n2; i+=2 ){
		Vec2f p1, p2;
		convert( points[links[i]],   p1 );  glVertex3f( p1.x, p1.y, z_layer );
		convert( points[links[i+1]], p2 );  glVertex3f( p2.x, p2.y, z_layer );
		//drawLine_d( points[links[i]], points[links[i+1]] );
		//printf ( " %i %i %i %f %f \n", i, links[i], links[i+1], points[links[i]].x, points[links[i+1]].x );
	}
	glEnd();
};

void Draw2D::drawConvexPolygon( int n, Vec2d * points, bool filled ){
	if( filled ){glBegin   ( GL_TRIANGLE_FAN );	}else{glBegin   ( GL_LINE_LOOP   );	}
	for( int i=0; i<n; i++ ){
		glVertex3f( (float)points[i].x, (float)points[i].y, z_layer );
	}
	glEnd();
};

void Draw2D::drawPolarFunc( double x0, double y0, double fscale, int n, double phi0, double * data ){
		double dphi = M_PI_2/n;
		glBegin(GL_LINE_STRIP);
		for( int i=-1; i<n; i++ ){
			int ii = i;	if( i<0 ) ii=n-1;
			double cd  = data[ii];
			double phi = dphi*i + phi0;
			double x   = cos( phi );
			double y   = sin( phi );
			glVertex3f( (float)( x0 + fscale*x ), (float)(y0 + fscale*y), z_layer );
		}
		glEnd();
};

void Draw2D::plot( int n, double * xs, double * ys ){
    glBegin(GL_LINE_STRIP);
    for( int i=0; i<n; i++ ){
        glVertex3f( (float)xs[i], (float)ys[i], z_layer );
    }
    glEnd();
};

void Draw2D::plot_cross( int n, double * xs, double * ys, double sz ){
    glBegin   (GL_LINES);
	for( int i=0; i<n; i++ ){
        float x = (float)xs[i]; float y = (float)ys[i];
		glVertex3f( x-sz, y   , z_layer );
		glVertex3f( x+sz, y   , z_layer );
        glVertex3f( x   , y-sz, z_layer );
		glVertex3f( x   , y+sz, z_layer );
	}
	glEnd();
};

void Draw2D::drawFunc( float xmin, float xmax, int n, Func1d func ){
    glBegin(GL_LINE_STRIP);
    float dx = (xmax-xmin)/n;
    for( float x=xmin; x<=xmax; x+=dx ){
        float y = (float) func( x );
        glVertex3f( x, y, z_layer );
    }
    glEnd();
};

void Draw2D::drawFuncDeriv( float xmin, float xmax, float d, int n, Func1d func ){
    glBegin(GL_LINE_STRIP);
    float dx = (xmax-xmin)/n;
    for( float x=xmin; x<=xmax; x+=dx ){
        float y  = (float) func( x );
        float y_ = (float) func( x + d );
        float dy = (y_ - y)/d;
        glVertex3f( x, dy, z_layer );
    }
    glEnd();
};

void Draw2D::drawGrid( float xmin, float ymin, float xmax, float ymax, float dx, float dy ){
    glBegin(GL_LINES);
    // X-grid
    int nmin,nmax;
    if( xmin>0 ){ nmin=(int)(xmin/dx)+1; }else{ nmin=(int)(xmin/dx);   }
    if( xmax>0 ){ nmax=(int)(xmax/dx);   }else{ nmax=(int)(xmax/dx)-1; }
    for( int i=nmin; i<=nmax; i++ ){
        float x = i*dx;
        glVertex3f( x, ymin, z_layer );
        glVertex3f( x, ymax, z_layer );
    }
    // Y-grid
    if( ymin>0 ){ nmin=(int)(ymin/dy)+1; }else{ nmin=(int)(ymin/dy); }
    if( ymax>0 ){ nmax=(int)(ymax/dy);   }else{ nmax=(int)(ymax/dy)-1; }
    for( int i=nmin; i<=nmax; i++ ){
        float y = i*dy;
        glVertex3f( xmin, y, z_layer );
        glVertex3f( xmax, y, z_layer );
    }
    glEnd();
};

void Draw2D::drawSimplex( float x, float y, bool s, float step ){
    glBegin   ( GL_TRIANGLES );
    if( s ){
        glVertex3f( (float)(x+0.5*step ), (float)(y+0.86602540378*step), 0.0f );
        glVertex3f( (float)(x+1.5*step ), (float)(y+0.86602540378*step), 0.0f );
        glVertex3f( (float)(x+1.0*step ), (float) y               , 0.0f );
    }else{
        glVertex3f( (float) x,           (float)y,                 0.0f );
        glVertex3f( (float)(x+step),     (float)y,                 0.0f );
        glVertex3f( (float)(x+0.5*step), (float)(y+0.86602540378*step), 0.0f );
    };
    glEnd( );
};

void Draw2D::drawSimplexGrid( int n, float step ){
    //glColor3f(0.1f,0.1f,0.1f);
    glBegin( GL_LINES );
    int n2 = 2 * n;
    float stepy = step*0.86602540378f;
    for( int i=-n; i<=n; i++ ){
        glVertex3f( -n*step,  i*stepy,  0.0f);
        glVertex3f(  n*step,  i*stepy,  0.0f);
    }
    for( int i=-n/2; i<=n/2; i++ ){
        glVertex3f( (i-n*0.5f)*step, -n*stepy,  0.0f);
        glVertex3f( (i+n*0.5f)*step,  n*stepy,  0.0f);

        glVertex3f( (i+n*0.5f)*step, -n*stepy,  0.0f);
        glVertex3f( (i-n*0.5f)*step,  n*stepy,  0.0f);
    }

    for( int i=-n/2; i<=n/2; i++ ){
        glVertex3f( -n*step,        -2*i*stepy,  0.0f);
        glVertex3f( (i-n*0.5f)*step,   n*stepy,  0.0f);

        glVertex3f(  n*step,        -2*i*stepy,  0.0f);
        glVertex3f( (i+n*0.5f)*step,  -n*stepy,  0.0f);

        glVertex3f(  n*step,         2*i*stepy,  0.0f);
        glVertex3f( (i+n*0.5f)*step,   n*stepy,  0.0f);

        glVertex3f( -n*step,         2*i*stepy,  0.0f);
        glVertex3f( (i-n*0.5f)*step,  -n*stepy,  0.0f);
    }

    glEnd();
}

void Draw2D::drawShape( const Vec2d& pos, const Vec2d& rot, int shape ){
	glPushMatrix();
	//glTranslatef( pos.x, pos.y , 0 );
	//glRotatef( phi*(180/M_PI), 0, 0, 1 );
	float glMat[16];
	toGLMat( pos, rot, glMat );
	glMultMatrixf( glMat );
	glCallList( shape );
	glPopMatrix();
};


// ===== image and sprite-text


void Draw2D::renderImage( GLuint itex, const Rect2d& rec ){
    glEnable( GL_TEXTURE_2D );
    glBindTexture( GL_TEXTURE_2D, itex );
    glColor3f(1.0f,1.0f,1.0f);
    //printf( " itex %i \n", itex );
    glBegin(GL_QUADS);
        glTexCoord2f( 0.0f, 1.0f ); glVertex3f( rec.a.x, rec.a.y, 3.0f );
        glTexCoord2f( 1.0f, 1.0f ); glVertex3f( rec.b.x, rec.a.y, 3.0f );
        glTexCoord2f( 1.0f, 0.0f ); glVertex3f( rec.b.x, rec.b.y, 3.0f );
        glTexCoord2f( 0.0f, 0.0f ); glVertex3f( rec.a.x, rec.b.y, 3.0f );
    glEnd();
};

void Draw2D::drawString( const char * str, int imin, int imax, float x, float y, float sz, int itex ){
    const int nchars = 95;
    float persprite = 1.0f/nchars;
    glEnable     ( GL_TEXTURE_2D );
    glBindTexture( GL_TEXTURE_2D, itex );
    //glColor4f(0.0f,0.5f,0.0f,1.0f);
    glEnable(GL_BLEND);
    glEnable(GL_ALPHA_TEST);
    //glBlendFunc( GL_ONE, GL_ZERO );
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    //glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    //glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
    //glBlendFunc(GL_ONE_MINUS_DST_ALPHA,GL_DST_ALPHA);
    glBegin(GL_QUADS);
    for(int i=imin; i<imax; i++){
        int isprite = str[i] - 33;
        float offset  = isprite*persprite+(persprite*0.57);
        float xi = i*sz + x;
        glTexCoord2f( offset          , 1.0f ); glVertex3f( xi,    y,      3.0f );
        glTexCoord2f( offset+persprite, 1.0f ); glVertex3f( xi+sz, y,      3.0f );
        glTexCoord2f( offset+persprite, 0.0f ); glVertex3f( xi+sz, y+sz*2, 3.0f );
        glTexCoord2f( offset          , 0.0f ); glVertex3f( xi,    y+sz*2, 3.0f );
    }
    glEnd();
    glDisable  ( GL_BLEND );
    glDisable  ( GL_ALPHA_TEST );
    glDisable  ( GL_TEXTURE_2D );
    glBlendFunc( GL_ONE, GL_ZERO );

    /*
    glColor4f(0.0f,1.0f,0.0f,1.0f);
    glEnable(GL_BLEND);
    glEnable(GL_ALPHA_TEST);
    //glBlendFunc( GL_ONE, GL_ZERO );
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    //glBlendFunc(GL_ONE_MINUS_DST_ALPHA,GL_DST_ALPHA);
    glBegin(GL_QUADS);
    for(int i=imin; i<imax; i++){
        int isprite = str[i] - 33;
        float offset  = isprite*persprite+(persprite*0.57);
        float xi = i*sz + x;
        glTexCoord2f( offset          , 1.0f ); glVertex3f( xi,    y,    3.0f );
        glTexCoord2f( offset+persprite, 1.0f ); glVertex3f( xi+sz, y,    3.0f );
        glTexCoord2f( offset+persprite, 0.0f ); glVertex3f( xi+sz, y+sz*2, 3.0f );
        glTexCoord2f( offset          , 0.0f ); glVertex3f( xi,    y+sz*2, 3.0f );
    }
    glEnd();
    glDisable  ( GL_BLEND );
    glDisable  ( GL_ALPHA_TEST );
    glDisable  ( GL_TEXTURE_2D );
    glBlendFunc( GL_ONE, GL_ZERO );
    */
}

void Draw2D::drawString( const  char * str, float x, float y, float sz, int itex ){
    const int nchars = 95;
    float persprite = 1.0f/nchars;
    glEnable( GL_TEXTURE_2D );
    glBindTexture( GL_TEXTURE_2D, itex );
    glEnable(GL_BLEND);
    glEnable(GL_ALPHA_TEST);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    //glColor4f(0.0f,1.0f,0.0f,1.0f);
    glBegin(GL_QUADS);
    for(int i=0; i<65536; i++){
        if( str[i] == 0 ) break; // 0-terminated string
        int isprite = str[i] - 33;
        float offset  = isprite*persprite+(persprite*0.57);
        float xi = i*sz + x;
        glTexCoord2f( offset          , 1.0f ); glVertex3f( xi,    y,    3.0f );
        glTexCoord2f( offset+persprite, 1.0f ); glVertex3f( xi+sz, y,    3.0f );
        glTexCoord2f( offset+persprite, 0.0f ); glVertex3f( xi+sz, y+sz*2, 3.0f );
        glTexCoord2f( offset          , 0.0f ); glVertex3f( xi,    y+sz*2, 3.0f );
    }
    glEnd();
    glDisable( GL_TEXTURE_2D );
    glDisable  ( GL_BLEND );
    glDisable  ( GL_ALPHA_TEST );
    glDisable  ( GL_TEXTURE_2D );
    glBlendFunc( GL_ONE, GL_ZERO );
}


void Draw2D::drawText( const char * str, const Vec2d& pos, int fontTex, float textSize, int istart, int iend ){
    glDisable    ( GL_LIGHTING   );
    glDisable    ( GL_DEPTH_TEST );
    glShadeModel ( GL_FLAT       );
    glPushMatrix();
        //glMatrixMode(GL_MODELVIEW);
        //glMatrixMode(GL_PROJECTION);
        glTranslatef( pos.x, pos.y, z_layer );
        //Draw::billboardCam( );
        //Draw::billboardCamProj( );
        //Draw2D::drawString( inputText.c_str(), 0, 0, textSize, fontTex );
        Draw::drawText( str, fontTex, textSize, istart, iend );
    glPopMatrix();
};




//}; // namespace Draw2D
