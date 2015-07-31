


void printVec( const Vec2d& v ){
	printf( " %f %f \n", v.x, v.y );
};


void drawPoint( const Vec2f& vec ){
	glDisable (GL_LIGHTING);
	glBegin   (GL_POINTS);	          	     
		glVertex3f( vec.x, vec.y, 0 );
	glEnd();
}

void drawPointCross( const Vec2f& vec, float d ){
	glDisable (GL_LIGHTING);
	glBegin   (GL_LINES);	          	     
		glVertex3f( vec.x-d, vec.y,   0   );  glVertex3f( vec.x+d, vec.y,   0 );
		glVertex3f( vec.x,   vec.y-d, 0   );  glVertex3f( vec.x,   vec.y+d, 0  );
	glEnd();
}

void drawVec( const Vec2f& vec ){
	glDisable (GL_LIGHTING);
	glBegin   (GL_LINES);	          	     
		glVertex3f( 0, 0, 0 ); glVertex3f( vec.x, vec.y, 0 );
	glEnd();
}

void drawVecInPos( const Vec2f& v, const Vec2f& pos ){
	glDisable (GL_LIGHTING);
	glBegin   (GL_LINES);	          	     
		glVertex3f( pos.x, pos.y, 0 ); glVertex3f( pos.x+v.x, pos.y+v.y, 0 );
	glEnd();
}

void drawLine( const Vec2f& p1, const Vec2f& p2 ){
	glDisable (GL_LIGHTING);
	glBegin   (GL_LINES);	          	     
		glVertex3f( p1.x, p1.y, 0 ); glVertex3f( p2.x, p2.y, 0 );
	glEnd();
}


void drawTriangle( const Vec2f& p1, const Vec2f& p2, const Vec2f& p3 ){
	glBegin   (GL_TRIANGLES);
		glNormal3f( 0, 0, 1 );	          	     
		glVertex3f( p1.x, p1.y, 0 ); 
		glVertex3f( p2.x, p2.y, 0 );
		glVertex3f( p3.x, p3.y, 0 );
	glEnd();
}


void drawPoint      ( const Vec2d& vec                   ){ Vec2f vec_;    convert( vec, vec_ ); drawPoint     ( vec_ );          }
void drawPointCross ( const Vec2d& vec, double d         ){ Vec2f vec_;    convert( vec, vec_ ); drawPointCross( vec_, (float)d); }
void drawVec        ( const Vec2d& vec                   ){ Vec2f vec_;    convert( vec, vec_ ); drawVec       ( vec_ );          }
void drawVecInPos   ( const Vec2d& v,   const Vec2d& pos ){ Vec2f v_,pos_; convert( v, v_ );     convert( pos, pos_ ); drawVecInPos( v_, pos_); }
void drawLine       ( const Vec2d& p1,  const Vec2d& p2  ){ Vec2f p1_,p2_; convert( p1, p1_ );   convert( p2, p2_ );   drawLine( p1_, p2_); }

void drawTriangle( const Vec2d& p1,  const Vec2d& p2, const Vec2d& p3 ){ Vec2f p1_,p2_,p3_;  convert( p1, p1_ );  convert( p2, p2_ ); convert( p3, p3_ ); drawTriangle( p1_, p2_, p3_ );  }

void drawLines( int nlinks, int * links, Vec2d * points ){
	int n2 = nlinks<<1;
	for( int i=0; i<n2; i+=2 ){
		drawLine( points[links[i]], points[links[i+1]] );
		//printf ( " %i %i %i %f %f \n", i, links[i], links[i+1], points[links[i]].x, points[links[i+1]].x );
	}
}






