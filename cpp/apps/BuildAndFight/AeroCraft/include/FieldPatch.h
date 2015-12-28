

inline float trinagle_area( float x12, float y12, float x13, float y13 ) {  return x12*y13 - x13*y12; };

class Rect{
	public:
	float value;
	float x1,y1,z1;
	float x2,y2,z2;
	float x3,y3,z3;
	float x4,y4,z4;
	//float x1,x2,x3,x4;
	//float y1,y2,y3,y4;
	//float z1,z2,z3,z4;

/*
	void draw() const {
		glBegin(GL_QUADS);
			glColor3f( value, 0.5, 0.2 );		          	     
			glNormal3f(0,0,-1); 
			glVertex3f( x1, z1, y1 ); 
			glVertex3f( x2, z2, y2 ); 
			glVertex3f( x4, z4, y4 ); 
			glVertex3f( x3, z3, y3 ); 
		glEnd();
	};
*/

	void draw() const {
		glBegin(GL_QUADS);	     
			glNormal3f( 0, 1, 0 ); 
/*
			glColor3f ( 0.25+z1*0.02, 0.5, 0.2 );	 glVertex3f( x1, z1, y1 ); 
			glColor3f ( 0.25+z2*0.02, 0.5, 0.2 );	 glVertex3f( x2, z2, y2 ); 
			glColor3f ( 0.25+z4*0.02, 0.5, 0.2 );	 glVertex3f( x4, z4, y4 ); 
			glColor3f ( 0.25+z3*0.02, 0.5, 0.2 );	 glVertex3f( x3, z3, y3 ); 
*/
			glColor3f ( 0.5+z1*0.01, 0.5, 0.2 );	 glVertex3f( x1, z1, y1 ); 
			glColor3f ( 0.5+z2*0.01, 0.5, 0.2 );	 glVertex3f( x2, z2, y2 ); 
			glColor3f ( 0.5+z4*0.01, 0.5, 0.2 );	 glVertex3f( x4, z4, y4 ); 
			glColor3f ( 0.5+z3*0.01, 0.5, 0.2 );	 glVertex3f( x3, z3, y3 ); 
		glEnd();
	};

	inline float area() const {
		float x12 = x1 - x2; float y12 = y1 - y2;
		float x13 = x1 - x3; float y13 = y1 - y3;
		float x42 = x4 - x2; float y42 = y4 - y2;
		float x43 = x4 - x3; float y43 = y4 - y3;
		return abs( trinagle_area( x12, y12, x13, y13 ) ) + abs(trinagle_area( x42, y42, x43, y43 ) );
	};
};

namespace FieldPatch {  
	//float thresh[] = { 0.5, 0.9 };
	float thresh[] = { 0.6 };
	int   startkill = 4;
	float dval      = 0.2;
	float dmax      = 0.3;
	float minarea   = 1000;
	float dheight    =  100; 
	float dmheight   =   50;

	int quadCount;
	
	void divide_1( int nlevels, int level, const Rect& rc   );
	void divide  ( int nlevels, int level, const Rect& rect );

	void divide_1( int nlevels, int level, const Rect& rc ){
		int dlevel = nlevels - level;
		//rc.paint();
		//float dmax= FIELD_PATCH_dmax;
		//float dval= FIELD_PATCH_dval;
		float fx = 0.5+randf(-dmax,dmax); float mfx = 1-fx;
		float fy = 0.5+randf(-dmax,dmax); float mfy = 1-fy;
		float topx = (rc.x1+rc.x2)*0.5;   float topy = (rc.y1+rc.y2)*0.5; float topz = (rc.z1+rc.z2)*0.5;
		float botx = (rc.x3+rc.x4)*0.5;   float boty = (rc.y3+rc.y4)*0.5; float botz = (rc.z3+rc.z4)*0.5;
		float lftx = (rc.x1+rc.x3)*0.5;   float lfty = (rc.y1+rc.y3)*0.5; float lftz = (rc.z1+rc.z3)*0.5;
		float rgtx = (rc.x2+rc.x4)*0.5;   float rgty = (rc.y2+rc.y4)*0.5; float rgtz = (rc.z2+rc.z4)*0.5;
		float ctrx = mfy*(mfx*rc.x1+fx*rc.x2) + fy*(mfx*rc.x3+fx*rc.x4);   
		float ctry = mfy*(mfx*rc.y1+fx*rc.y2) + fy*(mfx*rc.y3+fx*rc.y4);
		//float ctrz = (topz+botz+lftz+rgtz)*0.25 + randf( -dheight, dheight );
		float szscale = ( 4.0/( dlevel + 0.5 ) );
		float ctrz = (topz+botz+lftz+rgtz)*0.25 + randf( -dmheight, dheight )*szscale;
		//float ctrz = (topz+botz+lftz+rgtz)*0.25 + 0;
		//float ctrz = (topz+botz+lftz+rgtz)*0.25 + 0;
		divide( nlevels, level-1, { rc.value + randf(-dval,dval),   rc.x1, rc.y1, rc.z1,   topx,  topy,  topz,    lftx,  lfty,  lftz,    ctrx,  ctry,  ctrz  } );
		divide( nlevels, level-1, { rc.value + randf(-dval,dval),   topx,  topy,  topz,    rc.x2, rc.y2, rc.z2,   ctrx,  ctry,  ctrz,    rgtx,  rgty,  rgtz  } );
		divide( nlevels, level-1, { rc.value + randf(-dval,dval),   lftx,  lfty,  lftz,    ctrx,  ctry,  ctrz,    rc.x3, rc.y3, rc.z3,   botx,  boty,  botz  } );
		divide( nlevels, level-1, { rc.value + randf(-dval,dval),   ctrx,  ctry,  ctrz,    rgtx,  rgty,  rgtz,    botx,  boty,  botz,    rc.x4, rc.y4, rc.z4 } );
	};

	void divide( int nlevels, int level, const Rect& rect ){
		float rnd  = randf(); 
		float area = rect.area();
		if((level>0)&&( area>minarea )){
			if      ( rnd < thresh[0] ){ divide_1( nlevels, level, rect ); }
			//else if ( rnd < FIELD_PATCH_thresh[1] ){ divide_2( nlevels, level, rect ); }
			else{ 
				if( level>(nlevels-startkill) ){ divide_1( nlevels, level, rect );  }
				else                   { rect.draw(); quadCount++;    }; 
			}; 
		}else{ rect.draw(); quadCount++; };
	};

	int makeList( int nlevels, const Rect& rect ){
		int ilist=glGenLists(1);
		quadCount = 0;
		glNewList( ilist, GL_COMPILE );
			divide( nlevels, nlevels, rect );
		glEndList();
		printf(" FieldPatch made %i quads \n", quadCount );
		return( ilist );
	}

/*
	void divide_2( int nlevels, int level, Rect rc ){
	  //rc.paint();
	  float rnd = randf(1);
	  if(rnd>0.5){
		float f1 = 0.5+randf(-maxdev,maxdev); float mf1 = 1-f1;
		float f2 = 0.5+randf(-maxdev,maxdev); float mf2 = 1-f2;
		float topx = (mf1*rc.x1+f1*rc.x2);   float topy = (mf1*rc.y1+f1*rc.y2); float topz = (mf1*rc.z1+f1*rc.z2);
		float botx = (mf2*rc.x3+f2*rc.x4);   float boty = (mf2*rc.y3+f2*rc.y4); float botz = (mf2*rc.z3+f2*rc.z4);
		divide( nlevels, level-1, new Rect( rc.value + randf(-dval,dval),   rc.x1, rc.y1, rc.y2,   topx,  topy,  topz,    rc.x3, rc.y3, rc.z3,   botx,  boty,  botz  ) );
		divide( nlevels, level-1, new Rect( rc.value + randf(-dval,dval),   topx,  topy,  topz,    rc.x2, rc.y2, rc.z2,   botx,  boty,  botz,    rc.x4, rc.y4, rc.y4 ) );
	  }else{
		float f1 = 0.5+randf(-maxdev,maxdev); float mf1 = 1-f1;
		float f2 = 0.5+randf(-maxdev,maxdev); float mf2 = 1-f2;
		float lftx = (mf1*rc.x1+f1*rc.x3);   float lfty = (mf1*rc.y1+f1*rc.y3); float lftz = (mf1*rc.z1+f1*rc.z3);
		float rgtx = (mf2*rc.x2+f2*rc.x4);   float rgty = (mf2*rc.y2+f2*rc.y4); float rgtz = (mf2*rc.z2+f2*rc.z4);
		divide( nlevels, level-1, new Rect( rc.value + randf(-dval,dval),   rc.x1, rc.y1, rc.z1,   rc.x2, rc.y2, rc.z2,   lftx,  lfty,  lftz,     rgtx,  rgty,  rgtz   ) );
		divide( nlevels, level-1, new Rect( rc.value + randf(-dval,dval),   lftx,  lfty,  lftz,    rgtx, rgty, rgtz,      rc.x3, rc.y3, rc.z3,    rc.x4, rc.y4, rc.y4  ) );
	  }
	};
*/

}






