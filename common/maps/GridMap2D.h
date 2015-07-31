
#include <forward_list>


class Point2D{
	public:
	double x, y;
	bool equals( const Point2D& p ){ return ( x == p.x ) && ( y == p.y ); };
};

class Line2D{
	public:
	int id;
	Point2D *a,*b;
	
	Line2D( int id_, Point2D* a_, Point2D* b_ ){ id=id_; a=a_; b=b_; };

	void draw( ){ 
		glColor3f( randFuncf( id ), randFuncf( id+ 16874 ), randFuncf( id+ 98774 )  ); 
		glBegin( GL_LINES );
			glVertex3f( (float) a->x, (float) a->y, 0 );
			glVertex3f( (float) b->x, (float) b->y, 0 );

			glVertex3f( (float) a->x +0.1, (float) a->y    , 0 );
			glVertex3f( (float) a->x -0.1, (float) a->y    , 0 );
			glVertex3f( (float) a->x     , (float) a->y+0.1, 0 );
			glVertex3f( (float) a->x     , (float) a->y-0.1, 0 );
		glEnd();
	};
};

class Triangle2D{
	public:
	int id;
	Point2D *a,*b,*c;
	
	Triangle2D( int id_, Point2D* a_, Point2D* b_, Point2D* c_ ){ id=id_; a=a_; b=b_; c=c_; };

	void draw( ){ 
		glColor3f( randFuncf( id ), randFuncf( id+ 16874 ), randFuncf( id+ 98774 )  ); 
		glBegin( GL_LINE_LOOP );
			glVertex3f( (float) a->x, (float) a->y, 0 );
			glVertex3f( (float) b->x, (float) b->y, 0 );
			glVertex3f( (float) c->x, (float) c->y, 0 );
		glEnd();
	};
};



// ========== GridMap


template <class TYPE1 > class GridMap2D {
	public:
	double step, invStep;
	int nx, ny, nxy;
	bool isStatic = false;
	
	std::forward_list<TYPE1*>** store;
	TYPE1*** 					store_static;
	int*						ns;
	
	void init( int nx_, int ny_, double step_ ){
		step    = step_; 
		invStep = 1/step;
		nx = nx_; ny=ny_;
		nxy = nx*ny;
		ns     = new int[ nxy ];
		store  = new std::forward_list<TYPE1*>*[ nxy ];
		for (int i=0; i<nxy; i++){ store[i] = new std::forward_list<TYPE1*>(); ns[i] = 0; }
	};

	inline int getIx( double x ){ return (int)( invStep * x ); };
	inline int getIy( double y ){ return (int)( invStep * y ); };

	inline int getIndex( int ix, int iy     ){ return nx*iy + ix; }
	inline int getIndex( double x, double y ){ return getIndex( getIx(x), getIy(y) ); };

	void makeStatic( ){
		store_static = new TYPE1**[nxy];
		for (int i=0; i<nxy; i++){ 
			int ni = ns[i];
			if ( ni>0 ){ 
				//printf( " ni %i \n", ni );
				TYPE1** arr = new TYPE1*[ ni ];
				store_static[ i ] = arr;
				auto list = store [ i ];
				int j = 0;
				for ( auto it = list->cbegin(); it != list->cend(); ++it){
					arr[ j ] = *it;
					j++;
				}  
			}
		}
		isStatic = true;
		// dealocation
		for (int i=0; i<nxy; i++){ delete store[i]; }
		delete [] store;
	}


	inline void insert( TYPE1* p, int i ){
		store  [ i ] -> push_front( p );
		ns     [ i ] ++;
	};
	inline void insert( TYPE1* p, int ix, int iy ){ insert( p, getIndex( ix,iy ) ); };
	
	void plot( int ix, int iy ){ 
		glRect2D( ix*step, iy*step, (ix+1)*step, (iy+1)*step );
	}

	// look also
	// http://www.flipcode.com/archives/Raytracing_Topics_Techniques-Part_4_Spatial_Subdivisions.shtml
	// http://stackoverflow.com/questions/10350258/find-all-tiles-intersected-by-line-segment/10350503#10350503
	// http://www.cse.yorku.ca/~amana/research/grid.pdf

	//  http://stackoverflow.com/questions/27719906/triangle-rasterization-outer-bound-all-cells-which-contian-a-piece
/*
	int insertTriangle( Triangle2D* t ){
		Point2D* a = t->a;
		Point2D* b = t->b;
		Point2D* c = t->c;
		if( b->y > a->y ) { Point2D* tmp = a; a = b; b = tmp; }
		if( c->y > a->y ) { Point2D* tmp = a; a = c; c = tmp; }
		if( c->y > b->y ) { Point2D* tmp = b; b = c; c = tmp; }
		double xa  = a->x;
		double ya  = a->y;
		double xb  = b->x;
		double yb  = b->y;
		double xc  = c->x;
		double yc  = c->y;
	glColor3f( 0.2f*randFuncf( t->id ), 0.2f*randFuncf( t->id+ 16874 ), 0.2f*randFuncf( t->id+ 98774 )  ); 
		int   ixa  = getIx( xa );
		int   iya  = getIy( ya );
		int   ixb  = getIx( xb );
		int   iyb  = getIy( yb );
		int   ixc  = getIx( xc );
		int   iyc  = getIy( yc );
		double dxab = fabs( xb - xa );
		double dyab = fabs( yb - ya );
		double dxac = fabs( xc - xa );
		double dyac = fabs( yc - ya );

		int dixab = ( xa < xb ) ? 1 : -1;
		//int diyab = ( ya < yb ) ? 1 : -1;
		int dixac = ( xa < xc ) ? 1 : -1;
		//int diyac = ( ya < yc ) ? 1 : -1;

		double xab=0,yab=0,xac=0,yac=0; 
		int   ixab=ixa,ixac=ixa;
		int  iy=iya;
		int i=0;
		while ( iy != iyb ) {
			while ( xab < yab ) { xab += dyab;	ixab += dixab;	plot( ixab, iy ); } 
			while ( xac < yac ) { xac += dyac;	ixac += dixac;	plot( ixac, iy ); } 			
			iy--;
			yab += dxab;  plot( ixab, iy );
			yac += dxac;  plot( ixac, iy );
			//
			i++;
			if(i>100) break;
		}
	};
*/

	int insertTriangle( Triangle2D* t ){
		Point2D* a = t->a;
		Point2D* b = t->b;
		Point2D* c = t->c;
		if( b->y < a->y ) { Point2D* tmp = a; a = b; b = tmp; }
		if( c->y < a->y ) { Point2D* tmp = a; a = c; c = tmp; }
		if( c->y < b->y ) { Point2D* tmp = b; b = c; c = tmp; }
		double xa  = a->x;
		double ya  = a->y;
		double xb  = b->x;
		double yb  = b->y;
		double xc  = c->x;
		double yc  = c->y;
	glColor3f( 0.2f*randFuncf( t->id ), 0.2f*randFuncf( t->id+ 16874 ), 0.2f*randFuncf( t->id+ 98774 )  ); 
		int   ixa  = getIx( xa );
		int   iya  = getIy( ya );
		int   ixb  = getIx( xb );
		int   iyb  = getIy( yb );
		int   ixc  = getIx( xc );
		int   iyc  = getIy( yc );

		// up pass
		double dab = ( xb - xa )/( yb - ya );   double cab = xa - dab * ya;
		double dac = ( xc - xa )/( yc - ya );   double cac = xa - dac * ya;
		double y = iya * step;
		int  oixab = ixa, oixac = ixa, iy = iya; 
		while ( iy <= iyb ) {			
		//while ( true ) {
			y += step;
			double xab = dab * y + cab;  int ixab = getIx( xab ); if( iy == iyb ) ixab = ixb;
			double xac = dac * y + cac;  int ixac = getIx( xac );
			int ix1,ix2;
			if( dab < dac ){  
				ix1 = ( ixab < oixab ) ? ixab : oixab;
				ix2 = ( ixac > oixac ) ? ixac : oixac;
			}else{ 
				ix1 = ( ixac < oixac ) ? ixac : oixac;
				ix2 = ( ixab > oixab ) ? ixab : oixab;
			}	
			for ( int ix = ix1; ix <= ix2; ix++ ){	plot( ix, iy );  } // plot that line
			oixab = ixab; oixac = ixac;
			iy++;
		}
		// down pass
		double dca = ( xa - xc )/( ya - yc );   double cca = xa - dca * ya;
		double dcb = ( xb - xc )/( yb - yc );   double ccb = xb - dcb * yb;
		y   = iyc * step;
		int  oixca = ixc, oixcb = ixc; iy = iyc; 
		while ( iy > iyb ) {			
		//while ( true ) {
			double xca = dca * y + cca;  int ixca = getIx( xca );
			double xcb = dcb * y + ccb;  int ixcb = getIx( xcb );
			int ix1,ix2;
			if( dca > dcb ){  
				ix1 = ( ixca < oixca ) ? ixca : oixca;
				ix2 = ( ixcb > oixcb ) ? ixcb : oixcb;
			}else{ 
				ix1 = ( ixcb < oixcb ) ? ixcb : oixcb;
				ix2 = ( ixca > oixca ) ? ixca : oixca;
			}	
			printf( " ix1 ix2 %i %i \n ", ix1, ix2 );
			for ( int ix = ix1; ix <= ix2; ix++ ){	plot( ix, iy );  } // plot that line
			oixca = ixca; oixcb = ixcb; 
			y -= step;
			iy--;
		}
	};


	//  http://stackoverflow.com/questions/13542925/line-rasterization-4-connected-bresenham/27719652#27719652
	int insertLine( Line2D* l ){
		Point2D* a = l->a;
		Point2D* b = l->b;
		//if( b->x < a->x ) { Point2D* tmp = a; a = b; b = tmp;  }
		double ax  = a->x;
		double ay  = a->y;
		double bx  = b->x;
		double by  = b->y;

		double dx = fabs( bx - ax );
		double dy = fabs( by - ay );
		int dix = ( ax < bx ) ? 1 : -1;
		int diy = ( ay < by ) ? 1 : -1;
		int ix    = getIx( ax );
		int iy    = getIy( ay );
		int ixb   = getIx( bx );
		int iyb   = getIy( by );
		double x=0, y=0;
		int i=0;
		printf( " === dx dy %f %f \n", dx, dy );
		//glColor3f( 0.2f*randFuncf( l->id ), 0.2f*randFuncf( l->id+ 16874 ), 0.2f*randFuncf( l->id+ 98774 )  ); 
		//glRect2D( ix*step, iy*step, (ix+1)*step, (iy+1)*step );
		//glRect2D( ixb*step, iyb*step, (ixb+1)*step, (iyb+1)*step );
		insert( (TYPE1*)l, ix, iy   );
		insert( (TYPE1*)l, ixb, iyb );
		while ( ( ix != ixb ) && ( iy != iyb  ) ) {
			if ( x < y ) {
				x  += dy;
		        ix += dix;
 			} else {
				y  += dx;
		        iy += diy;
		    }
			insert( (TYPE1*)l, ix, iy );
			//glRect2D( ix*step, iy*step, (ix+1)*step, (iy+1)*step );
			//i++;
			//if(i>30) break;
		}

	};

// ======= just testing

	void paintGridLines(){
		glBegin( GL_LINES );
		glColor3f ( 0.1,0.1,0.1 );
		for (int i=0; i<=nx; i++){ glVertex3f( (float)(i*step), 0 , 0 ); glVertex3f( (float)(i*step),  (float)(ny*step), 0 ); };
		for (int i=0; i<=ny; i++){ glVertex3f( 0, (float)(i*step) , 0 ); glVertex3f( (float)(nx*step), (float)(i*step), 0 ); };
		glEnd( );
	};

	void paintStorage(){
		
		for (int iy=0; iy<ny; iy++){	for ( int ix=0; ix<nx; ix++  ){
			float red=0,green=0,blue=0;
			int  i = getIndex( ix, iy ); 
			auto list = store[i];
			for ( auto it = list->cbegin(); it != list->cend(); ++it ){
				TYPE1* p = *it;
				int   id = p->id;
				red    += 0.2f* randFuncf( id        ); 
				green  += 0.2f* randFuncf( id+ 16874 );
				blue   += 0.2f* randFuncf( id+ 98774 ); 
			}
			glColor3f( red, green, blue );
			glRect2D( ix*step, iy*step, (ix+1)*step, (iy+1)*step );
		} };

	};

};



