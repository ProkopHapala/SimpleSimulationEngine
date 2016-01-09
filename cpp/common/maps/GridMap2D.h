
#ifndef  GridMap2D_h
#define  GridMap2D_h

#include <forward_list>

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"

// ==============================================
// ========== CLASS TEMPLATE : GridMap2D
// ==============================================

template <class OBJECT, class TILE> 
class GridMap2D {
	public:
	double    step, invStep;
	int       nx, ny, nxy;
	bool      isStatic = false;
	TILE *    tiles;

	// ==== function declarations

	//void makeStatic    (               );
	//int  insertTriangle( Triangle2D* t );
	//int  insertLine    ( Line2D* l     );

	// ==== inline functions
	
	inline int getIx( double x ){ return (int)( invStep * x ); };
	inline int getIy( double y ){ return (int)( invStep * y ); };

	inline double getX( int ix ){ return Step * ix ); };
	inline double getY( int iy ){ return Step * iy ); };

	inline int getIndex( int ix,   int iy   ){ return nx*iy + ix; }
	inline int getIndex( double x, double y ){ return getIndex( getIx(x), getIy(y) ); };

/*
	inline void insert( OBJECT* p, int i ){
		tiles[ i ] -> push_back( p );
	}
	inline void insert( OBJECT* p, int ix, int iy ){ insert( p, getIndex( ix,iy ) ); };
*/

	void init( int nx_, int ny_, double step_, int tile_n0 ){
		step    = step_; 
		invStep = 1/step;
		nx = nx_; ny=ny_;
		nxy = nx*ny;
		ns     = new int[ nxy ];
		store  = new std::forward_list<TYPE1*>*[ nxy ];
		for (int i=0; i<nxy; i++){ 	
			store[i] = new TILE( ); 
			if ( tile_n0 != 0 ){
				store[i].reserve( tile_n0 );
			}
		}
	}
	
};


// ==================================================================================
// ========== CLASS SPECIALIZATION : GridMap2D< Segment2d, std::vector<Segment2d> >
// ==================================================================================

int GridMap2D::insertLine( Segment2d* l ){
	Vec2D* a = l->a;
	Vec2D* b = l->b;
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


// ==================================================================================
// ========== CLASS SPECIALIZATION : GridMap2D< Triangle2d, std::vector<Triangle2d> >
// ==================================================================================

template<> class GridMap2D< Triangle2d, std::vector<Triangle2d> >{
	
	inline void insert( Triangle2d * p, int i ){
		tiles[ i ] -> push_back( p );
	}
	inline void insert( Triangle2d * p, int ix, int iy ){ insert( p, getIndex( ix,iy ) ); };

	int insertTriangle( Triangle2d * t ){
		int n_inserts = 0;
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
		//glColor3f( 0.2f*randFuncf( t->id ), 0.2f*randFuncf( t->id+ 16874 ), 0.2f*randFuncf( t->id+ 98774 )  ); 
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
			for ( int ix = ix1; ix <= ix2; ix++ ){	
				n_inserts ++;
				insert( TYPE1* p, int ix, int iy );
				//plot( ix, iy );  
			} 
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
			//printf( " ix1 ix2 %i %i \n ", ix1, ix2 );
			for ( int ix = ix1; ix <= ix2; ix++ ){	
				n_inserts ++;
				insert( TYPE1* p, int ix, int iy );
				//plot( ix, iy );  
			}
			oixca = ixca; oixcb = ixcb; 
			y -= step;
			iy--;
		}
		return n_inserts;
	};

}





#endif 

