
#ifndef  GridMap2D_h
#define  GridMap2D_h

#include <vector>

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"


// ==============================================
// ========== CLASS STUB : GridMap2D_stub
// ==============================================

class GridMap2D_stub {
    public:
	double    step, invStep;
	int       nx, ny, nxy;

	inline int getIx( double x ){ return (int)( invStep * x ); };
	inline int getIy( double y ){ return (int)( invStep * y ); };

	inline double getX( int ix ){ return step * ix ; };
	inline double getY( int iy ){ return step * iy ; };

	inline int getIndex ( int ix,   int iy   ){ return nx*iy + ix;                     };
	inline int getIndex ( double x, double y ){ return getIndex( getIx(x), getIy(y) ); };

    inline void init_stub( int nx_, int ny_, double step_ ){
        step    = step_;
		invStep = 1/step;
		nx = nx_; ny=ny_;
		nxy = nx*ny;
    }

};


// ==============================================
// ========== CLASS TEMPLATE : GridMap2D
// ==============================================

template <class OBJECT, class TILE = std::vector<OBJECT*> >
class GridMap2D : public GridMap2D_stub {
	public:
	TILE *    tiles;

    inline TILE* getTile( double x, double y ){ return  tiles + getIndex( x, y );      };

	inline void insert( OBJECT* p, int i ){
		tiles[ i ] -> push_back( p );
	}
	inline void insert( OBJECT* p, int ix, int iy ){ insert( p, getIndex( ix,iy ) ); };

    // this is most general and safest method for inserting 3D objects
    inline int insert( OBJECT* p ){
        int n_inserts = 0;
        Rect2d bbox;
        p.boundingBox( &bbox );
        int ix0 = getIx( bbox.x0 ); // TODO: bound check ?
        int iy0 = getIy( bbox.y0 );
        int ix1 = getIx( bbox.x1 );
        int iy1 = getIy( bbox.y1 );
        for( int iy=iy0; iy<=iy1; iy++ ){
            for( int ix=ix0; ix<=ix1; ix++ ){
                n_inserts++;
                insert( p, ix, iy );
                //int i = getIndex( ix, iy );
                //tiles[i].push_back( p );
            }
        }
        return n_inserts;
    }

	void init( int nx_, int ny_, double step_, int tile_n0 ){
        init_stub( nx_, ny_, step_ );
		tiles  = new TILE[ nxy ];
		for (int i=0; i<nxy; i++){
			tiles[i] = new TILE( );
			if ( tile_n0 != 0 ){
				tiles[i].reserve( tile_n0 );
			}
		}
	}

};


// ==================================================================================
// ========== CLASS SPECIALIZATION : GridMap2D< Segment2d, std::vector<Segment2d> >
// ==================================================================================

template<> class GridMap2D< Segment2d, std::vector<Segment2d*> > : public GridMap2D_stub {
    public:
    std::vector<Segment2d*> *    tiles;

    inline void insert( Segment2d * p, int i ){
		tiles[ i ].push_back( p );
	}
	inline void insert( Segment2d* p, int ix, int iy ){ insert( p, getIndex( ix,iy ) ); };

    inline int insert( Segment2d* l ){
        int n_inserts = 0;
        Vec2d* a = l->a;
        Vec2d* b = l->b;
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
        //printf( " === dx dy %f %f \n", dx, dy );
        //glColor3f( 0.2f*randFuncf( l->id ), 0.2f*randFuncf( l->id+ 16874 ), 0.2f*randFuncf( l->id+ 98774 )  );
        //glRect2D( ix*step, iy*step, (ix+1)*step, (iy+1)*step );
        //glRect2D( ixb*step, iyb*step, (ixb+1)*step, (iyb+1)*step );
        insert( l, ix, iy   );
        insert( l, ixb, iyb );
        while ( ( ix != ixb ) && ( iy != iyb  ) ) {
            if ( x < y ) {
                x  += dy;
                ix += dix;
            } else {
                y  += dx;
                iy += diy;
            }
            insert( l, ix, iy );
            n_inserts++;
            //int index = getIndex( ix, iy );
            //tiles[ index ].push_back( p );
            //glRect2D( ix*step, iy*step, (ix+1)*step, (iy+1)*step );
            //i++;
            //if(i>30) break;
        }
        return n_inserts;
    }

    void init( int nx_, int ny_, double step_, int tile_n0 ){
        init_stub( nx_, ny_, step_ );
		tiles  = new std::vector<Segment2d*>[ nxy ];
		for (int i=0; i<nxy; i++){
			//tiles[i] = new std::vector<Segment2d*>( );
			if ( tile_n0 != 0 ){
				tiles[i].reserve( tile_n0 );
			}
		}
	}

};

// ==================================================================================
// ========== CLASS SPECIALIZATION : GridMap2D< Triangle2d, std::vector<Triangle2d> >
// ==================================================================================

template<> class GridMap2D< Triangle2d, std::vector<Triangle2d*> > : public GridMap2D_stub {

    public:
    std::vector<Triangle2d*> *    tiles;

    inline void insert( Triangle2d * p, int i ){
        tiles[ i ].push_back( p );
    }
    inline void insert( Triangle2d* p, int ix, int iy ){ insert( p, getIndex( ix,iy ) ); };

	int insert( Triangle2d * t ){
		int n_inserts = 0;
		Vec2d* a = t->a;
		Vec2d* b = t->b;
		Vec2d* c = t->c;
		if( b->y < a->y ) { Vec2d* tmp = a; a = b; b = tmp; }
		if( c->y < a->y ) { Vec2d* tmp = a; a = c; c = tmp; }
		if( c->y < b->y ) { Vec2d* tmp = b; b = c; c = tmp; }
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
                insert( t, ix, iy );
                n_inserts++;
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
                insert( t, ix, iy );
                n_inserts++;
				//plot( ix, iy );
			}
			oixca = ixca; oixcb = ixcb;
			y -= step;
			iy--;
		}
		return n_inserts;
	}

    void init( int nx_, int ny_, double step_, int tile_n0 ){
        init_stub( nx_, ny_, step_ );
		tiles  = new std::vector<Triangle2d*>[ nxy ];
		for (int i=0; i<nxy; i++){
			//tiles[i] = new std::vector<Segment2d*>( );
			if ( tile_n0 != 0 ){
				tiles[i].reserve( tile_n0 );
			}
		}
	}

};

#endif

