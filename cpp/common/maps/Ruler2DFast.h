
#ifndef  Ruler2DFast_h
#define  Ruler2DFast_h

#include "Vec2.h"

const static int Ruler2D_nEdges = 4;
const static int Ruler2D_nVerts = 4;

// ===== class CubicRuler

// TODO:
//  * insertion of  LineSegment from GridMap2D
//  * rayMarching from CubicRuler and/or SimplexRuler
//



class Ruler2DFast { public:
    // --- variables
    Vec2d pos0    = (Vec2d){0.0,0.0};
	Vec2d step    = (Vec2d){1.0,1.0};
	Vec2d invStep = (Vec2d){1.0,1.0};

    Vec2i n       = (Vec2i){0,0};
	int   ntot    = 0;

	// --- inline functions

    inline void  setN(Vec2i n_)       { n=n_; ntot=n.x*n.y;                }
    inline Vec2i i2ip(int   i ) const { return {i%n.x,i/n.x};   } // https://stackoverflow.com/questions/7070346/c-best-way-to-get-integer-division-and-remainder
    inline int   ip2i(Vec2i ip) const { return (n.x*ip.y+ip.x); }

	inline int    x2i( double x  ) const { return (int)( invStep.x*(x - pos0.x ) ); }
	inline int    y2i( double y  ) const { return (int)( invStep.y*(y - pos0.y ) ); }
	inline double i2x( double ix ) const { return      (step.x*ix)   + pos0.x; }
	inline double i2y( double iy ) const { return      (step.y*iy)   + pos0.y; }

	//inline int dx ( double ix, double x ) const { return (x-pos0.x) - (ix*step.x);  }
	//inline int dxi( double ix, double x ) const { return invStep.x*(x-pos0.x) - ix; }

	inline int x2id( double x, double& dix ) const { double x_ = invStep.x*(x-pos0.x); int ix = (int)x_; dix = x_ - ix; return ix; }
	inline int y2id( double y, double& diy ) const { double y_ = invStep.y*(y-pos0.y); int iy = (int)y_; diy = y_ - iy; return iy; }

	inline int x2idx( double x, double& dx ) const { double x_ = invStep.x*(x-pos0.x); int ix = (int)x_; dx = step.x*(x - ix); return ix; }
	inline int y2idy( double y, double& dy ) const { double y_ = invStep.y*(y-pos0.y); int iy = (int)y_; dy = step.y*(y - iy); return iy; }

	inline void pos2index( const Vec2d& pos, Vec2d& dipos, Vec2i& ipos ) const {
		ipos.x = x2id( pos.x, dipos.x );
		ipos.y = y2id( pos.y, dipos.y );
	}

    inline void index2pos( const Vec2i& ipos, const Vec2d& dipos, Vec2d& pos ) const {
        pos.x = i2x( ipos.x + dipos.x );
        pos.y = i2y( ipos.y + dipos.y );
    }

    inline void setup( const Vec2d& pos0_, const Vec2d& step_ ){
        pos0  .set( pos0_ );
        step  .set( step_  );
		invStep.set_inv( step );
    }


    int getOverlapingTiles( Vec2d pos, double r, int * results ){
        Vec2d dpos;
        Vec2i ipos;
        pos2index( pos, dpos, ipos );
        dpos.mul( step );
        results[0] = ip2i(ipos); int nret=1;
        int dix=0,diy=0;
        double dr2   = 0;
        double mr    = 1-r;
        if     (  dpos.x < r  ){ results[nret]=ip2i( {ipos.x-1  , ipos.y    }); nret++; dix=-1; dr2 += sq(  dpos.x); }
        else if(  dpos.x > mr ){ results[nret]=ip2i( {ipos.x+1  , ipos.y    }); nret++; dix=+1; dr2 += sq(1-dpos.x); }
        if     (  dpos.y < r  ){ results[nret]=ip2i( {ipos.x    , ipos.y-1  }); nret++; diy=-1; dr2 += sq(  dpos.y); }
        else if(  dpos.y > mr ){ results[nret]=ip2i( {ipos.x    , ipos.y+1  }); nret++; diy=+1; dr2 += sq(1-dpos.y); }
        if     ( dr2 < (r*r)  ){ results[nret]=ip2i( {ipos.x+dix, ipos.y+diy}); }
        //if( (dix!=0)&&(diy!=0) ){ insert( o, ipos.x+dix, ipos.y+diy ); }
        //printf( " %1.3f %1.3f  (%1.3f,%1.3f) (%i,%i) %1.3f \n", r, mr, dpos.x,dpos.y, dix, diy, dr2 );
        return nret;
    }

    inline int insertSegment( int * results, Vec2d* a, Vec2d* b ){
        int n_inserts = 0;
        //Vec2d* a = l->a;
        //Vec2d* b = l->b;
        //if( b->x < a->x ) { Point2D* tmp = a; a = b; b = tmp;  }
        double ax  = a->x;
        double ay  = a->y;
        double bx  = b->x;
        double by  = b->y;

        double dx = fabs( bx - ax );
        double dy = fabs( by - ay );
        int dix = ( ax < bx ) ? 1 : -1;
        int diy = ( ay < by ) ? 1 : -1;
        int ix    = x2i( ax );
        int iy    = y2i( ay );
        int ixb   = x2i( bx );
        int iyb   = y2i( by );
        double x=0, y=0;
        int i=0;
        //printf( " === dx dy %f %f \n", dx, dy );
        //glColor3f( 0.2f*randFuncf( l->id ), 0.2f*randFuncf( l->id+ 16874 ), 0.2f*randFuncf( l->id+ 98774 )  );
        //glRect2D( ix*step, iy*step, (ix+1)*step, (iy+1)*step );
        //glRect2D( ixb*step, iyb*step, (ixb+1)*step, (iyb+1)*step );
        //insert( l, ix, iy   );
        //insert( l, ixb, iyb );
        results[n_inserts]=ip2i({ix, iy });  n_inserts++;
        results[n_inserts]=ip2i({ixb,iyb});  n_inserts++;
        while ( ( ix != ixb ) && ( iy != iyb  ) ) {
            if ( x < y ) {
                x  += dy;
                ix += dix;
            } else {
                y  += dx;
                iy += diy;
            }
            //insert( l, ix, iy );
            results[n_inserts]=ip2i({ix, iy }); n_inserts++;
            //int index = getIndex( ix, iy );
            //tiles[ index ].push_back( p );
            //glRect2D( ix*step, iy*step, (ix+1)*step, (iy+1)*step );
            //i++;
            //if(i>30) break;
        }
        return n_inserts;
    }

    int insertTriangle( int * results, Vec2d* a, Vec2d* b, Vec2d* c ){
		int n_inserts = 0;
		//Vec2d* a = t->a;
		//Vec2d* b = t->b;
		//Vec2d* c = t->c;
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
		int   ixa  = x2i( xa );
		int   iya  = y2i( ya );
		int   ixb  = x2i( xb );
		int   iyb  = y2i( yb );
		int   ixc  = x2i( xc );
		int   iyc  = y2i( yc );

		// up pass
		double dab = ( xb - xa )/( yb - ya );   double cab = xa - dab * ya;
		double dac = ( xc - xa )/( yc - ya );   double cac = xa - dac * ya;
		double y = iya * step.y;
		int  oixab = ixa, oixac = ixa, iy = iya;
		while ( iy <= iyb ) {
		//while ( true ) {
			y += step.y;
			double xab = dab * y + cab;  int ixab = x2i( xab ); if( iy == iyb ) ixab = ixb;
			double xac = dac * y + cac;  int ixac = x2i( xac );
			int ix1,ix2;
			if( dab < dac ){
				ix1 = ( ixab < oixab ) ? ixab : oixab;
				ix2 = ( ixac > oixac ) ? ixac : oixac;
			}else{
				ix1 = ( ixac < oixac ) ? ixac : oixac;
				ix2 = ( ixab > oixab ) ? ixab : oixab;
			}
			for ( int ix = ix1; ix <= ix2; ix++ ){

                //insert( t, ix, iy );
                //results[(n_inserts<<1)  ];
                //results[(n_inserts<<1)+1];
                results[n_inserts]=ip2i({ix,iy});
                n_inserts++;
				//plot( ix, iy );
			}
			oixab = ixab; oixac = ixac;
			iy++;
		}
		// down pass
		double dca = ( xa - xc )/( ya - yc );   double cca = xa - dca * ya;
		double dcb = ( xb - xc )/( yb - yc );   double ccb = xb - dcb * yb;
		y   = iyc * step.y;
		int  oixca = ixc, oixcb = ixc; iy = iyc;
		while ( iy > iyb ) {
		//while ( true ) {
			double xca = dca * y + cca;  int ixca = x2i( xca );
			double xcb = dcb * y + ccb;  int ixcb = x2i( xcb );
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
                //insert( t, ix, iy );
                results[n_inserts]=ip2i({ix,iy});
                n_inserts++;
				//plot( ix, iy );
			}
			oixca = ixca; oixcb = ixcb;
			y -= step.y;
			iy--;
		}
		return n_inserts;
	}

};

#endif

