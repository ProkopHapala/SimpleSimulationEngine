
//#include "arrayAlgs.h"
#include "Convex2d.h" // THE HEADER

/*
bool Convex2d::lineIntersections( const Line2d& line, const Vec2d& p1, Vec2d& p2 ){
    int i0,i1;
    double sign0 = line.dist_unitary( corners[0] );
    double sign  = line.dist_unitary( corners[1] );
    for( i0=2; i0<n; i0++    ){  double sign  = line.dist_unitary( corners[i0] );  if( sign * sign0 ) < 0 ) break; }
    if( i0 >= n )return false;
    for( i1=i0+1; i1<n; i1++ ){  double sign  = line.dist_unitary( corners[i1] );  if( sign * sign0 ) > 0 ) break; }
    if( i1 >= n ) i1 = 0;

}
*/

bool Convex2d::cut( const Line2d& line, Convex2d& left, Convex2d& right ){
    //printf( " inside Convex2d::cut \n" );
    // ---- find index of corners left and right of line

    int i0,i1,i1_;
    double sign0 = line.dist_unitary( corners[0] );
    for( i0=1; i0<(n-1); i0++ ){ double sign  = line.dist_unitary( corners[i0] ); if( ( sign * sign0 ) < 0 ) break; }
    if ( i0 >= n )return false;
    for( i1=i0+1; i1<n;  i1++ ){ double sign  = line.dist_unitary( corners[i1] );  if( ( sign * sign0 ) > 0 ) break; }
    // ---- find intersection points
    if( i1 >= n ){ i1_ = 0; }else{ i1_ = i1; }
    Vec2d p1,p2;
    //printf( " %i %i   %i %i \n", i0-1, i0, i1, i1_  );
    line.intersectionPoint( corners[i0-1], corners[i0 ], p1 );
    line.intersectionPoint( corners[i1-1], corners[i1_], p2 );

    // ---- copy the points
    int nleft = i1-i0;
    int szleft  = 2 +       nleft  ;
    int szright = 2 + ( n - nleft );
    //printf( " %i %i   %i   %i %i %i \n", i0, i1, nleft,  n, szleft, szright  );
    //exit(0);
    left .reallocate( szleft  );
    right.reallocate( szright );
    //printf( " allocated \n" );
    right.corners[i0  ].set(p1);
    right.corners[i0+1].set(p2);
    left .corners[0].set(p2);
    left .corners[1].set(p1);
    int ileft  = 2;
    int iright = 0;
    for( int i=0; i<n; i++ ){
        if( ( i>=i0 ) && ( i<i1 ) ){
            left .corners[ileft ].set( corners[i] );
            ileft++;
        }else{
            if( iright == i0 ) iright+=2;
            right.corners[iright].set( corners[i] );
            iright++;
        }
    }
    //printf( " %i %i \n", ileft, iright );
    return false;
}


void Convex2d::make_corners( ){
    delete corners;
}


void Convex2d::fromPoints( int np, Vec2d * points, int * permut, double * vals ){
    // ---- find left-most pivot point ... to make sure that all angles are within (-pi/2,pi/2)
    int ipivot  = 0;
    double xmin = points[0].x;
    for( int i=1; i<np; i++ ){
        double xi = points[i].x;
        if( xmin > xi ){
            xmin   = xi;
            ipivot = i;
        }
    }
    // ---- order points by angle to x-axis
    int ii = 0;
    for( int i=0; i<np; i++ ){
        if( i != ipivot ){
            Vec2d d; d.set_sub( points[i], points[ipivot] );
            vals  [ii] = d.y / d.x;  // within (-pi/2,pi/2) if pivot is left-most point
            permut[ii] = ii;
            ii++;
        }
    }
    quickSort<double>( vals, permut, 0, np-1 );
    for( int i=0; i<(np-1); i++ ){ if( i<ipivot ) permut[i]++; }
	//permute  <double>( permut, A, out, 0, n );
	// ---- copy points, exclude inward kinks
	Vec2d *p1,*p2,*p3;
	p1 = points + ipivot;
	p2 = points + permut[ 0 ];
	p3 = points + permut[ 1 ];
	int nput = 0;
	int ip   = 0;
    for( int i=2; i<n; i++ ){
        double side  = line_side( *p2, *p1, *p3 );
        if( side > 0 ){
            permut[ nput ] = ip;
            nput++;
            p1 = p2;
            p2 = p3;
            p3 = points + permut[ i ];
        }else{
            p3 = points + permut[ i ];
        }
    }
    // finally we copy the points
    corners[ 0 ].set( ipivot );
    for( int i=1; i<nput; i++ ){
        corners[ i ].set( points[ permut[i-1] ] );
    }
    update_lines( );
}

void Convex2d::fromPoints( int np, Vec2d * points ){
	//int     * permut = new int   [ np ];
	//double  * vals   = new double[ np ];
	int     permut[np];
	double  vals  [np];
	fromPoints( np, points, permut, vals );
	//delete permut;
	//delete vals;
}

void Convex2d::projectToLine( const Vec2d& dir, double * xs, double * yLs, double * yRs ) const {   // FIXME move this to .cpp
    // -- order corners by pos along direction
    //printf( " DEBUG 0 \n" );
    double xs_   [ n ];
    int ibot,itop;
    double xbot=+1e+300,xtop=-1e+308;
    for( int i=0; i<n; i++ ){
        double xi = dir.dot( corners[i] );
        _minit( i, xi, ibot, xbot );   //if( xi<xbot ){ xbot=xi; ibot=i; }
        _maxit( i, xi, itop, xtop );   //if( x>xmax ){ xmax=x; imax=i; }
        xs_[i] = xi;
        //printf( " %i %f \n", i, xi );
    }
    //printf( " imin %i xmin %f \n", ibot, xbot );
    //printf( " imax %i xmax %f \n", itop, xtop );
    //printf( " DEBUG 1 \n" );
    // -- initialize left and right bonder from bottom point "ibot" to top point "itop"
    Vec2d oleft,oright,pleft,pright;
    oleft .set( xs_[ibot], dir.dot_perp( corners[ibot] ) );
    oright.set( oleft );
    xs [ 0 ]  = oleft.x;
    yRs[ 0 ]  = oleft.y;
    yLs[ 0 ]  = oleft.y;
    int index  = 0;
    int ileft  = ibot;  _circ_inc( ileft,  n );
    int iright = ibot;  _circ_dec( iright, n );
    pleft .set( xs_[ileft],  dir.dot_perp( corners[ileft ] ) );
    pright.set( xs_[iright], dir.dot_perp( corners[iright] ) );
    //printf( " DEBUG 3 \n" );
    // -- iterate over left and right border resolving points acording to its order along the direction
    do {
        index++;
        //printf( " DEBUG index %i %i %i %f %f \n", index, ileft, iright, pleft.x, pright.x );
        if( pleft.x < pright.x ){ // left is closer
            double yright = oright.y +  ( pleft.x - oright.x ) * ( pright.y - oright.y ) / ( pright.x - oright.x );
            yLs[ index ]  = pleft.y;
            yRs[ index ]  = yright;
            xs [ index ]  = pleft.x;
            oleft.set( pleft );
            _circ_inc( ileft,  n );
            pleft .set( xs_[ileft],  dir.dot_perp( corners[ileft ] ) );
            //printf( " left  %i %i \n", index, ileft );
            // FIXME : we should take care when it come to end ? probably not then ileft=itop
        }else{
            double yleft  = oleft.y +  ( pright.x - oleft.x ) * ( pleft.y - oleft.y ) / ( pleft.x - oleft.x );
            yLs[ index ]  = yleft;
            yRs[ index ]  = pright.y;
            xs [ index ]  = pright.x;
            oright.set( pright );
            _circ_dec( iright,  n );
            pright .set( xs_[iright],  dir.dot_perp( corners[iright ] ) );
            //printf( " right %i %i \n", index, iright );
        }
        if( index >= (n-1) ){ printf( " loop should end %i %i %i %i \n", index, n, ileft, iright ); break; } // FIXME DEBUG just to prevent infinite loop
    } while( !( ( itop == ileft ) && ( itop == iright ) ) );
    //printf( " index %i ileft %i iright %i itop %i \n", index, ileft, iright, itop );
    index = n-1;
    yLs[ index ] = dir.dot_perp( corners[ itop ] );
    yRs[ index ] = yLs[ index ];
    xs [ index ] = xs_[itop];
}













