#ifndef TerrainSimplex_h
#define TerrainSimplex_h

#include "Noise.h"

//typedef unsigned short  UHALF;
//const UHALF MAP_OFFSET = 0x7FFF;

const int MAP_OFFSET = 0x7FFF;

class TerrainSimplex{
public:

    static constexpr int     nneigh   = 8;
    static constexpr double  isqrt2  = 0.70710678118;
    static constexpr int     neighs     [nneigh*2] = { 0,-1, 0,1, -1,0, +1,0, -1,-1, -1,+1, +1,-1, +1+1 };
    static constexpr double  neigh_dists[nneigh  ] = { 1.0, 1.0, 1.0, 1.0, isqrt2, isqrt2, isqrt2, isqrt2 };

    double        step, invStep;

	int power;
	int mask;
	int nx, ny, ntot;
    int istep;

	double * ground;
	double * water;


	// countor step algorithm
	//int    nContourMax;
	int    nContour=0,nContour_old=0;
	bool   * known;
	int    * contour1;
	int    * contour2;
	double * water_;

    double c_errode    = 0.4;
    double c_rain      = 0.00;
    //double c_blur      = 0.0;
    double c_sediment  = 0.001;
    double f_sediment  = 0.1;
    double f_mix       = 0.1;
    double mf_mix      = 1.0-f_mix;

    int droplet_ix, droplet_iy;
    double droplet_size,droplet_h;

	// ==== function declaration

	int raster_line( Vec2d dirHat, Vec2d pos0, Vec2d pos1, Vec2d * hits, int * boundaries, int * edges );

    void genTerrainNoise( int n, double scale, double hscale, double fdown, double strength, int seed, const Vec2d& pos0 );
    void init_outflow();
    void outflow_step();
    void extend_path( float val, int oi, int i );

    void initErrosion( double w );
    void flow_errosion_step  ( );
    void rain_and_evaporation( );
    int  flow_errosion_step_noRain( );
    void initDroplet  ( double size_ );
    bool droplet_step ( );
    void errodeDroples( int n, int nStepMax, double size_ );

	// ==== inline functions


    inline void simplexIndexBare( double x, double y, int& ia, int& ib ) const {
        double a = ( invStep *    y * 1.15470053839       ) + MAP_OFFSET;
        double b = ( invStep * (  x - 0.57735026919*y   ) ) + MAP_OFFSET;
    	ia = (int)a;
        ib = (int)b;
    }

    inline bool simplexIndex( double x, double y, int& ia, int& ib, double& da, double& db ) const {
        double a = ( invStep *    y * 1.15470053839       ) + MAP_OFFSET;
        double b = ( invStep * (  x - 0.57735026919*y   ) ) + MAP_OFFSET;
    	ia = (int)a;
        ib = (int)b;
        da = a - ia;
        db = b - ib;
        return ( ( da + db ) > 1.0d );
    }


    inline void nodePoint( int ia, int ib, double& x, double& y ) const {
        int ia_ = ((int)ia) - MAP_OFFSET;
        y = step * ia_ * 0.86602540378;
        x = step * ( ( ((int)ib) - MAP_OFFSET ) + 0.5*ia_ );
    };

    inline void tilePoint( int ia, int ib, bool s, double& x, double& y ) const {
        nodePoint( ia, ia, x, y );
        if( s ){ x+=step; y+=0.57735026919*step; }else{ x+=0.5d*step; y+=0.28867513459*step;  }
    };

    //inline ULONG  getBucketInt    ( UHALF  ix, UHALF iy )const{ return ( iy << 16 ) + ix;                     };
	//inline void   unfoldBucketInt ( ULONG bucket, UHALF& ix, UHALF& iy  )const{ ix = bucket&0xFFFF; iy = (bucket>>16)&0xFFFF; }

	inline int   getBucketInt      (               int  ia, int ib   )const{ return ( ib << power ) + ia;                     };
	inline void  unfoldBucketInt ( int bucket, int& ia, int& ib  )const{ ia = bucket&mask; ib = (bucket>>power)&mask; };
    inline int   getBucket         ( double  x, double y )const{ int ia,ib; simplexIndexBare( x, y, ia, ib ); /* printf( " getBucket %3.3f %3.3f %i %i\n", x,y,ia,ib ); */ return getBucketInt( ia, ib );    };

    inline void initRuler( double step_, int power_ ){
        step    = step_;
		invStep = 1.0d/step;
		power   = power;
		nx      = 1<<power;
		ny      = nx;
		mask    = nx - 1;
    }

	inline int xy2i( int ix, int iy  ){ return iy*nx + ix; };
	//inline i2xy( int i, int ix, int iy ){     };

	void allocate( int nx_, int ny_ ){
		nx = nx_; ny = ny_; ntot = nx * ny;
		ground   = new double[ntot];
		water    = new double[ntot];
		water_   = new double[ntot]; // FIXME: this is just temporary array
        known    = new bool  [ntot]; // FIXME: this is just temporary array
	    contour1 = new int   [ntot]; // FIXME: just temporary, 1D => size < ntot
	    contour2 = new int   [ntot]; // FIXME: just temporary, 1D => size < ntot

	}

};



// by numlock
//   7-8-9
//  4-5-6
//   2-3


//inline int ipbc( int i, int n ){  if(i>n) return i-n; if(i<) }


#endif
