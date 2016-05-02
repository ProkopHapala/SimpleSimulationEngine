#ifndef TerrainHydraulics_h
#define TerrainHydraulics_h

#include "Noise.h"



class TerrainHydraulics{
public:

    static constexpr int     nmask   = 8;
    static constexpr double  isqrt2  = 0.70710678118;
    static constexpr int     mask [nmask*2] = { 0,-1, 0,1, -1,0, +1,0, -1,-1, -1,+1, +1,-1, +1+1 };
    static constexpr double  dists[nmask  ] = { 1.0, 1.0, 1.0, 1.0, isqrt2, isqrt2, isqrt2, isqrt2 };

    int istep;
	int    nx, ny, ntot;
	double * ground;
	double * water;


	// countor step algorithm
	//int    nContourMax;
	int    nContour=0,nContour_old=0;
	bool   * known;
	int    * contour1;
	int    * contour2;
	double * water_;

/*
    double c_errode    = 0.4;
    double c_rain      = 0.0001;
    double c_blur      = 0.0;
    double c_sediment  = 0.001;
    double f_sediment  = 0.1;
    double f_mix       = 0.1;
    double mf_mix      = 1.0-f_mix;
*/

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

#endif
