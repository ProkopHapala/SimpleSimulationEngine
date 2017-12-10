#ifndef TerrainHydraulics_h
#define TerrainHydraulics_h

#include <vector>
#include "Noise.h"
#include "arrayAlgs.h"

class River{ public:
    River* mouth = NULL;
    std::vector<int>    path;
    std::vector<double> flow;
};

class TerrainHydraulics{
public:

    static constexpr int     nmask   = 8;
    static constexpr double  isqrt2  = 0.70710678118;
    static constexpr int     mask [nmask*2] = { 0,-1, 0,1, -1,0, +1,0, -1,-1, -1,+1, +1,-1, +1+1 };
    static constexpr double  dists[nmask  ] = { 1.0, 1.0, 1.0, 1.0, isqrt2, isqrt2, isqrt2, isqrt2 };

    int      istep=0;
	int      nx=0, ny=0, ntot=0;

	double * ground  = NULL;
	double * water   = NULL;

	// countor step algorithm
	//int    nContourMax;
	bool   isOutflow = true;
	int    nContour=0,nContour_old=0;
	bool   * known     = NULL;
	int    * contour1  = NULL;
	int    * contour2  = NULL;
	double * water_    = NULL;

	std::vector<int>    sinks;
	std::vector<River*> rivers;

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
    double droplet_h;
    double droplet_w,droplet_disolve,droplet_sediment;

	// ==== function declaration

	void gatherRain( double minSinkFlow );
	int  traceDroplet( int ix, int iy, int nmax, int * trace );
	int  trackRiver( int sink, double minFlow, std::vector<int>& river, std::vector<int>& feeders );
	int  trackRiverRecursive( int sink, double minFlow, River * mouth );
	int  findAllRivers( double minFlow );

    void genTerrainNoise( int n, double scale, double hscale, double fdown, double strength, int seed, const Vec2d& pos0 );
    void init_outflow( double water_level );
    void outflow_step( );
    void extend_outflow( float val, int oi, int i );
    //void inflow_step( );
    //void extend_inflow( float val, int oi, int i );

    void initErrosion( double w );
    void flow_errosion_step  ( );
    void rain_and_evaporation( );
    int  flow_errosion_step_noRain( );
    void initDroplet  ( double w, double disolve, double sediment, int ix0, int iy0, int ix1, int iy1 );
    bool droplet_step ( );
    void errodeDroples( int n, int nStepMax, double w, double disolve, double sediment, int ix0, int iy0, int ix1, int iy1 );

	// ==== inline functions

	inline void setSize( int nx_, int ny_ ){ nx=nx_; ny=ny_; ntot=nx*ny; };
	inline int xy2i( int ix, int iy  ){ return iy*nx + ix; };
	//inline i2xy( int i, int ix, int iy ){     };

	void allocate_outflow(){
        known     = new bool[ntot];
        contour1  = new int [ntot];
        contour2  = new int [ntot];
    }

    void deallocate_outflow(){
        delete [] known; delete [] contour1; delete [] contour2;
    }

	void allocate( int nx_, int ny_ ){
		nx = nx_; ny = ny_; ntot = nx * ny;
		if(ground  ==NULL) ground   = new double[ntot];
		if(water   ==NULL) water    = new double[ntot];
		if(water_  ==NULL) water_   = new double[ntot]; // FIXME: this is just temporary array
        if(known   ==NULL) known    = new bool  [ntot]; // FIXME: this is just temporary array
	    if(contour1==NULL) contour1 = new int   [ntot]; // FIXME: just temporary, 1D => size < ntot
	    if(contour2==NULL) contour2 = new int   [ntot]; // FIXME: just temporary, 1D => size < ntot
	}

    void deallocate( bool all ){
		if(all){
            delete [] ground;
            delete [] water;
		}
		delete [] water_;
        delete [] known;
	    delete [] contour1;
	    delete [] contour2;
    };

};

#endif
