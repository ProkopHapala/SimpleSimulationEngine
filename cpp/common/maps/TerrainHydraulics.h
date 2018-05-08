#ifndef TerrainHydraulics_h
#define TerrainHydraulics_h

#include <vector>
#include "Noise.h"
#include "arrayAlgs.h"
#include "Grid2DAlgs.h"

class River{ public:
    River* mouth = NULL;
    std::vector<int>    path;
    std::vector<double> flow;
};

class HydraulicGrid2D :public Grid2DAlg { public:
    double * ground  = NULL;
	double * water   = NULL;

    Vec2i  droplet;
    double droplet_h;
    double droplet_w,droplet_disolve,droplet_sediment;

	// TODO: Inflow/outflow  +  Rivers
	bool   isOutflow = true;
	int    nContour=0,nContour_old=0;
	bool   * known     = NULL;
	int    * contour1  = NULL;
	int    * contour2  = NULL;
	double * water_    = NULL;
	std::vector<int>    sinks;
	std::vector<River*> rivers;

	// ==== function declaration

    // TODO: Inflow/outflow  +  Rivers
	double gatherRain( double minSinkFlow );
	int  traceDroplet( Vec2i ipd, int nmax, int * trace );
	int  trackRiver( int sink, double minFlow, std::vector<int>& river, std::vector<int>& feeders );
	int  trackRiverRecursive( int sink, double minFlow, River * mouth );
	int  findAllRivers( double minFlow );
    void init_outflow( double water_level );
    void outflow_step( );
    void extend_outflow( float val, int oi, int i );

    // TODO: Inflow/outflow  +  Rivers
	void allocate_outflow(){
        //known     = new bool[ntot];
        //contour1  = new int [ntot];
        //contour2  = new int [ntot];
        _realloc(known,ntot);
		_realloc(contour1,ntot);
		_realloc(contour2,ntot);
    }
    void deallocate_outflow(){
        delete [] known; delete [] contour1; delete [] contour2;
    }

    void initDroplet  ( double w, double disolve, double sediment, Vec2i ipmin, Vec2i ipmax );
    bool droplet_step ( );
    void errodeDroples( int n, int nStepMax, double w, double disolve, double sediment, Vec2i ipmin, Vec2i ipmax );

	void allocate( Vec2i n ){
		setN(n);
		_realloc(ground,ntot);
		_realloc(water,ntot);
		//n.x = nx_; n.y = ny_; ntot = n.x * n.y;
		//if(ground  ==NULL) ground   = new double[ntot];
		//if(water   ==NULL) water    = new double[ntot];
	}

    void deallocate( bool all ){
        delete [] ground;
        delete [] water;
    };
};

#endif
