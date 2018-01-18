#ifndef Grid2DAlgs_h
#define Grid2DAlgs_h

#include "Vec2.h"
#include "GridIndex2D.h"

static const int    SquareNeighs       []={-1,0, 1,0, 0,-1, 0,1, -1,-1,    -1,+1,     +1,-1,     +1,+1     };
static const double SquareNeighsDist   []={ 1.0, 1.0,  1.0, 1.0, M_SQRT2,   M_SQRT2,   M_SQRT2,   M_SQRT2  };
static const double SquareNeighsDistInv[]={ 1.0, 1.0,  1.0, 1.0, M_SQRT1_2, M_SQRT1_2, M_SQRT1_2, M_SQRT1_2};

void bisecNoise( Vec2i n, double * hs, double frndMin, double frndMax );

class Grid2DAlg : public GridIndex2D{ public:
    static const int nNeighMax = 8;
    int nneigh    = 0;
    Vec2i  neighs        [nNeighMax];
    double neigh_dist    [nNeighMax];
    double neigh_invDist [nNeighMax];

    void initNeighsSquareN( int n ){ // typically 4 or 8
        nneigh = n;
        memcpy(neighs,       SquareNeighs       ,2*nneigh*sizeof(int)   );
        memcpy(neigh_dist,   SquareNeighsDist   ,  nneigh*sizeof(double));
        memcpy(neigh_invDist,SquareNeighsDistInv,  nneigh*sizeof(double));
    }
    void initNeighs_6(bool switched){
        nneigh = 0;
        int  sela[]={0,1,2,3, 4,7};
        int  selb[]={0,1,2,3, 5,6};
        int* sel;
        if(switched){sel=sela;}else{sel=selb;}
        for(int i=0; i<nneigh; i++){
            int j=sel[i];
            neighs       [i]=((Vec2i*)SquareNeighs)[j];
            neigh_dist   [i]=1.0;
            neigh_invDist[i]=1.0;
        }
    }
    void initNeighsSquareMask( uint8_t mask ){
        uint8_t bit =1;
        nneigh=0;
        for(int i=0; i<8; i++){
            if(mask&bit){
                neighs       [nneigh]=((Vec2i*)SquareNeighs)[i];
                neigh_dist   [nneigh]=SquareNeighsDist[i];
                neigh_invDist[nneigh]=SquareNeighsDistInv[i];
                nneigh++;
            }
            bit<<=1;
        }
    }

};

class HydraulicGrid2D :public Grid2DAlg { public:
    double * ground  = NULL;
	double * water   = NULL;

    Vec2i  droplet;
    double droplet_h;
    double droplet_w,droplet_disolve,droplet_sediment;

	/*
	// TODO: Inflow/outflow  +  Rivers
	bool   isOutflow = true;
	int    nContour=0,nContour_old=0;
	bool   * known     = NULL;
	int    * contour1  = NULL;
	int    * contour2  = NULL;
	double * water_    = NULL;
	std::vector<int>    sinks;
	std::vector<River*> rivers;
	*/



	// ==== function declaration
    /*
    // TODO: Inflow/outflow  +  Rivers
	void gatherRain( double minSinkFlow );
	int  traceDroplet( int ix, int iy, int nmax, int * trace );
	int  trackRiver( int sink, double minFlow, std::vector<int>& river, std::vector<int>& feeders );
	int  trackRiverRecursive( int sink, double minFlow, River * mouth );
	int  findAllRivers( double minFlow );
    void genTerrainNoise( int n, double scale, double hscale, double fdown, double strength, int seed, const Vec2d& pos0 );
    void init_outflow( double water_level );
    void outflow_step( );
    void extend_outflow( float val, int oi, int i );
    */

    void initDroplet  ( double w, double disolve, double sediment, Vec2i ipmin, Vec2i ipmax );
    bool droplet_step ( );
    void errodeDroples( int n, int nStepMax, double w, double disolve, double sediment, Vec2i ipmin, Vec2i ipmax );

    /*
    // TODO: Inflow/outflow  +  Rivers
	void allocate_outflow(){
        known     = new bool[ntot];
        contour1  = new int [ntot];
        contour2  = new int [ntot];
    }
    void deallocate_outflow(){
        delete [] known; delete [] contour1; delete [] contour2;
    }
    */

	void allocate( int nx_, int ny_ ){
		n.x = nx_; n.y = ny_; ntot = n.x * n.y;
		if(ground  ==NULL) ground   = new double[ntot];
		if(water   ==NULL) water    = new double[ntot];
	}

    void deallocate( bool all ){
        delete [] ground;
        delete [] water;
    };
};

class PathFindingGrid2D :public Grid2DAlg{ public:

};

#endif
