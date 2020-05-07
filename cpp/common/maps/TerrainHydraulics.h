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

class Lake{ public:
    double level;
    double base;
};

class HydraulicGrid2D :public Grid2DAlg { public:
    double * ground  = NULL;
	double * water   = NULL;

	int  nlakes;
	int  * lakeMap = NULL;
	Lake * lakes   = NULL;

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


	// PathFinder
    //std::vector<int> rivers;
    double* moveCost = NULL;
    int* toBasin     = NULL;
    int* toTile      = NULL;




    inline void relaxWater( const Vec2i& ip0 ){
        double hs[nneigh+1]; // heights
        int    is[nneigh+1]; // indexes

        int i = ip2i(ip0);
        //-- load neighbors to work array
        double g = ground[i];
        is[nneigh]=i;
        hs[nneigh]=g;

        double wsum = water[i]-g;
        double wsum_DEBUG = wsum;

        for(int ing=0; ing<nneigh; ing++){
            Vec2i ip = wrap_index( ip0 + neighs[ing] );
            int j    = ip2i(ip);
            is[ing]=j;
            double wj = water [j];
            double gj = ground[j];
            wsum_DEBUG += wj-gj;
            if      (wj<g){
                //wsum += 0;
                hs[ing]=wj;
            }else if(gj<g){
                wsum += wj-g;
                hs[ing]=g;
            }else{
                wsum += wj-gj;
                hs[ing]=gj;
            }
            printf( "%i: %i(%i,%i) g %g w %g -> %g \n", ing, j, ip.x,ip.y,  ground[j], water[j], hs[ing] );
        }
        //printf( "%i: %i(%i,%i) g %g w %g -> %g \n", nneigh, i, ip0.x,ip0.y,  ground[i], water[i], hs[nneigh] );
        printf(" --- wsum %g wsum_DEBUG %g \n", wsum, wsum_DEBUG );
        if(wsum<0) return;
        //printf(" --- loaded \n");
        //for(int ing=0; ing<=nneigh; ing++){  printf( "%i : %i %g \n", ing, is[ing], hs[ing] ); }
        //-- sort work array by height ( bouble sort )
        int j=nneigh+1; bool s=true;
        while (s) {
            s = false;
            for (int i=1;i<j;i++) {
                int i_=i-1;
                //printf(  "bouble(%i,%i) [%i,%i] %g<%g ", j, i,   i,i_,  hs[i], hs[i_]  );
                if (hs[i]<hs[i_]){
                    SWAP(hs[i],hs[i_],double);
                    SWAP(is[i],is[i_],int);
                    s = true;
                }
                //printf(  "swap? %i \n", s );
            }
            j--;
        }
        //printf(" --- sorted \n");
        for(int ing=0; ing<=nneigh; ing++){  printf( "%i : %i %g \n", ing, is[ing], hs[ing] ); }
        // find common water level
        g = hs[0];
        printf( "--- g %g wsum %g \n", g, wsum );
        for(i=1; i<=nneigh; i++){
            double g_ = hs[i];
            double dg = g_-g;
            double dw = dg*i;
            // 0 = w -
            if(wsum>dw){
                g     = g_;
                wsum -= dw;
            }else{
                g += dg*(wsum/dw);
                //printf(  "[%i] dg %g dw %g | g %g wsum %g \n", i, dg, dw, g, wsum );
                break;
            }
            //printf(  "[%i] dg %g dw %g | g %g wsum %g \n", i, dg, dw, g, wsum );
        }
        printf(" --- watter level %g \n", g );
        // set water level
        for(i=0; i<=nneigh; i++){
            water[is[i]] = g;
        }
    }



    inline void relaxWater(){
        for(int iy=0; iy<n.y; iy++){
            for(int ix=0; ix<n.x; ix++){
                relaxWater( {ix, iy} );
                //double g = ground[i];
                //double w = water [i];
                /*
                // sum watter
                // 1) sum watter above ground of this pixel
                // 2) sum void   below ground level of this pixel
                double wsum = w-g;
                double vsum = 0;
                double gsum = 0;
                for(int ing=0; ing<nneigh; ing++){
                    Vec2i ip = wrap_index( ip0 + neighs[ing] );
                    int j    = ip2i(ip);
                    //double gj=ground[j];
                    double wj=water[j];
                    if(wj>g){
                        double gj = ground[j];
                        if(gj>g){ gsum += gj-g; wsum+=wj-gj; }else{  wsum+=wj-g; }
                    }else{
                        vsum += g-wj;
                    }
                }
                if()
                double w = (wsum+gsum-vsum)/(nneigh+1);
                for(int ing=0; ing<nneigh; ing++){
                    Vec2i ip = wrap_index( ip0 + neighs[ing] );
                    int j    = ip2i(ip);
                    double g=ground[j]; if(g<val){ imin=j; val=g; }
                }
                */
            }
        }



        /*
        double SAFETY = 1e-8;
        for(int ii=ntot; ii>0; ii--){
            int i = contour1[ii];
            Vec2i ip0 = i2ip(i);
            int imin = -1;
            double val = ground[i]-SAFETY;
            for(int ing=0; ing<nneigh; ing++){
                Vec2i ip = wrap_index( ip0 + neighs[ing] );
                int j    = ip2i(ip);
                double g=ground[j]; if(g<val){ imin=j; val=g; }
            }
            //printf("%i %i %i %f %i\n", ii, i, imin, val, nneigh );
            if(imin>=0){
                water[imin] += water[i];
            }else{
                // sink-less pixel
                known[i] = true;
                if( water[i] > minSinkFlow ) sinks.push_back(i);
            }
        }
        double wmax = 0.0;
        for(int i=0; i<ntot; i++){ wmax = fmax(wmax,water[i]); }
        //printf( "max water %f \n", wmax);
        return wmax;
        //for(int i=0; i<ntot; i++){ water[contour1[i]]=i*0.001; }
        */
    }





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
