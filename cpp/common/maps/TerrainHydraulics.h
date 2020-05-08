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
    // continuous patch of all pixels where
    // 1) ground is below g0
    // 2) water level is above w0=g0+dw
    double level;
    double base;
};

struct TerrainCell{
    double g;
    double w;
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

    inline void sumWater( double& wsum, double& E ){
        for(int i=0; i<ntot;i++){
            double g  = ground[i];
            double w  = water [i];
            double dw = w-g;
            wsum += dw;
            E    += (w+g)*dw*0.5;
        }
    }

    // TODO - later this relaxator can search best way (curved), e.g. along river

    inline void relaxWaterRasterX( int iy, int ix0, int ix1, double hmin ){
        int i0 = iy*n.x;
        int i1 = i0+ix1; i0+=ix0;
        int istart=-1;
        //int iend;
        double wsum=0;
        double wmin,gmax;
        //bool started=false;
        for(int i=i0; i<i1; i++){
            double g = ground[i];
            double w = water [i];
            double h = w-g;
            if(istart<0){
                if( h>hmin ){ // start new block of water ?
                    istart=i;
                    wsum  =w;
                    wmin  =w;
                    gmax  =g;
                }
            }else{
                wmin=fmin(w,wmin);
                gmax=fmax(g,gmax);
                if( wmin-gmax < hmin ){ // end block ?
                    double wlevel = wsum/(i-istart); // calculate average flat level
                    double DEBUG_wsum = 0.0;
                    for(int j=istart; j<i; j++){     // fill block of water
                        DEBUG_wsum += water[j];
                        water[j] = wlevel;
                    }
                    //printf( "wsum[%i|%i] %g %g |wl %g \n", istart, i-1, wsum, DEBUG_wsum, wlevel );
                    istart=-1;
                }else{
                    wsum += w;
                }
            }
        }
    }

    inline void relaxWaterRasterY( int ix0, int iy0, int iy1, double hmin ){
        int i0 = iy0*n.x+ix0;
        int i1 = iy1*n.x+ix0;
        int istart=-1;
        //int iend;
        double wsum=0;
        double wmin,gmax;
        //bool started=false;
        for(int i=i0; i<i1; i+=n.x){
            double g = ground[i];
            double w = water [i];
            double h = w-g;
            if(istart<0){
                if( h>hmin ){ // start new block of water ?
                    istart=i;
                    wsum  =w;
                    wmin  =w;
                    gmax  =g;
                }
            }else{
                wmin=fmin(w,wmin);
                gmax=fmax(g,gmax);
                if( wmin-gmax < hmin ){ // end block ?
                    double wlevel = (wsum*n.x)/(i-istart); // calculate average flat level
                    double DEBUG_wsum = 0.0;
                    for(int j=istart; j<i; j+=n.x){     // fill block of water
                        DEBUG_wsum += water[j];
                        water[j] = wlevel;
                    }
                    //printf( "wsum[%i|%i] %g %g |wl %g \n", istart, i-1, wsum, DEBUG_wsum, wlevel );
                    istart=-1;
                }else{
                    wsum += w;
                }
            }
        }
    }

    inline void relaxWater2cells( TerrainCell& Ci, TerrainCell& Cj ){
        double  w = Ci.w + Cj.w;
        double dg = Ci.g - Cj.g;
        if(dg>0){
            if(dg>w){
                Cj.w=Cj.g+w;
                Ci.w=Ci.g;
            }else{
                w -= dg;
                w *= 0.5;
                Cj.w=Cj.g+w;
                Ci.w=Ci.g+w;
            }
        }else{
            if(-dg>w){
                Cj.w=Cj.g+w;
                Ci.w=Ci.g;
            }else{
                w += dg;
                w *= 0.5;
                Cj.w=Cj.g+w;
                Ci.w=Ci.g+w;
            }
        }
    }

    inline void relaxWaterHexCells( const Vec2i& ip0 ){
        TerrainCell cells[nneigh];
        int         is   [nneigh];
        int i0 = ip2i(ip0);
        TerrainCell cell0 = (TerrainCell){ground[i0],water[i0]};
        for(int ing=0; ing<nneigh; ing++){
            Vec2i ip   = wrap_index( ip0 + neighs[ing] );
            int i      = ip2i(ip);
            is[ing]    = i;
            cells[ing] = (TerrainCell){ground[i],water[i]};
        }
        relaxWater2cells( cells[0], cell0           );
        relaxWater2cells( cells[0], cells[nneigh-1] );
        for(int ing=1; ing<nneigh; ing++){
            relaxWater2cells( cells[ing], cell0        );
            relaxWater2cells( cells[ing], cells[ing-1] );
        }
        for(int ing=0; ing<nneigh; ing++){
            int i    = is[ing];
            ground[i] = cells[ing].g;
            water[i] = cells[ing].w;
        }
    }

    inline void relaxWater( const Vec2i& ip0 ){
        double hs[nneigh+1]; // heights
        int    is[nneigh+1]; // indexes
        int i0 = ip2i(ip0);
        //-- load neighbors to work array
        double g0 = ground[i0];
        is[nneigh]=i0;
        hs[nneigh]=g0;

        double wsum = water[i0]-g0;
        //double wsum_DEBUG = wsum;
        for(int ing=0; ing<nneigh; ing++){
            Vec2i ip = wrap_index( ip0 + neighs[ing] );
            int j    = ip2i(ip);
            is[ing]=j;
            double wj = water [j];
            double gj = ground[j];
            //wsum_DEBUG += wj-gj;
            if      (wj<g0){
                //wsum += 0;
                hs[ing]=wj;
            }else if(gj<g0){
                wsum += wj-g0;
                hs[ing]=g0;
            }else{
                wsum += wj-gj;
                hs[ing]=gj;
            }
            //printf( "%i: %i(%i,%i) g %g w %g -> %g \n", ing, j, ip.x,ip.y,  ground[j], water[j], hs[ing] );
        }

        //water[i0] = ground[i0] + wsum;

        //printf( "%i: %i(%i,%i) g %g w %g -> %g \n", nneigh, i, ip0.x,ip0.y,  ground[i], water[i], hs[nneigh] );
        //printf(" --- wsum %g wsum_DEBUG %g \n", wsum, wsum_DEBUG );
        //if(wsum<0) return;

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
        //for(int ing=0; ing<=nneigh; ing++){  printf( "%i : %i %g \n", ing, is[ing], hs[ing] ); }

        //double DEBUG_wsum_bak = wsum;
        // find common water level
        double g = hs[0];
        //printf( "--- g %g wsum %g \n", g, wsum );
        {
            //double dg,dw;
            int ing;
            for(ing=1; ing<=nneigh; ing++){
                double g_ = hs[ing];
                double dg = g_-g;
                double dw = dg*ing;
                // 0 = w -
                if(wsum>dw){
                    g     = g_;
                    wsum -= dw;
                }else{
                    break;
                }
                //printf(  "[%i] dg %g dw %g | g %g wsum %g \n", ing, dg, dw, g, wsum );
            }
            //double dg_ = dg*(wsum/dw);
            //double dw_ = dg*jng;
            g += wsum/ing;
            //printf(  "DEBUG finish [%i] wsum %g dg %g g %g \n", ing, wsum,  wsum/ing, g );
        }


        // DEBUG check
        //double DEBUG_wsum_pref = 0.0;
        //double DEBUG_wsum_new  = 0.0;
        //double DEBUG_wsum_h    = 0.0;
        //double DEBUG_wsum_gh   = 0.0;
        //for(int ing=0; ing<=nneigh; ing++){
        //    int i     = is[ing];
        //    double gi = ground[i];
        //    double h  = hs[ing];
        //   //DEBUG_wsum_pref += water[i] - gi;
        //    DEBUG_wsum_h    += water[i] - h;
        //    //DEBUG_wsum_gh   += g - h;
        //    if( g>h ){ DEBUG_wsum_gh   += g - h; }
        //    //if( gi<g ){ DEBUG_wsum_new += g - gi; }
        //}
        //printf( "DEBUG[%i,%i] wsum  %g -> %g \n", ip0.x,ip0.y,  DEBUG_wsum_pref, DEBUG_wsum_new  );
        //printf( "DEBUG wsum  %g  %g %g \n", DEBUG_wsum_bak,  DEBUG_wsum_h, DEBUG_wsum_gh  );
        //printf( "DEBUG wsum  %g  %g | %g \n", DEBUG_wsum_bak,  DEBUG_wsum_gh+wsum, wsum  );


        //printf(" --- watter level %g \n", g );
        // set water level
        for(int ing=0; ing<=nneigh; ing++){
            int i     = is[ing];
            double h  = hs[ing];
            //double gi = ground[i];
            //if( gi<g ){ gi = g; }
            if( g>h ){ h=g; }
            water[i] = h;
        }

        //double DEBUG_whsum = 0.0;
        //for(int ing=0; ing<=nneigh; ing++){
        //    DEBUG_whsum += water[is[ing]] - hs[ing];
        //}
        //printf( " DEBUG check wsum %g %g \n",  DEBUG_whsum, DEBUG_wsum_bak );
        //if( fabs(DEBUG_wsum_bak-DEBUG_whsum) >1e-6 ){
        //    printf( "DEBUG pix[%i,%i] check wsum %g %g \n", ip0.x,ip0.y,  DEBUG_whsum, DEBUG_wsum_bak );
        //    exit(0);
        //}

    }

    inline void relaxWater(){
        for(int iy=0; iy<n.y; iy++){
            for(int ix=0; ix<n.x; ix++){
                relaxWater( {ix, iy} );
            }
        }
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
