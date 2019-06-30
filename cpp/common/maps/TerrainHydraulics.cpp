
#include "TerrainHydraulics.h"  // THE HEADER

// ====================    HydraulicGrid2D
// ====================    DropletErrosion

void HydraulicGrid2D::initDroplet ( double w, double disolve, double sediment, Vec2i ipmin, Vec2i ipmax ){
    droplet_w        = w;
    droplet_disolve  = disolve;
    droplet_sediment = sediment;
    droplet.x = ipmin.x + rand()%(ipmax.x-ipmin.x);
    droplet.y = ipmin.y + rand()%(ipmax.y-ipmin.y);
    droplet_h  = ground[ip2i( droplet )];
    //printf( "initDroplet  %i %i    %f  \n", droplet_ix, droplet_iy, droplet_h  );
}

bool HydraulicGrid2D::droplet_step( ){
    if( !innerIndex(droplet) ) return true;
    Vec2i ipmin;
    double h,hmin=droplet_h+droplet_w;
    bool found=false;
    for(int ing=0; ing<nneigh; ing++){
        // Vec2i ip = wrap_index( droplet + neighs[ing] ); // This would be slightly slower
        Vec2i ip = droplet + neighs[ing];
        //printf( "ineigh %i (%i,%i) \n", ing, ip.x, ip.y );
        int i    = ip2i(ip);
        double h = ground[i];
        if( h<hmin ){ ipmin = ip; hmin=h; found=true; }
    }
    //exit(0);
    if( found ){
        int    i  = ip2i( droplet );
        int    j  = ip2i( ipmin   );
        double hi = droplet_h;
        double hj = ground[j];
        double dh = (hi-hj);
        //printf( " %f %f %f \n", hi, hj, dh  );
        if( dh > 0 ){
            dh *= droplet_disolve;
            ground[i] = hi - dh;
            droplet_h = hj + dh*droplet_sediment;
            ground[j] = droplet_h;
            //ground[jj] = hj + dh;
        }
        droplet = ipmin;
        //droplet_h   = hmin;
        //ground[iimin] += mud*f_sediment;
        return false;
    }
    return true;
}

void HydraulicGrid2D::errodeDroples( int n, int nStepMax, double w, double disolve, double sediment, Vec2i ipmin, Vec2i ipmax ){
    for(int i=0; i<n; i++ ){
        initDroplet( w, disolve, sediment, ipmin, ipmax );
        for(int j=0; j<nStepMax; j++ ){
            if( droplet_step( ) ) break;
        }
    }
};

// ====================    HydraulicGrid2D
// ====================    DropletErrosion

double HydraulicGrid2D::gatherRain( double minSinkFlow ){
    double SAFETY = 1e-8;
    sinks.clear();
    for(int i=0; i<ntot; i++){ contour1[i]=i; water[i]=1.0; known[i]=false; }
    quickSort( ground, contour1, 0, ntot);
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
}

int HydraulicGrid2D::traceDroplet( Vec2i ipd, int nmax, int * trace ){
    double SAFETY = 1e-8;
    //double SAFETY = -0.5;
    //int i=ip2i(ip);
    int ii=0;
    for(ii=0; ii<nmax; ii++){
        //ip = i2ip(i);
        Vec2i ipmin;
        int   imin=-1;
        bool found = false;
        double val = ground[ ip2i(ipd) ]-SAFETY;
        for(int ing=0; ing<nneigh; ing++){
            Vec2i ip = wrap_index( ipd + neighs[ing] );
            int i    = ip2i(ip);
            double g=ground[i]; if(g<val){ ipmin=ip; imin=i; val=g; }
        }
        if(imin){
            ipd       = ipmin;
            trace[ii] = imin;
        }else{
            break;
        }
    }
    return ii;
}

int HydraulicGrid2D::trackRiver( int sink, double minFlow, std::vector<int>& river, std::vector<int>& feeders ){
    river  .clear();
    feeders.clear();
    double SAFETY = 1e-8;
    //printf( " TerrainHydraulics::trackRiver( %i, %f ) \n", sink, minFlow );
    int ii;
    Vec2i ipd = i2ip(sink);
    for(ii=0; ii<10000; ii++){
        Vec2i ipmax;
        int imax  = -1;
        int i     = ip2i(ipd);
        double vallim = water[i]-SAFETY;
        double valmax = 0;
        int nhigh=0;
        for(int ing=0; ing<nneigh; ing++){
            Vec2i ip = wrap_index( ipd + neighs[ing] );
            int j    = ip2i(ip);
            double w=water[j]; if( (!known[j])&&(w>minFlow) ){ nhigh++; if( (w<vallim)&&(w>valmax) ){ imax=j; ipmax=ip; valmax=w; } }
        }
        //printf( "trackRiver %i : (%i,%i) %i %f \n", ii, ipmax.x,ipmax.y,  imax, valmax );
        if(imax>=0){
            //if(imax2>=0) feeders.push_back(i);
            if(nhigh>1){
                bool cannot = (feeders.size()>0) && ( feeders[feeders.size()-1] == river[river.size()-1] ); // check if included in previous step
                if( !cannot ){  feeders.push_back(i); }
            }
            river.push_back(i);
            known[i]=true;
            //i = imax;
            ipd = ipmax;
        }else{
            // sink-less pixel
            //known[i] = true;
            //if( water[i] > minSinkFlow ) sinks.push_back(i);
            break;
        }
    }
    //double wmax = 0.0;
    //for(int i=0; i<ntot; i++){ wmax = fmax(wmax,water[i]); }
    //printf( "max water %f \n", wmax);
    //for(int i=0; i<ntot; i++){ water[contour1[i]]=i*0.001; }
    return ii;
}

int HydraulicGrid2D::trackRiverRecursive( int sink, double minFlow, River * mouth ){
    //printf( "trackRiverRecursive  %i \n", sink );
    std::vector<int> feeders;
    River* river = new River();
    river->mouth=mouth;
    int n=trackRiver( sink, minFlow, river->path, feeders );
    //printf( "path size  %i %i \n", river->path.size(), n );
    int nriv=1;
    if( river->path.size() > 5 ){
        rivers.push_back(river);

        // save flow
        river->flow.reserve(river->path.size());
        for( int i=0; i<river->path.size(); i++ ){
            river->flow[i] = water[river->path[i]];
        }

        for(int ifeeder : feeders ){
            nriv+=trackRiverRecursive( ifeeder, minFlow, river );
        }
    }else{
        delete river;
    }
    return nriv++;
}

int HydraulicGrid2D::findAllRivers( double minFlow ){
    for(River * river : rivers) delete river;
    rivers.clear();
    for(int i=0; i<ntot; i++){ known[i]=false; }
    int n=0;
    for( int sink : sinks ){
        n+=trackRiverRecursive( sink, minFlow, NULL );
    }
    printf( "findAllRivers DONE \n" );
    return n;
}

void HydraulicGrid2D::init_outflow( double water_level ){
    for (int i=0; i<ntot; i++){
        //water[i] = 1e+300;
        water[i] = water_level;
        known[i] = false;
    }
}

void HydraulicGrid2D::outflow_step(){
    // clean known
    // swap old and new endpoints
    nContour_old = nContour;  nContour = 0;
    int * tmp = contour1; contour1 = contour2; contour2 = tmp;
    for ( int ii=0; ii<nContour_old; ii++ ){ known[contour1[ii]] = false; }
    // check all countor point neighbors for possible path extension
    //if( nContour_old>0 ) printf( "nContour_old %i nContour %i \n", nContour_old, nContour );
    for ( int ii=0; ii<nContour_old; ii++ ){
        int    i   = contour1[ii];
        Vec2i ip0  = i2ip(i);
        double val = water[i];
        //val +=  dval/( val - ground[i] + 0.1 ); // this is just to simulate non-zero viscosity of watter
        //printf( " %i %i (%i,%i)\n", ii, i, ix, iy );
        // extend in four directions, check boundary overflow
        for(int ing=0; ing<nneigh; ing++){
            Vec2i ip = wrap_index( ip0 + neighs[ing] );
            if( validIndex( ip ) ){ extend_outflow( val, i, ip2i(ip) ); };
        }
    }
}

void HydraulicGrid2D::extend_outflow( float val, int oi, int i ){
    // evaluate objective function for proposed path
    //val += dval*ground[i];
    //printf( " %i val %f g %f w %f k %i \n", i, val, ground[i], water[i], known[i] );
    val = fmax( val, ground[i] );
    // check if new path is advantegeneous
    bool doit;
    if( isOutflow ){ doit = val < water[i]; }  // ouflow
    else           { doit = val > water[i]; }; // inflow
    if( doit ){
    //if( (val < water[i]) == isOutflow ){
        water[i] = val;       //  extend path
        //printf( "mod %i  \n", i );
        //gofrom[i] = oi;     //  ... optionaly we may store oi to reconstruct path backwards later
        if( !known[i] ){
                //printf( "added %i to contour \n", i );
                known[i]=true; contour2[nContour]=i; nContour++;
        } // add to list only if not already known
    }
}


