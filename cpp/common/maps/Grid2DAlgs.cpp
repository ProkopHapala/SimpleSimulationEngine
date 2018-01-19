#include "Grid2DAlgs.h"

void bisecNoise( int npow, double * hs, double frndMin, double frndMax ){
    int n    = 1<<npow;
    int mask = n-1;
    #define ip2i(ix,iy) ((iy&mask)<<npow)+(ix&mask)
    //for(int iy=0;iy<n;iy++){ for(int ix=0;ix<n;ix++){ hs[ip2i(ix,iy)]=0.5; } };
    for(int ipow=0;ipow<npow;ipow++){
        int p  = 1<<(ipow       );
        int q2 = 1<<(npow-ipow  );
        int q  = 1<<(npow-ipow-1);
        //printf( "==== %i %i %i  %f %f\n", ipow, p, q, frndMin*q2,frndMax*q2 );
        for(int iy=0;iy<n;iy+=q2){ for(int ix=0;ix<n;ix+=q2){
            //printf( "%i %i\n", ix,iy );
            double h00 = hs[ip2i(ix   ,iy   )];
            double h01 = hs[ip2i(ix+q2,iy   )];
            double h10 = hs[ip2i(ix   ,iy+q2)];
            double h11 = hs[ip2i(ix+q2,iy+q2)];
            hs[ip2i(ix+q,iy+q)]=0.25*(h00+h01+h10+h11) + randf(frndMin,frndMax)*q2;
        } }
        for(int iy=0;iy<n;iy+=q2){ for(int ix=q;ix<n;ix+=q2){
            double h00 = hs[ip2i(ix-q,iy   )];
            double h01 = hs[ip2i(ix+q,iy   )];
            double h10 = hs[ip2i(ix  ,iy-q)];
            double h11 = hs[ip2i(ix  ,iy+q)];
            hs[ip2i(ix,iy)] = 0.25*(h00+h01+h10+h11) + randf(frndMin,frndMax)*q2*M_SQRT1_2;
        } }
        for(int iy=q;iy<n;iy+=q2){ for(int ix=0;ix<n;ix+=q2){
            double h00 = hs[ip2i(ix-q,iy   )];
            double h01 = hs[ip2i(ix+q,iy   )];
            double h10 = hs[ip2i(ix  ,iy-q)];
            double h11 = hs[ip2i(ix  ,iy+q)];
            hs[ip2i(ix,iy)] = 0.25*(h00+h01+h10+h11) + randf(frndMin,frndMax)*q2*M_SQRT1_2;
        } }
        //if (ipow>=2) break;
    }
    #undef ip2i
};

void bisecPaternNoise( int npow, double * hs, double frndMin, double frndMax ){
    int n    = 1<<npow;
    int mask = n-1;
    #define ip2i(ix,iy) ((iy&mask)<<npow)+(ix&mask)
    //for(int iy=0;iy<n;iy++){ for(int ix=0;ix<n;ix++){ hs[ip2i(ix,iy)]=0.5; } };
    for(int ipow=0;ipow<npow;ipow++){
        int p  = 1<<(ipow       );
        int q2 = 1<<(npow-ipow  );
        int q  = 1<<(npow-ipow-1);
        for(int iy=0;iy<n;iy+=q2){ for(int ix=0;ix<n;ix+=q2){
            double h00 = hs[ip2i(ix   ,iy   )];
            double h01 = hs[ip2i(ix+q2,iy   )];
            double h10 = hs[ip2i(ix   ,iy+q2)];
            double h11 = hs[ip2i(ix+q2,iy+q2)];
            hs[ip2i(ix+q,iy  )]=0.5*(h00+h01) + randf(frndMin,frndMax)*q2;
            hs[ip2i(ix  ,iy+q)]=0.5*(h00+h10) + randf(frndMin,frndMax)*q2;
        } }
        for(int iy=q;iy<n;iy+=q2){ for(int ix=q;ix<n;ix+=q2){
            double h00 = hs[ip2i(ix-q,iy   )];
            double h01 = hs[ip2i(ix+q,iy   )];
            double h10 = hs[ip2i(ix  ,iy-q2)];
            double h11 = hs[ip2i(ix  ,iy+q2)];
            hs[ip2i(ix,iy)] = 0.25*(h00+h01+h10+h11) + randf(frndMin,frndMax)*q;
        } }
    }
    #undef ip2i
};

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
    if( (droplet.x<0)&(droplet.x>(n.x-1))&(droplet.y<0)&(droplet.y>(n.y-1)) ) return true;
    Vec2i ipmin;
    double h,hmin=droplet_h+droplet_w;
    bool found=false;
    for(int ing=0; ing<nneigh; ing++){
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

void HydraulicGrid2D::gatherRain( double minSinkFlow ){
    double SAFETY = 1e-8;
    sinks.clear();
    for(int i=0; i<ntot; i++){ contour1[i]=i; water[i]=1.0; known[i]=false; }
    quickSort( ground, contour1, 0, ntot);
    for(int ii=ntot; ii>0; ii--){
        int i = contour1[ii];
        Vec2i ip0 = i2ip(i);
        //int iy = i/nx;
        //int ix = i%nx;
        int imin = -1;
        double val = ground[i]-SAFETY;
        for(int ing=0; ing<nneigh; ing++){
            Vec2i ip = wrap_index( ip0 + neighs[ing] );
            int j    = ip2i(ip);
            double g=ground[j]; if(g<val){ imin=j; val=g; }
        }
        /*
        #define MINSTEP double g=ground[i_]; if(g<val){ imin=i_; val=g; }
        if( ix>0                ){ int i_=i-1;    MINSTEP }
        if( ix<(nx-1)           ){ int i_=i+1;    MINSTEP }
        if( iy>0                ){ int i_=i-nx;   MINSTEP }
        if( iy<(ny-1)           ){ int i_=i+nx;   MINSTEP }
        if( (iy>0)&&(ix<(nx-1)) ){ int i_=i-nx+1; MINSTEP }
        if( (iy<(ny-1))&&(ix>0) ){ int i_=i+nx-1; MINSTEP }
        #undef MINSTEP
        */
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
    printf( "max water %f \n", wmax);
    //for(int i=0; i<ntot; i++){ water[contour1[i]]=i*0.001; }
}

/*
class River{ public:
    River* mouth = NULL;
    std::vector<int>    path;
    std::vector<double> flow;
};
*/

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
        /*
        #define MINSTEP double g=ground[i_]; if(g<val){ imin=i_; val=g; }
        if( ix>0                ){ int i_=i-1;    MINSTEP }
        if( ix<(nx-1)           ){ int i_=i+1;    MINSTEP }
        if( iy>0                ){ int i_=i-nx;   MINSTEP }
        if( iy<(ny-1)           ){ int i_=i+nx;   MINSTEP }
        if( (iy>0)&&(ix<(nx-1)) ){ int i_=i-nx+1; MINSTEP }
        if( (iy<(ny-1))&&(ix>0) ){ int i_=i+nx-1; MINSTEP }
        //printf( "trace[%i] (%i,%i) : %i %g\n", ii, ix, iy,   imin, val );
        #undef MINSTEP
        */
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
        /*
        int iy = i/nx;
        int ix = i%nx;
        int imax  = -1;
        //int imax2 = -1;
        int nhigh=0;
        double vallim = water[i]-SAFETY;
        double valmax = 0;
        #define MAXSTEP double w=water[i_]; if( (!known[i_])&&(w>minFlow) ){ nhigh++; if( (w<vallim)&&(w>valmax) ){ imax=i_; valmax=w; } }
        if( ix>0                ){ int i_=i-1;    MAXSTEP }
        if( ix<(nx-1)           ){ int i_=i+1;    MAXSTEP }
        if( iy>0                ){ int i_=i-nx;   MAXSTEP }
        if( iy<(ny-1)           ){ int i_=i+nx;   MAXSTEP }
        if( (iy>0)&&(ix<(nx-1)) ){ int i_=i-nx+1; MAXSTEP }
        if( (iy<(ny-1))&&(ix>0) ){ int i_=i+nx-1; MAXSTEP }
        #undef MAXSTEP
        */
        //printf( " (%i,%i) %i %f \n", ix,iy,  imax, valmax );
        if(imax>=0){
            //if(imax2>=0) feeders.push_back(i);
            if(nhigh>1){
                bool cannot = (feeders.size()>0) && ( feeders[feeders.size()-1] == river[river.size()-1] ); // check if included in previous step
                if( !cannot ){  feeders.push_back(i); }
            }
            river.push_back(i);
            known[i]=true;
            i = imax;
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
        for( ifeeder : feeders ){
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
    if( nContour_old>0 ) printf( "nContour_old %i nContour %i \n", nContour_old, nContour );
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
        /*
        if( ix>0      ){  extend_outflow( val, i, i-1   ); }
        if( ix<(nx-1) ){  extend_outflow( val, i, i+1   ); }
        if( iy>0      ){  extend_outflow( val, i, i-nx  ); }
        if( iy<(ny-1) ){  extend_outflow( val, i, i+nx  ); }
        if( (iy>0)&&(ix<(nx-1)) ){  extend_outflow( val, i, i-nx+1); }
        if( (iy<(ny-1))&&(ix>0) ){  extend_outflow( val, i, i+nx-1); }
        */
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
