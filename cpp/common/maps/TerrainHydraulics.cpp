
#include "TerrainHydraulics.h"  // THE HEADER

void TerrainHydraulics::genTerrainNoise( int n, double scale,  double hscale,  double fdown, double strength, int seed, const Vec2d& pos0 ){
    Vec2d pos,dpos,rot,a,b;
    rot.fromAngle( seed*0.1 );
    a.set( 1.0d, 0.0d           ); a.mul(scale);
    b.set( 0.5d, 0.86602540378d ); b.mul(scale);
    double vmin=+1e+300,vmax=-1e+300;
    int ii = 0;
    for (int iy=0; iy<ny; iy++){
        for (int ix=0; ix<nx; ix++){
            pos.set( ix*a.x+iy*b.x +pos0.x, ix*a.y+iy*b.y+pos0.y );
            Noise::warpNoise3R( pos, rot, fdown, strength, n, dpos );
            //ground[ii] = 0.5*sin( 1000*dpos.x * dpos.y )+0.5;
            double val = fabs(dpos.x * dpos.y);
            ground[ii] = val;
            if(val<vmin){vmin=val;}
            if(val>vmax){vmax=val;}
            ii++;
        }
    }
    printf( " vmin %e vmax %e \n", vmin, vmax );
    double renorm = hscale/(vmax-vmin);
    for(int i=0; i<ntot; i++){ ground[i] = renorm*(ground[i]-vmin); }
}

void TerrainHydraulics::init_outflow( double water_level ){
    for (int i=0; i<ntot; i++){
        //water[i] = 1e+300;
        water[i] = water_level;
        known[i] = false;
    }
}

void TerrainHydraulics::outflow_step(){
    // clean known
    // swap old and new endpoints
    nContour_old = nContour;  nContour = 0;
    int * tmp = contour1; contour1 = contour2; contour2 = tmp;
    for ( int ii=0; ii<nContour_old; ii++ ){ known[contour1[ii]] = false; }
    // check all countor point neighbors for possible path extension
    if( nContour_old>0 ) printf( "nContour_old %i nContour %i \n", nContour_old, nContour );
    for ( int ii=0; ii<nContour_old; ii++ ){
        int    i   = contour1[ii];
        double val = water[i];
        //val +=  dval/( val - ground[i] + 0.1 ); // this is just to simulate non-zero viscosity of watter
        int iy = i/nx;
        int ix = i%nx;
        //printf( " %i %i (%i,%i)\n", ii, i, ix, iy );
        // extend in four directions, check boundary overflow
        if( ix>0      ){  extend_outflow( val, i, i-1   ); }
        if( ix<(nx-1) ){  extend_outflow( val, i, i+1   ); }
        if( iy>0      ){  extend_outflow( val, i, i-nx  ); }
        if( iy<(ny-1) ){  extend_outflow( val, i, i+nx  ); }
        if( (iy>0)&&(ix<(nx-1)) ){  extend_outflow( val, i, i-nx+1); }
        if( (iy<(ny-1))&&(ix>0) ){  extend_outflow( val, i, i+nx-1); }
    }
}

void TerrainHydraulics::extend_outflow( float val, int oi, int i ){
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

void TerrainHydraulics::initErrosion( double w ){
    for (int i=0; i<ntot; i++){
            if(ground[i]>1.0) ground[i] = 1.0;
            if(ground[i]<0.0) ground[i] = 0.0;
            water[i] = w; water_[i] = 0.0;
    }
}

int TerrainHydraulics::flow_errosion_step_noRain( ){
// gradient evaluation
    int ii = 0;
    int npix=0;
    for(int iy=1;iy<ny-1;iy++){
        for(int ix=1;ix<ny-1;ix++){
            //int    ii   = xy2i( ix, iy );
            ii++;
            double hij  = ground[ii];
            double dh,dhmin=0,dhmin2=0;
            int jj,iimin;
            jj = ii-nx  ; dh = ground[jj]-hij; if( dh<dhmin ){ iimin = jj; dhmin = dh; }
            jj = ii+nx  ; dh = ground[jj]-hij; if( dh<dhmin ){ iimin = jj; dhmin = dh; }
            jj = ii-1   ; dh = ground[jj]-hij; if( dh<dhmin ){ iimin = jj; dhmin = dh; }
            jj = ii+1   ; dh = ground[jj]-hij; if( dh<dhmin ){ iimin = jj; dhmin = dh; }
            jj = ii-nx+1; dh = ground[jj]-hij; if( dh<dhmin ){ iimin = jj; dhmin = dh; }
            jj = ii+nx-1; dh = ground[jj]-hij; if( dh<dhmin ){ iimin = jj; dhmin = dh; }
            if( dhmin<-0.00001 ){
                double wij       = water[ii];
                if(wij>0.0001) npix++;
                double sediment  = wij*c_sediment/(1-dhmin);
                wij             -= sediment;
                water_[iimin]   += wij;
                water[ii]        = 0;
                double mud =  c_errode * sqrt(wij)*(-dhmin);
                ground[ii   ]    = hij - mud + sediment;
                ground[iimin]   += mud*f_sediment;
            }
        }
    }
    ii =0;
    for(int iy=0;iy<ny;iy++){
        for(int ix=0;ix<nx;ix++){
            double wij = water_[ii];
            water [ii] = wij;
            water_[ii] = 0.0;
            ii++;
        }
    }
    //double * tmp = water; water=water_; water = tmp;
    return npix;
}

void TerrainHydraulics::flow_errosion_step( ){
// gradient evaluation
//int ii = 0;

    int ii =0;
    for(int iy=1;iy<ny-1;iy++){
        for(int ix=1;ix<ny-1;ix++){
            //int    ii   = xy2i( ix, iy );
            double hij  = ground[ii];
            double dh,dhmin=0,dhmin2=0;
            int jj,iimin;

            jj = ii-nx  ; dh = ground[jj]-hij; if( dh<dhmin ){ iimin = jj; dhmin = dh; }
            jj = ii+nx  ; dh = ground[jj]-hij; if( dh<dhmin ){ iimin = jj; dhmin = dh; }
            jj = ii-1   ; dh = ground[jj]-hij; if( dh<dhmin ){ iimin = jj; dhmin = dh; }
            jj = ii+1   ; dh = ground[jj]-hij; if( dh<dhmin ){ iimin = jj; dhmin = dh; }
            jj = ii-nx+1; dh = ground[jj]-hij; if( dh<dhmin ){ iimin = jj; dhmin = dh; }
            jj = ii+nx-1; dh = ground[jj]-hij; if( dh<dhmin ){ iimin = jj; dhmin = dh; }
            //println( i+" "+j+" "+dhmin+" "+dhmin2 );
            if( dhmin<-0.00001 ){
                /*
                double wij       = water[ii];
                water_[iimin]   += wij;
                double mud       = wij*10;
                ground[ii   ]   -= mud;
                ground[iimin]   += mud;
                */

                double wij       = water[ii];

                double sediment  = wij*c_sediment/(1-dhmin);
                wij             -= sediment;
                water_[iimin]   += wij;
                //float mud =  c_errode * wij*(-dhmin);
                double mud =  c_errode * sqrt(wij)*(-dhmin);
                // double mud =  c_errode * wij*(-dhmin);
                //mud = fmin( mud, (dhmin2 - dhmin) * 0.1  );   // WTF IS THIS ????
                ground[ii   ]  = hij - mud + sediment;
                //h[i][j]            =    f_orig*h_orig[i][j] + mf_orig*hij  - mud + sediment;
                ground[iimin] += mud*f_sediment;
            }
            ii++;
        }
    }
    //printf( "npix %i \n", npix);
}

void TerrainHydraulics::rain_and_evaporation( ){
    int ii =0;
    for(int iy=0;iy<ny;iy++){
        for(int ix=0;ix<nx;ix++){
            //int ii     = xy2i( ix, iy );
            double wij = water_[ii]        +  c_rain;
            water [ii] = water [ii]*mf_mix +  wij*f_mix;
            water_[ii] = 0.0;
            ii++;
        }
    }
}

void TerrainHydraulics::initDroplet ( double w, double disolve, double sediment, int ix0, int iy0, int ix1, int iy1 ){
    droplet_w        = w;
    droplet_disolve  = disolve;
    droplet_sediment = sediment;
    //droplet_ix = rand()%nx;
    //droplet_iy = rand()%ny;
    droplet_ix = ix0 + rand()%(ix1-ix0);
    droplet_iy = iy0 + rand()%(iy1-iy0);
    droplet_h  = ground[ xy2i( droplet_ix, droplet_iy ) ];
    //printf( "initDroplet  %i %i    %f  \n", droplet_ix, droplet_iy, droplet_h  );
}

bool TerrainHydraulics::droplet_step( ){
    if( (droplet_ix > 0)&(droplet_ix < nx)&(droplet_iy > 0)&(droplet_iy < ny) ){
        int    ii   = xy2i( droplet_ix, droplet_iy );
        int xmin,ymin;
        double h,hmin=droplet_h+droplet_w;
        bool found=false;
        h = ground[ii-nx  ]; if( h<hmin ){ ymin = -1; xmin =  0; hmin = h; found=true; }
        h = ground[ii+nx  ]; if( h<hmin ){ ymin = +1; xmin =  0; hmin = h; found=true; }
        h = ground[ii   -1]; if( h<hmin ){ ymin =  0; xmin = -1; hmin = h; found=true; }
        h = ground[ii   +1]; if( h<hmin ){ ymin =  0; xmin = +1; hmin = h; found=true; }
        h = ground[ii-nx+1]; if( h<hmin ){ ymin = -1; xmin = +1; hmin = h; found=true; }
        h = ground[ii+nx-1]; if( h<hmin ){ ymin = +1; xmin = -1; hmin = h; found=true; }
        //if( (hmin-droplet_h)<-0.00001 ){
        //found=true;
        //printf( " %i %i   %i %i   %f %f %f  \n", droplet_ix, droplet_iy,  xmin, ymin, h, hmin, droplet_h  );
        if( found ){
                //printf( " found \n" );
                //int    jj   = xy2i( droplet_ix, droplet_iy );

                //ground[ii   ]  -= droplet_dh;
                droplet_ix += xmin;
                droplet_iy += ymin;
                int    jj       = xy2i( droplet_ix, droplet_iy );
                //double hi = ground[ii];
                double hi = droplet_h;
                double hj = ground[jj];
                double dh = (hi-hj);
                //printf( " %f %f %f \n", hi, hj, dh  );
                if( dh > 0 ){
                    dh *= droplet_disolve;
                    ground[ii] = hi - dh;
                    droplet_h  = hj + dh*droplet_sediment;
                    ground[jj] = droplet_h;
                    //ground[jj] = hj + dh;
                }
                //droplet_h   = hmin;

                //ground[iimin] += mud*f_sediment;
                return false;
        }
    }
    return true;
}

void TerrainHydraulics::errodeDroples( int n, int nStepMax, double w, double disolve, double sediment, int ix0, int iy0, int ix1, int iy1 ){
    for(int i=0; i<n; i++ ){
        initDroplet( w, disolve, sediment, ix0, iy0, ix1, iy1 );
        for(int j=0; j<nStepMax; j++ ){
            if( droplet_step( ) ) break;
        }
    }
};



