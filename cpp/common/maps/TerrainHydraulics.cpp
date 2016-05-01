
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

void TerrainHydraulics::init_outflow(){
    for (int i=0; i<ntot; i++){
        water[i] = 1e+300;
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
    printf( "nContour_old %i nContour %i \n", nContour_old, nContour );
    for ( int ii=0; ii<nContour_old; ii++ ){
        int    i   = contour1[ii];
        double val = water[i];
        //val +=  dval/( val - ground[i] + 0.1 ); // this is just to simulate non-zero viscosity of watter
        int iy = i/nx;
        int ix = i-(iy*nx);
        //printf( " %i %i (%i,%i)\n", ii, i, ix, iy );
        // extend in four directions, check boundary overflow
        if( ix>0      ){  extend_path( val, i, i - 1  ); }
        if( ix<(nx-1) ){  extend_path( val, i, i + 1  ); }
        if( iy>0      ){  extend_path( val, i, i - nx ); }
        if( iy<(ny-1) ){  extend_path( val, i, i + nx ); }
    }
}

void TerrainHydraulics::extend_path( float val, int oi, int i ){
    // evaluate objective function for proposed path
    //val += dval*ground[i];
    //printf( " %i val %f g %f w %f k %i \n", i, val, ground[i], water[i], known[i] );
    val = fmax( val, ground[i] );
    // check if new path is advantegeneous
    if( val < water[i] ){
        water[i] = val;       //  extend path
        //printf( "mod %i  \n", i );
        //gofrom[i] = oi;     //  ... optionaly we may store oi to reconstruct path backwards later
        if( !known[i] ){
                //printf( "added %i to contour \n", i );
                known[i]=true; contour2[nContour]=i; nContour++;
        } // add to list only if not already known
    }
}


/*
void flow_errosion_step( int i ){
    // gradient evaluation
    //int ii = 0;
    int which_mask = istep%mask;
    int    imaskx   = mask [which_mask<<1  ];
    int    imasky   = mask [which_mask<<1+1];
    double cdist    = dists[which_mask];
    for(int iy=1; iy<ny; iy++){
        for(int iy=1; j<nx; ix++){
            int ii =  xy2i( ix, iy );
            double hij   = ground[ii];
            int    iimin = imins [ii];
            double dhmin = ground[iimin]-hij;
            double dh,dhmin2=0;
            int itest = ii + imask;
            dh = (ground[itest]-hij)*cdist;   if( dh<dhmin ){  iimin = itest;  imins[ii] = itest; dhmin2=dhmin; dhmin = dh;   }
            //println( i+" "+j+" "+dhmin+" "+dhmin2 );
            if( dhmin<-0.00001 ){
                double wij   = watter[ii];
                double sediment = wij*c_sediment/(1-dhmin); wij-=sediment;
                watter_[iimin] += wij;
                //float mud =  c_errode * wij*(-dhmin);
                //float mud =  c_errode * 5.0* wij*(-dhmin);
                double mud =  c_errode * sqrt(wij)*(-dhmin);
                //mud = min( mud, (dhmin2 - dhmin) * 0.1  );  // WTF IS THIS ????
                h[ii]            =   hij  - mud + sediment;
                //h[i][j]            =    f_orig*h_orig[i][j] + mf_orig*hij  - mud + sediment;
                h[iimin] += mud*f_sediment;
            }
        }
    }
}
*/


void TerrainHydraulics::flow_errosion_step( ){
// gradient evaluation
//int ii = 0;
    for(int iy=1;iy<ny-1;iy++){
        for(int ix=1;ix<ny-1;ix++){
            int    ii   = xy2i( ix, iy );
            double hij  = ground[ii];
            double dh,dhmin=0,dhmin2=0;
            int jj,iimin;

            jj = ii+nx  ; dh = ground[jj]-hij; if( dh<dhmin ){ iimin = jj; dhmin = dh; }
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
                //mud = fmin( mud, (dhmin2 - dhmin) * 0.1  );   // WTF IS THIS ????
                ground[ii   ]  = hij - mud + sediment;
                //h[i][j]            =    f_orig*h_orig[i][j] + mf_orig*hij  - mud + sediment;
                ground[iimin] += mud*f_sediment;

            }
        }
    }
}

void TerrainHydraulics::rain_and_evaporation( ){
    for(int iy=0;iy<ny;iy++){
        for(int ix=0;ix<nx;ix++){
            int ii     = xy2i( ix, iy );
            double wij = water_[ii]        +  c_rain;
            water [ii] = water [ii]*mf_mix +  wij*f_mix;
            water_[ii] = 0.0;
        }
    }
}


