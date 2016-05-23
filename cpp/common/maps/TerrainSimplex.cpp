
#include "TerrainSimplex.h"  // THE HEADER


int TerrainSimplex::raster_line( Vec2d dirHat, Vec2d pos0, Vec2d pos1, Vec2d * hits, int * boundaries, int * edges ){
    double t0    = dirHat.dot( pos0 );
    double t1    = dirHat.dot( pos1 );
    double tspan = t1-t0;
    double pa,pb,pc, invPa,invPb,invPc;
    int    ia,ib,ic,i;
    printf( " %f %f \n", step, invStep );
    double mda,mdb,mdc, mta,mtb,mtc;
    pa  = dirHat.dot( { 0.0d        ,1.15470053839*invStep} );
    pb  = dirHat.dot( { 1.0d*invStep,0.57735026919*invStep} );
    pc  = dirHat.dot( {-1.0d*invStep,0.57735026919*invStep} );
    mda = pos0.dot  ( { 0.0d        ,1.15470053839*invStep} );
    mdb = pos0.dot  ( { 1.0d*invStep,0.57735026919*invStep} );
    mdc = pos0.dot  ( {-1.0d*invStep,0.57735026919*invStep} );
    if( pa < 0 ){ pa=-pa; mda = 1-mda; };
    if( pb < 0 ){ pb=-pb; mdb = 1-mdb; };
    if( pc < 0 ){ pc=-pc; mdc = 1-mdc; };
    ia=(int)(mda + MAP_OFFSET);   mda = 1-(mda - (ia - MAP_OFFSET) );
    ib=(int)(mdb + MAP_OFFSET);   mdb = 1-(mdb - (ib - MAP_OFFSET) );
    ic=(int)(mdc + MAP_OFFSET);   mdc = 1-(mdc - (ic - MAP_OFFSET) );
    invPa = 1/pa; invPb = 1/pb; invPc = 1/pc;
    //printf( " t_1,2  %f %f   p_a,b,c %f %f %f  \n", t0, t1, pa, pb, pc );
    printf( " pa invPa \n", pa, invPa );
    double t = 0;
    i=0;
    int ia_,ib_;
    simplexIndexBare( pos0.x, pos0.y, ia_, ib_ );
    while( t<tspan ){
        double tma = mda * invPa;
        double tmb = mdb * invPb;
        double tmc = mdc * invPc;
        //t += tma; boundaries[i] = 0;  mda = 1;
        //t += tmb; boundaries[i] = 1;  mdb = 1;
        //t += tmc; boundaries[i] = 2;  mdc = 1;
        int ii = i<<2;
        if( tma < tmb ){
           if( tma < tmc ){  // a min
                t    += tma;
                mda   = 1; mdb -= pb*tma; mdc -= pc*tma;
                boundaries[i] = 0; ia_++;
                edges[ii  ] = ia_; edges[ii+1] = ib_;
                edges[ii+2] = ia_; edges[ii+3] = ib_+1;
           }else{            // c min
                t    += tmc;
                mda  -= pa*tmc; mdb -= pb*tmc; mdc = 1;
                boundaries[i] = 2; ib_++;
                edges[ii  ] = ia_; edges[ii+1] = ib_;
                edges[ii+2] = ia_+1; edges[ii+3] = ib_;
           }
        }else{
           if( tmb < tmc ){  // b min
                t    += tmb;
                mda  -= pa*tmb; mdb = 1; mdc -= pc*tmb;
                boundaries[i] = 1;
                edges[ii  ] = ia_;   edges[ii+1] = ib_+1;
                edges[ii+2] = ia_+1; edges[ii+3] = ib_;
           }else{            // c min
                t    += tmc;
                mda  -= pa*tmc; mdb -= pb*tmc; mdc = 1;
                boundaries[i] = 2; ib_++;
                edges[ii  ] = ia_; edges[ii+1] = ib_;
                edges[ii+2] = ia_+1; edges[ii+3] = ib_;
           }
        }
        hits[i].set( pos0 );
        hits[i].add_mul( dirHat, t );
        //hits[i].set_mul( dirHat, t );
        //printf( "%i %i  (%f,%f,%f)     %f (%f,%f) \n", i, boundaries[i], tma, tmb, tmc,       t, hits[i].x, hits[i].y );
        printf( "%i %i  (%f,%f,%f)     %f (%f,%f) \n", i, boundaries[i], mda, mdb, mdc,       t, hits[i].x, hits[i].y );
        i++;
    }
    return i;
}


void TerrainSimplex::genTerrainNoise( int n, double scale,  double hscale,  double fdown, double strength, int seed, const Vec2d& pos0 ){
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

void TerrainSimplex::init_outflow(){
    for (int i=0; i<ntot; i++){
        water[i] = 1e+300;
        known[i] = false;
    }
}

void TerrainSimplex::outflow_step(){
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

void TerrainSimplex::extend_path( float val, int oi, int i ){
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

void TerrainSimplex::initErrosion( double w ){
    for (int i=0; i<ntot; i++){
            if(ground[i]>1.0) ground[i] = 1.0;
            if(ground[i]<0.0) ground[i] = 0.0;
            water[i] = w; water_[i] = 0.0;
    }
}

int TerrainSimplex::flow_errosion_step_noRain( ){
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

void TerrainSimplex::flow_errosion_step( ){
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

void TerrainSimplex::rain_and_evaporation( ){
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

void TerrainSimplex::initDroplet ( double size_ ){
    droplet_size = size_;
    droplet_ix = rand()%nx;
    droplet_iy = rand()%ny;
    droplet_h = ground[ xy2i( droplet_ix, droplet_iy ) ];

    //printf( "initDroplet  %i %i    %f  \n", droplet_ix, droplet_iy, droplet_h  );
}

bool TerrainSimplex::droplet_step( ){
    if( (droplet_ix > 0)&(droplet_ix < nx)&(droplet_iy > 0)&(droplet_iy < ny) ){
        int    ii   = xy2i( droplet_ix, droplet_iy );
        int xmin,ymin;
        double h,hmin=droplet_h;
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
                ground[ii   ]  -= droplet_size;
                droplet_h   = hmin;
                droplet_ix += xmin;
                droplet_iy += ymin;
                //ground[iimin] += mud*f_sediment;
                return false;
        }
    }
    return true;
}

void TerrainSimplex::errodeDroples( int n, int nStepMax, double size_ ){
    for(int i=0; i<n; i++ ){
        initDroplet( size_ );
        for(int j=0; j<nStepMax; j++ ){
            if( droplet_step( ) ) break;
        }
    }
};



