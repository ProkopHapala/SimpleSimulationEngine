#include "Grid2DAlgs.h"

void bisecNoise( int npow, double * hs, double frndMin, double frndMax ){
    int n    = 1<<npow;
    int mask = n-1;
    #define ip2i(ix,iy) ((iy&mask)<<npow)+(ix&mask)
    for(int iy=0;iy<n;iy++){ for(int ix=0;ix<n;ix++){ hs[ip2i(ix,iy)]=0.5; } };
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
            //printf( "%i %i  %f %f %f %f \n", ix+q,iy+q,  h00, h01, h10, h11 );
            //double hx          =0.5*(h00+h01) + randf(frndMin,frndMax)*fabs(h00-h01);
            //double hy          =0.5*(h00+h10) + randf(frndMin,frndMax)*fabs(h00-h10);
            //hs[ip2i(ix+q,iy+q)]=0.5*(hx+hy) + randf(frndMin,frndMax)*fabs(hx-hy);
            //hs[ip2i(ix+q,iy+q)]=0.5*(hx +hy ) + randf(frndMin,frndMax)*q;
            hs[ip2i(ix+q,iy  )]=0.5*(h00+h01) + randf(frndMin,frndMax)*q2;
            hs[ip2i(ix  ,iy+q)]=0.5*(h00+h10) + randf(frndMin,frndMax)*q2;
            //printf( "%i %i  %f %f\n", ix+q,iy+q,  hx, hy );
        } }
        for(int iy=q;iy<n;iy+=q2){ for(int ix=q;ix<n;ix+=q2){
            double h00 = hs[ip2i(ix-q,iy   )];
            double h01 = hs[ip2i(ix+q,iy   )];
            double h10 = hs[ip2i(ix  ,iy-q2)];
            double h11 = hs[ip2i(ix  ,iy+q2)];
            hs[ip2i(ix,iy)] = 0.25*( h00+h01+h00+h10 ) + randf(frndMin,frndMax)*q;
        } }
        //if (ipow>=2) break;
    }
    #undef ip2i
};

void HydraulicGrid2D::initDroplet ( double w, double disolve, double sediment, Vec2i ipmin, Vec2i ipmax ){
    droplet_w        = w;
    droplet_disolve  = disolve;
    droplet_sediment = sediment;
    droplet.x = ipmin.x + rand()%(ipmax.x-ipmin.x);
    droplet.y = ipmin.x + rand()%(ipmax.x-ipmin.x);
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
        int i    = ip2i(ip);
        double h = ground[i];
        if( h<hmin ){ ipmin = ip; hmin=h; found=true; }
    }
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
