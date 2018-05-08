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


