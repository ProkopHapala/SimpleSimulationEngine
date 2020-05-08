
#ifndef  dataprocess1D_h
#define  dataprocess1D_h

#include "macroUtils.h"

void bisecNoise1D( int npow, double * hs, double frndMin, double frndMax ){
    int n    = 1<<npow;
    int mask = n-1;
    for(int ipow=0;ipow<npow;ipow++){
        int p  = 1<<(ipow       );
        int q2 = 1<<(npow-ipow  );
        int q  = 1<<(npow-ipow-1);
        for(int i=0;i<n;i+=q2){
            double h0 = hs[i   ];
            double h1 = hs[i+q2];
            hs[i+q]=0.5*(h0+h1) +  randf(frndMin,frndMax)*q2;
        }
    }
}


void findMax(int i0, int i1, double* hs, int& imax, double& hmax){
    for(int i=i0; i<i1; i++){
        double h = hs[i];
        if(h>hmax){ hmax=h; imax=i; }
    }
}

void runningMax( int n, int m, double* hs, double* hos){
    int    imax = 0;
    double hmax=-1e+300;
    findMax( 0, m, hs, imax, hmax);
    for(int i=m; i<n; i++){
        if(imax<(i-m)){ // max outside interval - we need new max
            hmax=-1e+300;
            findMax( i-m, i-1, hs, imax, hmax);
        }
        double h = hs[i];
        if( h>hmax){ hmax=h; imax=i; };
        hos[i] = hmax;
    }
}


inline void findMax2( int im, int im2, double& hm, double& hm2, int n, double* hs ){
    while( (hm<hm2)&&(im2<n) ){
        hm=hm2;
        im=im2;
        im2=im+1;
        hm2=hs[im2];
    }
}


void runningMax2( int n, int m, double* hs, double* hos){
    int    im=0,im2=1;
    double hm =hs[0];
    double hm2=hs[1];
    //if(hm2>hm){ _swap(hm,hm2); _swap(im,im2);  };
    findMax2( im, im2, hm, hm2, n, hs );
    for(int i=im2; i<n; i++){
        double h = hs[i];
        if(h>hm){
            hm=h;
            im=i;
            findMax2( im, im2, hm, hm2, n, hs );
        }
    }
}




int convexApprox( int n, int m, double* hs, int* is, double* his){
    int     i0= 0;
    int     i1= 1;
    double  h0=hs[i0];
    double  dh=hs[i1]-h0;
    is [0]=i0;
    his[0]=h0;
    int j=1;
    for(int i=2; i<n; i++){
        int di = i-i0;
        if(di>m){
            i0=i1;
            i1=i;
            h0=hs[i0];
            dh=(hs[i1]-h0)/(i1-i0);
            // store
            is [j]=i0;
            his[j]=h0;
            //printf( " \n", j, is[j], his[j] );
            j++;
        }else{
            double hi = hs[i];
            double h  = h0 + dh*di;
            if(h<hi){
                i1=i;
                dh=(hs[i1]-h0)/(i1-i0);
            }
        }
    }
    return j;
}


void slopeSweep( int n, double dh, double* hs, double* hos){
    int     i0= 0;
    double  h0=hs[i0];
    hos[0]=h0;
    for(int i=1; i<n; i++){
        int    di = i-i0;
        double h  = hs[i];
        double h_ = h0+dh*di;
        if(h>h_){
            i0=i;
            h0=h;
            h_=h;
        }
        hos[i]=h_;
    }
}

void slopeSweepBack( int n, double dh, double* hs, double* hos){
    int     i0= n-1;
    double  h0=hs[i0];
    hos[0]=h0;
    for(int i=n-2; i>=0; i--){
        int    di = i0-i;
        double h  = hs[i];
        double h_ = h0+dh*di;
        if(h>h_){
            i0=i;
            h0=h;
            h_=h;
        }
        hos[i]=h_;
    }
}




#endif
