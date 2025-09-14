#ifndef hydraulics1D_h
#define hydraulics1D_h

void bisectNoise1D(int npow, double* hs, double frndMin, double frndMax ){
    int n    = 1<<npow;
    int mask = n-1;
    for(int ipow=0;ipow<npow;ipow++){
        int p  = 1<<(ipow       );
        int q2 = 1<<(npow-ipow  );
        int q  = 1<<(npow-ipow-1);
        for(int i=0;i<n;i+=q2){
            // Use periodic wrap to avoid out-of-bounds at i+q2==n
            double h0 = hs[i];
            double h1 = hs[(i+q2) & mask];
            hs[(i+q) & mask] = 0.5*(h0+h1) + randf(frndMin,frndMax)*q2;
        }
    }
}


class Hydraulics1D{ public:
    int n = 0;
    double* ground = 0;
    double* water  = 0;

    void step( double rain, double evapor ){
        for(int i=1; i<n; i++){
            double dg   = ground[i] - ground[i-1];
            double w0   = fmax( 0, water[i-1]-evapor );
            double w1   = fmax( 0, water[i  ]-evapor );
            double wtot = w0 + w1 + 2*rain;
            double w=0;
            if(dg>0){
                if(dg>wtot){
                    water [i-1] = wtot;
                    water [i  ] = 0;
                }else{
                    w = 0.5*(wtot-dg);
                    water [i-1] = w + dg;
                    water [i  ] = w;
                }
            }else{
                if(-dg>wtot){
                    water [i-1] = 0;
                    water [i  ] = wtot;
                }else{
                    w = 0.5*(wtot+dg);
                    water [i-1] = w;
                    water [i  ] = w - dg;
                }
            }
            //printf( "%i : (dg,wtot,w): %g %g %g \n", i, dg, wtot, w );
        }
    }

    void deepAccel(double depth){
        for(int i=0; i<n; i++){
            double h1,h2;
            double wsum = 0;
            int j =i;
            while(true){
                double g=ground[j];
                double w=water [j];
                double h=w+g;
                if(i==j){
                    //double h0 = h+w*0.5;
                    h1=g+1e-8;
                    h2=h1+depth;
                }
                if( (j<n)&&(g<h1)&&(h>h2) ){
                    //printf( "h,h1 %g %g \n", h, h1 );
                    wsum+=h-h1;
                    j++;
                }else{
                    break;
                }
            }
            int nj = j-i;
            //printf( "%i(%i) : %g \n", i, nj, wsum );
            double hav = h1 + wsum/nj;
            double wsum1=0, wsum2=0;
            while(i<j){
                wsum1    += water[i];
                water[i]  = hav - ground[i];
                wsum2    += water[i];
                i++;
            }
            //printf( "%i(%i) : %g %g %g   %g \n", i, nj, wsum, wsum1, wsum2, hav );
        }
        //exit(0);
    }



    void clear(){
        for(int i=0; i<n; i++){ ground[i]=0; water[i]=0; }
    }

    void realloc(int n_){
        n=n_;
        _realloc(ground,n);
        _realloc(water ,n);
    }

};

#endif
