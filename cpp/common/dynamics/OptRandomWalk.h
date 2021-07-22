
#ifndef OptRandomWalk_h
#define OptRandomWalk_h

//#include <cstddef>
#include <math.h>
#include "macroUtils.h"
#include "VecN.h"

//typedef void (*ForceFunction)( int n, double * xs, double * dfs );

typedef double (*EnergyFunction)( int n, double * Xs );

class OptRandomWalk{ public:

    // ====  Variables

    int n=0;
    double * X     = 0;
    double * Xbest = 0; // currently best solution
    double * Xdir  = 0; // decent solution at reasonable distance from Xbest (not too far, not too close) which approximates soft-mode direction
    double * scales = 0;

    int nTryDir = 0;

    //double biasDir;
    //double biasSpread;

    EnergyFunction getEnergy = 0;

    double stepSize=0.1;
    double Ebest,E;
    double Edir,Rdir; // dir bias

    // ====  Functions

    void realloc( int n_, double* X_ = 0 ){
        printf( "OptRandomWalk::realloc    *X %li %li \n", X_, X_+(n_-1) );
        n=n_;
        if(X_){ X=X_; }else{ _realloc(X,n); };
        _realloc(Xbest ,n);
        _realloc(Xdir  ,n);
        _realloc(scales,n);
        for(int i=0; i<n; i++){
            //Xbest[i]=
            Xbest [i]=X[i];
            scales[i]=1;
        }
    }

    void dealloc( double bX=true ){
        _dealloc(X     );
        _dealloc(Xbest );
        _dealloc(Xdir  );
        _dealloc(scales);
    }

    void mutate( double step ){
        double rndDir = randf();
        double cfw    = 0.5 + 1.0/nTryDir;
        double cbak   = 0.5;
        double c      = cfw+cbak;
        for( int i=0; i<n; i++){
            //X[i] = X[i] + (randf()-0.5);
            //X[i] = Xbest[i] + (randf()-0.5)*scales[i]*step;
            //X[i] = Xbest[i] + (randf()-0.5)*scales[i]*step  +  (Xdir[i]-Xbest[i])*rndDir;
            //X[i] = Xbest[i] + (randf()-0.5)*scales[i]*step  +  (Xdir[i]-Xbest[i])*(rndDir*3.0-1.0);
            //X[i] = Xbest[i] + (randf()-0.5)*scales[i]*step  +  (Xdir[i]-Xbest[i])*(rndDir*0.75-0.25);
            X[i] = Xbest[i] + (randf()-0.5)*scales[i]*step  +  (Xdir[i]-Xbest[i])*(rndDir*1.2-0.2);
            //X[i] = Xbest[i] + (randf()-0.5)*scales[i]*step  +  (Xdir[i]-Xbest[i])*(rndDir*c-cbak);
        }
    }

    void run( int nstep ){
        for(int i=0; i<nstep;i++){
            mutate( stepSize );
            E = getEnergy( n, X );
            if(E<Ebest){
                printf( "Energy improved %g -> %g \n", Ebest, E );
                VecN::set(n,X,Xbest);
                //Rdir = VecN::err2( n, X, Xbest );
                Edir = 1e+300;
                Ebest=E;
            }else{
                double R = VecN::err2( n, X   , Xbest );
                Rdir     = VecN::err2( n, Xdir, Xbest );
                double curv    = (E   -Ebest)/(R*R);
                double curvDir = (Edir-Ebest)/(Rdir*Rdir);
                if( curv<curvDir ){
                    //printf( "dir update: (%g|dE %g R %g ) -> ( %g|dE %g R %g ) \n", curvDir,(Edir-Ebest), Rdir,    curv,(E   -Ebest), R    );
                    VecN::set(n,X,Xdir);
                    Rdir=R;
                    Edir=E;
                    nTryDir=0;
                }else{
                    nTryDir++;
                }
            }
        }
    }

    void start(){
         Ebest = getEnergy( n, X );
    }

};

#endif
