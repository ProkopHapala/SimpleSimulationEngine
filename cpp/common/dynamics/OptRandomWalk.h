
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

    double biasDir;
    double biasSpread;

    EnergyFunction getEnergy = 0;

    double stepSize=0.1;
    double Ebest,E;

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
        for( int i=0; i<n; i++){
            //X[i] = X[i] + (randf()-0.5);
            X[i] = Xbest[i] + (randf()-0.5)*scales[i]*step;
            //X[i] = Xbest[i] + (randf()-0.5)* ( scales[i] + biasSpread*(Xbest[i]-Xdir) );
        }
    }

    void run( int nstep ){
        for(int i=0; i<nstep;i++){
            mutate( stepSize );
            E = getEnergy( n, X );
            if(E<Ebest){
                printf( "Energy improved %g -> %g \n", Ebest, E );
                VecN::set(n,X,Xbest);
                Ebest=E;
            }
        }
    }

};

#endif
