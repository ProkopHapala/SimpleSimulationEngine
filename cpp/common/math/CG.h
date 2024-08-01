
#ifndef  Lingebra_h
#define  Lingebra_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "fastmath.h"
#include "VecN.h"

template<typename T>
inline void dot_ax(int n,int m, const T* a, const T* b, T* sum ){ 
    //double sum[m];
    for (int k=0; k<m; k++ ){ sum[k]=0.0; } 
    for (int i=0; i<n; i++ ){
        for (int k=0; k<m; k++ ){ const int ii=i*m+k; sum[k]+=a[ii]*b[ii]; };        
    }
}

template<typename T>
inline void fma_ax(int n,int m, const T* f, const T* a, const T* b, T* out ){ 
    //double sum[m]; 
    for (int i=0; i<n; i++ ){
        for (int k=0; k<m; k++ ){ const int ii=i*m+k; out[ii] += a[ii] + b[ii]*f[k]; };        
    }
}


//template<void dotFunc(int nx, int ny, double* x, double* y)>
class LinSolver{ public:
    int n;  // dimension of vectro 
    int m;  // number of rhs
    // temp
    int    istep=0;
    double *  r=0;   // [n,m] 
    double *  r2=0;  // [n,m] 
    double *  p =0;  // [n,m]   
    double *  Ap=0;  // [n,m] 
    double *  rho=0; // [m] 
    //double rho = 0;
    //double alpha = 0;

    realloc(int n_, int m_ ){
        n=n_;m=m_;
        _realloc(p ,n*m);
        _realloc(Ap,n*m);
        _realloc(r);
        _realloc(r2);
    }

    void setLinearProblem(int n_, int m_, double* x_, double* b_, bool bRealoc=true ){
        n=n_;m=m_;
        if(bRealoc)realloc(n_,m_);
        x=x_; b=b_;
    }

    template<typename Func>
    double step0_CG(  Func dotFunc ){
        // f  = f0 - np.dot(K,x)
        // d  = f.copy()
        // f2 = np.dot(f,f)
        dot_ax(n,m, p, Ap, r );
        dotFunc  ( n,m, x, r );      // r = A*x
        VecN::sub( n*m, b, r, r ); // r = b - A*x
        VecN::set( n*m, r, p    ); // p = r
        dot_ax(n,m, r,r, err2 );
        for(int k=0;k<m;k++){ rho[k]=err2[k]; };
        return err2;
    }

    template<typename Func>
    double step_CG_simple( Func dotFunc ){

        double err2 [m];
        double alpha[m];
        double beta [m];
        double rho2 [m];

        dotFunc( n,m, p, Ap);   // Kd = K*d

        //alpha = VecN::dot(n, r, p) / VecN::dot(n, p, Ap);
        //alpha = rho / VecN::dot(n, p, Ap);   // dt  = dot(d,f) / dot(d,Kd)    # step length such that f is orthogonal to d
        
        dot_ax(n,m, p, Ap, rho2 );

        for(int k=0;k<m;k++){  alpha[k]=rho[k]/rho2[k]; }

        //printf( "### CG_step %i dt=%g rho=%g \n", istep, alpha, rho );
        fma_ax( n,m, x, p ,  alpha,   x );  // x = x + d  * dt;
        fma_ax( n,m, r, Ap, -alpha,   r );  // f = f - Kd * dt;
        
        dot_ax(n,m, r, r, err2 );
        double err2 = VecN::dot(n, r,r);

        dot_ax(n,m, p, Ap, rho2 );
        double err2tot=0;
        for(int k=0;k<m;k++){ beta[k]=err2[k]/rho[k]; err2tot+=err2[k]*err2[k]; }
        
        fma_ax(( n,m, r, p, beta,   p );    // d   = f + d*(f2/f2old)
        
        for(int k=0;k<m;k++){ rho[k]=err2[k]; };
        istep++;
        return err2tot;
    }

    template<typename Func>
    int solve_CG( Func dotFunc, int maxIters=100, double maxErr2=1e-6 ){
        istep = 0;
        step0_CG();
        int i=0;
        for ( i =0; i<maxIters; i++) {
            //printf( "[%i]========\n", istep );
            double err2 = step_CG_simple( dotFunc );
            //printf("CG[%i]x :",istep  ); VecN::print_vector(n, x);
            //printf( "[%i]err2 %g maxErr2 %g \n", istep, err2, maxErr2 );
            if ( err2< maxErr2 ){ 
                //printf("### CG converged at step %i \n", istep );
                break;
            }
        }
        return i;
    }

};


#endif

