
#ifndef  CG_h
#define  CG_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>
#include <functional> 

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
inline void dot_ax_(int n,int m, const T* A, const T* b, T* sum ){ 
    for (int k=0; k<m; k++ ){ sum[k]=0.0; } 
    for (int i=0; i<n; i++ ){
        //printf( "dot_ax_[%i] nm[%i,%i]\n", i, n,m );
        double  ai = A[i];
        for (int k=0; k<m; k++ ){ sum[k]+=b[i*m+k]*ai; };        
    }
}

void dotM_ax( int n, int m, const double* A, const double* x, double* out ){
    for (int i=0; i<n; i++ ){
        //printf( "dotM_ax[%i] nm[%i,%i] \n", i, n,m );
        dot_ax_(n,m, A+i*n, x, out+i*m );
    }
}

template<typename T>
inline void fma_ax(int n,int m, const T* a, const T* b, const T* f, T* out ){ 
    //double sum[m]; 
    for (int i=0; i<n; i++ ){
        for (int k=0; k<m; k++ ){ const int ii=i*m+k; out[ii] = a[ii] + b[ii]*f[k]; };        
    }
}

template<typename T>
inline void mul_ax(int n,int m, const T* f, const T* a,  T* out ){ 
    //double sum[m]; 
    for (int i=0; i<n; i++ ){
        for (int k=0; k<m; k++ ){ const int ii=i*m+k; out[ii] = a[ii]*f[k]; };        
    }
}

template<typename T>
inline void mul_ax_(int n,int m, const T* f, const T* a,  T* out ){ 
    //double sum[m]; 
    for (int i=0; i<n; i++ ){
        const T fi=f[i];
        for (int k=0; k<m; k++ ){  const int ii=i*m+k; out[ii] = a[ii]*fi; };        
    }
}

//template<void dotFunc(int nx, int ny, double* x, double* y)>
class CGsolver{ public:
    int n;  // dimension of vectro 
    int m;  // number of rhs
    // temp
    int    istep=0;

    double *  invD =0; // [n] Jacobi diagonal preconditoroner

    double *  x   =0;  // [n,m] residual // force
    double *  b   =0;  // [n,m] residual // force

    double *  res =0;  // [n,m] residual // force
    double *  d   =0;  // [n,m]   
    double *  Ad  =0;  // [n,m] 
    double *  z  =0;  // [n,m] 
    
    double* alpha =0;
    double* beta  =0;
    double* err   =0;
    double* err2  =0;
    
    bool bConverged=false;

    std::function<void(int,double*,double*)> dotFunc;
    
    //double rho = 0;
    //double alpha = 0;

    void realloc(int n, int m ){
        _realloc0(invD, n  , 0.0 );
        _realloc0(res,  n*m, 0.0 );
        _realloc0(d,    n*m, 0.0 );
        _realloc0(Ad,   n*m, 0.0 );
        _realloc0(z,    n*m, 0.0 );
    }

    void initDiagPrecond( double* A, double safe=0.0 ){
        for (int i=0; i<n; i++ ) invD[i] = 1.0/( A[i*n+i] );
    }

    void setLinearProblem(int n_, int m_, double* x_, double* b_, bool bRealloc=true ){
        n=n_; m=m_;
        if(bRealloc) realloc(n,m);
        x=x_; b=b_;
    }


    int solve_m1( double tol=1e-6, int maxIters=100 ){
        //printf( "CGsolver::solve_m1()\n" );
        double alpha;
        double alpha_;
        double beta;
        double err;
        double err2;
        double tol2 = tol*tol;

        istep = 0;
        //dotFunc( n, d, Ad); 
        dotFunc( n, x, Ad);    
        VecN::sub( n, b,     Ad,    res );   
        VecN::mul( n, invD,  res,    d  ); 
        err=VecN::dot( n, d,     res ); 

        // printf("invD: "); printArray(n,invD);  
        // printf("x:    "); printArray(n,x);  
        // printf("Ad:   "); printArray(n,Ad);  
        // printf("b:    "); printArray(n,b); 
        // printf("res:  "); printArray(n,res);  
        // printf("d:    "); printArray(n,d);  
        // printf("err:  %g\n", err );

        //exit(0);

        int iter=0;
        for ( iter=0; iter<maxIters; iter++) {
            dotFunc( n, d, Ad );   // Kd = K*d            
            err2=VecN::dot(n, d, Ad );
            alpha=err/err2;
            //printf( "# --- iter=%i err2(%g) alpha(%g)\n", iter, err2, alpha );
            VecN::fma( n,   x,   d,  alpha, x   );  // x = x + d  * dt;
            VecN::fma( n,  res, Ad, -alpha, res );  // f = f - Kd * dt;
            VecN::mul( n,  invD  , res,     z   );  // z = f * invD;

            err2 = VecN::dot(n, z, res );

            beta = err2/err; 
            printf( "CGsolver: iter=%i err=%g \n", iter, err2 );
            if( err2<tol2){
                bConverged=true;
                return err2;
            }
            VecN::fma( n, z, d, beta, d );    // d   = f + d*(f2/f2old)

            // printf("Ad:   "); printArray(n,Ad);  
            // printf("x:    "); printArray(n,x);  
            // printf("res:  "); printArray(n,res);  
            // printf("z:    "); printArray(n,z); 
            // printf("d:    "); printArray(n,d); 
            // printf( "alpha %g \n", alpha );
            // printf( "beta  %g \n", beta  );
            
            err = err2;
            istep++;
        }
        exit(0);
        return iter;
    }


    int solve( double tol=1e-6, int maxIters=2 ){
        //printf( "CGsolver::solve()\n" );
        double alpha [m];
        double alpha_[m];
        double beta  [m];
        double err   [m];
        double err2  [m];
        double tol2 = tol*tol;

        istep = 0;
        dotFunc( n, x, Ad);  
        VecN::sub( n*m, b,     Ad,    res );  
        mul_ax_  ( n,m, invD,  res,    d  ); 
        dot_ax   ( n,m, d,     res,   err ); 
        double err2tot=0;                     
        for(int k=0;k<m;k++){ err2tot=fmax(err2[k],err2tot); }

        // printf("invD: "); printArray(n, invD );  
        // printf("x:    "); printArray(n,3,0, x    );  
        // printf("Ad:   "); printArray(n,3,0, Ad   );  
        // printf("b:    "); printArray(n,3,0, b    ); 
        // printf("res:  "); printArray(n,3,0, res  );  
        // printf("d:    "); printArray(n,3,0, d    );  
        // printf("err:  %g %g %g\n", err[0], err[1], err[2] );

        int iter=0;
        for ( iter=0; iter<maxIters; iter++) {
            dotFunc( n, d, Ad );   // Kd = K*d            
            dot_ax(n,m, d, Ad, err2 );

            for(int k=0;k<m;k++){  
                if(err2[k]<=1e-32){ alpha[k]=0.0;  }else{  alpha[k]=err[k]/err2[k]; }; 
                alpha_[k]=-alpha[k]; 
            }

            //printf( "# --- iter=%i err2(%g,%g,%g) alpha(%g,%g,%g)\n", iter, err2[0],err2[1],err2[2], alpha[0],alpha[1],alpha[2] );

            // VecN::fma( n,   x,   d,  alpha, x   );  // x = x + d  * dt;
            // VecN::fma( n,  res, Ad, -alpha, res );  // f = f - Kd * dt;
            // VecN::mul( n,  invD  , res,     z   );  // z = f * invD;

            //printf( "### CG_step %i dt=%g rho=%g \n", istep, alpha, rho );
            fma_ax ( n,m,  x,    d, alpha,   x   );  // x = x + d  * dt;
            fma_ax ( n,m,  res, Ad, alpha_,  res );  // f = f - Kd * dt;
            mul_ax_( n,m,  invD,   res,      z   );  // z = f * invD;

            dot_ax(n,m, z, res, err2 );

            err2tot=0;         
            for(int k=0;k<m;k++){ 
                beta[k] = err2[k]/err[k]; 
                err2tot = fmax(err2[k],err2tot); 
            }
            printf( "CGsolver: iter=%i err=%g \n", iter, err2tot );
            if( err2tot<tol2){
                bConverged=true;
                return err2tot;
            }
            fma_ax( n,m, z, d, beta, d );    // d   = f + d*(f2/f2old)

            // printf("Ad:   "); printArray(n,3,0,Ad);  
            // printf("x:    "); printArray(n,3,0,x);  
            // printf("res:  "); printArray(n,3,0,res);  
            // printf("z:    "); printArray(n,3,0,z); 
            // printf("d:    "); printArray(n,3,0,d); 
            // printf( "alpha %g %g %g\n", alpha[0], alpha[1], alpha[2] );
            // printf( "beta  %g %g %g\n", beta[0], beta[1], beta[2] );
            
            for(int k=0;k<m;k++){ err[k]=err2[k]; };
            istep++;
        }
        exit(0);
        //DEBUG
        return iter;
    }

    int solve_l( double tol=1e-6, int niter=100 ){
        double tol2 = tol*tol;
        //invD = 1.0/D;
        dotFunc( n, d, Ad);
        
        double err = 0.0;
        for(int i=0; i<n; i++){
            // res = b - A * x
            // z = M * r  
            // p = copy(z)
            // r_norm_sq = dot(r, z)
            double ri  = b[i] - Ad[i];
            double zi  = ri * invD[i];
            err       += ri * zi;
            res[i]     =  ri;
            d[i]       =  zi;
        }
        int iter=0;
        for (iter=0; iter<niter; iter++ ){
            dotFunc( n, d, Ad);

            double dAd =  VecN::dot(n,d,Ad);
            double alpha = err / dAd;
            for(int i=0; i<n; i++){
                // x .+= alpha .* p
                // r .-= alpha .* Ap
                // z = M * r  
                double xi = x  [i] + d [i]*alpha;
                double ri = res[i] - Ad[i]*alpha;
                double zi =          ri   *invD[i];
                x[i]      = xi;
                res[i]    = ri;
                d[i]      = zi;
            }

            double err2 = VecN::dot(n,res,z);
        
            if ( err2 < tol2 ){ break; }
            double beta = err2 / err;
            for(int i=0; i<n; i++){
                d[i] = d[i]*beta + z[i];
            }
            err = err2;
        }
        return iter;
    }


/*
    double step0(){
        // res = b - A * x
        // z = M * r  
        // p = copy(z)
        // r_norm_sq = dot(r, z)
        dotFunc( n, d, Ad);   // Kd = K*d
        VecN::sub( n*m, b, Ad,    res ); 
        mul_ax( n,m, res, invD,   d   );
        dot_ax( n,m, res,  d,     err );
        double err2tot=0;
        for(int k=0;k<m;k++){ err2tot += err[k]*err[k]; }
        return err2tot;
    }

    double step(){
        dotFunc( n, d, Ad);   // Kd = K*d

        //alpha = VecN::dot(n, r, p) / VecN::dot(n, p, Ap);
        //alpha = rho / VecN::dot(n, p, Ap);   // dt  = dot(d,f) / dot(d,Kd)    # step length such that f is orthogonal to d
        
        dot_ax(n,m, d, Ad, err2 );

        for(int k=0;k<m;k++){  alpha[k]=err[k]/err2[k]; }

        //printf( "### CG_step %i dt=%g rho=%g \n", istep, alpha, rho );
        fma_ax( n,m,  alpha, x,    d,     x  );  // x = x + d  * dt;
        fma_ax( n,m, -alpha, res, Ad,    res );  // f = f - Kd * dt;
        mul_ax( n,m, invD, res,          z   );  // z = f * invD;

        double err2 = VecN::dot(n, z,r);

        double err2tot=0;
        for(int k=0;k<m;k++){ 
            beta[k]  = err2[k]/err[k]; 
            err2tot += err2[k]*err2[k]; 
        }
        // if( err2tot<tol){
        //     bConverged=true;
        //     return err2tot;
        // }
        fma_ax(( n,m,  beta,  res, d,  d );    // d   = f + d*(f2/f2old)
        
        for(int k=0;k<m;k++){ err[k]=err2[k]; };
        istep++;
        return err2tot;
    }

    //void makeInvDiag(){}

    int solve_steps( int maxIters=100, double tol=1e-6 ){
        istep = 0;
         double tol2 = tol*tol;
        step0( dotFunc );
        int i=0;
        for ( i =0; i<maxIters; i++) {
            //printf( "[%i]========\n", istep );
            double err2 = step( dotFunc );
            //printf("CG[%i]x :",istep  ); VecN::print_vector(n, x);
            //printf( "[%i]err2 %g maxErr2 %g \n", istep, err2, maxErr2 );
            if ( err2< tol2 ){ 
                //printf("### CG converged at step %i \n", istep );
                break;
            }
        }
        return i;
    }
*/
};


#endif

