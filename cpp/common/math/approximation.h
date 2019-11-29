
#ifndef  approximation_h
#define  approximation_h

#include "Lingebra.h"

/*

Numerical Approximation Techniques & Tricks

lets approximate function y=y(x)

1) Taylor Series - provide very good approximation around a point, but accuracy quickly decrease with distance from the point
2) Linear Fitting
    y = sum_i{ a_i b_i(x) } where a_i is coefficient and b_i(x) is basis function non-lineralily transforming x. For computational convenience and speed are relevant only following functions:
    2) Polynominal bais: b_i = x^i
    3) Hramonic series:  b_i = 1/x^i
3) Non-linear fitting

1) Non-lineary Transform of the function followed by fitting:
    1) Invert the function ( 1/y ), then fit using basis:
        1/y = sum_i{ a_i b_i(x) }     =>     y = 1/(sum_i{ a_i b_i(x) })
    2) n-th Suare Root of the function
        y^(1/n) = sum_i{ a_i b_i(x) }     =>     y = ( sum_i{ a_i b_i(x) } )^(n)
2) Non-lineary Transform of the argument followed by fitting:
    1) Invert the argument:
        y = sum_i{ a_i b_i( 1/x ) }
    2)  Polynominal of argument
        p = sum_i{ c_i u_i( x ) }
        y = sum_i{ a_i b_i( p ) }

*/


namespace Approx {

void polyProject( int n, int m, double* xs, double* ys, double* BB, double* By ){
    double bas[m];  // eventually we can use any function instead of xns
    for(int ip=0; ip<n; ip++ ){
        double x  = xs[ip];
        double y  = ys[ip];
        //printf( "polyProject[%i] %g %g \n", ip, x, y );
        double xn = 1; bas[0]=xn;
        for(int i=1;i<m;i++){ xn*=x; bas[i]=xn; }
        for(int i=0;i<m;i++){
            double bi  = bas[i];
            By[i]      += bi*y;
            for(int j=0;j<=i;j++){ BB[i*m+j] += bi*bas[j]; }
        }
    }
    Lingebra::symCopy( m, BB, false );
}

void polyFit( int n, int m, double* xs, double* ys, double* coefs ){
    // --- build linear system
    const int mm = m*m;
    double BB[mm];
    double By [m];
    VecN::set( mm, 0., BB );
    VecN::set(  m, 0., By );
    polyProject( n, m, xs, ys, BB, By );
    // --- solve
    double*BB_[m];
    int  index[m];
    Lingebra::from_continuous( m, m, BB, BB_ );
    Lingebra::linSolve_gauss ( m, BB_, By, index, coefs );
}

void polyeval( int n, int m, double* xs, double* ys, double* coefs ){
    for(int ip=0; ip<n; ip++){
        double x  = xs   [ip];
        double y  = coefs[0];
        double xn = x;
        for(int i=1; i<m; i++){
            y += coefs[i]*xn;
            xn*=x;
        }
        ys[ip]=y;
    }
}

void ypowsApprox( int npows, int n, int m, double* xs, double* ys, double* coefs, double* ypows, double* errs ){
    double* yis   = new double[n];
    double* ytest = new double[n];
    for(int ipow=0; ipow<npows; ipow++ ){
        double alpha = ypows[ipow];
        for(int i=0; i<n; i++){
            yis[i]=pow(ys[i], alpha );
            //printf( "[%i]y^%g [%i] %g | %g %g \n", ipow, alpha, i, yis[i], ys[i], xs[i] );
        };
        polyFit ( n, m, xs, yis,   coefs );
        polyeval( n, m, xs, ytest, coefs );
        double err2=0,errmax=0;
        for(int i=0; i<n; i++){
            double derr =  yis[i]-ytest[i];
            errmax  = fmax( errmax, fabs(derr) );
            err2   += derr*derr;
        }
        //printf("ypowsApprox[%i] %g %g \n", ipow, errmax, err2 );
        errs[ipow*2+0]=errmax;
        errs[ipow*2+1]=err2;
        //printf("ypowsApprox[%i] %g %g \n", ipow, errs[ipow*2+0], errs[ipow*2+1] );
        //errs+=2;
        coefs+=m;
    }
    delete [] yis;
    delete [] ytest;
}


}; // namespace Approx

#endif

