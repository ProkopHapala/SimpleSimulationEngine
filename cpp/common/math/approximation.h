
#ifndef  approximation_h
#define  approximation_h

#include <vector>

#include "Lingebra.h"
#include "macroUtils.h"

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

void polyProject( int n, int m, int* pows, double* xs, double* ys, double* BB, double* By, double* ws=0 ){
    double bas[m];  // eventually we can use any function instead of xns
    int maxOrder = pows[m-1];

    //printf("polyProject: n %i m %i maxOrder %i \n", n, m, maxOrder );
    for(int ip=0; ip<n; ip++ ){
        double x  = xs[ip];
        double y  = ys[ip];
        double w  = 1;
        if(ws){ w = ws[ip]; y*=w; }
        double xn = 1;
        int ipow  = 0;
        //printf( "xyw %g %g %g \n",  x, y, w );
        for(int i=0;i<=maxOrder;i++){
            if(i==pows[ipow]){
                bas[ipow]=xn*w;
                ipow++;
            }
            xn*=x;
            //if(ip==0) printf( "i %i ipow %i xn %g x %g \n", i, ipow, xn, x );
            //if(ip==0) printf( "bas[%i] %g \n", i, bas[i] );
        }
        for(int i=0;i<m;i++){
            double bi  = bas[i];
            By[i]      += bi*y;
            //printf( "By[%i] %g %g %g \n", i, bi, y, By[i] );
            for(int j=0;j<=i;j++){
                //if(ip==0) printf( "i %i j %i i*m+j %i \n", i, j, i*m+j );
                BB[i*m+j] += bi*bas[j];
            }
        }
        //printf( "polyProject[%i] %g %g \n", ip, x, y );
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
    polyProject( n, m, xs, ys, BB, By );  // ToDo: This can be optimized: we don't need to recalculate whole BB, By matrix for each polynominal order (just the highest)
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

void polyeval( int n, int m,  int* pows, double* xs, double* ys, double* coefs ){
    int maxOrder = pows[m-1];
    //for(int i=0;i<=maxOrder;i++){ printf( "coefs[%i] %g \n ", i, coefs[i] ); }
    for(int ip=0; ip<n; ip++){
        double x  = xs[ip];
        double y  = 0;
        double xn = 1;
        int ipow  = 0;
        for(int i=0;i<=maxOrder;i++){
            if(i==pows[ipow]){
                y += coefs[ipow]*xn;
                ipow++;
            }
            xn*=x;
        }
        //printf( "polyeval[%i] %g %g %g \n", ip, x, y, xn );
        ys[ip]=y;
    }
}

/*
void polyFitMultiOrder( int n, int m0, int m1, int dm, double* xs, double* ys, double* coefs ){
    // --- build linear system
    const int mm = m1*m1;
    double BB[mm];
    double By[m1];
    VecN::set( mm, 0., BB );
    VecN::set( m1, 0., By );
    polyProject( n, m1, xs, ys, BB, By );  // ToDo: This can be optimized: we don't need to recalculate whole BB, By matrix for each polynominal order (just the highest)
    // --- solve
    double*BB_[m1];
    int  index[m1];
    for( int m=m0; m<=m1; m+=dm ){
        Lingebra::from_continuous( m1, m1, BB, BB_ );
        //Lingebra::submat( m1,m1, 0,0, BB, BB_ );
        Lingebra::linSolve_gauss ( m, BB_, By, index, coefs );
        coefs+=m;
    }
}
*/

/*
int multiFitAndCheck( int m0, int m1, int dm, double errTol, double* xs, double* ys_ref, double* ys_test, const double alpha=1, const bool relative=false ){
    const int mm = m1*m1;
    double BB[mm];
    double By[m1];
    VecN::set( mm, 0., BB );
    VecN::set( m1, 0., By );
    polyProject( n, m, xs, ys, BB, By );  // ToDo: This can be optimized: we don't need to recalculate whole BB, By matrix for each polynominal order (just the highest)
    // --- solve
    double*BB_[m1];
    int  index[m1];
    double Cs [m1];
    int m_conv=-1;
    for( int m=m0; m<=m1; m+=dm ){
        Lingebra::from_continuous( m1, m1, BB, BB_ );
        //Lingebra::submat( m1,m1, 0,0, BB, BB_ );
        Lingebra::linSolve_gauss ( m, BB_, By, index, Cs );
        //coefs+=m;
        polyeval( n, m, xs, ys, Cs );
        double errmax =0;
        double err0   =0;
        for(int i=0; i<np; i++){
            //double derr = pow( ys_test[i], alpha ) - yn_ref[i0+i]; // ?????? Premature optimization ?????
            double derr
            if(alpha==1){ derr =      ys_test[i]          -      ys_ref[i0+i]          ; }
            else        { derr = pow( ys_test[i], alpha ) - pow( ys_ref[i0+i] , alpha ); }
            if(relative) derr/=(ys[i]+errTol);
            errmax  = fmax( errmax, fabs(derr) );
            err2   += derr*derr;
        }
        printf( "multiFitAndCheck order[%i] %g %g \n", m, errmax, err2 );
        if(errmax<errTol){ m_conv=m; break;}
    }
    return m_conv;
}
*/

int ascendingPolyFit( double& err, int ipow0, int npows, int* pows, int n, double* xs, double* ys_ref, double* ys_test, double* coefs, double* ws=0, const double alpha=1, bool allCoefs=false, const bool relative=false ){
    double errTol = err;
    //const int m1 = pows[npows-1];
    const int m1 = npows;
    const int mm = m1*m1;
    //printf("multiFitAndCheck: m1 %i mm %i \n", m1, mm );
    //printf("multiFitAndCheck: ipow0 %i npows %i m1 %i mm %i |  pows %i %i %i %i \n", ipow0, npows, m1, mm,     pows[0], pows[1], pows[2], pows[npows-1] );
    double BB[mm];
    double By[m1];
    VecN::set( mm, 0., BB );
    VecN::set( m1, 0., By );
    //polyProject( n, m, xs, ys, BB, By );
    //printf("multiFitAndCheck: m1 %i mm %i \n", m1, mm );
    polyProject( n, npows, pows, xs, ys_ref, BB, By, ws );
    //printf("multiFitAndCheck: m1 %i mm %i \n", m1, mm );

    // --- solve
    double BBcp[mm];
    double Bycp[m1];
    double*BB_ [m1];
    int  index[m1];
    //double Cs [m1];
    int m_conv=-1;
    Lingebra::from_continuous( m1, m1, BBcp, BB_ );

    //for(int i=0;i<m1;i++){
    //    for(int j=0;j<m1;j++){ printf("%3.3f ", BB[i*m1+j] ); };
    //    printf(" BB \n");
    //}
    //for(int i=0;i<m1;i++){
    //    for(int j=0;j<m1;j++){ printf("%3.3f ", BB_[i][j] ); };
    //    printf(" BB_ \n");
    //}
    //return 0;

    for( int ipw=ipow0; ipw<npows; ipw++ ){
        //int pwi = pows[ipw];
        //printf("m1 %i mm %i ipw %i \n", m1, mm, ipw );
        for(int i=0;i<ipw;i++){ Bycp[i] = By[i];  int ioff=i*m1; for(int j=0;j<ipw;j++){ BBcp[ioff+j] = BB[ioff+j]; }; }  // must copy otherwise linSolve_gauss destroy it
        //printf("By:  "); for(int j=0;j<ipw;j++){ printf("%3.3f ", By[j] ); };
        //for(int i=0;i<ipw;i++){
        //    for(int j=0;j<ipw;j++){ printf("%3.3f ", BB_[i][j] ); };
        //    printf("\n");
        //}
        Lingebra::linSolve_gauss ( ipw, BB_, Bycp, index, coefs );

        polyeval( n, ipw, pows, xs, ys_test, coefs );
        double errmax =0;
        double err2   =0;
        for(int i=0; i<n; i++){
            //double derr = pow( ys_test[i], alpha ) - yn_ref[i0+i]; // ?????? Premature optimization ?????
            double derr;
            //printf( "[%i][%i] %g %g \n", ipw, i, ys_test[i], ys_ref[i] );
            if(alpha==1){ derr =      ys_test[i]          -      ys_ref[i]          ; }
            else        { derr = pow( ys_test[i], alpha ) - pow( ys_ref[i] , alpha ); }
            if(relative) derr/=(ys_ref[i]+errTol);
            errmax  = fmax( errmax, fabs(derr) );
            err2   += derr*derr;
        }
        double errRMS = sqrt(err2/n);
        printf( "multiFitAndCheck[%i] m1[%i] %g %g \n", ipw, pows[ipw], errmax, errRMS );
        err=errmax;
        if(errmax<errTol){  m_conv=ipw;  break;}
        if(allCoefs) coefs+=ipw;
    }
    return m_conv;
}

/*
    int ypowsApprox( int npows, int n, int m, double* xs, double* ys, double* coefs, double* ypows, double* errs, bool relative ){
        double* yis   = new double[n];
        double* ytest = new double[n];
        int    ipow_best = 0;
        double errBest=1e+300;
        for(int ipow=0; ipow<npows; ipow++ ){
            // --- function transform
            double alpha    = ypows[ipow];
            double invAlpha = 1/alpha;
            for(int i=0; i<n; i++){
                yis[i]=pow(ys[i], invAlpha );
                //printf( "[%i]y^%g [%i] %g | %g %g \n", ipow, alpha, i, yis[i], ys[i], xs[i] );
            };
            // --- poly fit
            polyFit ( n, m, xs, yis,   coefs );
            polyeval( n, m, xs, ytest, coefs );
            // --- error evaluation
            double err2=0,errmax=0;
            for(int i=0; i<n; i++){
                //double derr =  yis[i]-[i];
                double derr = pow( ytest[i], alpha ) - ys[i];
                if(relative) derr/=ys[i];
                errmax  = fmax( errmax, fabs(derr) );
                err2   += derr*derr;
            }
            //printf("ypowsApprox[%i] %g %g \n", ipow, errmax, err2 );
            if(errmax<errBest){ errBest=errmax; ipow_best=ipow; };
            errs[ipow*2+0]=errmax;
            errs[ipow*2+1]=err2;
            //printf("ypowsApprox[%i] %g %g \n", ipow, errs[ipow*2+0], errs[ipow*2+1] );
            //errs+=2;
            coefs+=m;
        }
        delete [] yis;
        delete [] ytest;
        return ipow_best;
    }

    void fitOnIntervals( int npows, int n, int m, double* xs, double* ys, double* coefs, double* ypows, double* errs, bool relative, int nint, int* i0s, int* i1s ){
        for(int ii=0; ii<nint; ii++){
            int i0        = i0s[ii];
            int ni        = i1s[ii]-i0s[ii];
            int ipow_best = ypowsApprox( npows, ni, m, xs+i0, ys+i0, coefs, ypows, errs, relative );
        }
    }
*/


inline double pow2m(double x, int m){
    for(int i=0;i<m;i++){ x*=x; }
    return x;
};


struct ApproxVariant{
    double  cost;    // computational cost estimate (heuristics)
    double  errMax;  //
    double  errRMS;  //
    int     i0,i1;   // initial end final starting point
    int     ipow;    // power of y_ref
    int     ncoefs;  // order of polynominal
    double* coefs;   // pointer to array of cooefs

    /*
    //double evalError(double* xs_ref, double* yn_ref, bool relative, double* ys_test=0, double* pows=0 ){
    double evalError(double* xs, double* ys_ref, powered=false, bool relative=false, double* ys_test=0, double* pows=0 ){
        int np = i1-i0;
        bool bAlloc = (ys_test==0);
        if(bAlloc) ys_test = new double[np];
        double alpha;
        if(pows){ alpha=pows[ipow]; }else{ alpha=(1<<ipow); };
        polyeval( np, ncoefs, xs+i0, ys_test, coefs );
        errMax=0;
        errRMS=0;
        for(int i=0; i<np; i++){
            //double derr = pow( ys_test[i], alpha ) - yn_ref[i0+i]; // ?????? Premature optimization ?????
            double derr = pow( ys_test[i], alpha ) - pow( ys_ref[i0+i] , alpha );
            if(relative) derr/=ys[i];
            errmax  = fmax( errmax, fabs(derr) );
            err2   += derr*derr;
        }
        if(bAlloc) delete [] ys_test;
    }

    double polyfit( double* xs, double* ys_ref, bool powered=false, double* pows=0 ){
        int np = i1-i0;
        double* ys = ys_ref+i0;
        if(!powered){
            double invAlpha;
            ys = new double[np];
            if(pows){ invAlpha=1/pows[ipow]; }else{ invAlpha=1./(1<<ipow); };
            for(int i=0; i<np; i++){
                ys[i]=pow(ys_ref[i], invAlpha );
                //printf( "[%i]y^%g [%i] %g | %g %g \n", ipow, alpha, i, yis[i], ys[i], xs[i] );
            }
        }
        polyFit ( n, m, xs+i0, ys, coefs ); // ToDo: This can be optimized: we don't need to recalculate whole BB, By matrix for each polynominal order (just the highest)
        if(!powered) delete [] ys;
    }
    */

};


class AutoApprox{ public:

    double errTol         = 1e-6;
    double errmax =0;  // ToDo make this class variable
    double err2   =0;
    bool bRelativeError = false;
    bool bAllCoefs      = true;

    double err;

    int np=0;
    double* xs       = 0;  // x-argument values
    double* ys_ref   = 0;  // reference function values
    double* ys_test  = 0;
    double* ys_trans = 0;  // transformed function

    //double* ws       = 0;
    Buf<double> ync;  // fixed component of function to be added befor exponentiation
    Buf<double> ws;   // weights to be used for sampling points in least square fit

    int     ipoly0 = 3;
    int     npoly  = 0;
    int*    ipolys = 0;
    double* coefs  = 0;

    int     npows = 0;
    double* pows  = 0;


    double* BB = 0;
    double* By = 0;

    int*    index = 0;
    double* BBcp = 0;
    double* Bycp = 0;
    double**BB_  = 0;

    std::vector<ApproxVariant> variants;

    // ===== Functions

    void bindOrRealloc(  int npoly_, int npows_, int np_, int* ipolys_=0, double* pows_=0, double* xs_=0, double* ys_ref_=0 ){
        np   =np_;
        npoly=npoly_;
        npows=npows_;
        _bindOrRealloc( np   , xs_    , xs     );
        _bindOrRealloc( np   , ys_ref_, ys_ref );
        ys_test  = new double[np      ];
        ys_trans = new double[np*npows];
        if(bAllCoefs){ coefs    = new double[npoly*npoly]; }
        else         { coefs    = new double[npoly      ]; };

        _bindOrRealloc( npoly, ipolys_, ipolys );
        _bindOrRealloc( npows, pows_  , pows   );
        //printf( "npoly %i ipolys[1,2,3] %i %i %i \n", npoly, ipolys[0], ipolys[1], ipolys[2] );
    }
    // ToDo :: deallocate or destructor (?)

    void preparePowers(){
        for(int ipow=0;ipow<npows;ipow++){
            double invAlpha = 1./pows[ipow];
            //printf( "ipow %i invAlpha %g \n", ipow, invAlpha );
            double *ys_i = ys_trans + ipow*np;
            for(int i=0; i<np; i++){
                double       y  = pow( ys_ref[i], invAlpha );
                if(ync.data) y -= ync[i];
                ys_i[i] = y;
            }
        }
    }

    inline double* getYrefPow(int ipow ){ return ys_trans + ipow*np; }


    int tryVariant( int i0, int i1, int ipow ){
        //double evalError( xs+i0, ys_trans + ipow*np + i0, powered=false, bool relative=false, double* ys_test=0, double* pows=0 ){
        int     ni   = i1-i0;
        err =  errTol;
        int i = ascendingPolyFit( err, ipoly0, npoly, ipolys, ni, xs+i0, getYrefPow(ipow)+i0, ys_test, coefs, ws.data, pows[ipow], bRelativeError );
        return i;
    }

    // ======= ascendingPolyFit

    void reallocAux(){
        int mm = npoly*npoly;
        _realloc(BB  ,mm);
        _realloc(BBcp,mm);
        _realloc(By  ,npoly);
        _realloc(Bycp,npoly);
        _realloc(BB_ ,npoly);
        _realloc(index,npoly);
        Lingebra::from_continuous( npoly, npoly, BBcp, BB_ );
    }

    double checkError( double alpha, int i0, int n ){
        //printf( "checkError \n" );
        errmax =0;  // ToDo make this class variable
        err2   =0;
        for(int i=0; i<n; i++){
            double derr;
            double y    = ys_test[i0+i];
            double yref = ys_ref [i0+i];
            //if(ynf)   y*=ynf;
            if(ync.data)   y+=ync.data[i0+i];
            if(alpha!=1){ y = pow( y, alpha ); }
            derr =      y - yref;
            if(bRelativeError) derr/=(yref+errTol);
            errmax  = fmax( errmax, fabs(derr) );
            err2   += derr*derr;
            //printf( "[%i] y %g ys_test %g yref %g derr %g errmax %g \n", i, y, ys_test[i0+i], yref, derr, errmax );
        }
        //errRMS = sqrt(err2/n);
        return errmax;
    }

    int ascendingPolyFit_(int i0, int i1, int ipow ){
        int n = i1-i0;
        //printf("i0 %i i1 %i ipow %i n %i \n", i0, i1, ipow, n ); exit(0);

        double* xs_ = xs               + i0;
        double* ys_ = getYrefPow(ipow) + i0;
        double mm = npoly*npoly;
        VecN::set( mm,    0., BB );
        VecN::set( npoly, 0., By );
        DEBUG
        polyProject( n, npoly, ipolys, xs_, ys_, BB, By, ws.data );
        int m_conv=-1;
        DEBUG
        Lingebra::from_continuous( npoly, npoly, BBcp, BB_ );
        DEBUG

        /*
        //for(int i=0;i<npoly;i++){ Bycp[i] = By[i];  int ioff=i*npoly; for(int j=0;j<npoly;j++){ BBcp[ioff+j] = BB[ioff+j]; }; }
        for(int i=0;i<npoly;i++){
            for(int j=0;j<npoly;j++){ printf("%3.3f ", BB[i*npoly+j] ); };
            printf(" BB By %g \n", By[i] );
        }

        for(int i=0;i<npoly;i++){
            for(int j=0;j<npoly;j++){ printf("%3.3f ", BB_[i][j] ); };
            printf(" BB_ \n");
        }
        */
        //return 0;

        for( int ipw=ipoly0; ipw<npoly; ipw++ ){
            //printf( "=================== ipw %i \n", ipw );
            for(int i=0;i<ipw;i++){ Bycp[i] = By[i];  int ioff=i*npoly; for(int j=0;j<ipw;j++){ BBcp[ioff+j] = BB[ioff+j]; }; }  // must copy otherwise linSolve_gauss destroy it

            //for(int i=0;i<ipw;i++){
            //    for(int j=0;j<ipw;j++){ printf("%3.3f ", BB_[i][j] ); };
            //    //printf(" BB_ Bycp %g \n", Bycp[i] );
            //}
            Lingebra::linSolve_gauss ( ipw, BB_, Bycp, index, coefs );
            //for(int i=0;i<ipw;i++){ printf("coefs[%i] %g \n", i, coefs[i] ); }
            polyeval  ( n, ipw, ipolys, xs, ys_test+i0, coefs );
            checkError( pows[ipow], i0, n );

            printf( "multiFitAndCheck[%i] errmax %g rmse %g \n", ipw, errmax, sqrt(err2/n) );
            if(errmax<errTol){ m_conv=ipw; break; }
            if(bAllCoefs) coefs+=ipw;

        }
        return m_conv;
    }

};




/*
class AutoApprox{ public:

    int np=0;             // number of data points
    double* xs       =0;  // x-argument values
    double* ys_ref   =0;  // reference function values
    double* ys_trans =0;  // transformed function
    double* ytest    =0;

    int     ncoefMin=0;   // starting number of coefs
    int     ncoef=0;      // number for fit coefficients (degree of polynominal, or number of orther basis function)
    double* coefs=0;      //

    double* errMaxs = 0;
    double* errRMSs = 0;

    int nint=0; //  number of intervals
    int* i0s=0; //
    int* i1s=0; //

    void realloc(int np, ){
        double* yis   = new double[n];
        double* ytest = new double[n];
    }

    int exhaustive( ){
        int    ipow_best = 0;
        double errBest=1e+300;
        for(int ipow=0; ipow<npows; ipow++ ){ // loop over y-powers
            // --- function transform
            double alpha    = ypows[ipow];
            double invAlpha = 1/alpha;
            for(int i=0; i<np; i++){
                yis[i]=pow(ys[i], invAlpha );
                //printf( "[%i]y^%g [%i] %g | %g %g \n", ipow, alpha, i, yis[i], ys[i], xs[i] );
            };
            // --- poly fit
            for(int ii=0; ii<nint; ii++){ // loop over intervals
                int i0        = i0s[ii];
                int ni        = i1s[ii]-i0s[ii];
                for(int ii=0; ii<nint; ii++){ // loop over polynom-order
                    polyFit ( n, m, xs, yis,   coefs );
                    polyeval( n, m, xs, ytest, coefs );
                    // --- error evaluation
                    double err2=0,errmax=0;
                    for(int i=0; i<n; i++){
                        //double derr =  yis[i]-[i];
                        double derr = pow( ytest[i], alpha ) - ys[i];
                        if(relative) derr/=ys[i];
                        errmax  = fmax( errmax, fabs(derr) );
                        err2   += derr*derr;
                    }
                    //printf("ypowsApprox[%i] %g %g \n", ipow, errmax, err2 );
                    if(errmax<errBest){ errBest=errmax; ipow_best=ipow; };
                    errs[ipow*2+0]=errmax;
                    errs[ipow*2+1]=err2;
                }
            }
        }
        return ipow_best;
    }

    void fitOnIntervals( int npows, int n, int m, double* xs, double* ys, double* coefs, double* ypows, double* errs, bool relative, int nint, int* i0s, int* i1s ){
        for(int ii=0; ii<nint; ii++){
            int i0        = i0s[ii];
            int ni        = i1s[ii]-i0s[ii];
            int ipow_best = ypowsApprox( npows, ni, m, xs+i0, ys+i0, coefs, ypows, errs, relative );
        }
    }

};
*/

}; // namespace Approx

#endif

