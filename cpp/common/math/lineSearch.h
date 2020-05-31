#ifndef  lineSearch_h
#define  lineSearch_h

/*

Application : RayMarching  (and Curve Marching):
 - finish Auto-extend (non-bound intersection)
 - utilize distance-function estimate ( like in RayMarching )

*/

#include "fastmath.h"

class LineSearch{ public:
    double xmin,xmax;
    double ymin,ymax;
    int side=0;

    bool   bound=false;
    double xExtendMax =  0.3;
    double xErrConv   = 1e-6;

    //Func1d func;

    inline bool checkConv(){
        return fabs(xmax-xmin)<xErrConv;
    }

    inline double guess_RegulaFalsa(){
        return (ymin*xmax - ymax*xmin)/(ymin-ymax);
    }

    inline double step_RegulaFalsa( double x, double y ){
        //double x =
        //if (fabs(xmax-xmin) < errMax*fabs(xmax+xmin)) break; // Convergence reached !
        //double y = func(x);
        if (y * ymax > 0){
            xmax = x;
            ymax = y;
            if (side==-1) ymin /= 2; //   Illinois Algorithm  correction
            side = -1;
        }else if (y * ymin > 0){
            xmin = x;
            ymin = y;
            if (side==+1) ymax /= 2; //   Illinois Algorithm  correction
            side = +1;
        }//else{ break; } // fr * f_ very small (looks like zero)
        return guess_RegulaFalsa();
    }

    double linear_extend( double x, double y,  double x0, double y0 ){
        double dxdy = (x-x0)/(y-y0);
        return x - clamp( dxdy*y, -xExtendMax, xExtendMax );
    }

    double init(double xmin_, double ymin_, double xmax_, double ymax_){
        bound = (ymin_*ymax_<0);
        if(ymin_>ymax_){
            xmin=xmax_;  ymin=ymax_;
            xmax=xmin_;  ymax=ymin_;
        }else{
            xmin=xmin_;  ymin=ymin_;
            xmax=xmax_;  ymax=ymax_;
        }
        //return step( );
        return 0;
    }

    double step( double x, double y ){
        // Automatically choose interpolation vs. extrapolation step (bound and un-bound solution)
        // ToDo : consider ray-arching with distance estimator
        if(bound){
            return step_RegulaFalsa( x, y );
        }else{
            return linear_extend(x,y,xmin,ymin);
        }
    }

    double solve( Func1d func, double xmin, double xmax, int nMaxIter=1000, double xErrConv=1e-6 ){
        double ymin = func(xmin);
        double ymax = func(xmax);
        init( xmin, ymin, xmax, ymax );
        // ToDo : extend step later !!!
        double x =  guess_RegulaFalsa();
        for(int i=0; i<nMaxIter; i++){
            double y = func(x);
            x = step_RegulaFalsa( x, y );
            if( checkConv() ) break;
        }
        return x;
    }

};

template <typename Func>
double lineSearch( LineSearch& lsearch, Func func, double xmin, double xmax, int nMaxIter=1000, double xErrConv=1e-6 ){
    double ymin = func(xmin);
    double ymax = func(xmax);
    lsearch.init( xmin, ymin, xmax, ymax );
    // ToDo : extend step later !!!
    double x =  lsearch.guess_RegulaFalsa();
    for(int i=0; i<nMaxIter; i++){
        double y = func(x);
        x = lsearch.step_RegulaFalsa( x, y );
        //printf( "lineSearch([%i] x %g | <%g..%g> |span| %g \n", i, x, lsearch.xmin,lsearch.xmax, lsearch.xmax-lsearch.xmin );
        if( lsearch.checkConv() ) break;
    }
    return x;
}


/*

from here:
http://www1.fpl.fs.fed.us/Fmin.java

Here is a copy of the Netlib documentation:


      An approximation x to the point where f attains a minimum on
  the interval (ax,bx) is determined.

  input..

  ax    left endpoint of initial interval
  bx    right endpoint of initial interval
  f     function subprogram which evaluates f(x) for any x
        in the interval (ax,bx)
  tol   desired length of the interval of uncertainty of the final
        result (.ge.0.)

  output..

  fmin  abcissa approximating the point where  f  attains a
        minimum

     The method used is a combination of golden section search and
  successive parabolic interpolation.  Convergence is never much slower
  than that for a Fibonacci search.  If f has a continuous second
  derivative which is positive at the minimum (which is not at ax or
  bx), then convergence is superlinear, and usually of the order of
  about 1.324.
     The function f is never evaluated at two points closer together
  than eps*abs(fmin)+(tol/3), where eps is approximately the square
  root of the relative machine precision.  If f is a unimodal
  function and the computed values of f are always unimodal when
  separated by at least eps*abs(x)+(tol/3), then fmin approximates
  the abcissa of the global minimum of f on the interval (ax,bx) with
  an error less than 3*eps*abs(fmin)+tol.  If f is not unimodal,
  then fmin may approximate a local, but perhaps non-global, minimum to
  the same accuracy.

*/

// template see here :  Ben Supnik  http://stackoverflow.com/questions/1174169/function-passed-as-template-argument

//typedef double( *scalar_func )(double);
//template<scalar_func pot_func>

template<class ScalarFunction1D>
double lineSearch_Brent (ScalarFunction1D pot_func, double x0, double a, double b, double tol ) {
    double  x,  u, v, w;
    double fx, fu,fv,fw;
    double p,q,r;
    double tol1,t2;

    const double golden = 0.5d*( 3.0d - sqrt(5.0d) );
    double d = 0.0d;

    double eps = 1.2e-16;
    tol1 = eps + 1.0d;
    eps = sqrt(eps);

    //v = a + golden*(b-a);
    v = x0;
    w = v;
    x = v;
    double e = 0.0d;
    fx = pot_func(x); //ellipse( sX(x),  sY(fx), 2,2 );
    fv = fx;
    fw = fx;
    const double tol3 = tol/3.0d;

    double xm = 0.5d*(a + b);
    tol1 = eps* fabs(x) + tol3;
    t2 = 2.0d*tol1;

// main loop

    while ( fabs(x-xm) > (t2 - 0.5d*(b-a))) {
       p = q = r = 0.0d;
       if ( fabs(e) > tol1) {   // fit the parabola
          r = (x-w)*(fx-fv);
          q = (x-v)*(fx-fw);
          p = (x-v)*q - (x-w)*r;
          q = 2.0d*(q-r);
          if (q > 0.0d)  p = -p;
          else          q = -q;
          r = e;
          e = d;
       }

       if (( fabs(p) < fabs(0.5d*q*r)) && (p > q*(a-x)) && (p < q*(b-x))) {  // a parabolic interpolation step
          d = p/q;
          u = x+d;
          if (((u-a) < t2) || ((b-u) < t2)) {   // f must not be evaluated too close to a or b
             d = tol1;             // offset from left
             if (x >= xm) d = -d;  // offset from right
          }
          //print( " parabolic "  );
       } else {   // a golden-section step
          if (x < xm) e = b-x;
          else        e = a-x;
          d = golden*e;
          //print( " golden    "  );
       }

       if ( fabs(d) >= tol1) {    // f must not be evaluated too close to x
          u = x+d;
       } else {
          if (d > 0.0d) u = x + tol1;
          else         u = x - tol1;
       }

       fu = pot_func(u);         //ellipse( sX(u),  sY(fu), 2,2 );

       if (fx <= fu) {
          if (u < x) a = u;
          else       b = u;
       }
       if (fu <= fx) {
          if( fabs(u-x)<tol) break;
          if (u < x)  b = x;
          else        a = x;
          v = w; fv = fw;
          w = x; fw = fx;
          x = u; fx = fu;
       } else {
          if ( (fu <= fw) || (w == x) ) {
             v = w; fv = fw;
             w = u; fw = fu;
          } else {
            if (!((fu > fv) && (v != x) && (v != w))) {
               v = u;  fv = fu;
            }
          }
       }
       xm   = 0.5d*(a + b);
       tol1 = eps * fabs(x) + tol3;
       t2   = 2.0d*tol1;

       //println( " t2 "+t2 +"  tol1   "+ tol1 + "    x   "+ x );
    }
    return x;
 };

 #endif
