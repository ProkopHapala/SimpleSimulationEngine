
#ifndef  CubicBSpline_h
#define  CubicBSpline_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"

//========================
//   Hermite splines
//========================
// https://math.stackexchange.com/questions/699113/relation-of-cubic-b-splines-with-cubic-splines
//      x   1  x  x2  x3
//  ----------------------
// c0   1   -3    3  -1     /6.0
// c1   4    0   -6   3     /6.0
// c2   1   -6    3  -3     /6.0
// c3   0    0    0   1     /6.0


// ============ optimized

namespace CubicBSpline{

const static double B[4][4] = {
    {  0.16666666666d, -0.5d, 0.5d, -0.16666666666d },
    {  0.66666666666d, 0.0d, -1.0d,  0.5d           },
    {  0.16666666666d, 0.5d,  0.5d, -0.5d           },
    {  0.0d,           0.0d,  0.0d,  0.16666666666d }
};

    inline double val( double u, double c0, double c1, double c2, double c3 ){
      return    (         c2 + 4*c1  +   c0
            + u*(      3*(c2         -   c0)
            + u*(      3*(c2 - 2*c1  +   c0)
            + u*(  c3 -3*(c2 -  c1)  -   c0   )))) * 0.16666666666d;
    /*
    return      (  B[0][0]*c0 + B[1][0]*c1 + B[2][0]*c2 + B[3][0]*c3
            + u*(  B[0][1]*c0 + B[1][1]*c1 + B[2][1]*c2 + B[3][1]*c3
            + u*(  B[0][2]*c0 + B[1][2]*c1 + B[2][2]*c2 + B[3][2]*c3
            + u*(  B[0][3]*c0 + B[1][3]*c1 + B[2][3]*c2 + B[3][3]*c3     ))));
    */
    }

    inline double dval( double u, double c0, double c1, double c2, double c3 ){
      return    (            c2         -  c0
            + u*(         2*(c2 - 2*c1  +  c0)
            + u*(    c3 - 3*(c2 -   c1) -  c0   ))) * 0.5;
    }

    inline double ddval( double u, double c0, double c1, double c2, double c3 ){
      return    (         c2 - 2*c1  +  c0
            + u*(  c3 -3*(c2 -   c1) -  c0  ));
    }

    inline void basis( double u, double& b0, double& b1, double& b2, double& b3 ){
        /*
        double i6  = 0.16666666666d;
        double u2  = u*u;
        double u3  = u*u2;
        b0  = i6*( 1 - 3*u + 3*u2 -   u3 );
        b1  = i6*( 4       - 6*u2 + 3*u3 );
        b2  = i6*( 1 + 3*u + 3*u2 - 3*u3 );
        b3  = i6*(                    u3 );
        */
        double u2  = u*u;
        double u3  = u*u2;
        b0  =  0.16666666666d*(1-u3) - 0.5d*(u-u2);
        b1  =  0.66666666666d        - u2 + 0.5d*u3;
        b2  =  0.16666666666d        + 0.5d*(u+u2-u3);
        b3  =  0.16666666666d*u3;
    }

    inline void dbasis( double u, double& b0, double& b1, double& b2, double& b3 ){
        /*
        double u2  = u*u*0.5d;
        b0  =  -0.5d +   u -   u2;
        b1  =        - 2*u + 3*u2;
        b2  =   0.5d +   u - 3*u2;
        b3  =                  u2;
        */
        double u2 = u*u*0.5d;
        u  -=  u2;
        b0  =  -0.5d +   u;
        b1  =        - 2*u +   u2;
        b2  =   0.5d +   u - 2*u2;
        b3  =                  u2;
    }

    inline void ddbasis( double u, double& b0, double& b1, double& b2, double& b3 ){
        b0 =   1 -   u;
        b1 =  -2 + 3*u;
        b2 =   1 - 3*u;
        b3 =         u;
    }

    /*
    inline void basis( double u, double& b0, double& b1, double& b2, double& b3 ){
        double i6  = 0.16666666666d;
        double u2  = u*u;
        double u3  = u*u2;
        b0  = i6*( 1 - 3*u + 3*u2 -   u3 );
        b1  = i6*( 4       - 6*u2 + 3*u3 );
        b2  = i6*( 1 + 3*u + 3*u2 - 3*u3 );
        b3  = i6*(                    u3 );
    }

    inline void dbasis( double u, double& b0, double& b1, double& b2, double& b3 ){
        double i6  = 0.16666666666d;
        double u2  = u*u;
        b0  = i6*( -3 + 6*u - 3*u2 );
        b1  = i6*(    -12*u + 9*u2 );
        b2  = i6*(  3 + 6*u - 9*u2 );
        b3  = i6*(           3*u2 );
    }

    inline void ddbasis( double u, double& b0, double& b1, double& b2, double& b3 ){
        double i6  = 0.16666666666d;
        b0 = i6*(  6 -  6*u);
        b1 = i6*(-12 + 18*u);
        b2 = i6*(  6 - 18*u);
        b3 = i6*(       6*u);
    }
    */


/*
    // optimized version
    // coefficients in reverse order - for compactibility
    inline double val( double u, double c0, double c1, double c2, double c3 ){
      return    (         c1 + 4*c2  +   c3
            + u*(      3*(c1         -   c3)
            + u*(      3*(c1 - 2*c2  +   c3)
            + u*(  c0 -3*(c1 -   c2) -   c3   )))) * 0.16666666666d;
    }

    inline double dval( double u, double c0, double c1, double c2, double c3 ){
      return    (            c1         -  c3
            + u*(         2*(c1 - 2*c2  +  c3)
            + u*(    c0 - 3*(c1 -   c2) -  c3   ))) * 0.5;
    }

    inline double ddval( double u, double c0, double c1, double c2, double c3 ){
      return    (         c1 -  2*c2 +  c3
            + u*(  c0 -3*(c1 -   c2) -  c3  ));
    }
*/


/*
// From Java boundary conditions (val, deriv )

void boundaryLeft( double [] Cs, double y, double dy ){
  final double  B0 = 2.0d/3.0d;
  final double  B1 = 1.0d/6.0d;
  double c2 = Cs[2];
  double c0 = -( dy*2 - c2 );
  //double c1 =  ( y -  B1 * (c0+c2) )/B0;
  double c1 =  ( y -  2*B1*( c2 - dy ) )/B0;
  Cs[1] = c1; Cs[0] = c0;
};

void boundaryRight( double [] Cs, double y, double dy ){
  final double  B0 = 2.0d/3.0d;
  final double  B1 = 1.0d/6.0d;
  double c2 = Cs[Cs.length - 3];
  double c0 = ( dy*2 + c2 );
  //double c1 = ( y -  B1 * (c0+c2) )/B0;
  double c1 = ( y -  2*B1*( c2 + dy ) )/B0;
  Cs[Cs.length-2] = c1; Cs[Cs.length-1] = c0;
};

*/

};

#endif



