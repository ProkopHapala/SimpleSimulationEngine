
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
//       x3   x2  x   1
//  ----------------------
// c0   1    0    0   0     /6.0
// c1  -3    3    3   1
// c2  +3   -6    0   4
// c3  -1    3    3   1


// ============ optimized

namespace CubicBSpline{

const static double B[4][4] = {
    {  1,  0,  0,  0 },
    {  0,  1,  0,  0 },
    { -3, -2,  3, -1 },
    {  2,  1, -2,  1 }
};


    inline double val( double u, double c0, double c1, double c2, double c3 ){
      return    (         c2 + 4*c1  +   c0
            + u*(      3*(c2         -   c0)
            + u*(      3*(c2 - 2*c1  +   c0)
            + u*(  c3 -3*(c2 -  c1) -   c0   )))) * 0.16666666666d;
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



