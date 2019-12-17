
#ifndef  Besier_h
#define  Besier_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include <fastmath.h>
#include <Vec3.h>

//========================
//   Cubic B-spline
//========================
// http://www2.cs.uregina.ca/~anima/408/Notes/Interpolation/UniformBSpline.htm
// http://alex.vlachos.com/graphics/CurvedPNTriangles.pdf
// // see http://www.gamasutra.com/view/feature/131389/b%C3%A9zier_triangles_and_npatches.php?print=1
//
//       x3   x2   x   1
//  ----------------------
//  p0  -1   +3   -3    1  /6
//  p1  +3   -6         4  /6
//  p2  -3   +3    +3   1  /6
//  p3  +1                 /6
//
//  dval
//       x2   x   1
//  -----------------
//  p0  -3   +6   -3    /0
//  p1  +9   -12        /6
//  p2  -9   +6   +3    /6
//  p3  +3              /6
//       x2   x   1
//  -----------------
//  p0  -1   +2   -1    /2
//  p1  +3   -4         /2
//  p2  -3   +2   +1    /2
//  p3  +1              /2
//
//  ddval
//       x   1
//  ---------------
//  p0  -1   +1
//  p1  +3   -2
//  p2  -3   +1
//  p3  +1


// ============ optimized

namespace Besier{

    inline void triangleBasisO3( double a, double b, double * bas ){
        double c  = 1-a-b;
        double a2 = a*a;
        bas[0]    = a2*a;      // a3
        bas[1]    = a2*b*3;    // a2b
        bas[2]    = a2*c*3;    // ca2
        double b2 = b*b;
        bas[3]    = b2*b;      // b3
        bas[4]    = b2*a*3;    // ab2
        bas[5]    = b2*c*3;    // b2c
        double c2 = c*c;
        bas[6]    = c2*c;      // c3
        bas[7]    = c2*a*3;    // c2a
        bas[8]    = c2*b*3;    // bc2
        bas[9]    = 6*a*b*c;   // abc
    }

    inline void triangleBasisO2( double a, double b, double * bas ){
        double c  = 1-a-b;
        bas[0]    = a*a;    // a3
        bas[1]    = b*b;    // b3
        bas[2]    = c*c;    // c3
        bas[3]    = a*b*2;  // a2b
        bas[4]    = a*c*2;  // c2a
        bas[5]    = b*c*2;  // ab2
    }

};

const int BesierTriangle_borderHull[9] = {0,1,4,3,5,8,6,7,2};
const int BesierTriangle_innerHull [9] = {1,9,8,7,9,4,5,9,2};

class BesierTriangle{
    public:
    Vec3d ps[10];
    Vec3d ns[6];

    // see http://www.gamasutra.com/view/feature/131389/b%C3%A9zier_triangles_and_npatches.php?print=1
    void fromPNs(
        const Vec3d& A,  const Vec3d& B,  const Vec3d& C,
        const Vec3d& nA, const Vec3d& nB, const Vec3d& nC
    ){
        constexpr double c16 = 1.0d/6.0d;
        constexpr double c13 = 1.0d/3.0d;
        Vec3d p;
        ps[0].set(A);             p.set_lincomb(-c13, A, c13, B ); p.add_mul( nA, -p.dot(nA) );
        ps[1].set_add( A, p);     p.set_lincomb(-c13, A, c13, C ); p.add_mul( nA, -p.dot(nA) );
        ps[2].set_add( A, p);

        ps[3].set(B);             p.set_lincomb(-c13, B, c13, A ); p.add_mul( nB, -p.dot(nB) );
        ps[4].set_add( B, p);     p.set_lincomb(-c13, B, c13, C ); p.add_mul( nB, -p.dot(nB) );
        ps[5].set_add( B, p);

        ps[6].set(C);             p.set_lincomb(-c13, C, c13, A ); p.add_mul( nC, -p.dot(nC) );
        ps[7].set_add( C, p);     p.set_lincomb(-c13, C, c13, B ); p.add_mul( nC, -p.dot(nC) );
        ps[8].set_add( C, p);

        ps[9].set_lincomb(0.25d, ps[8]+ps[7]+ps[5]+ps[4]+ps[2]+ps[1], c16, A+B+C );  // 6/4- 3/6 = 1
        //ps[9].set_lincomb(0.5d, ps[8]+ps[7]+ps[5]+ps[4]+ps[2]+ps[1], -c13, A+B+C );  // 6/2- 3/1 = 1
        //ps[9].set_lincomb(1.0, ps[8]+ps[7]+ps[5]+ps[4]+ps[2]+ps[1], -5.0/3.0, A+B+C );    // 6 - 3 = 1
        //ps[9].set_lincomb(5.0/3.0, ps[8]+ps[7]+ps[5]+ps[4]+ps[2]+ps[1], -3.0, A+B+C );
        //ps[9].set_lincomb(0.25d, ps[8]+ps[7]+ps[5]+ps[4]+ps[2]+ps[1], 1, A+B+C );    // this makes it smooth - probably some bug

        ns[0].set(nA);
        ns[1].set(nB);
        ns[2].set(nC);

        Vec3d n;
        constexpr double cn = -1.5d;
        //constexpr double cn = -2.0d;
        n.set_add(nA,nB);   p.set_sub(A,B);  n.add_mul(p, cn*p.dot(n) / p.norm2() );  n.normalize();
        ns[3].set(n);
        n.set_add(nA,nC);   p.set_sub(A,C);  n.add_mul(p, cn*p.dot(n) / p.norm2() );  n.normalize();
        ns[4].set(n);
        n.set_add(nB,nC);   p.set_sub(B,C);  n.add_mul(p, cn*p.dot(n) / p.norm2() );  n.normalize();
        ns[5].set(n);



    }

    Vec3d getPoint( double a, double b ) const {
        double bas[10];
        Besier::triangleBasisO3( a, b, bas );
        Vec3d p; p.set(0.0d);
        for(int i=0; i<10; i++){
            //printf( " %i (%3.3f) (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) \n",  i, bas[i], ps[i].x, ps[i].y, ps[i].z );
            p.add_mul(ps[i],bas[i]);
        }
        return p;
    }

    Vec3d getNormal( double a, double b ) const {
        double bas[6];
        Besier::triangleBasisO2( a, b, bas );
        Vec3d p; p.set(0.0d);
        for(int i=0; i<6; i++){
            //printf( " %i (%3.3f) (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) \n",  i, bas[i], ps[i].x, ps[i].y, ps[i].z );
            p.add_mul(ns[i],bas[i]);
        }
        return p;
    }

};

#endif



