
// according to :
//   C1- and C2-continuous spline-interpolation of a regulartriangular net of points
//   Albert Wiltsche,Computers & Graphics 27 (2003) 917-930
//   https://www.sciencedirect.com/science/article/pii/S0097849303001705

// see also python : /Dropbox/MyDevSW/Python/_math/triangular_spline/C1tri.py

// tau(u,v,w) :=  v*v * w*w
// pi(u,v,w)  :=  u**4   + 4*u**3*(v+w)  +  u**2*( 2*(v*v+w*w) + 8*v*w ) + 2*u*(v**2*w + v*w**2) + v**2*w**2
//  A*fA + B*fB + C*fC   + D*fD + E*fE + F*fF

//====================
//--- upper trinalge
//            //
//       E    //
//            //
//     A-C    //
//     |/     //
//   F B D    //
//            //
//--- Lower triangle
//            //
//     D B F  //
//      /|    //
//     C-A    //
//            //
//     E      //
//            //
//===================

//fA = pi(w,u,v)
//fB = pi(u,v,w)
//fC = pi(v,w,u)
//fD = tau(w,u,v)
//fE = tau(u,v,w)
//fF = tau(v,w,u)


// matrix form
//    uuuu uuuv uuuw   uuvv uuvw uuww   uvvv uvvw uvww    vvvv vvvw vvww   vwww wwww
// 15*6 matrix ? is it worh it ?

#ifndef  spline_triC1_h
#define  spline_triC1_h

//#include <math.h>
//#include <cstdlib>
//#include <stdio.h>

// ============ optimized

namespace Spline_triC1{

template <class T>
inline T pi(T u, T v, T w){
    //return  u*u*u*u   + 4*u*u*u*(v+w)  +  u*u*( 2*(v*v+w*w) + 8*v*w ) + 2*u*(v*v*w + v*w*w) + v*v*w*w;
    //return ( ( u        +       4*(v+w) )*u  +      ( 2*(v*v+w*w) + 8*v*w ) )*u + 2*(v*v*w + v*w*w) )*u + v*v*w*w;
    //return v*v*w*w +  2*u*(v*v*w + v*w*w)   +  u*u*(   2*(v*v+w*w) + 8*v*w ) + 4*u*u*u*(v+w)   u*u*u*u;
    //return v*v*w*w +  u*( 2*(v*v*w + v*w*w)   +  u*( ( 2*(v*v+w*w) + 8*v*w ) +   u*( 4*(v+w) + u ) ) );
    //T v2 = v*v; T w2 = w*w;
    //return v2*w2 +  u*( 2*(v2*w + v*w2)   +  u*( ( 2*(v2+w2) + 8*v*w ) +   u*( 4*(v+w) + u ) ) );
    T vw = v*w; T vAw = v+w;
    //return vw*vw +  u*( 2*( vw*(v+w) +  u*( ( 2*(v*v+w*w) + 8*vw ) +   u*( 4*(v+w) + u ) ) );
    //return vw*vw +  u*( 2*( vw*vAw )+  u*( ( 2*vAw*vAw + 4*vw ) +   u*( 4*vAw + u ) ) );
    return vw*vw +    u*(  2*vw*vAw +  u*( 2*( vAw*vAw + vw+vw ) +   u*( 4*vAw + u ) ) );
}

template <class T>
inline T tau(T u, T v, T w){
    return v*v*w*w;
}

template <class T>
inline T pitau(T u, T v, T w, T A, T D ){
    T vw = v*w; T vAw = v+w;
    return vw*vw*(A+D) +  A*u*(  2*vw*vAw +  u*( 2*( vAw*vAw + vw+vw ) +   u*( 4*vAw + u ) ) );
}

template <class T>
inline T tau_(T u, T v, T w){
    return v*v*w*w;
}

template <class T>
inline T val( T u,T v,T w,   T A, T B, T C,   T D, T E, T F ){
    //return A*pi (w,u,v) + B*pi (u,v,w) + C*pi (v,w,u)
    //+      D*tau(w,u,v) + E*tau(u,v,w) + F*tau(v,w,u);
    return pitau(w,u,v, A,D) + pitau(u,v,w, B,E) + pitau(v,w,u, C,F);
}

};

#endif



