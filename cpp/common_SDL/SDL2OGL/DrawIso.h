
#ifndef  DrawIso_h
#define  DrawIso_h

#include "Draw3D.h"

namespace Draw3D{

//template< double (*func)(const Vec3d& p)>
template< typename Func >
void drawIso( Func func, Vec3i n, Vec3d pmin, Mat3d dcell, double iso ){
    const int np = 4+6+1;
    Vec3f  ps  [np];
    double vals[np];
    const Vec3d p0s[np]{
        0.0,0.0,0.0,    // 0
        1.0,0.0,0.0,    // 1
        0.0,1.0,0.0,    // 2
        0.0,0.0,1.0,    // 3
        +0.5,+0.5,+0.5, // 4
        +0.5,-0.5,-0.5, // 5
        -0.5,+0.5,-0.5, // 6
        -0.5,-0.5,+0.5, // 7
        -0.5,+0.5,+0.5, // 8
        +0.5,-0.5,+0.5, // 9
        +0.5,+0.5,-0.5, // 10
    };
    const int nt = 12;
    const Quat4i tis[nt]{
        0,1, 4,9,
        0,1, 9,5,
        0,1, 5,10,
        0,1, 10,4,

        0,2, 4,10,
        0,2, 10,6,
        0,2, 6,8,
        0,2, 8,4,

        0,3, 4,8,
        0,3, 8,7,
        0,3, 7,9,
        0,3, 9,4
    };
    //glBegin(GL_LINES);
    glBegin(GL_TRIANGLES);
    for(int ix=0; ix<n.x; ix++){
        for(int iy=0; iy<n.y; iy++){
            for(int iz=0; iz<n.z; iz++){
                Vec3d p =  pmin + dcell.a*ix + dcell.b*iy + dcell.c*iz ;
                for(int i=0; i<np; i++){
                    Vec3d pi = p+dcell.dot( p0s[i] );
                    ps   [i] = (Vec3f)pi;
                    vals [i] = func(pi) - iso;
                    //printf( "(%i,%i,%i|%i) %g \n", ix,iy,iz, i, vals[i] );
                }
                Vec3f* pps[4];
                for(int i=0; i<nt; i++){
                    Quat4i ti = tis[i];
                    pps[0]=ps+ti.x; pps[1]=ps+ti.y; pps[2]=ps+ti.z; pps[3]=ps+ti.w;
                    Draw3D::drawTetraIso( pps, (Quat4d){vals[ti.x],vals[ti.y],vals[ti.z],vals[ti.w]} );

                }
            }
        }
    }
    glEnd();
}

template< typename Func >
void drawIso( Func func, double h, Vec3d pmin, Vec3d pmax, double  iso ){
    Vec3d dp = pmax-pmin;
    double invh = 1/h;
    Vec3i n = { round(dp.x*invh), round(dp.y*invh), round(dp.z*invh) };
    Mat3d dcell = Mat3dZero;
    dcell.xx = dp.x/n.x;
    dcell.yy = dp.y/n.x;
    dcell.zz = dp.z/n.x;
    drawIso( func, n, pmin, dcell, iso );
}

}; // namespace Draw3D

#endif

