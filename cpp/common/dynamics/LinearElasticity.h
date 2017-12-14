
#ifndef LinearElasticity_h
#define LinearElasticity_h

#include "Vec3.h"
#include "Mat3.h"

// see
// /home/prokop/git/SimpleSimulationEngine/cpp/common/dynamics/Notes
// https://en.wikipedia.org/wiki/Flexibility_method
// https://en.wikipedia.org/wiki/Direct_stiffness_method

// A x_i =  g_i
// g_i = Sum s_i (x_i - x_j) (dk)^2


void build_stiffness_matrix_3D( int n, int m,   int * links, double * stiffness, Vec3d * points, double * M ){
    int n3 = n*3;
    for(int k=0; k<m; k++){
        int i = links[0];
        int j = links[1];
        links+=2;

        Vec3d d   = points[j] - points[i];
        //double l2 = d.norm2();
        Mat3d D;
        D.set_outer(d,d);
        double s  = stiffness[k];
        D.mul(s/d.norm2());

        int i3 = i*3;
        int j3 = j*3;

        // diagonal

        double* M_ = M + n3*i3 + i3;
        M_[0] -= D.xx;
        M_[1] -= D.xy;
        M_[2] -= D.xz;
        M_    += n3;
        M_[0] -= D.yx;
        M_[1] -= D.yy;
        M_[2] -= D.yz;
        M_    += n3;
        M_[0] -= D.zx;
        M_[1] -= D.zy;
        M_[2] -= D.zz;

        M_ = M + n3*j3 + j3;
        M_[0] -= D.xx;
        M_[1] -= D.yx;
        M_[2] -= D.zx;
        M_    += n3;
        M_[0] -= D.xy;
        M_[1] -= D.yy;
        M_[2] -= D.zy;
        M_    += n3;
        M_[0] -= D.xz;
        M_[1] -= D.yz;
        M_[2] -= D.zz;

        // cross coupling

        M_ = M + n3*i3 + j3;
        M_[0] += D.xx;
        M_[1] += D.xy;
        M_[2] += D.xz;
        M_    += n3;
        M_[0] += D.yx;
        M_[1] += D.yy;
        M_[2] += D.yz;
        M_    += n3;
        M_[0] += D.zx;
        M_[1] += D.zy;
        M_[2] += D.zz;

        M_ = M + n3*j3 + i3;
        M_[0] += D.xx;
        M_[1] += D.yx;
        M_[2] += D.zx;
        M_    += n3;
        M_[0] += D.xy;
        M_[1] += D.yy;
        M_[2] += D.zy;
        M_    += n3;
        M_[0] += D.xz;
        M_[1] += D.yz;
        M_[2] += D.zz;

    }
}


#endif
