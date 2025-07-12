
#ifndef  SDfuncs_h
#define  SDfuncs_h

#include "Vec3.h"

struct SDF_Sphere {
    Vec3d center;
    double radius;
    double operator()(const Vec3d& p) const { return (p - center).norm() - radius; }
};

struct SDF_AABB {
    Vec3d center;
    Vec3d halfSpan;

    double operator()(const Vec3d& p) const {
        Vec3d  q = (p-center).abs() - halfSpan;
        double l = q.maxComponent(); if(l< 0.0) return l;
        return q.max(Vec3dZero).norm() + l;
    }
};

#endif


