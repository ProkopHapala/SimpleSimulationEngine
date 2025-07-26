
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

struct SDF_Cylinder {
    Vec3d p0;
    Vec3d hdir;
    double l;
    double r;
    bool capped;

    SDF_Cylinder( Vec3d p0_, Vec3d p1, double r_, bool capped_ ){
        from_ends(p0_,p1); r = r_; capped = capped_; 
        //printf("SDF_Cylinder: p0(%8.4f,%8.4f,%8.4f) p1(%8.4f,%8.4f,%8.4f) hdir(%8.4f,%8.4f,%8.4f) l: %8.4f r: %8.4f \n", p0.x,p0.y,p0.z, p1.x,p1.y,p1.z, hdir.x,hdir.y,hdir.z, l, r );
    }
    void from_ends( Vec3d p0_, Vec3d p1 ){ p0 = p0_; hdir = p1-p0_; l = hdir.normalize(); }
    double operator()(const Vec3d& p) const {
        Vec3d q = p - p0;
        double h = q.dot(hdir);
        double dist;
        if     (h<0.0){ if(capped){ dist = (p- p0        ).norm() - r; } else dist = -h ; }
        else if(h>l  ){ if(capped){ dist = (p-(p0+hdir*h)).norm() - r; } else dist = h-l; }
        else          { 
            Vec3d proj = p0 + hdir * h;
            dist = (p - proj).norm() - r;
        }
        printf("SDF_Cylinder: p(%8.4f,%8.4f,%8.4f) h: %8.4f dist: %8.4f \n", p.x,p.y,p.z, h, dist );
        return dist;
    }
};


#endif


