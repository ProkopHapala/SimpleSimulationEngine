#ifndef  BucketUtils_h
#define  BucketUtils_h

/// @file BucketUtils.h
/// @brief Ray-picking helpers that use Buckets AABBs to pre-filter candidates — works with any bucketed point/edge arrays.

#include "Buckets.h"
#include "raytrace.h"

/// @brief Pick nearest point to ray, using bucket AABBs to skip non-intersecting buckets.
/// @param buckets  CSR layout mapping cells → point indices
/// @param BBs      per-cell bounding boxes (AABB), indexed same as buckets.ncell
/// @param points   point array (must support .f returning Vec3d)
/// @return index of nearest point within Rmax, or -1
template<typename PointT>
inline int pick_point_bucket( const Buckets& buckets, const Quat8d* BBs, const PointT* points, const Vec3d& ray0, const Vec3d& hray, double Rmax ){
    if( buckets.ncell <= 0 || !BBs ) return -1;
    double r2min = Rmax*Rmax;
    int    imin  = -1;
    Vec3d  tmpPos, tmpNrm;
    for(int ib=0; ib<buckets.ncell; ib++){
        int n = buckets.cellNs[ib];
        if(n<=0) continue;
        const Quat8d& bb = BBs[ib];
        double t = rayBox( ray0, hray, (Vec3d)bb.lo.f, (Vec3d)bb.hi.f, tmpPos, tmpNrm );
        if( t >= 1e30 ) continue;
        int i0 = buckets.cellI0s[ib];
        for(int j=0; j<n; j++){
            int i = buckets.cell2obj[i0+j];
            double tt;
            double r2 = rayPointDistance2( ray0, hray, points[i].f, tt );
            if(r2<r2min){ imin=i; r2min=r2; }
        }
    }
    return imin;
}

/// @brief Pick nearest bond (edge) to ray, using bucket AABBs to skip non-intersecting buckets.
/// @param buckets  CSR layout mapping cells → edge indices
/// @param BBs      per-cell bounding boxes (AABB), indexed same as buckets.ncell
/// @param bonds    edge index pairs (int2)
/// @param points   point array (must support .f returning Vec3d)
/// @return index of nearest bond within Rmax, or -1
template<typename PointT>
inline int pick_bond_bucket( const Buckets& buckets, const Quat8d* BBs, const int2* bonds, const PointT* points, const Vec3d& ray0, const Vec3d& hRay, double Rmax ){
    if( buckets.ncell <= 0 || !BBs ) return -1;
    double dist_min = Rmax;
    int    imin     = -1;
    Vec3d  tmpPos, tmpNrm;
    for(int ib=0; ib<buckets.ncell; ib++){
        int n = buckets.cellNs[ib];
        if(n<=0) continue;
        const Quat8d& bb = BBs[ib];
        double t = rayBox( ray0, hRay, (Vec3d)bb.lo.f, (Vec3d)bb.hi.f, tmpPos, tmpNrm );
        if( t >= 1e30 ) continue;
        int i0 = buckets.cellI0s[ib];
        for(int j=0; j<n; j++){
            int ie = buckets.cell2obj[i0+j];
            int2 b = bonds[ie];
            double t1,t2;
            Vec3d p0 = (Vec3d)points[b.x].f;
            Vec3d d  = (Vec3d)points[b.y].f - p0;
            double l = d.normalize();
            double dist = rayLine( ray0, hRay, p0, d, t1, t2 );
            if( (dist<dist_min) && (t2>0) && (t2<l) ){ imin=ie; dist_min=dist; }
        }
    }
    return imin;
}

#endif
