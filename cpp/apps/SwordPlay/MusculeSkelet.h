#ifndef MusculeSkelet_h
#define MusculeSkelet_h


#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "geom3D.h"
#include "SoftBody.h"


struct Muscule{
    Vec2i bones;    //
    Vec2d anchors;  // anchor points along bones
};


class MusculeSkelet : public SoftBody{ public:
    Sphere3d head;
    int      nmuscules = 0;
    Muscule* muscules = 0;

    inline void interpolateBond(int i, double t, Vec3d& p)const{
        const Bond& b=bonds[i];
        p.set_lincomb( (1-t), points[b.i], t, points[b.j] );
    }

    inline void getMusculeAnchors(int i, Vec3d& p1,Vec3d& p2)const{
        const Muscule& mus =muscules[i];
        interpolateBond(mus.bones.a, mus.anchors.a, p1);
        interpolateBond(mus.bones.b, mus.anchors.b, p2);
    }

    void allocate( int npoints_, int nbonds_, int nfix_, int nmuscules_ ){
        nmuscules=nmuscules_;
        _realloc( muscules, nmuscules );
        SoftBody::allocate( npoints_, nbonds_, nfix_ );
    }

};

#endif
