#ifndef CombatMoves_h
#define CombatMoves_h

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "geom3D.h"
#include "SoftPolyLine3D.h"

class CombatMove{ public:

    Vec3d*   pose;
    double*  times;      // time samples
    double** controlSeq; // commads to actuators int time-slots

    double compare( const SoftPolyLine3d& arm, const Vec3d& p0, const Mat3d& rot ){
        Vec3d* posei = pose;
        double r2 = 0.0;
        for(int i=0; i<arm.n; i++){
            Vec3d p,v;
            rot.dot_to( arm.ps[i] - p0, p );
            rot.dot_to( arm.vs[i]     , v );
            r2+=posei->dist2(p);
            r2+=posei->dist2(v);
            posei+=2;
        }
        return r2;
    }

};

#endif
