#ifndef SwordArm_h
#define SwordArm_h


#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "geom3D.h"
#include "SoftPolyLine3D.h"

class SwordArm : public SoftPolyLine3d{ public:
    Sphere3d head;

};

#endif
