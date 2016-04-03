
#ifndef Warrior3D_h
#define Warrior3D_h

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "Body.h"

class Warrior3D : public RigidBody {
	public:
    bool    landed = false;
    bool    trigger = false;
    Vec3d   surf;
    Quat4d  look;
    Vec3d   gun_rot;
    double until_reaload = 0;

    int shape;

};

#endif
