#ifndef Object3D_h
#define Object3D_h

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "geom3D.h"
#include "Collisions.h"
#include "Body.h"

class Object3D{
	public:
	int id, kind, shape;
	Ellipsoide       bounds;
	CollisionShape * coll      = NULL;
	RigidBody      * controler = NULL;


    virtual bool pointIn( const Vec3d& point ){
        // distance from bounding ellipsoide
        return bounds.pointIn( point );
    };

    virtual double ray( const Vec3d& ray0, const Vec3d& hRay, Vec3d * normal ){
        // distance from bounding volume
        // transform space
        return bounds.ray( hRay, ray0, normal );
    };

};

#endif
