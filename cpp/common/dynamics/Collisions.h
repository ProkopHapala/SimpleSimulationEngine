#ifndef Collisions_h
#define Collisions_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "Vec3.h"

class CollisionShape{
	public:
	double collision_radius;
	int    displayList;
	virtual bool colideWithLineSegment( const Vec3d& p1, const Vec3d& p2, const Vec3d& center, Vec3d * where, Vec3d * normal );
};

class CollisionObject{
	public:
	CollisionShape * collisionShape;
	virtual bool colideWithLineSegment( const Vec3d& p1, const Vec3d& p2, Vec3d * where, Vec3d * normal ) = 0;
};

#endif
