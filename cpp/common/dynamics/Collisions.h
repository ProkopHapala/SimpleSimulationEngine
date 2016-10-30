#ifndef Collisions_h
#define Collisions_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "Vec3.h"
#include "Mesh.h"
#include "Terrain25D.h"
#include "Body.h"

/*
Task:
separate collision shape from body dynamics
body reaction on collision should be specified by user, but rutines for evaluation of collision should be prefabricated


*/

class CollisionShape{
	public:
	double collision_radius=1.0d;
	int    displayList;
	virtual bool colideWithLineSegment( const Vec3d& p1, const Vec3d& p2, const Vec3d& center, Vec3d * where, Vec3d * normal );
	virtual bool colideWithTerrain    ( Terrain25D * terrain, RigidBody& rb  );
};

class MeshCollisionShape : public CollisionShape{
	public:
	Mesh * mesh;
	//virtual bool colideWithTerrain( Terrain25D * terrain, const Mat3d& orientation, const Vec3d& cog, Vec3d& force, Vec3d& torq );
	virtual bool colideWithTerrain( Terrain25D * terrain, RigidBody& rb );
};


class CollisionObject{
	public:
	CollisionShape * collisionShape;
	virtual bool colideWithLineSegment( const Vec3d& p1, const Vec3d& p2, Vec3d * where, Vec3d * normal ) = 0;
};


#endif
