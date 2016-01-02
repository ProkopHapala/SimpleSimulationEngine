
#include "fastmath.h"
#include "raytrace.h"

#include "Collisions.h" // THE HEADER

// ===================
// ==== CollisionShape 
// ===================

bool CollisionShape::colideWithLineSegment( const Vec3d& p1, const Vec3d& p2, const Vec3d& center, Vec3d * where, Vec3d * normal ){
	Vec3d hRay; 
	hRay.set_sub( p2, p1 );
	double r = hRay.norm();
	hRay.mul( 1.0d/r );
	double t = raySphere( p1, hRay, collision_radius, center );
	if( ( t > 0 ) && ( t < r ) ){
		if( where != NULL ){
			where->set_mul( hRay , t );
			where->add    ( p1       );
		}
		if( normal != NULL ){
			normal->set_sub( p1, center );
			normal->add_mul( hRay, t    );
		}
		return true;
	}else{
		return false;
	}
};

// =====================
// ==== CollisionObject
// =====================
