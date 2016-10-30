
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

bool CollisionShape::colideWithTerrain( Terrain25D * terrain, const Mat3d& orientation, const Vec3d& cog, Vec3d& force, Vec3d& torq ){
    //torq.set(0.0);
    Vec2d dv;
    double h  = terrain->eval( {cog.x,cog.z}, dv);
    double dh = h+collision_radius-cog.y;
    if(dh>0){
        force.add(-dv.x*dh,dh,-dv.y*dh);
    };
};

// ========================
// ==== MeshCollisionShape
// ========================

bool MeshCollisionShape::colideWithTerrain( Terrain25D * terrain, const Mat3d& orientation, const Vec3d& cog, Vec3d& force, Vec3d& torq ){
    //torq.set(0.0);
    //force.set(0.0);
    for(int i=0; i<mesh->points.size(); i++){
        Vec3d p = mesh->points[i];
        Vec2d dv;
        double h = terrain->eval( {cog.x+p.x,cog.z+p.z}, dv);
        double dh = h-cog.y;
        if(dh>0){
            Vec3d f;
            f.set(-dv.x*dh,dh,-dv.y*dh);
            torq .add_cross( p, f );
            force.add(f);
        }
    }
};

// =====================
// ==== CollisionObject
// =====================
