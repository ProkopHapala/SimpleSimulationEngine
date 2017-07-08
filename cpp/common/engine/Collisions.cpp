
#include "fastmath.h"
#include "raytrace.h"

#include "Collisions.h" // THE HEADER

// just for debug
//#include <SDL2/SDL.h>
//#include <SDL2/SDL_opengl.h>
//#include "Draw3D.h"


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

/*
bool CollisionShape::colideWithTerrain( Terrain25D * terrain, const Mat3d& orientation, const Vec3d& cog, Vec3d& force, Vec3d& torq ){
    //torq.set(0.0);
    Vec2d dv;
    double h  = terrain->eval( {cog.x,cog.z}, dv);
    double dh = h+collision_radius-cog.y;
    if(dh>0){
        force.add(-dv.x*dh,dh,-dv.y*dh);
    };
};
*/


bool CollisionShape::colideWithTerrain( Terrain25D * terrain, RigidBody& rb ){
    //torq.set(0.0);
    Vec2d dv;
    double h  = terrain->eval( {rb.pos.x,rb.pos.z}, dv);
    double dh = h+collision_radius-rb.pos.y;
    if(dh>0){
        rb.force.add(-dv.x*dh,dh,-dv.y*dh);
    };
};

// ========================
// ==== MeshCollisionShape
// ========================

bool MeshCollisionShape::colideWithTerrain( Terrain25D * terrain, RigidBody& rb ){
    //torq.set(0.0);
    //force.set(0.0);
    for(int i=0; i<mesh->points.size(); i++){
        Vec3d gdp,gv,gp;
        rb.velOfPoint(mesh->points[i], gv, gdp);
        gp.set_add(gdp, rb.pos);
        //glColor3f( 0.0f,0.0f,0.8f); Draw3D::drawVecInPos( v, p+rb.pos );

        Vec2d dv;
        double h = terrain->eval( {gp.x,gp.z}, dv);
        double dh = h-gp.y;
        if(dh>0){
            Vec3d f;
            //double clat  = 30.0*(v.x*dv.x + v.z*dv.y)/dv.norm2(); if(clat<0.0) clat*=0.5f;
            double clat  = 30.0;
            double cvert = (gv.y<0)?30.0:10.0;
            f.set(-dv.x*dh*clat, dh*cvert, -dv.y*dh*clat);
            rb.apply_force(f,gdp);
            //glColor3f( 0.8f,0.0f,0.0f); Draw3D::drawVecInPos( f, p+rb.pos );
        }

    }
};

// =====================
// ==== CollisionObject
// =====================
