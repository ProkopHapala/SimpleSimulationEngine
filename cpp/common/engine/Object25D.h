#ifndef Object25D_h
#define Object25D_h

#include "fastmath.h"
#include "Vec2.h"
//#include "Mat3.h"
//#include "geom3D.h"
#include "Collisions.h"
#include "Body.h"

class Object25D {
	public:
	int id, kind, shape;
	// geometry
    Vec2d lpos,span;
    Vec2d lrot;
    // axuliaries
    Vec2d gpos;
    Vec2d grot;

	//Ellipsoide       bounds;
	CollisionShape * coll      = NULL;
	RigidBody      * controler = NULL;

	inline void updateTransforms_Object25D( const Vec2d& pos0, const Vec2d& rot0 ){
        //rot0.dot_to( lpos, gpos ); gpos.add( pos0 );
        //grot.set_mmul( lrot, rot0 );
	}

	virtual void updateTransforms( const Vec2d& pos0, const Vec2d& rot0 ){
	    updateTransforms_Object25D( pos0, rot0 );
	}

    virtual bool pointIn( const Vec2d& point ){
        return false; // place holder
        // distance from bounding ellipsoide
        //return bounds.pointIn( point );
        //return 0>dist_Ellipsoide( point, gpos, grot, span );
    }

    virtual double ray( const Vec2d& ray0, const Vec2d& hRay, Vec2d * normal ){
        return 1e+300; // place holder
        // distance from bounding volume
        // transform space
        //return bounds.ray( hRay, ray0, normal );
        //normal_Ellipsoide( const Vec3d& phit, const Vec3d& pos, const Mat3d& orientation, const Vec3d& invspan );
        //return ray_Ellipsoide( ray0, hRay, normal, gpos, grot, span );
    }


    inline void initOne(){
    //    lpos.set(0.0,0.0,0.0);
    //    span.set(1.0,1.0,1.0);
    //    lrot.a.set(1.0,0.0,0.0);
    //    lrot.b.set(0.0,1.0,0.0);
    //    lrot.c.set(0.0,0.0,1.0);
    //    gpos=lpos; grot=lrot;
    }

    //inline bool fromVecs( const Mat3d& mat ){
    //    lrot.set(mat);
    //    span.a = lrot.a.normalize();
    //    span.b = lrot.b.normalize();
    //    span.c = lrot.c.normalize();
    //    //gpos=lpos;
    //    grot=lrot;
    //}

    // place-holder functions
    //virtual void draw();

};

#endif
