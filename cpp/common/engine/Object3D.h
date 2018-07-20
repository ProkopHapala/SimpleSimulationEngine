#ifndef Object3D_h
#define Object3D_h

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "geom3D.h"
#include "Collisions.h"
#include "Body.h"

#include "ShooterCommon.h"

class BodyControler{ public:
};

//class CollisionShape{};

class VisualShape{ public:
};

class VehicleType{ public:
};

class Object3d{ public:
    int id,kind;
    Vec3d pos;                      // global pos  (we do not need local pos - that is eventually done by controler)
    Mat3d rot;                      // global rot
    double R;                       // radisus for spherical bounding box
    BodyControler  * control     = 0; // this can modify/override state of object
    CollisionShape * colShape    = 0;
    VisualShape    * visualShape = 0;

    inline bool pointIn_Sphere( const Vec3d& p ){
        return sq(R)>pos.dist2(p);
    };

    inline bool ray_Sphere( const Vec3d& ray0, const Vec3d& hRay, Vec3d * normal ){
        double t; return sq(R)>rayPointDistance2( ray0, hRay, pos, t );
    };

    inline bool getShot_Sphere( const Vec3d& p0, const Vec3d& p1 ){
        Vec3d hdir; hdir.set_sub(p1,p0);
        double l = hdir.normalize();
        double t; return sq(R)>linePointDistance2( p0, hdir, pos, l );
    };

    virtual bool pointIn( const Vec3d& p ){ return pointIn_Sphere( p ); };
    virtual double ray( const Vec3d& ray0, const Vec3d& hRay, Vec3d * normal ){ return ray_Sphere( ray0, hRay, normal ); };

    virtual bool getShot( const Vec3d& p0, const Vec3d& p1, const ProjectileType& prjType, double dt ){
        return getShot_Sphere(p0,p1);
    };

};

class Vehicle3d : public Object3d {  public:
    Vec3d  vel;
    Vec3d  angMomentum;
    double mass;
    Vec3d  invIbody;     // This may be in type, but 1) this is faster 2) mass distribution may change during simulation
    VehicleType * type;
};




class Object3D_Interface{ public:
	//int id, kind, shape;
	virtual void   updateTransforms( const Vec3d& pos0, const Mat3d& rot0 ) = 0;
	virtual bool   pointIn         ( const Vec3d& point ) = 0;
	virtual double ray             ( const Vec3d& ray0, const Vec3d& hRay, Vec3d * normal ) = 0;

	// virtual getAABB();
	// virtual getBoundingSphere();
};

class Object3D : public Object3D_Interface { public:
	int id, kind, shape;
	// geometry
    Vec3d lpos,span;
    Mat3d lrot;
    // axuliaries
    Vec3d gpos;
    Mat3d grot;

	//Ellipsoide       bounds;
	CollisionShape * coll      = NULL;
	RigidBody      * controler = NULL;

	inline void updateTransforms_Object3D( const Vec3d& pos0, const Mat3d& rot0 ){
        rot0.dot_to( lpos, gpos ); gpos.add( pos0 );
        grot.set_mmul( lrot, rot0 );
	}

	virtual void updateTransforms( const Vec3d& pos0, const Mat3d& rot0 ){
        updateTransforms_Object3D( pos0, rot0 );
	}

    virtual bool pointIn( const Vec3d& point ){
        // distance from bounding ellipsoide
        //return bounds.pointIn( point );
        return 0>dist_Ellipsoide( point, gpos, grot, span );
    };

    virtual double ray( const Vec3d& ray0, const Vec3d& hRay, Vec3d * normal ){
        // distance from bounding volume
        // transform space
        //return bounds.ray( hRay, ray0, normal );
        //normal_Ellipsoide( const Vec3d& phit, const Vec3d& pos, const Mat3d& orientation, const Vec3d& invspan );
        return ray_Ellipsoide( ray0, hRay, normal, gpos, grot, span );
    };


    inline bool initOne(){
        lpos.set(0.0,0.0,0.0);
        span.set(1.0,1.0,1.0);
        lrot.a.set(1.0,0.0,0.0);
        lrot.b.set(0.0,1.0,0.0);
        lrot.c.set(0.0,0.0,1.0);
        gpos=lpos; grot=lrot;
    }

    inline bool fromVecs( const Mat3d& mat ){
        lrot.set(mat);
        span.a = lrot.a.normalize();
        span.b = lrot.b.normalize();
        span.c = lrot.c.normalize();
        //gpos=lpos;
        grot=lrot;
    }

    // place-holder functions
    //virtual void draw();

};


#endif
