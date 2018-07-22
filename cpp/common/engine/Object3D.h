#ifndef Object3D_h
#define Object3D_h

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "raytrace.h"
#include "geom3D.h"
#include "Collisions.h"
#include "Body.h"

#include "ShooterCommon.h"

class BodyControler{ public:
};

//class CollisionShape{};

class VisualShape{ public:
};

class ObjectType{ public:
    Vec3d span     = Vec3dZero;
    CollisionShape * shape = 0;         // in case can be also used to store visual shape ?
    //VisualShape    * visualShape = 0; // DEPRECATED: for drawing complex shapes, can be added in derived types
    int ogl=0;                          // for drawing simple shapes usually sufficient
};

class VehicleType : ObjectType { public:
};

class Object3d{ public:
    int id=-1;
    //int kind=-1;                      // DEPRECATED : we use type instead
    Vec3d pos=Vec3dZero;                // global pos  (we do not need local pos - that is eventually done by controler)
    Mat3d rot=Mat3dIdentity;            // global rot
    double R=0;                         // radisus for spherical bounding box
    BodyControler  * control = 0;       // this can modify/override state of object
    //CollisionShape * colShape    = 0; // DEPRECATED: inside type
    //VisualShape    * visualShape = 0; // DEPRECATED: inside type
    ObjectType * type = 0;              //

    inline bool hitPoint_Object3d( const Vec3d& p, double r ){
        return sq(R+r)>pos.dist2(p);
    };

    inline bool hitLine_Object3d( const Vec3d& p0, const Vec3d& p1, double r ){
        Vec3d hdir; hdir.set_sub(p1,p0);
        double l = hdir.normalize();
        double t; return sq(R+r)>linePointDistance2( p0, hdir, pos, l );
    };

    inline bool hitRay_Object3d( const Vec3d& ray0, const Vec3d& hRay ){
        double t; return sq(R)>rayPointDistance2( ray0, hRay, pos, t );
    };

    inline double ray_Object3d( const Vec3d& ray0, const Vec3d& hRay, Vec3d* normal ){
        double t = raySphere( ray0, hRay, R, pos );
        if(normal)sphereNormal( t, ray0, hRay, pos, *normal );
    };

    //inline bool collide_Object3d( Object3d* obj ){ return hitPoint(obj->pos, obj->R); }
    inline bool collide_Object3d( Object3d* obj ){ return sq(R+obj->R)>pos.dist2(obj->pos); }

    virtual bool   collide ( Object3d* obj )                                                             { return collide_Object3d(obj);                         };
    virtual bool   hitPoint( const Vec3d& p, double r=0.0 )                                              { return hitPoint_Object3d( p    , r );                 };
    virtual bool   hitLine ( const Vec3d& p0, const Vec3d& p1, double r=0.0 )                            { return hitLine_Object3d( p0, p1, r );                 };
    virtual double ray     ( const Vec3d& ray0, const Vec3d& hRay, Vec3d* normal )                       { return ray_Object3d( ray0, hRay, normal );            };
    virtual bool   getShot ( const Vec3d& p0, const Vec3d& p1, const ProjectileType& prjType, double dt ){ return hitLine_Object3d(p0,p1, prjType.caliber*0.5 ); };

    Object3d(){};
    Object3d( double R_, Vec3d pos_, ObjectType* type_, int id_ ){R=R_;pos=pos_;type=type_;id=id_;  };
    //Object3d( double R, Vec3d pos ){};
};

class Vehicle3d : public Object3d {  public:
    Vec3d  vel=Vec3dZero;
    Vec3d  angMomentum=Vec3dZero;
    double mass=1.0;
    Vec3d  invIbody=Vec3dOne;  // This may be in type, but 1) this is faster 2) mass distribution may change during simulation
    //VehicleType * type;      // DEPRECATED: we rather cast Object3d::type
};

////////////////////////////////////////
//        Object3D ( DEPRECATED )
//////////////////////////////////
// - currently used only in "Tanks" demo


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
