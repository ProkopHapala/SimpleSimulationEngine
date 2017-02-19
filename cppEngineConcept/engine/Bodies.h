    #ifndef Objects_h
    #define Objects_h

    #include "fastmath.h"
    #include "Vec3.h"
    #include "Mat3.h"
    #include "quaternion.h"

    // see 11.11 Action State Machines

    // ==== interfaces
    class Collidable{
        public:
    /*
        virtual double point_in      ( const Vec3d& point,  InsideResult*    result );    // returns distance estimate
        virtual double raycast       ( const Ray&  ray, RayHit*          result );  // return "t" of ray
        virtual void   project1D     ( const Ray&  ray, Vec2d*           result );
        virtual void   project2D     ( const Mat3& cam, Convex2d*        result );
        virtual void   getAABB       ( AABB*   box );
        virtual void   getBoundSphere( Sphere* sph );
    */
    };

    class Updateable{
        public:
        virtual void clean_temp(){};
        virtual void update( double t, double dt ){};
    };

    // ==== Prefabricates

    class _identifiers {
        public:
        int id;            // instance identifier
        int kind;          // type identifier - used also to identify GL objects
    };
    class _pos         {
        public:
        Vec3d pos;        // global position
    };
    class _rot         {
        public:
        Quat4d rot;
    };
    class _size        {
        public:
        double Rbound;  // bounding sphere radius
        Vec3d  span;    // size of bounding box or ellipsoide
    };
    class _attached  {
        public:
        Vec3d  lpos;   //  local translation with respect to parrent scene graph node
        Quat4d lrot;   //  local rotation    with respect to parrent scene graph node
        virtual void update_glob( const Vec3d& pos0, const Quat4d& rot0 ){};  // update global position with respect to
    };

    // ==== Classes

    class KinematicBody : public virtual _pos, public _rot{
    //class KinematicBody : public _pos, public _rot{
        void toGLmat( float * mat );
    };

    class PointBody : public virtual _pos, public Updateable {
    //class PointBody : public  _pos, public Updateable {
        public:
        Vec3d   vel,force;
        double  mass;

        inline  void update_(double dt){ vel.add_mul(force, dt); pos.add_mul(vel,(dt/mass)); };  // can be called as obj->PointBody::update_(dt);
        virtual void update (double t, double dt ){ update_(dt); };
        inline  void clean_temp_(  ){ force.set(0.0);       };
        virtual void clean_temp (  ){ clean_temp_(); };
    };

    class RigidBody : public PointBody, public KinematicBody {  // WARRNING : both PointBody and KinematicBody inheriate "pos"
    //class RigidBody : public PointBody, _rot {  // WARRNING : both PointBody and KinematicBody inheriate "pos"
        public:
        Mat3d  momentOfInertia;
        Vec3d  torque;
        Vec3d  angular_velocity;
        // pos, vel, force, mass, rot, span,  id, kind  ... should be inheriate
    };


    class RigidBody_flat {  // WARRNING : both PointBody and KinematicBody inheriate "pos"
        public:
        Vec3d   pos;
        Quat4d  rot;
        Vec3d   vel,force;
        double  mass;
        Mat3d   momentOfInertia;
        Vec3d   torque;
        Vec3d   angular_velocity;
        // pos, vel, force, mass, rot, span,  id, kind  ... should be inheriate
    };

    class EngineObject : public _identifiers, public _size, public Collidable, public Updateable {
    };

    class PassiveObject : public EngineObject, public KinematicBody {
    };

    class ActiveObject : public EngineObject, public RigidBody {  // WARRNING : both EngineObject and RigidBody are "Updateable"
    };

    class AttachedObject : public PassiveObject, public _attached {
    };

    /*
    class Projectile : public PointBody, _identifiers {
        ProjectileType * type;
    };
    */

    #endif
