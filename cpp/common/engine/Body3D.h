#ifndef Body3D_h
#define Body3D_h

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "geom3D.h"
#include "Collisions.h"
#include "Body.h"

class Body3D_Interface{ public:
	//int id, kind, shape;
	virtual void   updateTransforms( const Vec3d& pos0, const Mat3d& rot0 ) = 0;
	virtual bool   pointIn         ( const Vec3d& point ) = 0;
	virtual double ray             ( const Vec3d& ray0, const Vec3d& hRay, Vec3d * normal ) = 0;
};

class Body3DType{
    Vec3d span = Vec3dOne;
    double mass=1.0,invMass=1.0;
    Mat3d invInertia;
    Mat3d Inertia;
}

class Body3D{
    Body3DType* type = 0;
    Vec3d vel;
    Vec3d angMoment;
    //Vec3d omega;     // redudant but usefull  _angVelocity = invI * _angMomentum;
    Vec3d pos;
    Mat3d rot;
    
    void moveRigidBody(double dt, const Vec3d& force, const Vec3d& torq ){
      vel.add_mul( force,  dt * invMass );
      pos.add_mul( vel, dt );
      Vec3d angMomentum; 
      _angMomentum += torq * dt;
      
      Matrix3 invI = rot * invIbody * rot.InverseRotation();
      _angVelocity = invI * _angMomentum;
      float r2Omega = _angVelocity.SquareSize();
      if (r2Omega>1e-8)
      {
        float rOmega = sqrt(r2Omega);
        Matrix3 drot = M3Identity;
        drot.SetRotationAxis(_angVelocity / rOmega, dt*rOmega);
        _orientation = drot * _orientation;
        _orientation.Orthogonalize(); // we probably do not need this each iteration
      }
    }
}

#endif
