
#ifndef AeroCraftWarrior_h
#define AeroCraftWarrior_h

#include "AeroSurf.h"
#include "AeroCraft.h"
#include "Warrior3D.h"

class AeroCraftWarrior : public AeroCraft, public Warrior3D { public:

    virtual RigidBody* asRigidBody( ){ return static_cast<RigidBody*>(this); };

    virtual void move_warrior( double dt, Vec3d& wind_speed, Vec3d& gravity, Terrain25D * terrain ){
        clean_temp();
		force.add ( gravity*mass );
		applyAeroForces( wind_speed );
		move(dt);
    };

};

#endif

