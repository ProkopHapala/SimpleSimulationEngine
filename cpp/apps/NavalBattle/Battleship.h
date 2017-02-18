
#ifndef Battleship_h
#define Battleship_h



#include <vector>

#include "Vec3.h"
#include "Ship2D.h"
#include "Projectile3D.h"
#include "Gun3D.h"
#include "Turret.h"
#include "Collisions.h"

//#include <SDL2/SDL.h>
//#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "Draw3D.h"

class Battleship : public Ship2D, public CollisionObject {
	public:
    char * name;
	int nguns;
	int nsalvo = 1;
	//Turret * turrets;
	std::vector<TurretType*> turretTypes;
	std::vector<Turret*>     turrets;

	double life_max               = 1.0d;
	double life                   = 1.0d;
	double life_regeneration_rate = 2.9d;

	// ==== function declarations

	Turret ** initGuns( int n, Vec3d pos1, Vec3d pos2, Vec3d ldir, double muzzle_velocity );
	//void initAllGuns( int n );
	//void drawGun( Gun * gun );
	//virtual void draw( );
	//virtual void drawHitBox( );
	virtual bool colideWithLineSegment( const Vec3d& p1, const Vec3d& p2, Vec3d * where, Vec3d * normal );
	virtual bool loadFromFile( char const* filename );
    virtual void update( double dt, const Vec3d& wind_speed, const Vec3d& watter_speed );   // override Warrior25D,Ship2D

    virtual void draw();   // overtide RigidBody2D::draw()
    virtual void render();

    void aim( const Vec3d& target, double G ){
        for( Turret* tur : turrets ){
            tur->aim( target, G, phi );
        }
    }

    virtual int shoot( std::vector<Projectile3D*>& projectiles ); // override Warrior25D,Ship2D

};

#endif  // #ifndef Battleship_h

