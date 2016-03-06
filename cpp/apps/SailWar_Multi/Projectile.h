
#ifndef Projectile_h
#define Projectile_h

#include <vector>

#include "Vec3.h"
#include "Body.h"
#include "Collisions.h"

class SailWarWorld;

class Projectile : public PointBody {
	public:
	SailWarWorld * world;
	double dragCoef;
	double penetration;
	double damage;
	Vec3d  old_pos;

	void ground_hit( );
	bool check_hit_ground(  );

	template <typename T> bool check_hit_vector( std::vector<T*>& objlist );
	//template <typename T> bool check_hit( T* obj );
	void update_old_pos();
	void addDragForce( const Vec3d& vwind, Vec3d& aeroForce );
	virtual void evalForce( );
	//virtual void move( double dt );
	virtual void draw();
};

// ============== template implementation


/*
template <typename T>
bool Projectile::check_hit( T* obj ){
    obj->colideWithLineSegment( old_pos, pos, NULL, NULL );
}
*/

template <typename T>
bool Projectile::check_hit_vector( std::vector<T*>& objlist ){
    Vec3d hit_pos;
    bool hitted = false;
    for( auto obj : objlist ){
        hitted = obj->colideWithLineSegment( old_pos, pos, &hit_pos, NULL );
        if( hitted ) break;
    }
    return hitted;
}

#endif
