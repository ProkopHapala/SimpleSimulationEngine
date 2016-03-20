
#include <SDL2/SDL_opengl.h>

#include <vector>
#include <algorithm>

#include "fastmath.h"

#include "Draw2D.h"

#include "SailWarWorld.h" // THE HEADER

const int     npts = 4;
static double poss[npts*2] = { -1.0, 0.0,   0.0, -0.1,   0.0, +0.1,   +1.0, 0.0  };
static double mass[npts  ] = {  10.0, 50.0, 50.0, 10.0  };


void SailWarWorld::update_world( ){
	for( int i=0; i<perFrame; i++ ){

        auto it_ship = ships.begin();
		while( it_ship != ships.end() ) {
			Frigate2D * ship = *it_ship;
			ship->clean_temp( );
            ship->applySailForces(  *(Vec2d*)&wind_speed,  watter_speed );
            ship->move( dt );
            ship->update( dt );
            ++it_ship;
		}

		auto it_proj = projectiles.begin();
		while( it_proj != projectiles.end() ) {
			Projectile * proj = *it_proj;
			proj ->update_old_pos( );
			proj ->evalForce(    );
			proj ->move     ( dt );
			bool hitted = false;
			hitted |= proj->check_hit_ground( );
			hitted |= proj->check_hit_vector<Frigate2D>( ships );
			if( hitted ){
                //printf( " removing projectile \n" );
                it_proj = projectiles.erase( it_proj );
                delete proj;
            }else{
                ++it_proj;
            }
		}

	}
};

void SailWarWorld::makeShip( const Vec2d& pos, double angle, char * filename, int shape, CollisionShape * collisionShape ){
    int ith = ships.size();
	printf( " >>> Setup  ship %i \n", ith );
	Frigate2D* ship = new Frigate2D();
	ship->loadFromFile( filename );
	ship->from_mass_points( 2, mass, (Vec2d*)poss );
	//printf( " I invI  %f %f \n", ship1->I, ship1->invI );
	ship->setDefaults();
	ship->setAngle( angle );
	ship->pos.set ( pos );
	ship->omega = 0.0;
	ship->shape = shape;
	printf( "DEBUG 1 \n" );
	ship->initAllGuns( 6 );
	printf( "DEBUG 2 \n" );
	ship->world = this;
    ship->collisionShape = collisionShape;
    printf( "DEBUG 3 \n" );
    ship->name = new char[7];
    sprintf( ship->name, "Ship_%02i", ith );
    printf( "DEBUG 4 \n" );
    ships.push_back( ship );
    printf( "DEBUG 5 \n" );
}

void SailWarWorld::init_world( ){

    // ---- misc.

    printf( " SailWarWorld::init_world: misc. \n" );

	ground_level = 0.0d;
	watter_speed.set(   0.0, 0.0     );
	wind_speed  .set( -10.0, 0.0, 0.0 );

	projectiles.reserve(100);

    // ---- ships

    int ifree,igl,nvert,ndiv;

    printf( " SailWarWorld::init_world: hitBoxShape \n" );

    int hitBoxShape = glGenLists(1);
	glNewList( hitBoxShape , GL_COMPILE );
        Draw2D::drawCircle_d( {0.0f,0.0f}, 1.0f, 32, false );
	glEndList();

    printf( " SailWarWorld::init_world: defaultShipShape \n" );

    defaultShipShape = glGenLists(1);
	glNewList( defaultShipShape , GL_COMPILE );
	glBegin   (GL_TRIANGLE_FAN);
		glNormal3f( 0.0f, 0.0f, 1.0f );
		glVertex3f( +1.5,  0.0, 0 );
 		glVertex3f( +0.5,  0.2, 0 );
		glVertex3f( -1.0,  0.2, 0 );
 		glVertex3f( -1.0, -0.2, 0 );
		glVertex3f( +0.5, -0.2, 0 );
		glVertex3f( +1.5,  0.0, 0 );
	glEnd();
	glEndList();

    printf( " SailWarWorld::init_world: defaultCollisionShape \n" );

    defaultCollisionShape  = new CollisionShape();
    defaultCollisionShape->collision_radius = 1.0;
    defaultCollisionShape->displayList      = hitBoxShape;

    //SailWarWorld::makeShip( { 3.0, -3.0}, M_PI*0.6, "data/FrigateType.txt", defaultShipShape, defaultCollisionShape );
    //SailWarWorld::makeShip( {-3.0, -3.0}, M_PI*0.6, "data/FrigateType.txt", defaultShipShape, defaultCollisionShape );

/*
	printf( " >>> Setup  ship1: \n" );
	Frigate2D* ship1 = new Frigate2D();
	ship1->loadFromFile( "data/FrigateType.txt" );
	ship1->from_mass_points( 2, mass, (Vec2d*)poss );
	//printf( " I invI  %f %f \n", ship1->I, ship1->invI );
	ship1->setDefaults();
	ship1->setAngle( M_PI*0.6   );
	ship1->pos.set ( {3.0, -3.0} );
	ship1->omega = 0.0;
	ship1->shape = FigateShape;
	ship1->initAllGuns( 6 );
	ship1->world = this;
    ship1->collisionShape = collisionShape;
    ship1->name = "ship1";
    ships.push_back( ship1 );

	printf( " >>> Setup  ship2: \n" );
	Frigate2D* ship2 = new Frigate2D();
	ship2->loadFromFile( "data/FrigateType.txt" );
	ship2->from_mass_points( 2, mass, (Vec2d*)poss );
	//printf( " I invI  %f %f \n", ship1->I, ship1->invI );
	ship2->setDefaults();
	ship2->setAngle( M_PI*0.6   );
	ship2->pos.set ( {-3.0, -3.0} );
	ship2->omega = 0.0;
	ship2->shape = FigateShape;
	ship1->initAllGuns( 6 );
    ship2->world = this;
    ship2->collisionShape = collisionShape;
    ship2->name = "ship2";
    ships.push_back( ship2 );
*/

    // ---- isles

    printf( " SailWarWorld::init_world: isles \n" );

    int ncorners = 5;
    int nisles   = 30;
    std::vector<double> phi_buf(ncorners);
    for( int i=0; i<nisles; i++ ){
        for( double& phi : phi_buf ){ phi = randf() * M_PI * 2; }
        std::sort( phi_buf.begin(), phi_buf.end() );
        Convex2d * isle = new Convex2d( ncorners );
        double x0 = randf( -10.0, 10.0 );
        double y0 = randf( -8.0, 8.0 );
        for( int i=0; i<ncorners; i++ ){
            double phi = phi_buf[i];
            isle->corners[i].set( x0 + cos( phi ), y0 + sin(phi) );
        }
        isle->update_lines();
        isles.push_back( isle );
    }

};

/*
void SailWarWorld::projectile_collisions(){
    Vec3d hit_pos;
    for( auto proj: projectiles ){
        for( auto ship: ships ){
            ship->colideWithLineSegment( proj->old_pos, proj->pos, hit_pos, NULL );
        }
    }
};
*/

