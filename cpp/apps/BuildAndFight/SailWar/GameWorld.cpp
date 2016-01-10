
#include <SDL2/SDL_opengl.h>

#include <vector>
#include <algorithm>

#include "fastmath.h"
#include "drawMath2D.h"

#include "GameWorld.h" // THE HEADER

const int npts = 4;
static double poss[npts*2] = { -1.0, 0.0,   0.0, -0.1,   0.0, +0.1,   +1.0, 0.0  };
static double mass[npts  ] = {  10.0, 50.0, 50.0, 10.0  };


void GameWorld::update( ){
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


void GameWorld::init( ){
	int ifree,igl,nvert,ndiv;

	ground_level = 0.0d;
	watter_speed.set(   0.0, 0.0     );
	wind_speed  .set( -10.0, 0.0, 0.0 );

    int hitBoxShape = glGenLists(1);
	glNewList( hitBoxShape , GL_COMPILE );
        Draw2D::drawCircle_d( {0.0f,0.0f}, 1.0f, 32 );
	glEndList();

    int FigateShape = glGenLists(1);
	glNewList( FigateShape , GL_COMPILE );
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

    CollisionShape * collisionShape = new CollisionShape();
    collisionShape->collision_radius = 1.0;
    collisionShape->displayList = hitBoxShape;

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

	printf( " >>> Setup  ship1 DONE \n" );

	projectiles.reserve(100);

/*
	Convex2d * isle1 = new Convex2d( 5 );
	isle1->corners[0].set( -1.0, -1.0 );
    isle1->corners[1].set( +1.0, -1.0 );
    isle1->corners[2].set( +2.0,  0.0 );
	isle1->corners[3].set( +1.0, +1.0 );
	isle1->corners[4].set( -1.0, +1.0 );
    isle1->update_lines();
	isles.push_back( isle1 );
*/

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
void GameWorld::projectile_collisions(){
    Vec3d hit_pos;
    for( auto proj: projectiles ){
        for( auto ship: ships ){
            ship->colideWithLineSegment( proj->old_pos, proj->pos, hit_pos, NULL );
        }
    }
};
*/

