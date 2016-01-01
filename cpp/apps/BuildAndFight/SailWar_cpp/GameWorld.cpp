
#include <SDL2/SDL_opengl.h>

#include "GameWorld.h" // THE HEADER



const int npts = 4;
static double poss[npts*2] = { -1.0, 0.0,   0.0, -0.1,   0.0, +0.1,   +1.0, 0.0  };
static double mass[npts  ] = {  10.0, 50.0, 50.0, 10.0  };


void GameWorld::update( ){
	for( int i=0; i<perFrame; i++ ){
		ship1.clean_temp( );
		ship1.applySailForces(  *(Vec2d*)&wind_speed,  watter_speed );
		ship1.move( dt );

		ship2.clean_temp( );
		ship2.applySailForces(  *(Vec2d*)&wind_speed,  watter_speed );
		ship2.move( dt );

		std::vector<Projectile*>::iterator it = projectiles.begin();
		while( it != projectiles.end() ) {
			Projectile * p = *it; 
			p -> evalForce(    );
			p -> move     ( dt );
			if( p -> check_hit( ) ){ it = projectiles.erase( it ); }
			else                   { ++it;                  }
		}

	}
};


void GameWorld::init( ){
	int ifree,igl,nvert,ndiv;

	ground_level = 0.0d;
	watter_speed.set(   0.0, 0.0     );
	wind_speed  .set( -10.0, 0.0, 0.0 );

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

	printf( " >>> Setup  ship1: \n" );
	ship1.loadFromFile( "data/FrigateType.txt" );
	ship1.from_mass_points( 2, mass, (Vec2d*)poss );  printf( " I invI  %f %f \n", ship1.I, ship1.invI );
	ship1.setDefaults();
	ship1.setAngle( M_PI*0.6   );
	ship1.pos.set ( {0.0, 0.0} );
	ship1.omega = 0.0;
	ship1.shape = FigateShape;
	ship1.initAllGuns( 1 );

	printf( " >>> Setup  ship2: \n" );
	ship2.loadFromFile( "data/FrigateType.txt" );
	ship2.from_mass_points( 2, mass, (Vec2d*)poss );  printf( " I invI  %f %f \n", ship1.I, ship1.invI );
	ship2.setDefaults();
	ship2.setAngle( M_PI*0.6   );
	ship2.pos.set ( {1., 1.0} );
	ship2.omega = 0.0;
	ship2.shape = FigateShape;

	printf( " >>> Setup  ship1 DONE \n" );

	projectiles.reserve(100);

};


