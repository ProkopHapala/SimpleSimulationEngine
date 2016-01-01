
#include <math.h>
#include <cstdlib>
#include <stdio.h>
#include <vector>

#include <SDL2/SDL_opengl.h>

#include "drawMath2D.h"

#include "Frigate2D.h" // THE HEADER

Gun ** Frigate2D::initGuns( int n, Vec3d pos1, Vec3d pos2, Vec3d ldir, double muzzle_velocity ){
	Vec3d pos;
	double d = 1.0d / n;
	Gun ** guns = new Gun*[ n ];
	for( int i=0; i<n; i++){
		double t = ( i + 0.5d ) * d;
		Gun * gun = new Gun( );
		gun->lpos.set_lincomb( 1-t,  pos1,  t, pos2 );
		gun->set_direction( ldir );
		gun->muzzle_velocity = 	muzzle_velocity;
		guns[i] = gun;
		printf( " gun %i \n", i );
	}
	return guns;
}

void Frigate2D::initAllGuns( int n ){
	double dy = 0.15; 
	double dz = 0.1;
	double muzzle_vel = 1.0;
	nguns = n;
	Vec3d pos1,pos2,ldir;
	pos1.set(  0.5, dy, dz );
	pos2.set( -1.0, dy, dz );
	ldir.set(  0.0, 1.0, 0.2 ); ldir.normalize();
	left_guns = initGuns( nguns, pos1, pos2, ldir, muzzle_vel );
	pos1.set(  0.5, -dy, dz );
	pos2.set( -1.0, -dy, dz );
	ldir.set( 0.0, -1.0, 0.2 ); ldir.normalize();
	right_guns = initGuns( nguns, pos1, pos2, ldir, muzzle_vel );
}

void Frigate2D::fire_gun_row( int n, Gun ** guns, std::vector<Projectile*> * projectiles ){
	printf( " fire_left !!! \n" );
	Vec3d vel3D,pos3D;
	Mat3d rot3D;
	vel3D.set( vel.x, vel.y, 0 );
	pos3D.set( pos.x, pos.y, 0 );
	rot3D.a.set( rot.x,  rot.y, 0.0d );
	rot3D.b.set( rot.y, -rot.x, 0.0d );
	rot3D.c.set(  0.0d,   0.0d, 1.0d ); 
	for( int i=0; i<n; i++ ){ 
		Projectile * p = guns[i]->fireProjectile( pos3D, rot3D, vel3D ); 
		projectiles->push_back( p );
	}
	printf( " %i projectiles in air !!! \n", projectiles->size() );
}
void Frigate2D::fire_left ( std::vector<Projectile*> * projectiles ){ fire_gun_row( nguns, left_guns , projectiles ); }
void Frigate2D::fire_right( std::vector<Projectile*> * projectiles ){ fire_gun_row( nguns, right_guns, projectiles ); }


void Frigate2D::drawGun( Gun * gun ){
	Vec2d lpos, lrot; 
	lpos.set(  gun->lpos.x,   gun->lpos.y   );
	lrot.set( -gun->lrot.c.y, gun->lrot.c.x );
	Vec2d gpos, grot;  		
	grot  .set_mul_cmplx( rot, lrot );
	gpos  .set_mul_cmplx( rot, lpos );
	gpos.add( pos );
	float lperp = 0.1;
	float llong = 0.5;
	glBegin(GL_LINES);
		glVertex3f( (float)( gpos.x-grot.x*lperp), (float)(gpos.y-grot.y*lperp), 1 );   glVertex3f( (float)(gpos.x+grot.x*lperp), (gpos.y+grot.y*lperp), 1 );
		//glVertex3f( (float)( gpos.x-grot.y*llong), (float)(gpos.y+grot.x*llong), 1 );   glVertex3f( (float)(gpos.x+grot.y*llong), (gpos.y-grot.x*llong), 1 );
	glEnd();
}

void Frigate2D::draw( ){ 
	keel  .draw  ( *this );
	rudder.draw( *this );
	mast  .draw  ( *this );
	if( left_guns != NULL ){
		//printf( " plotting guns \n" );
		for( int i=0; i<nguns; i++ ){ 
			drawGun( left_guns [i] ); 
			drawGun( right_guns[i] ); 
		};
	}
}

bool Frigate2D::loadFromFile( char const* filename ){
	//using namespace std;
	printf(" filename: %s \n", filename );
	FILE * pFile;
	pFile = fopen ( filename, "r" );
	const int nbuf = 1000; 
	char line [ nbuf ];
	fgets( line, nbuf, pFile );   keel  .fromString( line );  //printf( "%s \n", keel  .toString( ) );
	fgets( line, nbuf, pFile );   rudder.fromString( line );  //printf( "%s \n", rudder.toString( ) );
	fgets( line, nbuf, pFile );   mast  .fromString( line );  //printf( "%s \n", mast  .toString( ) );
	fclose( pFile );
	//exit(0);
	return false;
}

