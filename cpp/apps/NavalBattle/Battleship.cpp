
#include <math.h>
#include <cstdlib>
#include <stdio.h>
#include <vector>

//#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "IO_utils.h"
#include "Draw2D.h"

#include "Battleship.h" // THE HEADER

bool Battleship::loadFromFile( char const* filename ){
    printf(" filename: %s \n", filename );
    FILE * pFile;
    constexpr int nbuf = 1024;
    char line[nbuf];
    pFile = fopen (filename,"r");
    fgets_comment( line, nbuf, pFile );  fromString(line);
    fgets_comment( line, nbuf, pFile );  keel  .fromString( line );
    fgets_comment( line, nbuf, pFile );  rudder.fromString( line );
    return false;
};


/*
Turret ** Battleship::initGuns( int n, Vec3d pos1, Vec3d pos2, Vec3d ldir, double muzzle_velocity ){
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

void Battleship::initAllGuns( int n ){
	double dy         = 0.15;
	double dz         = 0.2;
	double muzzle_vel = 50.0;
    double elevation  = 0.00;
    double sa = sin(elevation);
	double ca = cos(elevation);
	nguns = n;
	Vec3d pos1,pos2,ldir;
	pos1.set(  0.5,  dy, dz );
	pos2.set( -1.0,  dy, dz );
	ldir.set(  0.0,  ca, sa );
	left_guns = initGuns( nguns, pos1, pos2, ldir, muzzle_vel );
	pos1.set(  0.5, -dy, dz );
	pos2.set( -1.0, -dy, dz );
	ldir.set(  0.0, -ca, sa );
	right_guns = initGuns( nguns, pos1, pos2, ldir, muzzle_vel );
}

void Battleship::fire_gun_row( int n, Gun ** guns, std::vector<Projectile*> * projectiles ){
	Vec3d vel3D,pos3D;
	Mat3d rot3D;
	vel3D.set( vel.x, vel.y, 0 );
	pos3D.set( pos.x, pos.y, 0 );
	rot3D.a.set( rot.x,  rot.y, 0.0d );
	rot3D.b.set( rot.y, -rot.x, 0.0d );
	rot3D.c.set(  0.0d,   0.0d, 1.0d );
	for( int i=0; i<n; i++ ){
		Projectile * p = guns[i]->fireProjectile( pos3D, rot3D, vel3D );
		if( p == NULL ) { printf( " p is NULL !! \n" ); }
		// FIXME: make sure ship does not hit itself
		projectiles->push_back( p );
		p->world = world;
	}
	printf( " %i projectiles in air !!! \n", projectiles->size() );
};




void Battleship::drawGun( Gun * gun ){
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

void Battleship::draw( ){
	keel  .draw( *this );
	rudder.draw( *this );
	mast  .draw( *this );
	if( left_guns != NULL ){
		//printf( " plotting guns \n" );
		for( int i=0; i<nguns; i++ ){
			drawGun( left_guns [i] );
			drawGun( right_guns[i] );
		};
	}
}


void Battleship::drawHitBox( ){
    float clife = (float)( life / life_max );
    glColor4f( 1.0f, clife, clife, 1.0f );
    Draw2D::drawShape( pos, rot, collisionShape->displayList );
}
*/

bool Battleship::colideWithLineSegment( const Vec3d& p1, const Vec3d& p2, Vec3d * where, Vec3d * normal ){
    bool hitted = false;
    Vec3d pos3d; pos3d.set( pos.x, pos.y , 0 );
    //printf( " Battleship::colideWithLineSegment collisionShape:  %s %s \n",  collisionShape );
    hitted = collisionShape->colideWithLineSegment( p1, p2, pos3d, where, normal );

    /*
    Vec3d dp;
    dp.set_sub( pos3d, p2 );
    double r2   = dp.norm2();
    double rmax = collisionShape->collision_radius;
    if( r2 < (rmax*rmax) ) hitted = true;
    */

    //printf( " %f %f  %f %f   %d \n",  p1.x,p1.y,    p2.x,p2.y   , hitted );
    if( hitted ){
        printf( " Figate %s hitted \n", name );
        life = 0.0;
    }
    return hitted;
};


void Battleship::update( double dt ){
    if( life < life_max ) {
        life += life_regeneration_rate * dt;
        //printf( " %s %f \n", name, life );
    }
    if( gunload_left  < 1.0d ) { gunload_left  += reload_rate * dt; }
    if( gunload_right < 1.0d ) { gunload_right += reload_rate * dt; }
}
