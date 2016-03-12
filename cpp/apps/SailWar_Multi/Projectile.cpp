
#include <math.h>
#include <cstdlib>
#include <stdio.h>
#include <vector>

#include <SDL2/SDL_opengl.h>

#include "SailWarWorld.h"

#include "Projectile.h" // THE HEADER


void Projectile::ground_hit( ){
	printf( " Ground hit %f %f \n", pos.x, pos.y );
}

bool Projectile::check_hit_ground(  ){
    bool hitted = pos.z < world->ground_level;
    if( hitted ){ ground_hit( ); };
	return hitted;
}

void Projectile::addDragForce( const Vec3d& vwind, Vec3d& aeroForce ){
	Vec3d vair;
	vair.set_sub( vel, world->wind_speed );
	double vr2   = vair.norm2();
	double vr    = sqrt( vr2 );
	aeroForce.add_mul( vair, dragCoef * vr );
}

void Projectile::evalForce( ){
    //printf( " Projectile::evalForce \n" );
	force.set( 0.0, 0.0, -9.81f );
	addDragForce( world->wind_speed, force );
}

/*
void Projectile::move( double dt ){
	old_pos.set( pos );
	PointBody::move( dt );
}
*/

void Projectile::update_old_pos(){
    old_pos.set( pos );
}

void Projectile::draw(){
    double vscale = 0.01d;
	glBegin(GL_LINES);
		//glVertex3f( (float)( old_pos.x ), (float)( old_pos.y ), 0 );
		//glVertex3f( (float)( 0 ), (float)( 0 ), 0 );
		glVertex3f( (float)(     pos.x ), (float)(     pos.y     ), 0 );
		glVertex3f( (float)(     pos.x + vel.x * vscale), (float)(     pos.y + vel.y * vscale ), 0 );
	glEnd();
	//printf( " I'm projectile \n");
	//printf( " projectile   pos  %10.5f %10.5f %10.5f old_pos %10.5f %10.5f %10.5f \n",   pos.x, pos.y, pos.z,   old_pos.x, old_pos.y, old_pos.z );
}

char * Projectile::toBytes  ( char * buff ){
    (*((Vec3d* )buff)).set( pos );    buff += sizeof( Vec3d );
    (*((Vec3d* )buff)).set( vel );    buff += sizeof( Vec3d );
    //printf( "toBytes (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) \n", pos.x, pos.y, pos.z, vel.x, vel.y, vel.z );
    //double* bf = (double*)buff;
    //printf( "toBytes (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) \n", bf[-6], bf[-5], bf[-4], bf[-3], bf[-2], bf[-1] );
    return buff;
}

char * Projectile::fromBytes( char * buff ){
    pos.set ( *((Vec3d* )buff)  ); buff += sizeof(Vec3d );
    vel.set ( *((Vec3d* )buff)  ); buff += sizeof(Vec3d );
    //double* bf = (double*)buff;
    //printf( "fromBytes (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) \n", bf[-6], bf[-5], bf[-4], bf[-3], bf[-2], bf[-1] );
    //printf( "fromBytes (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) \n", pos.x, pos.y, pos.z, vel.x, vel.y, vel.z );
    return buff;
}

