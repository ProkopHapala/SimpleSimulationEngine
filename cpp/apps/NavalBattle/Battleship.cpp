
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

    int n;
    //---- turret types
    printf("//---- turret types \n");
    fgets_comment( line, nbuf, pFile ); sscanf(line,"%i\n", &n );
    for(int i=0; i<n; i++ ){
        if ( fgets_comment( line, nbuf, pFile ) ){
            TurretType * TT = new TurretType();
            TT->fromString(line);
            turretTypes.push_back(TT);
        };
    };
    //---- turrets
    printf("//---- turrets \n");
    fgets_comment( line, nbuf, pFile ); sscanf(line,"%i\n", &n );
    for(int i=0; i<n; i++ ){
        if ( fgets_comment( line, nbuf, pFile ) ){
            Turret * T = new Turret();
            T->fromString(line);
            T->type = turretTypes[T->kind];
            //printf("T->type %i %i \n", T->kind, T->type );
            turrets.push_back(T);
        };
    };
    return false;
};

void Battleship::render( ){
    float sc = 10.0;
    shape = glGenLists(1);
    glNewList( shape, GL_COMPILE );
		glScalef( sc, sc, sc ); Draw3D::drawMesh( *mesh );   // the ship model in obj is 1:10 smaller
		//Draw3D::drawSphere_oct( 3, 1.0, {0.0,0.0,0.0} );
	glEndList();


	for( TurretType* TT : turretTypes ){
        //glScalef( 1.0, 1.0, 1.0 );
        TT->shape = glGenLists(1);
        glNewList( TT->shape, GL_COMPILE );
            sc = TT->Rsize*0.3; glScalef( sc, sc, sc ); Draw3D::drawMesh( *(TT->mesh) );
            //printf( ">>>>> TT->shape %i \n",  TT->shape);
        glEndList();
	}
	for( Turret* tur : turrets ){
        //printf( ">>>>> %i %i \n",  tur->shape, tur->type );
        //printf( ">>>>> %i %i %f \n",  tur->shape, tur->type, tur->type->mass );
        //printf( ">>>>> %i %i %i\n",  tur->shape, tur->type, tur->type->shape);
        //exit(0);
        tur->shape = tur->type->shape;
    }

}

void Battleship::draw( ){

    glColor3f(0.8,0.8,0.8);
    Draw3D::drawShape( shape, pos3d, rot3d);
    glColor3f(0.8,0.0,0.0);
    for( Turret* tur : turrets ){
        printf( "(%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f)\n", tur->gpos.x, tur->gpos.y, tur->gpos.z, tur->gun_dir.x, tur->gun_dir.y, tur->gun_dir.z );
        //Draw2D::drawVecInPos( {tur->gpos.x,tur->gpos.z}, {tur->grot.c.x,tur->grot.c.z} );

        Draw3D::drawVecInPos( tur->gun_dir*50.0, tur->gpos );

        //Draw3D::drawShape( tur->shape, {0.0,0.0,0.0}, {1.0,0.0,0.0,   0.0,1.0,0.0,  1.0,0.0,0.0} );
        //Draw3D::drawShape( tur->shape, tur->gpos, {1.0,0.0,0.0,   0.0,1.0,0.0,  0.0,0.0,1.0} );
        Draw3D::drawShape( tur->shape, tur->gpos, tur->grot );




        //Draw2D::drawShape( shape, pos, rot );
        //Draw2D::drawShape( tur->shape, pos, rot );
    }

}

int Battleship::shoot( std::vector<Projectile3D*>& projectiles ){
    int nshots = 0;
    for( Turret* tur : turrets ){
        if( tur->reload >= 1.0 ){
            printf("reload rate %f\n", tur->type->reload_rate );
            nshots += tur->shoot( projectiles, {vel.x, 0.0, vel.y} );
        }
        if(nshots>=nsalvo) break;
    }
    return nshots;
};

void Battleship::update( double dt, const Vec3d& wind_speed, const Vec3d& watter_speed ){
    printf(" Battleship::update \n");

    //if( life < life_max ) { life += life_regeneration_rate * dt; }

    //Mat3d rotmat; rotmat.set( {rot.x,0,rot.y}, {0.0,1.0,0.0}, {-rot.y,0,rot.x} );

    for( Turret* tur : turrets ){
        tur->updateTransforms( pos3d, rot3d );
        tur->update( dt );
    }

    // from Ship2D
    clean_temp( );
    applyHydroForces( {0.0,0.0} );
    move_RigidBody2D(dt);
    sync3D();
}

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
    Draw2D::drawShape( collisionShape->displayList, pos, rot );
}
*/


