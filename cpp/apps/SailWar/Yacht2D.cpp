
#include <SDL2/SDL_opengl.h>

#include "Draw2D.h"

#include "Yacht2D.h" // THE HEADER

void Yacht2D::draw( ){
    /*
	keel.draw  ( *this );
	rudder.draw( *this );
	mast.draw  ( *this );
	*/
	Draw2D::draw_attached_vec( keel.pos,   keel.rot,   pos, rot, {0.1,0.5} );
	Draw2D::draw_attached_vec( rudder.pos, rudder.rot, pos, rot, {0.1,0.5} );
	Draw2D::draw_attached_vec( mast.pos,   mast.rot,   pos, rot, {0.1,0.5} );
};

void Yacht2D::applySailForces( const Vec2d& windSpeed, const Vec2d& watterSpeed ){
	Vec2d lvel;
	lvel.set_sub( watterSpeed, vel );
	//printf( " keel   ");
	keel  .assertAeroForce( *this, lvel, 1000.0 );
	//printf( " rudder ");
	rudder.assertAeroForce( *this, lvel, 1000.0 );
	lvel.set_sub( windSpeed, vel   );
	//printf( " mast   ");
	mast.assertAeroForce  ( *this, lvel, 1.225  );
	//glColor3f( 0.5f, 0.1f, 0.1f ); Draw2D::drawVecInPos_d( force*0.01, pos );
};

void Yacht2D::testSail( int n, int nsub, double dt, double windSpeed, Vec2d * poss, Vec2d * vels, Vec2d * rots ){
    //pos.set(0.0);
    //vel.set(0.0);
    Vec2d wind_speed;   wind_speed  .set( windSpeed, 0.0d );
    Vec2d watter_speed; watter_speed.set(            0.0d );
    for( int i=0; i<n; i++ ){
        for ( int j=0; j<nsub; j++ ){
            clean_temp( );
            applySailForces( wind_speed, watter_speed );
            move( dt );
        }
        //Draw2D::drawBody2d_d( rot, pos, 0.5, 0.2 );
        //Draw2D::drawPointCross_d( pos, 0.2 );
        if( poss != NULL ) poss[i].set( pos );
        if( vels != NULL ) vels[i].set( vel );
        if( rots != NULL ) rots[i].set( rot );
        //printf( " (%3.3f,%3.3f) (%3.3f,%3.3f) (%3.3f,%3.3f) \n", poss[i], vels[i], rots[i] );
    }
}

int Yacht2D::convergeSail( int nmax, int nsub, double dt, double windSpeed, double vtol, double rtol, Vec2d& vel_conv, Vec2d& rot_conv ){
    double vtol2 = vtol*vtol;
    double rtol2 = rtol*rtol;
    Vec2d wind_speed;   wind_speed  .set( windSpeed, 0.0d );
    Vec2d watter_speed; watter_speed.set(            0.0d );
    Vec2d old_vel; old_vel.set(0.0d);
    Vec2d old_rot; old_rot.set(0.0d);
    for( int i=0; i<nmax; i++ ){
        for ( int j=0; j<nsub; j++ ){
            clean_temp( );
            applySailForces( wind_speed, watter_speed );
            move( dt );
            glVertex3f( (float)pos.x, (float)pos.y, 0.0 );
        }
        Vec2d err;
        err.set_sub( vel, old_vel ); double vel_err2 = err.norm2();
        err.set_sub( rot, old_rot ); double rot_err2 = err.norm2();
        //printf( "pos(%3.3f,%3.3f)vel(%3.3f,%3.3f)(%3.3f,%3.3f)rot(%3.3f,%3.3f)(%3.3f,%3.3f)\n", pos, vel, old_vel,  rot, old_rot,  );
        //printf( "vel(%3.3f,%3.3f)(%3.3f,%3.3f)rot(%3.3f,%3.3f)(%3.3f,%3.3f) %1.1e %1.1e\n", vel, old_vel, rot, old_rot, vel_err2, rot_err2  );
        if( ( vel_err2 < vtol2 ) && ( rot_err2 < rtol2 ) ){
            vel_conv.set( vel );
            rot_conv.set( rot );
            return i;
        }
        old_rot.set( rot );
        old_vel.set( vel );
    }
    return nmax;
}

void Yacht2D::evalPolar( int n, double dt, double vtol, double rtol, double * phi_rudder, double * phi_mast, double * wind_speed, Vec2d * vels, Vec2d * rots, bool reseting ){
    for( int i =0; i<n; i++ ){
        pos.set(0.0d);
        if( reseting ){
            vel.set( vels[i] ); // it is probably better start from velocity of previous step as estimate
            rot.set( rots[i] );
        }
        constexpr double vtol = 1.0e-8;
        constexpr double rtol = 1.0e-8;
        Vec2d vel_conv, rot_conv;
        mast  .setAngle( phi_mast[i]   );
        rudder.setAngle( phi_rudder[i] );
        int nmax  = 10000;
        glBegin( GL_LINE_STRIP );
        int nconv = convergeSail( nmax, 20, dt, wind_speed[i], vtol, rtol, vel_conv, rot_conv );
        glEnd();
        if ( nconv < nmax ){
            vels[i].set( vel_conv );
            rots[i].set( rot_conv );
        }
    }
};
