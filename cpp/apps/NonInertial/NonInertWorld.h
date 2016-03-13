

#ifndef NonInertWorld_h
#define NonInertWorld_h

#include <vector>
//#include <unordered_set>
//#include <algorithm>

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"
#include "Body2D.h"

//#include "NBodyWorld2D.h"       // TODO: We would like to add this later

#include "Warrior2D.h"
#include "Projectile2D.h"

// class NonInertWorld : public NBodyWorld2D {       //  TODO:  We would like to add this later
class NonInertWorld {
	public:

    int perFrame = 10;
    double dt    = 0.005d;

    double restitution = -0.8d;
    double airDrag     = -0.05d;
    double landDrag    = -0.5d;

    double rmax  = 15.0d;
    Vec2d  center;       // center of rotation
    double omega = 0.2d; // angular velocity of the space station
    double phi   = 0.0d;
    Vec2d  rot;

    int defaultWarriorShape;

	std::vector<Warrior2D*>    warriors;
	std::vector<Projectile2D*> projectiles;  // see http://stackoverflow.com/questions/11457571/how-to-set-initial-size-of-stl-vector

    // ==== function declarations

    virtual void update_world( );
    virtual void init_world  ( );

    void makeWarrior( const Vec2d& pos, double angle, char * filename, int shape );

    // ==== inline functions

    inline void addEnviroForces( const Vec2d& pos, const Vec2d& vel, Vec2d& fout, bool landed ){
        /*
        // nonInertial - not very stable numerically
        Vec2d rvec,ovec;
        rvec.set_sub( pos, center );
        double r     = rvec.norm();
        double inv_r = 1/( r + 1.0e-2 );  // we could use approximation https://en.wikipedia.org/wiki/Fast_inverse_square_root
        //double inv_r = 1/r;
        rvec.mul     ( inv_r );
        ovec.set_perp( rvec  );
        double vo = ovec.dot( vel ) - ( omega * r );
        double vr = rvec.dot( vel );
        fout.add_mul( ovec, -2*vr*vo*inv_r );
        fout.add_mul( ovec, 1 );
        fout.add_mul( rvec,    vo*vo*inv_r );
        //fout.add_mul( vel,     airDrag     );
        */

        Vec2d rvec,ovec,venv;
        rvec.set_sub( pos, center );
        double r     = rvec.norm();
        ovec.set_perp( rvec  );
        venv.set    ( vel );
        venv.add_mul( ovec, -omega  );
        if( landed ) {
            fout.set_mul( venv, landDrag );
        }else{
            fout.set_mul( venv, airDrag );
        };
    }

    inline bool collideWithWorld( const Vec2d& pos, Vec2d& vel, Vec2d& normal ){
        double r2 = pos.norm2( );
        if( r2 > sq(rmax) ){
            double r     = sqrt(r2);
            double inv_r = 1/( r + 1.0e-8 );
            normal.set_mul( pos, inv_r );
            if( r > (rmax+0.5) ){
                double vnormal = normal.dot( vel );
                //printf( "vel (%3.3f,%3.3f,%3.3f) normal (%3.3f,%3.3f,%3.3f) vnormal %3.3f \n" );
                if ( vnormal > 0 ){
                    vel.add_mul( normal, (restitution-1)*vnormal );
                    //vel.add_mul( normal, -vnormal );
                    //vel.mul( restitution );
                }
            }
            //if( norm.dot_perp( pos ) > 1.0e-4 ) normal.set_mul( pos, pos.norm() ); // we set normal only if it changed more than 1e-4 to save computation of sqrt
            return true;
        }
        return false;
    }

};

#endif
