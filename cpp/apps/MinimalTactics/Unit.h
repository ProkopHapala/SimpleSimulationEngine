

#ifndef Unit_h
#define Unit_h

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"

#include "fastmath.h"
#include "tresholdFunctions.h"
#include "Vec2.h"
#include "Vec3.h"
#include "geom2D.h"
#include "Body2D.h"

#include "MinimalTacticsCommon.h"
#include "UnitType.h"

const int UNIT_JOB_mask_move = 0b00010000;
const int UNIT_JOB_mask_pos  = 0b00100000;
const int UNIT_JOB_mask_unit = 0b01000000;
const int UNIT_JOB_mask_auto = 0b10000000;

const int UNIT_JOB_IDLE         = 0; // do whatever you want, regenerate at maximum speed
const int UNIT_JOB_COVER        = 1; // duck, cover, hide, do not move

const int UNIT_JOB_GOTO         = 2; // go to target with standard readyness, oppurtunistic shot at any target
const int UNIT_JOB_ATTACK       = 3; // go to target while shooting; prefer cover, suppresive fire and oppurtinity to destroy enemy at target aread
const int UNIT_JOB_HUNT         = 4; // chase unit in order to kill it (move and fire)
const int UNIT_JOB_SNEAK        = 5; // go to target carefully, and avoid any shooting
const int UNIT_JOB_ASSAULT      = 6; // go to target at maximum speed; prefer movement over shooting

const int UNIT_JOB_FIRE_AT_GOAL = 7;  // fire at position in terrain (do not move )
const int UNIT_JOB_FIRE_AT_UNIT = 8;  // fire at particular unit ( do not move )
const int UNIT_JOB_FIRE_AT_ARK  = 9;  // fire at any enemy which appears in the target; choose priority tragets automatically
const int UNIT_JOB_FIRE_AT_WILL = 10; // fire at any target in range;  choose priority tragets automatically


class Unit : public RigidBody2D {
	public:

    int n, ntot;
    double suppressed = 0.0; // between one and zero

    Faction  * faction;
    UnitType * type      = NULL;

    Vec2d   goal;
    Unit  * opponent;
    double  opponent_score = 0.0;

    Vec2d  attentionDir;
    double maxwf = 1.0;

    double radius = 1.0;

    int  job         = UNIT_JOB_IDLE;
    int  default_job = UNIT_JOB_IDLE;

    double reload = 0.0;
    int    ammo   = 100;

    // ===== inline functions

    void move_to_goal( double dt ){
        Vec2d d;
        d.set_sub( goal, pos );
        double dist = d.norm();
        double vdt = type->speed * dt;
        if( dist > vdt ){
            pos.add_mul( d, vdt/dist );
        }else{
            pos.set(goal);
            job = default_job;
        }
    }

    void fire_at_unit( Unit * target ){
        reload = 0.0;
        Vec2d d;
        d.set_sub( target->pos, pos );
        double dist = d.normalize();
        // this does not care about time-of-flight neither raytracing we may later emit projectile instead of this
        Draw2D::drawLine_d(pos,target->pos);
        target->getShot( d, type->fire_nburst, sq(type->fire_range * dist), type->fire_damage_const, type->fire_damage_kinetic, type->fire_penetration );
    }

    void update( double dt ){

        suppressed -= type->recovery_rate * dt; if( suppressed<0.0 ) suppressed=0.0;
        dt         *= (1.0-suppressed);
        reload     += dt/type->fire_period;

        switch(job){
            case UNIT_JOB_GOTO:         move_to_goal( dt ); break;
            case UNIT_JOB_FIRE_AT_UNIT: if(opponent==NULL){ job=default_job; break;} if( reload > 1.0 ) fire_at_unit( opponent ); break;
        }

    }

    double damage_ramp( double att, double def ){ return Treshold::r2( (att-def)/(att+def) ); }

    void getShot( const Vec2d& from, int nburst, double area, double damage_const, double damage_kinetic, double penetration ){

        // evaluate expected damage
        //double  ca  = from.dot( dir );
        // double armor = 0.5* (1+ca) * (type->armorFw - type->armorBg )  + type->armorBg;
        double pass = ( penetration - type->armorFw )/penetration;
        if( pass < 0 ) return;
        double damage = pass*damage_kinetic + damage_const;

        // evaluate hit probability
        double hit_prob  =            type->hit_area/(type->hit_area + area);
        double kill_prob = hit_prob * damage_ramp( damage, type->damage_tolerance );

        suppressed -= kill_prob;
        // TODO : carefullness to modulate suppression and damage ??
        for( int i = 0; i<nburst; i++ ){
            if( randf() > kill_prob ) continue;
            n--;  // soft kill
            if( randf() > type->heal_prob ){ ntot--; } // hard kill
        }
    }

    void setOpponent( Unit * opponent_ ){
        job = UNIT_JOB_FIRE_AT_UNIT;
        opponent = opponent_;
    };

    void setGoal( const Vec2d& goal_ ){
        job = UNIT_JOB_GOTO;
        goal.set( goal_ );
    };

    void setType( UnitType* type_ ){
        type = type_;
    }

    void render( Vec3f& c ){
        glColor3f( c.x, c.y, c.z );
        Draw2D::drawCircle_d( pos, radius, 16, false );
        Draw2D::drawVecInPos_d( attentionDir, pos );
        //printf( "render (%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) \n", pos.x, pos.y, pos.x, c.x, c.y, c.z );
    }

    void renderJob( Vec3f& c ){
        if(job == UNIT_JOB_GOTO         ) Draw2D::drawLine_d( pos, goal );
        if(job == UNIT_JOB_FIRE_AT_UNIT ) Draw2D::drawLine_d( pos, opponent->pos );
        //printf( "render (%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) \n", pos.x, pos.y, pos.x, c.x, c.y, c.z );
    }

    Unit(){};
    Unit( UnitType* type_, Faction* faction_, const Vec2d& pos_ ){
        pos.set(pos_);
        setType( type_ );
        faction   = faction_ ;
    }

};

#endif
