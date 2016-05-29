

#ifndef Soldier_h
#define Soldier_h

#include "fastmath.h"
#include "tresholdFunctions.h"
#include "Vec2.h"
#include "Vec3.h"
#include "geom2D.h"
#include "Body2D.h"


#include "FormationTacticsCommon.h"
#include "SoldierType.h"

static const double body_push = 10.0;

static const float injury_chances [8] = {  // MUST BE NORMALIZE !!!
 0.2,  // 1 disturbed
 0.15, // 2 frightened
 0.2,  // 3 minor wound
 0.05, // 4 deserted broken moral
 0.1,  // 5 cannot fight
 0.1,  // 6 cannot move
 0.1,  // 7 coma
 0.1   // 8 die
};

static const float impair_colors [8][3] = {
//  health  ability moral
   { 1.0f,  1.0f,   0.66f }, // 1 disturbed
   { 1.0f,  1.0f,   0.33f }, // 2 frightened
   { 0.5f,  1.0f,   1.0f },  // 3 minor wound
   { 1.0f,  1.0f,   0.0f },  // 4 deserted broken moral
   { 0.5f,  0.33f,  0.0f },  // 5 cannot fight
   { 0.5f,  0.66f,  1.0f },  // 6 cannot move
   { 1.0f,  0.0f,   0.0f },  // 7 coma        // fall to coma (head injury and trauma)
   { 0.0f,  0.0f,   0.0f }   // 8 die         // totaly dead
};

class Soldier : public RigidBody2D{
	public:
	Formation   * formation = NULL;
    SoldierType * type      = NULL;
    Soldier     * opponent;
    double opponent_score = 0.0;
    //bool alive   = true;
    //bool capable = true;

    Vec2d  attentionDir;
    Vec2d  willForce;
    double maxwf = 1.0;

    bool    fire = false;
    uint8_t impair_mask  = 0;   //  0=fit, 1==incapable, 4=deserted, 128==dead    probably bit mask
    double  moral        = 1.0; // pshycical wellness
    double  stamina      = 1.0; // pshycical wellness
    double  time_buf     = 1.0; // action point for attack ( both meele and ranged )
    double  shield       = 1.0;

    // ===== inline functions

    inline void   clean_temp   ( ){ force.set(0.0); willForce.set(0.0); attentionDir.set(0.0); opponent=NULL; opponent_score=0; }
    inline double action_period( ){ return fire ? type->fire_period : type->melee_period; }

    Vec3f impair2color(  ){
        Vec3f c; c.set(1.0,1.0,1.0);
        for( int i=0; i<8; i++ ){
            if( impair_mask  & ( 1<<i ) ){
                c.x = fmin( c.x, impair_colors[i][0] );
                c.y = fmin( c.y, impair_colors[i][1] );
                c.z = fmin( c.z, impair_colors[i][2] );
            }
        }
        if( impair_mask != 0 ){ printf( " injury %i color (%3.3f,%3.3f,%3.3f) \n", impair_mask, c.x, c.y, c.z ); };
        return c;
    }

    void update( double dt, double tAttack ){
        /*
        double  rwf2 = willForce.norm2();
        if( rwf2 > (maxwf*maxwf) ){
           willForce.mul( maxwf / sqrt( rwf2 ) );
        }
        */

        double dstamina = type->stamina_regain;
        if( time_buf > tAttack ){
            if( opponent ) attack_melee( opponent );
        }else{
            double ds  = dt*stamina;
            time_buf  += ds;
            stamina   *= ( 1 - 10*ds*dstamina );
            //stamina  -= ds;
        }
        if( stamina  < 1.0 ) stamina += dt*dstamina;


        //stamina = clip( stamina + dt, 0.0, type->stamina   );
        //charge  = clip( charge  + dt, 0.0, action_period() );
        //stamina = fmin( stamina + dt, type->stamina  );
        //charge  = fmin( charge  + dt, fire ? type->fire_period : type->melee_period );
        //moral   = fmax( moral +=dt, type->moral  ); // this should be done differently

        attentionDir.normalize();
        rot.set( attentionDir );
        vel.mul( 0.5 );

        force.add( willForce );
        move_PointBody2D( dt );
    }


    double skill_ramp( double att, double def ){ return Treshold::r2( (att-def)/(att+def) ); }
    double armor_ramp( double att, double def ){ return Treshold::r2( (att-def)/(att+def) ); }

    void getDamage( double damage ){
        //moral       -= damage;
        //stamina     -= damage;
        double xdmg = damage/type->damage_tolerance;
        float dice = 1-pow( randf(), xdmg );
        float s = 0;
        int i;
        for( i=0; i<8; i++ ){
            s+= injury_chances[i];
            if( s > dice ){
                impair_mask |=  1<<i;
                break;
            }
        }
        printf( "damage %f xdmg %f dice %f i %i impair_mask %i \n", damage, xdmg, dice, i, impair_mask );
    }

    bool getHit( double bare_damage, double penetration, double skill, const Vec2d& dir ){
        //if( rot.dot(dir) )
        //type->
        /*
        // FIXME :  this require some DEBUG
        double c_dir  = rot.dot( dir );
        if( time_buf > 0 ){  // defence action
            double skill_tresh = skill_ramp( skill, c_dir*type->defence_skill );
            double dice        = randf();
            if ( dice > skill_tresh ){
                time_buf -= type->defence_period;
                return false;
            };
        }
        double c_armor = armor_ramp( penetration, type->armorFw );
        double damage  = c_armor * bare_damage;
        getDamage( damage );
        */
        getDamage( bare_damage );
        return true;
    }

    bool getShot( double bare_damage, double penetration, const Vec2d& dir ){
        //if( rot.dot(dir) )
        //type->
        double c_dir        = rot.dot( dir );
        if( c_dir > type->shield_cos ){
            shield -= bare_damage;
            return false;
        }
        double c_armor = armor_ramp( penetration, type->armorFw );
        double damage  = c_armor * bare_damage;
        return true;
    }

    void attack_melee( Soldier * enemy ){
        time_buf -= type->melee_period;
        stamina  *= type->melee_fStamina;
        enemy->getHit( type->melee_damage, type->melee_penetration, type->attack_skill, rot );
    }

/*
    void check_opponent( Soldier * opponent_, const Vec2d& d, double r2 ){
        double score = ( 1+rot.dot( d ) )/( 1 + r2 );
        if( score > opponent_score ) opponent = opponent_;
    }
*/

    bool friend_interaction( Soldier * sj ){
        Vec2d d;
        d.set_sub( pos, sj->pos );
        double r2 = d.norm2( );
        double R  = type->radius + sj->type->radius;
        //double R = 1.0;
        double R2 = R*R;
        if( r2 > R2 ) return false;
        d.mul( (1-(r2/R2))*body_push );
        force.add( d );
        sj->force.sub( d );
    }

    bool enemy_interaction( Soldier * sj, bool melee ){
        Vec2d d;
        d.set_sub( pos, sj->pos );
        double r2  = d.norm2( );
        double Ri  = type->radius;
        double Rj  = sj->type->radius;
        if( melee ){
            double Rm  = type->melee_range;
            if( r2 > sq(Rm+Rj) ) return false;
            double r     = sqrt( r2 );
            double invr  = 1/(r+1e-8);
            d.mul( invr );
            double c_dot = rot.dot( d );
            double c_dir = c_dot * c_dot;

            // check enemy
            //attentionDir.add_mul( d, -1/(1+r2) );
            double score = c_dir * invr;
            if( score > opponent_score ){
                opponent = sj;
                attentionDir.set( d*-1 );
            }

            // push force
            double R  = Rj + c_dir*Rm + (1-c_dir)*Ri;
            if( r > (R*R) ) return false;
            d.mul( type->melee_push );
        }else{ // not melee
            double R     = Ri+Rj;
            //double R     = 1.0;
            double R2 = R*R;
            if( r2 > R2 ) return false;
            d.mul( (1-(r2/R2))*body_push );
        }
        // apply force
        force.add    ( d );
        sj->force.sub( d );
        return true;
    }

    void setType( SoldierType* type_ ){
        type = type_;
        //stamina = type->stamina;
        //moral   = type->moral;
        stamina  = 1.0;
        moral    = 1.0;
        shield   = type->shield_endurance;
        time_buf = 0.0;
    }

    Soldier(){};
    Soldier( SoldierType* type_, Formation* formation_  ){
        setType( type_ );
        formation_ = formation_;
    }

};

#endif
