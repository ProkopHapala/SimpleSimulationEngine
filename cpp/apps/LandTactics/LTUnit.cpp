
#include "LTUnit.h" // THE HEADER

//class LTUnit : public RigidBody2D {

void LTUnit::move_to_goal( double dt ){
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

void LTUnit::fire_at_unit( LTUnit * target ){
    reload = 0.0;
    Vec2d d;
    d.set_sub( target->pos, pos );
    double dist = d.normalize();
    // this does not care about time-of-flight neither raytracing we may later emit projectile instead of this
    Draw2D::drawLine_d(pos,target->pos);
    target->getShot( d, type->fire_nburst, sq(type->fire_spread_aim * dist), type->fire_damage_const, type->fire_damage_kinetic, type->fire_penetration );
}

void LTUnit::update( double dt ){
    suppressed -= type->recovery_rate * dt; if( suppressed<0.0 ) suppressed=0.0;
    dt         *= (1.0-suppressed);
    reload     += dt/type->fire_period;
    switch(job){
        case Unit_JOB_GOTO:         move_to_goal( dt ); break;
        case Unit_JOB_FIRE_AT_UNIT: if(opponent==NULL){ job=default_job; break;} if( reload > 1.0 ) fire_at_unit( opponent ); break;
    }
}

double LTUnit::damage_ramp( double att, double def ){ return 0.5*(1+Treshold::r2( (att-def)/(att+def) )); }

void LTUnit::getShot( const Vec2d& from, int nburst, double area, double damage_const, double damage_kinetic, double penetration ){
    // evaluate expected damage
    //double  ca  = from.dot( dir );
    // double armor = 0.5* (1+ca) * (type->armorFw - type->armorBg )  + type->armorBg;
    double pass = ( penetration - type->armorFw )/penetration;
    printf( "pass %f \n", pass );
    if( pass < 0 ) return;
    double damage = pass*damage_kinetic + damage_const;

    // evaluate hit probability
    double hit_prob  =            type->hit_area/(type->hit_area + area);
    double kill_prob = hit_prob * damage_ramp( damage, type->damage_tolerance );
    printf( "kill_prob %f   %f  (%f %f)   (%f %f) \n", kill_prob, hit_prob, damage, type->damage_tolerance, area, type->hit_area );
    suppressed -= kill_prob;
    // TODO : carefullness to modulate suppression and damage ??
    for( int i = 0; i<nburst; i++ ){
        if( randf() > kill_prob ) continue;
        n--;  // soft kill
        if( randf() > type->heal_prob ){ ntot--; } // hard kill
    }
}

void LTUnit::setOpponent( LTUnit * opponent_ ){
    job = Unit_JOB_FIRE_AT_UNIT;
    opponent = opponent_;
};

void LTUnit::setGoal( const Vec2d& goal_ ){
    job = Unit_JOB_GOTO;
    goal.set( goal_ );
};

void LTUnit::setType( LTUnitType* type_ ){
    type = type_;
}

void LTUnit::render( uint32_t color ){
    //glColor3f( c.x, c.y, c.z );
    Draw::setRGBA(color);
    Draw2D::drawCircle_d( pos, radius, 16, false );
    Draw2D::drawVecInPos_d( attentionDir, pos );
    //printf( " %f %f \n", attentionDir.x, attentionDir.y );
    //printf( "render (%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) \n", pos.x, pos.y, pos.x, c.x, c.y, c.z );
    char str[8];
    sprintf(str,"%4i",n);
    //Draw2D::drawString( str, (float)pos.x, (float)pos.y, 0.4f, default_font_texture );
    Draw2D::drawText( str, 0, {pos.x,pos.y}, 0.0, default_font_texture, 2.0 );

}

void LTUnit::renderJob( uint32_t c){
    if(job == Unit_JOB_GOTO         ) Draw2D::drawLine_d( pos, goal );
    if(job == Unit_JOB_FIRE_AT_UNIT ) Draw2D::drawLine_d( pos, opponent->pos );
    //printf( "render (%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) \n", pos.x, pos.y, pos.x, c.x, c.y, c.z );
}

LTUnit::LTUnit( LTUnitType* type_, LTFaction* faction_, const Vec2d& pos_ ){
    pos.set(pos_);
    setType( type_ );
    faction   = faction_ ;
}

//}
