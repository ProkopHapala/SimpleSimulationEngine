#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw2D.h"

#include "LTSquad.h" // THE HEADER

//class LTSquad : public RigidBody2D {

void LTSquad::move_to_goal( double dt ){
    Vec2d d;
    d.set_sub( goal, pos );
    double dist = d.norm();
    double vdt = type->maxSpeed * dt;  // TODO actual speed is not max speed
    if( dist > vdt ){
        pos.add_mul( d, vdt/dist );
    }else{
        pos.set(goal);
        job = default_job;
    }
}

void LTSquad::fire_at_squad( LTSquad * target ){
    reload = 0.0;
    Vec2d d;
    d.set_sub( target->pos, pos );
    double dist = d.normalize();
    // this does not care about time-of-flight neither raytracing we may later emit projectile instead of this
    Draw2D::drawLine_d(pos,target->pos);
    LTGunType& gunType = *type->guns[0]; // TODO: which gun ?
    double crossection = 1.0; // TODO: get target crossection
    double dHeight     = 0.0; // TODO: height difference
    //target->getShot( d, gunType.nburst, crossection, gunType.dmg, gunType.getKineticDamage(dist,dHeight), gunType.AP );
    // TODO
    //  - distribute over units inside
}

void LTSquad::update( double dt ){
    //suppressed -= type->recovery_rate * dt; if( suppressed<0.0 ) suppressed=0.0;
    //dt         *= (1.0-suppressed);
    //reload     += dt/type->fire_period;
    switch(job){
        case Unit_JOB_GOTO:         move_to_goal( dt ); break;
        case Unit_JOB_FIRE_AT_UNIT: if(opponent==NULL){ job=default_job; break;} if( reload > 1.0 ) fire_at_squad( opponent ); break;
    }
    for( LTUnit u: units ){
        u.update( dt);
    };
}

/*
void LTSquad::getShot( const Vec2d& from, int nburst, double area, double damage_const, double damage_kinetic, double penetration ){
    // evaluate expected damage
    //double  ca  = from.dot( dir );
    // double armor = 0.5* (1+ca) * (type->armorFw - type->armorBg )  + type->armorBg;

    for( LTUnit u:units){ u.getShot(); }

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
*/

void LTSquad::setOpponent( LTSquad * opponent_ ){
    job = Unit_JOB_FIRE_AT_UNIT;
    opponent = opponent_;
};

void LTSquad::setGoal( const Vec2d& goal_ ){
    job = Unit_JOB_GOTO;
    goal.set( goal_ );
};

void LTSquad::setType( LTUnitType* type_ ){
    type = type_;
}

void LTSquad::makeUnits( int n ){
    units.reserve(n);
    for(int i=0; i<n; i++){
        units.push_back( LTUnit(type) );
    }
};

void LTSquad::render( uint32_t color, int iLOD ){
    //printf( "squad \n" );
    //printf( "squad pos (%f,%f) \n", pos.x, pos.y );
    //glColor3f( c.x, c.y, c.z );
    Draw::setRGBA(color);
    Draw2D::drawCircle_d( pos, radius, 16, false );
    Draw2D::drawVecInPos_d( attentionDir*radius*0.5, pos );
    Draw2D::drawVecInPos_d( rot*radius,              pos );
    //printf( " %f %f \n", attentionDir.x, attentionDir.y );
    //printf( "render (%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) \n", pos.x, pos.y, pos.x, c.x, c.y, c.z );
    char str[8];
    sprintf(str,"%4i",n);
    //Draw2D::drawString( str, (float)pos.x, (float)pos.y, 0.4f, default_font_texture );
    Draw2D::drawText( str, 0, {pos.x,pos.y}, 0.0, default_font_texture, 2.0 );

    if(iLOD>0){
        for( LTUnit u: units ){
            //printf( "unit \n"  );
            u.render( color, iLOD );
        }
    }

}

void LTSquad::populate(int n){
    units.reserve(n);
    double rs = radius*0.7;
    for(int i=0; i<n; i++){
        LTUnit u(type);
        u.pos.set( pos.x+randf(-rs,rs), pos.y+randf(-rs,rs) );
        u.rot = rot;
        units.push_back(u);
    }
}

void LTSquad::fromString(const char * str_, const UnitTypeDict& dct ){
    //printf("1 \n");
    char *str    = strdup(str_);
    //printf( "%s\n", str );
    char * token = strtok(str, ";"); //printf( ">>%s<<\n", token );

    std::string s=stripWhite(token);
    auto it = dct.find( s );
    if(it!=dct.end()){
        printf("found >>%s<<\n", s.c_str() );
        type = it->second;
        //printf( "type.name >>%s<< \n", type->name.c_str() );
    }else{
        printf("not found >>%s<<\n", s.c_str() );
        exit(-1);
    };

    token = strtok(NULL, ";"); //printf( ">>%s<<\n", token );
    sscanf( token, "%li", &ntot ); n=ntot;
    //printf( "n \n", n );

    token = strtok(NULL, ";"); //printf( ">>%s<<\n", token );
    sscanf( token, "%lf %lf %lf", &pos.x, &pos.y ); n=ntot;
    //printf( "pos (%f,%f,%f) \n",   pos.x,  pos.y );

    rot.set(1.0,0.0);
    //populate(n);
};

void LTSquad::renderJob( uint32_t c){
    if(job == Unit_JOB_GOTO         ) Draw2D::drawLine_d( pos, goal );
    if(job == Unit_JOB_FIRE_AT_UNIT ) Draw2D::drawLine_d( pos, opponent->pos );
    //printf( "render (%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) \n", pos.x, pos.y, pos.x, c.x, c.y, c.z );
}

LTSquad::LTSquad( LTUnitType* type_, LTFaction* faction_, const Vec2d& pos_ ){
    pos.set(pos_);
    rot.set(1.0,0.0);
    setType( type_ );
    faction = faction_ ;
}

//}
