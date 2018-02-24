
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw2D.h"

#include "LTUnit.h" // THE HEADER

void LTUnit::setType      ( LTUnitType* type_ ){};

void LTUnit::update       ( double dt ){
    suppressed -= type->recovery_rate * dt; if( suppressed<0.0 ) suppressed=0.0;
    dt         *= (1.0-suppressed);
    for(int i=0; i<type->nGun; i++ ){
        guns[i].timer  += dt*type->guns[i].rps;
    }
    //reload     += dt/type->fire_period;
};

void LTUnit::getShot( const Vec3d& from, int nburst, double area, double damage_const, double damage_kinetic, double penetration ){
    // evaluate expected damage
    //double  ca  = from.dot( dir );
    // double armor = 0.5* (1+ca) * (type->armorFw - type->armorBg )  + type->armorBg;
    double pass = ( penetration - type->armorFront )/penetration; // TODO: evaluate directionality
    printf( "pass %f \n", pass );
    if( pass < 0 ) return;
    double damage = pass*damage_kinetic + damage_const;

    // evaluate hit probability
    double hit_area  = type->sz.a*type->sz.b;  // TODO: evaluate directionality
    double hit_prob  =            hit_area/(hit_area+ area);
    double kill_prob = hit_prob * damage_ramp( damage, type->HP );
    //printf( "kill_prob %f   %f  (%f %f)   (%f %f) \n", kill_prob, hit_prob, damage, type->HP, area, hit_area );
    suppressed -= kill_prob;
    // TODO : carefullness to modulate suppression and damage ??
    for( int i = 0; i<nburst; i++ ){
        if( randf() > kill_prob ) continue;
        wound = 1;
        //n--;  // soft kill
        //if( randf() > type->heal_prob ){ alive = false; } // hard kill
    }
}


void LTUnit::render( uint32_t color, int iLOD ){
    Draw::setRGBA(color);
    switch(type->kind){
        case(UnitKind::inf):
            Draw2D::drawRotT      (pos, rot, {type->sz.a, type->sz.b} );
            break;
        case(UnitKind::gun):
            Draw2D::drawRotTriangle(pos, rot, {type->sz.a, type->sz.b});
            Draw2D::drawVecInPos_d (turret_dir*type->sz.a,pos);
            break;
        case(UnitKind::tank):
            Draw2D::drawRotRect   (pos, rot, {type->sz.a, type->sz.b});
            Draw2D::drawCircle_d  (pos,type->sz.b,16,false);
            Draw2D::drawVecInPos_d(turret_dir*type->sz.a,pos);
            break;
        case(UnitKind::stug):
            Draw2D::drawRotRect   (pos, rot, {type->sz.a, type->sz.b});
            Draw2D::drawVecInPos_d(rot*type->sz.a,pos);
            break;
    }
};
