
#include <SDL2/SDL_opengl.h>
#include "Draw3D.h"
#include "Solids.h"

#include "Shooter.h" // THE HEADER

Projectile3D* Shooter::fireProjectile( Warrior3D * w, double speed, int kind ){
    Projectile3D * p = new Projectile3D();
    p->kind = kind;
    p->id   = shotsCount; shotsCount++;
    p->vel.set_mul ( w->gun_rot, speed );
    p->vel.add     ( w->vel );
    //p->vel.add( { randf(-0.1,0.1), randf(-0.1,0.1), randf(-0.1,0.1) } );
    p->pos.set     ( w->pos );
    p->pos.add_mul ( w->gun_rot, 5.0 );
    projectiles.push_back( p );
    return p;
};

Warrior3D* Shooter::makeWarrior( const Vec3d& pos, const Vec3d& dir, const Vec3d& up, int kind ){
    //int ith = warriors.size();
    //printf( " >>> Setup  ship %i \n", ith );
    Warrior3D* w = new Warrior3D();
    w->kind = kind; w->id = warriorCount; warriorCount++;
    //w->loadFromFile( filename );
    //w->from_mass_points( 2, mass, (Vec2d*)poss );
    w->initOne();

    //printf( " I invI  %f %f \n", ship1->I, ship1->invI );
    //w->setDefaults();
    //w->setAngle( angle );
    w->pos.set ( pos  );
    //w->omega = 0.0;

    w->rotMat.a.set      ( dir     );
    w->rotMat.b.set      ( up      );
    w->rotMat.c.set_cross( dir, up );
    w->qrot.fromMatrix   ( w->rotMat );

    //w->initAllGuns( 6 );
    //printf( "DEBUG 2 \n" );
    //ship->world = this;
    //ship->collisionShape = collisionShape;
    //ship->name = new char[7];
    //sprintf( ship->name, "Ship_%02i", ith );
    //printf( "DEBUG 4 \n" );
    warriors.push_back( w );
    //printf( "DEBUG 5 \n" );
    return w;
}



void Shooter::update_warriors3D(){
    auto itw = warriors.begin();
    while( itw != warriors.end() ) {
        Warrior3D * w = *itw;

        w->clean_temp( );
        //addEnviroForces              ( w->pos, w->vel, w->force,  w->landed );
        w->force.add( {0.0,gravity,0.0} );
        if(terrain){
            Vec2d dv;
            double h  = terrain->eval( {w->pos.x,w->pos.z}, dv );
            double dh = w->pos.y - w->hground - h;
            if( dh  < 0.0 ){
                //w->force.add( {0.0,gravity,0.0} );
                //w->force.add( {dv.x, dh*(-1-0.5*w->vel.y), dv.y} );
                w->force.add( { -dv.x, dh*(-100+2.5*w->vel.y), -dv.y } );
                w->force.add( {w->vel.x*landDrag,0.0,w->vel.z*landDrag} );
                //printf( " dv (%3.3f,%3.3f) (%3.3f,%3.3f)  \n", dv.x, dv.y, w->pos.x,w->pos.y );
                //w->force.add( {0, dh*(-100+2.5*w->vel.y), 0} );

            }
        }
        //w->landed = collideWithWorld ( w->pos, w->vel, w->surf );
        w->move( dt );

        //w->gun_rot.set_mul_cmplx( rot, w->rot );
        //w->update( dt );

        //printf( " warriro %i pos (%3.3f,%3.3f) vel (%3.3f,%3.3f) force (%3.3f,%3.3f) \n", itw, w->pos.x, w->pos.y, w->vel.x, w->vel.y, w->force.x, w->force.y );

        //printf( " trigger %i until_reaload %f \n ", w->trigger, w->until_reaload );

        ++itw;
    }
}

void Shooter::update_warriors25D(){
    auto itw = warriors25D.begin();
    while( itw != warriors25D.end() ) {
        Warrior25D * w = *itw;

        //w->clean_temp( );
        //w->force.add( {0.0,gravity,0.0} );
        //w->move( dt );

        w->update( dt, wind_speed, watter_speed );

        if(w->trigger) w->shoot(projectiles);

        ++itw;
    }

}

void Shooter::update_projectiles3D(){
    printf( "update_projectiles3D %i\n", projectiles.size() );
        auto it_proj = projectiles.begin();
        while( it_proj != projectiles.end() ) {
            Projectile3D * proj = *it_proj;

            //proj ->update_old_pos(    );
            //proj ->evalForce     (    );
            proj ->clean_temp();
            proj ->force.add( {0.0,gravity,0.0});
            proj ->update_Projectile3D( dt );

            printf( "prj (%3.3f,%3.3f,%3.3f)  (%3.3f,%3.3f,%3.3f) %f\n", proj->pos.x, proj->pos.y, proj->pos.z,    proj->old_pos.x, proj->old_pos.y, proj->old_pos.z,   proj->time );

            Vec3d hRay,normal;
            hRay.set_sub( proj->pos, proj->old_pos );
            double tmax = hRay.normalize();

            bool hitted = false;
            //hitted |= proj->check_hit_ground( );
            //hitted |= proj->check_hit_vector<Frigate2D>( warriors );
            for( o : objects ){
                o->ray( proj->old_pos, hRay, &normal );
                if (hitted) break;
            }

            if( proj->pos.y < ground_height ){
                hitted = true;
            }

            //if( hitted ){
            if( hitted || (proj->time > projLifetime) ){
                it_proj = projectiles.erase( it_proj );
                delete proj;
            }else{
                ++it_proj;
            }
        }
}

void Shooter::update_world( ){
	for( int i=0; i<perFrame; i++ ){
//        phi += omega * dt;
//        rot.fromAngle( phi );
        update_warriors3D();
        update_warriors25D();
        update_projectiles3D();
    }
};


void Shooter::init_world  ( ){ debug_shit=23; };
