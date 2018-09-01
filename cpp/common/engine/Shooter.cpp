
#include <SDL2/SDL_opengl.h>
//#include "Draw3D.h"
//#include "Solids.h"

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "raytrace.h"
#include "geom3D.h"

#include "Object3D.h"
#include "Terrain25D.h"
#include "Warrior3D.h"
#include "Warrior25D.h"
#include "Projectile3D.h"

#include "testUtils.h"

#include "Shooter.h" // THE HEADER

Projectile3D* Shooter::fireProjectile( Warrior3D * w, double speed, int kind ){
    Projectile3D * p = new Projectile3D();
    //RigidBody *b = (RigidBody*)p; // FIXME: This is very bad but get arround diamond problem
    RigidBody *b = w->asRigidBody();
    p->kind = kind;
    p->id   = shotsCount; shotsCount++;
    p->vel.set_mul ( w->gun_rot, speed );
    p->vel.add     ( b->vel );
    //p->vel.add( { randf(-0.1,0.1), randf(-0.1,0.1), randf(-0.1,0.1) } );
    p->pos.set     ( b->pos );
    p->pos.add_mul ( w->gun_rot, 5.0 );
    projectiles.push_back( p );
    return p;
};

int Shooter::registrWarrior( Warrior3D * w ){
    w->id = warriorCount;
    warriorCount++;
    warriors.push_back( w );
    return warriorCount;
}

/*
Warrior3D* Shooter::makeWarrior( const Vec3d& pos, const Vec3d& dir, const Vec3d& up, int kind ){
    //int ith = warriors.size();
    //printf( " >>> Setup  ship %i \n", ith );
    Warrior3D* w = new Warrior3D();
    w->kind = kind;
    //w->id = warriorCount; warriorCount++;
    //w->loadFromFile( filename );
    //w->from_mass_points( 2, mass, (Vec2d*)poss );
    w->setPose( pos, dir, up );
    //w->initAllGuns( 6 );
    //printf( "DEBUG 2 \n" );
    //ship->world = this;
    //ship->collisionShape = collisionShape;
    //ship->name = new char[7];
    //sprintf( ship->name, "Ship_%02i", ith );
    //printf( "DEBUG 4 \n" );
    //warriors.push_back( w );
    //printf( "DEBUG 5 \n" );
    registrWarrior( w );
    return w;
}
*/

void Shooter::update_warriors3D(){
    auto itw = warriors.begin();
    while( itw != warriors.end() ) {
        Warrior3D * w = *itw;
        /*

        w->clean_temp( );
        //addEnviroForces              ( w->pos, w->vel, w->force,  w->landed );
        w->force.add( {0.0,gravity*w->mass,0.0} );
        if(terrain){ w->interact(terrain); }
        //w->landed = collideWithWorld ( w->pos, w->vel, w->surf );
        w->move( dt );
        */

        Vec3d G = (Vec3d){0.0,gravity,0.0};
        w->move_warrior( dt, wind_speed, G, terrain );

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
    //printf( "update_projectiles3D %i\n", projectiles.size() );
    auto it_proj = projectiles.begin();
    while( it_proj != projectiles.end() ) {
        Projectile3D * proj = *it_proj;

        //proj ->update_old_pos(    );
        //proj ->evalForce     (    );
        proj ->clean_temp();
        proj ->force.add( {0.0,gravity,0.0});
        proj ->update_Projectile3D( dt );

        //printf( "prj (%3.3f,%3.3f,%3.3f)  (%3.3f,%3.3f,%3.3f) %f\n", proj->pos.x, proj->pos.y, proj->pos.z,    proj->old_pos.x, proj->old_pos.y, proj->old_pos.z,   proj->time );

        Vec3d hRay,normal;
        hRay.set_sub( proj->pos, proj->old_pos );
        double tmax = hRay.normalize();

        bool hitted = false;
        //hitted |= proj->check_hit_ground( );
        //hitted |= proj->check_hit_vector<Frigate2D>( warriors );
        for(Object3d* o : objects ){
            o->ray( proj->old_pos, hRay, &normal );
            if (hitted) break;
        }

        if( proj->pos.y < ground_height ){ // TODO - real terrain height
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

void Shooter::burstObjectCollision( Object3d& obj, Burst3d& burst ){
    //printf( "burstObjectCollision b: %i o: %i ", burst.id, obj.id );
    if( sq( obj.R)> burst.bbox.dist2_Cilinder( obj.pos ) ){
        //int np = burst->shots.size();
        //for( int i=0; i<np; i++ ){
        //printf( " hit " );
        //printf( " hit (b:%i,o:%i) ", burst.id, obj.id );
        int j = 0;
        auto ip = burst.shots.begin();
        while( ip != burst.shots.end() ){
        // FIXME : to erase particles we need to use iterator
            //Particle3d& p = burst->shots[i];
            opCount_ShotObject++;
            Particle3d& p = *ip;
            Vec3d op; p.getOldPos(dt,op);
            if( obj.getShot( op, p.pos, *burst.type, dt ) ){
                //printf( "-" );
                burst.hit(j);
                //ip = burst.shots.erase( ip );
            };
            ++ip;
            j++;
        }
        //printf( "\n" );
    };
    //printf( "\n" );
}

void Shooter::update_bursts(){
    //printf( "update_bursts n=%i\n", bursts.size() );
    Vec3d vGrav = (Vec3d){0.0,gravity,0.0};
    auto it = bursts.begin();
    stopWatch.start();

    opCount_ShotObject = 0;
    while( it != bursts.end() ) {
        Burst3d * burst = *it;

        double airDensity = 1.27; // TODO : read from atmosphere
        burst->move( dt, vGrav, airDensity );

        for( Object3d* o : objects ){ // TODO : this is BruteForce - make optimized BroadPhase
            burstObjectCollision( *o, *burst );
        }

        /*
        for( Vehicle3d* vh : vehicles ){ // TODO : this is BruteForce - make optimized BroadPhase
            burstObjectCollision( *vh, *burst );
        }
        */
        ++it;
    }
    stopWatch.stop();
    if(debugFlags&DEBUG_FLAG_SHOTS)printf( " update_bursts  %g [Mticks] nbursts %i objs %i nbboxOps %i nShotOps %i \n", stopWatch.T*1e-6, bursts.size(), objects.size(), bursts.size()*objects.size(), opCount_ShotObject );

}

void Shooter::update_world( ){
	for( int i=0; i<perFrame; i++ ){
//        phi += omega * dt;
//        rot.fromAngle( phi );
        for( AnyControler* c : controlers ) c->update(dt);
        update_warriors3D ();
        update_warriors25D();
        update_projectiles3D();
        update_bursts();
        time += dt;
    }
};


void Shooter::init_world  ( ){ debug_shit=23; };
