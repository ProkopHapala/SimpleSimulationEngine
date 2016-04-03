
#include <SDL2/SDL_opengl.h>
#include "Draw3D.h"
#include "Solids.h"

#include "MultiFight3DWorld.h" // THE HEADER

void MultiFight3DWorld::fireProjectile( Warrior3D * w ){
    Projectile3D * p = new Projectile3D();
    p->vel.set_mul( w->gun_rot, 10.0 );
    p->vel.add( w->vel );
    p->vel.add( { randf(-0.1,0.1), randf(-0.1,0.1), randf(-0.1,0.1) } );
    p->pos.set( w->pos );
    p->pos.add_mul( w->gun_rot, 5.0 );
    projectiles.push_back( p );
};

Warrior3D* MultiFight3DWorld::makeWarrior( const Vec3d& pos, const Vec3d& dir, const Vec3d& up, int shape ){
    int ith = warriors.size();
    //printf( " >>> Setup  ship %i \n", ith );
    Warrior3D* w = new Warrior3D();
    //w->loadFromFile( filename );
    //w->from_mass_points( 2, mass, (Vec2d*)poss );
    w->setMass( 1.0 );
    w->Ibody.a.set(1,0,0);
    w->Ibody.b.set(0,1,0);
    w->Ibody.c.set(0,0,1);
    w->Ibody.invert_to( w->invIbody );

    //printf( " I invI  %f %f \n", ship1->I, ship1->invI );
    //w->setDefaults();
    //w->setAngle( angle );
    w->pos.set ( pos  );
    //w->omega = 0.0;

    w->rotMat.a.set      ( dir     );
    w->rotMat.b.set      ( up      );
    w->rotMat.c.set_cross( dir, up );
    w->qrot.fromMatrix   ( w->rotMat );

    w->shape = shape;

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


void MultiFight3DWorld::update_world( ){
	for( int i=0; i<perFrame; i++ ){
//        phi += omega * dt;
//        rot.fromAngle( phi );

        auto itw = warriors.begin();
        while( itw != warriors.end() ) {
            Warrior3D * w = *itw;

            w->clean_temp( );
            //addEnviroForces              ( w->pos, w->vel, w->force,  w->landed );
            //w->landed = collideWithWorld ( w->pos, w->vel, w->surf );
            w->move( dt );

            //w->gun_rot.set_mul_cmplx( rot, w->rot );
            //w->update( dt );

            //printf( " warriro %i pos (%3.3f,%3.3f) vel (%3.3f,%3.3f) force (%3.3f,%3.3f) \n", itw, w->pos.x, w->pos.y, w->vel.x, w->vel.y, w->force.x, w->force.y );

            //printf( " trigger %i until_reaload %f \n ", w->trigger, w->until_reaload );
            if( w->trigger && ( w->until_reaload <= 0.0 ) ){
                fireProjectile( w );
                w->until_reaload  = reload_time;
                printf( " fire !!! %i-th projectile % \n ", projectiles.size() );
            }
            w->until_reaload -= dt;

            ++itw;
        }

        auto it_proj = projectiles.begin();
        while( it_proj != projectiles.end() ) {
            Projectile3D * proj = *it_proj;

            //proj ->update_old_pos(    );
            //proj ->evalForce     (    );
            proj ->move          ( dt );

            bool hitted = false;
            //hitted |= proj->check_hit_ground( );
            //hitted |= proj->check_hit_vector<Frigate2D>( warriors );
            Vec3d normal;
            //hitted = collideWithWorld ( proj->pos, proj->vel, normal );
            for( o : objects ){
                hitted = collideWithObject( o->lpos, objR, proj->pos, proj->vel, normal );
                if (hitted) break;
            }
            //if( hitted ){
            if( proj->time > 8.0 ){
                it_proj = projectiles.erase( it_proj );
                delete proj;
            }else{
                ++it_proj;
            }
        }

    }
};


void MultiFight3DWorld::init_world  ( ){

    defaultWarriorShape = glGenLists(1);
    glNewList( defaultWarriorShape , GL_COMPILE );
    glBegin   (GL_TRIANGLE_FAN);
        glNormal3f( 0.0f, 0.0f, 1.0f );
        glVertex3f( +1.5,  0.0, 0 );
        glVertex3f( +0.5,  0.2, 0 );
        glVertex3f( -1.0,  0.2, 0 );
        glVertex3f( -1.0, -0.2, 0 );
        glVertex3f( +0.5, -0.2, 0 );
        glVertex3f( +1.5,  0.0, 0 );
    glEnd();
    //Draw2D::drawCircle_d( {0.0, 0.0}, 0.5, 16, true );
    glEndList();


    defaultObjectShape = glGenLists(1);
    glNewList( defaultObjectShape , GL_COMPILE );
        //glPushMatrix();

        //glDisable ( GL_LIGHTING );
        //Draw3D::drawAxis ( 3.0f );
        //glColor3f( 1.0f, 0.0f, 1.0f ); Draw3D::drawLines   ( Solids::Icosahedron_nedges, Solids::Icosahedron_edges, Solids::Icosahedron_verts                             );

        glEnable( GL_LIGHTING );
        glColor3f( 0.8f, 0.8f, 0.8f ); Draw3D::drawPolygons( Solids::Icosahedron_nfaces, Solids::Icosahedron_ngons, Solids::Icosahedron_faces,  Solids::Icosahedron_verts );

        //glEnable( GL_LIGHTING ); glColor3f( 0.8f, 0.8f, 0.8f );   Draw3D::drawSphere_oct( 3, 1.0, {0.0,0.0,0.0} );
        //glPopMatrix();
    glEndList();
    //objects = new KinematicBody();

    defaultObjectHitShape = glGenLists(1);
    glNewList( defaultObjectHitShape , GL_COMPILE );
        //glPushMatrix();
        glDisable ( GL_LIGHTING );
        Draw3D::drawAxis ( 3.0f );
        //glColor3f( 1.0f, 0.0f, 1.0f ); Draw3D::drawLines   ( Solids::Icosahedron_nedges, Solids::Icosahedron_edges, Solids::Icosahedron_verts                             );
        glColor3f( 0.8f, 0.0f, 0.8f ); Draw3D::drawSphereOctLines( 16, 2.0, {0.0,0.0,0.0} );
        //glPopMatrix();
    glEndList();
    //objects = new KinematicBody();

    defaultProjectileShape = glGenLists(1);
    glNewList( defaultProjectileShape , GL_COMPILE );
        //glPushMatrix();
        glDisable ( GL_LIGHTING );
        //glColor3f( 1.0f, 0.0f, 1.0f ); Draw3D::drawLines ( Solids::Tetrahedron_nedges, Solids::Tetrahedron_edges, Solids::Tetrahedron_verts );
        glScalef( 0.2f, 0.2f, 0.2f );
        glColor3f( 1.0f, 0.5f, 0.0f ); Draw3D::drawPolygons( Solids::Octahedron_nfaces, Solids::Octahedron_ngons, Solids::Octahedron_faces,  Solids::Octahedron_verts );
    glEndList();
    //objects = new KinematicBody();

    int nobjects = 100;
    float Lspan = 50.0;
    objects.reserve( nobjects );
    for( int i=0; i<nobjects; i++ ){
        KinematicBody* obj = new KinematicBody();
        obj->lpos.set( randf(-Lspan,Lspan), randf(-Lspan,Lspan), randf(-Lspan,Lspan) );
        obj->lrot.setOne();
        objects.push_back( obj );
    }

    //center.set(0.0d);

};
