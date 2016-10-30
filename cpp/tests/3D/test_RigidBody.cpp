
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw3D.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"
#include "Mesh.h"
#include "Body.h"
#include "Object3D.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

// ============= Application

class TestAppRigidBody : public AppSDL2OGL_3D {
	public:

    std::vector<Object3D*>         objects;
    std::vector<SpringConstrain*>  springs;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppRigidBody( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppRigidBody::TestAppRigidBody( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    Mesh * mesh = new Mesh();
    mesh->fromFileOBJ( "common_resources/turret.obj" );
    mesh->rendered_shape = glGenLists(1);
    glNewList( mesh->rendered_shape , GL_COMPILE );
        glEnable( GL_LIGHTING );
        glColor3f( 0.8f, 0.8f, 0.8f );
        Draw3D::drawMesh( *mesh );
    glEndList();

    Object3D * o;

    o = new Object3D();
    o->id = 1;
    o->bounds.initOne();
    //o->bounds.pos.set(-1.0,4.0,-2.0);
    o->controler = new RigidBody();
    o->controler->initOne();
    o->controler->pos.set(-2.0,4.0,1.0);
    MeshCollisionShape * coll = new MeshCollisionShape();
    coll->mesh = mesh;
    o->coll = coll;
    o->shape = mesh->rendered_shape;
    objects.push_back(o);


}

void TestAppRigidBody::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	int perFrame = 10;
	double dt    = 0.001;
	for(int itr=0; itr<perFrame; itr++){
        for( Object3D * o : objects ){
            RigidBody * rb = o->controler;
            if(rb){
                rb->clean_temp();
                rb->vel.mul( 1-dt );
                rb->L.  mul( 1-dt );
                rb->apply_force({0,-9.81,0},{0.0,0.0,0.0});
                Vec3d force, torq; force.set(0.0); torq.set(0.0);
                //o->coll->colideWithTerrain(terrain, rb->rotMat, rb->pos, force, torq );
                rb->force.add_mul( force,10.0 );
                rb->torq .add_mul( torq ,10.0 );
                rb->move_RigidBody(dt);
            }
        }
	}

	for( Object3D * o : objects ){
        if(o->controler){
            o->bounds.pos         = o->controler->pos;
            o->bounds.orientation = o->controler->rotMat;
        };
	}

	// ======== draw
    glEnable( GL_LIGHTING );
    glEnable(GL_DEPTH_TEST);
    //glCallList(terrain->shape);

    for( Object3D * o : objects ){
        if (o->shape){
            float glMat[16];
            glPushMatrix();
            //printf("------ %i \n", o->id );
            //printf("(%3.3f,%3.3f,%3.3f)\n", o->bounds.pos.x, o->bounds.pos.y, o->bounds.pos.z );
            //printf("(%3.3f,%3.3f,%3.3f)\n", o->bounds.orientation.a.x, o->bounds.orientation.a.y, o->bounds.orientation.a.z );
            //printf("(%3.3f,%3.3f,%3.3f)\n", o->bounds.orientation.b.x, o->bounds.orientation.b.y, o->bounds.orientation.b.z );
            //printf("(%3.3f,%3.3f,%3.3f)\n", o->bounds.orientation.c.x, o->bounds.orientation.c.y, o->bounds.orientation.c.z );
            Draw3D::toGLMat( o->bounds.pos, o->bounds.orientation, o->bounds.span, glMat );
            glMultMatrixf( glMat );
            glCallList( o->shape );
            glPopMatrix();
        }
    }

    glDisable( GL_LIGHTING );
    Draw3D::drawAxis(1.0);

    glDisable( GL_DEPTH_TEST );

    //glColor3f(0,0,0); Draw3D::drawPointCross( mesh.points[ipicked], 0.2 );

};


void TestAppRigidBody::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}

void TestAppRigidBody::drawHUD(){
    glDisable ( GL_LIGHTING );

}

// ===================== MAIN

TestAppRigidBody * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppRigidBody( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















