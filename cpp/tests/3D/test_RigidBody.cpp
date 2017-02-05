
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
#include "SoftBody.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"


// ============= Application

void drawTruss( SoftBody * truss, float vsc, float fsc, bool DEBUG ){
    glBegin( GL_LINES );
    glColor3f(0.0f,0.0f,0.0f);
    for( int i=0; i<truss->nbonds; i++ ){
        Bond& bond = truss->bonds[i];
        Vec3d& pi  = truss->points[bond.i];
        Vec3d& pj  = truss->points[bond.j];
        glVertex3f( (float)pi.x, (float)pi.y, (float)pi.z );
        glVertex3f( (float)pj.x, (float)pj.y, (float)pj.z );
        if( DEBUG ){
            printf( " %i  %i %i   (%3.3f,%3.3f,%3.3f)  (%3.3f,%3.3f,%3.3f)\n",  i,   bond.i, bond.j,  pi.x, pi.y, pi.z,  pj.x, pj.y, pj.z  );
        }
    }
    if(vsc>0.0){
        glColor3f(0.0f,0.0f,1.0f);
        for( int i=0; i<truss->npoints; i++ ){
            Vec3d p,v;
            p = truss->points[i]; v = truss->velocities[i];
            glVertex3f( (float)p.x, (float)p.y, (float)p.z );
            glVertex3f( (float)(p.x+vsc*v.x), (float)(p.y+vsc*v.y), (float)(p.z+vsc*v.z) );
        }
    }
    if(fsc>0.0){
        glColor3f(1.0f,0.0f,0.0f);
        for( int i=0; i<truss->npoints; i++ ){
            Vec3d p,f;
            p = truss->points[i]; f = truss->forces[i];
            glVertex3f( (float)p.x, (float)p.y, (float)p.z );
            glVertex3f( (float)(p.x+fsc*f.x), (float)(p.y+fsc*f.y), (float)(p.z+fsc*f.z) );
        }
    }
    for( int i=0; i<truss->nfix; i++ ){
        Vec3d p = truss->points[truss->fix[i]];
        Draw3D::drawPointCross( p, 0.3 );
        if( DEBUG ){
            printf( " fix %i %i (%3.3f,%3.3f,%3.3f) \n",  i,   truss->fix[i], p.x, p.y, p.z  );
        }
    }
    glEnd();
}


void drawRigidBody( const RigidBody& rb, int npoints, Vec3d* points ){
    glColor3f(0.0f,0.0f,0.0f);

    glPushMatrix();
    float glMat[16];
    Mat3d rotT;
    rotT.setT(rb.rotMat);
    Draw3D::toGLMat( rb.pos, rotT, glMat );
    glMultMatrixf( glMat );
    if( points ) for(int i=0; i<npoints; i++) Draw3D::drawPointCross(points[i],0.3);
    glPopMatrix();

    glColor3f(1.0f,0.0f,1.0f); Draw3D::drawVecInPos( rb.L,     {0.0,0.0,0.0} );
    glColor3f(0.0f,0.0f,1.0f); Draw3D::drawVecInPos( rb.omega, {0.0,0.0,0.0} );

    if( points ){
        for(int i=0; i<npoints; i++){
            Vec3d p,v;
            //rb.rotMat.dot_to_T(points[i], p);
            //v.set_cross( p, rb.omega );        // omega is in global coordinates
            rb.velOfPoint( points[i], v, p );
            p.add(rb.pos);
            Draw3D::drawVecInPos  ( v, p );
            Draw3D::drawPointCross( p, 0.2 );
        }
    };
};

// ==========================
// TestAppRigidBody
// ==========================

class TestAppRigidBody : public AppSDL2OGL_3D {
	public:

    std::vector<Object3D*>         objects;
    std::vector<SpringConstrain*>  springs;

    SoftBody  truss;

    RigidBody rb;

    static const int npoints =7;
    static const int nbonds  =18;
    static const int nfix    =1;
    double points[npoints*3] = {-8,0,0, 8,0,0, 0,-4,0, 0,4,0, 0,0,-6, 0,0,6, 0,0,0};
    double masses[npoints*3] = { 1,     1,     1,      1,     1,      1,     1    };
    int    fixes [nfix]      = {1};
    int    bips  [nbonds*2]  = {0,2, 0,3, 0,4, 0,5,  1,2, 1,3, 1,4, 1,5,  2,4, 2,5,  3,4, 3,5,    6,0, 6,1, 6,2, 6,3, 6,4, 6,5  };
    int    bits  [nbonds  ]  = {0,   0,   0,   0,    0,   0,   0,   0,    0,   0,    0,   0,      0,   0,   0,   0,   0,   0    };

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppRigidBody( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppRigidBody::TestAppRigidBody( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    // ------------ SoftBody Test
    BondType * btyp = new BondType();
    btyp->id     = 1;
    btyp->kPress = 100;
    btyp->kTens  = 100;
    btyp->linearDensity = 0.0;
    btyp->sPress = 1e+8;
    btyp->sTens  = 1e+8;

    truss.setPoints    ( npoints, (Vec3d*)points, masses, NULL );
    truss.setConstrains( nfix, fixes                           );
    btyp = new BondType(); *btyp = default_BondType;
    truss.setBonds     ( nbonds, bips, bits, btyp              );
    truss.prepareBonds (  true );
    truss.preparePoints(  true, 0.0d, 1.0d );

    //rb.from_mass_points( npoints, masses, (Vec3d*)points );
    //rb.qrot.setOne();
    //rb.L.set(100.0,0.0,0.0);
    //rb.qrot.set(0.1,0.2,0.3,0.2); rb.qrot.normalize();

    // TODO : test fly-wheel  https://en.wikipedia.org/wiki/Rigid_body_dynamics#/media/File:Gyroscope_precession.gif

    rb.initOne();
    //rb.L.set(1.0,0.0,0.0);
    //rb.L.set(1.0,1.8,1.0);
    srand(10454);
    rb.L.set(randf(-1.0,1.0),randf(-1.0,1.0),randf(-1.0,1.0));
    printf("L (%3.3f,%3.3f,%3.3f)\n",rb.L.x,rb.L.y,rb.L.z);

    //camMat.c.set(rb.L);                     camMat.c.normalize();
    //camMat.a.set_cross(camMat.b,camMat.c);  camMat.a.normalize();
    //camMat.b.set_cross(camMat.a,camMat.c);  camMat.b.normalize();
    //qCamera.fromMatrix(camMat); // this not work ?


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
    o->initOne();
    //o->bounds.pos.set(-1.0,4.0,-2.0);
    o->controler = new RigidBody();
    o->controler->initOne();
    //o->controler->initSpherical( 1.0, 2.0 );
    o->controler->pos.set(5.0,0.0,0.0);
    MeshCollisionShape * coll = new MeshCollisionShape();
    coll->mesh = mesh;
    o->coll = coll;
    o->shape = mesh->rendered_shape;
    objects.push_back(o);

    SpringConstrain * sp;
    sp = new SpringConstrain( 40.0, 0.0, 5.0, o->controler, NULL, mesh->points[0], {4.0,3.0,0.0} ); springs.push_back(sp);
    sp = new SpringConstrain( 40.0, 0.0, 5.0, o->controler, NULL, mesh->points[2], {6.0,3.0,0.0} ); springs.push_back(sp);
    //sp = new SpringConstrain( 10.0, 0.0, 5.0, o->controler, NULL, mesh->points[3], o->controler->pos ); springs.push_back(sp);

}

void TestAppRigidBody::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    int perFrame = 1;
	double dt    = 0.01;

	// ---- soft body test
    //truss.dt   = dt;
    //truss.damp = 0.0;
    //for(int itr=0; itr<perFrame; itr++){ truss.step(); }
    //drawTruss( &truss, 1.0, 1.0, false );

    // ---- RidigBody basics - drag
    Vec3d gdp,gv,gf;
    rb.clean_temp();
    rb.velOfPoint ( {2.0,5.0,6.0}, gv, gdp );
    gf.set_mul(gv, -0.01 );
    rb.apply_force( gf, gdp );
    rb.move(dt);
    rb.pos.set(0.0);
    rb.vel.set(0.0);

    glPushMatrix();
    glTranslatef(-5.0,0.0,0.0);
    glScalef(0.6,0.6,0.6);
    drawRigidBody( rb, npoints, (Vec3d*)points );
    glColor3f(1.0f,1.0f,1.0f);
    Draw3D::drawPointCross( gdp, 0.2 );
    Draw3D::drawVecInPos  ( gf*1000.0, gdp );
    glPopMatrix();


    // ---- RigidBody springs
	for(int itr=0; itr<perFrame; itr++){
        for(SpringConstrain * sp : springs){
            sp->apply();
        }
        for( Object3D * o : objects ){
            RigidBody * rb = o->controler;
            if(rb){
                rb->vel.mul( 1-dt*0.2 );
                rb->L.  mul( 1-dt*0.2 );
                rb->apply_force({0,-9.81,0},{0.0,0.0,0.0});
                rb->move_RigidBody(dt);
                rb->clean_temp();
            }
        }
	}
	for( Object3D * o : objects ){
        if(o->controler){
            o->gpos         = o->controler->pos;
            //o->bounds.orientation = o->controler->rotMat;
            o->grot.setT(o->controler->rotMat);
        };
	}
    glEnable( GL_LIGHTING );
    glEnable(GL_DEPTH_TEST);
    for( Object3D * o : objects ){
        if (o->shape){
            float glMat[16];
            glPushMatrix();
            Draw3D::toGLMat( o->gpos, o->grot, o->span, glMat );
            glMultMatrixf( glMat );
            glCallList( o->shape );
            glPopMatrix();
        }
    }
    glDisable( GL_LIGHTING );
    for(SpringConstrain * sp : springs ){
        Vec3d gp1,gp2,f;
        sp->getPoints( gp1, gp2 );
        f = sp->getForce( gp1, gp2 );
        glColor3f(1.0f,1.0f,1.0f); Draw3D::drawLine( gp1, gp2 );
        //glColor3f(1.0f,0.0f,0.0f); Draw3D::drawVecInPos( f, gp1 );
    }
    Draw3D::drawAxis(1.0);


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
















