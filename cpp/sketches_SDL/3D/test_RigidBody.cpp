
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

void testApprox_Vec3Normalize(){
    Vec3d v;
    v.fromRandomSphereSample();
    double d = 1e-12;
    for(int i=0; i<12; i++){
        Vec3d v_;
        v_ = v*(1+d);
        printf("%e", v_.norm()-1 );
        v_.normalize_taylor3();
        printf("-> %e \n", v_.norm()-1 );
        d*=10;
    }
}

void testApprox_Vec3Rotate(){
    Vec3d v;
    Vec3d omg;
    v.  fromRandomSphereSample();
    omg.fromRandomSphereSample();
    double d = 1e-12;
    glColor3f(1.0,1.0,1.0); Draw3D::drawVec( omg  );
    glColor3f(0.0,0.0,0.0); Draw3D::drawVec( v  );
    printf("=== testApprox_Vec3Rotate\n");
    for(int i=0; i<20; i++){
        d = 0.05*i;    
        Vec3d omg_;
        Vec3d v1=v,v2=v;
        double angle = d;
        omg_ = omg*angle;
        v1.rotate       (angle,omg);
        //v2.drotate_omega(omg_);
        //v2.drotate_omega2(omg_);
        v2.drotate_omega6(omg_);
        printf("alfa %g err %e \n", angle, (v1-v2).norm() );
        glColor3f(0.0,0.0,1.0); Draw3D::drawVec( v1 );
        glColor3f(1.0,0.0,0.0); Draw3D::drawVec( v2 );
        //d*=10;
    }
    exit(0);
}

void testApprox_Mat3Orthonormalize(){
    Quat4d  qrot; qrot.fromUniformS3( {randf(),randf(),randf()} );
    printf( "qrot %g %g %g %g\n", qrot.x, qrot.y, qrot.z, qrot.w );
    Mat3d   mrot; qrot.toMatrix(mrot);
    printf("mrot: "); mrot.printOrthoErr();
    double d = 1e-6;    
    for(int i=0; i<30; i++){
        Mat3d   m = mrot;
        for(int j=0;j<9;j++){ m.array[j]+=randf(-d,d); }
        m.orthogonalize_taylor3(2,1,0);
        printf("m(%e): ", d); m.printOrthoErr();
        d*=2;
    }
    exit(0);
}

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
	virtual void keyStateHandling( const Uint8 *keys );

	TestAppRigidBody( int& id, int WIDTH_, int HEIGHT_ );

};

void TestAppRigidBody::keyStateHandling( const Uint8 *keys ){
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

    //rb.initOne();
    //rb.L.set(1.0,0.0,0.0);
    //rb.L.set(1.0,1.8,1.0);
    //srand(10454);
    //rb.L.set(randf(-1.0,1.0),randf(-1.0,1.0),randf(-1.0,1.0));
    //rb.L.set(0.0);
    //printf("L (%3.3f,%3.3f,%3.3f)\n",rb.L.x,rb.L.y,rb.L.z);

    //rb.invIbody.a.set(1.0,0.0,0.0 );
    //rb.invIbody.b.set(0.0,0.5,0.0 );
    //rb.invIbody.c.set(0.0,0.0,0.25);
    
    rb.setInertia_box( 1.0, {1.0,0.5,0.25} ); 
    rb.rotMat.setOne();

    //camMat.c.set(rb.L);                     camMat.c.normalize();
    //camMat.a.set_cross(camMat.b,camMat.c);  camMat.a.normalize();
    //camMat.b.set_cross(camMat.a,camMat.c);  camMat.b.normalize();
    //qCamera.fromMatrix(camMat); // this not work ?


    OMesh * mesh = new OMesh();
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
    o->controler->setInertia_box( 1.0, {1.0,1.0,1.0} );
    //o->controler->initSpherical( 1.0, 2.0 );
    o->controler->pos.set(5.0,0.0,0.0);
    MeshCollisionShape * coll = new MeshCollisionShape();
    coll->mesh = mesh;
    o->coll = coll;
    o->shape = mesh->rendered_shape;
    objects.push_back(o);

    SpringConstrain * sp;
    sp = new SpringConstrain( 40.0, 0.0, 5.0, o->controler, NULL, mesh->points[0], {4.0,3.0,0.0} ); springs.push_back(sp);
    sp = new SpringConstrain( 40.0, 0.0, 5.0, o->controler, NULL, mesh->points[1], {6.0,3.0,0.0} ); springs.push_back(sp);
    //sp = new SpringConstrain( 10.0, 0.0, 5.0, o->controler, NULL, mesh->points[3], o->controler->pos ); springs.push_back(sp);

}

void TestAppRigidBody::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    int perFrame = 1;
	double dt    = 0.003;

    srand(54564564);    
    //testApprox_Vec3Normalize( ); exit(0);
    //testApprox_Vec3Rotate(); return;
    //testApprox_Mat3Orthonormalize();
    
	// ---- soft body test
    //truss.dt   = dt;
    //truss.damp = 0.0;
    //for(int itr=0; itr<perFrame; itr++){ truss.step(); }
    //drawTruss( &truss, 1.0, 1.0, false );

    // ---- RidigBody basics - drag
    /*
    Vec3d gdp,gv,gf;
    rb.clean_temp();
    rb.velOfPoint ( {2.0,5.0,6.0}, gv, gdp );
    gf.set_mul(gv, -0.01 );
    rb.apply_force( gf, gdp );
    rb.move(dt);
    rb.pos.set(0.0);
    rb.vel.set(0.0);
    */

    const Uint8 *keys = SDL_GetKeyboardState(NULL);
    rb.clean_temp();
    float step = 10.2;
    //rb.L.mul( 1-0.01 );
    if( keys[ SDL_SCANCODE_W  ] ){ rb.torq.add_mul(rb.rotMat.a,  step ); }
	if( keys[ SDL_SCANCODE_S  ] ){ rb.torq.add_mul(rb.rotMat.a, -step ); }
	if( keys[ SDL_SCANCODE_A  ] ){ rb.torq.add_mul(rb.rotMat.b,  step ); }
	if( keys[ SDL_SCANCODE_D  ] ){ rb.torq.add_mul(rb.rotMat.b, -step ); }
    if( keys[ SDL_SCANCODE_Q  ] ){ rb.torq.add_mul(rb.rotMat.c,  step ); }
	if( keys[ SDL_SCANCODE_E  ] ){ rb.torq.add_mul(rb.rotMat.c, -step ); }
    rb.move(dt);
    rb.pos.set(0.0);
    rb.vel.set(0.0);
    glColor3f(0.0,0.0,0.0);
    Draw3D::drawTriclinicBoxT( rb.rotMat, {-1.0,-0.5,-0.25}, {1.0,0.5,0.25} );
    Draw3D::drawMatInPos     ( rb.rotMat, rb.pos );
    glColor3f(1.0,1.0,1.0); Draw3D::drawVec( rb.L );
    
    //printf( " omega %g L %g I %g \n", rb.omega.norm(), rb.L.norm(),   rb.L.norm()/rb.omega.norm()  );


    /*
    glPushMatrix();
        glTranslatef(-5.0,0.0,0.0);
        glScalef(0.6,0.6,0.6);
        drawRigidBody( rb, npoints, (Vec3d*)points );
        glColor3f(1.0f,1.0f,1.0f);
        Draw3D::drawPointCross( gdp, 0.2 );
        Draw3D::drawVecInPos  ( gf*1000.0, gdp );
    glPopMatrix();
    */


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
            //o->grot.setT(o->controler->rotMat);
            o->grot.set(o->controler->rotMat);
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


    //Draw3D::drawAxis(1.0);


};


void TestAppRigidBody::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
                case SDLK_r:  rb.rotMat = Mat3dIdentity; rb.L=Vec3dZero;
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
















