
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
//#include "Body.h"


#include "Shooter.h"

#include "broadPhaseCollision.h" // Move to Sooter later

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"


/*
## TODO:
 - using Sooter
 - moving porejectiles along predefined paths (Lissajous curves?)
 - Projectiles ejected from user camera (fixed position or simple movement)
 - NxM problem (N bodies, M projectiles) each frame
## HOW TO:
 #### Naieve NxM
 #### Hashmap
   - large boxes, insert & remove only when leave the box
   - Test my version of hashmap vs STL
 #### Projectile Bunches
   - projectile fired from automatic canones are grouped to sommon bounding object
   - when projectiles are scattered ( e.g. when they hit something ) they are released from bunch

*/

double projLifetime = 10.0;
double airDensity = 1.27;

class Target : public Object3d{ public:
    int nhits=0;
    int nboxhits=0;
    //virtual bool getShot( const Vec3d& p0, const Vec3d& p1, const ProjectileType& prjType, double dt ){
    virtual bool getShot( const Vec3d& p0, const Vec3d& p1, const ProjectileType& prjType, double dt ) override {
        //printf("getShot \n");
        bool bHit = Object3d::getShot(p0,p1,prjType,dt); // TODO: does it doe virtual dispatch ?
        //bool bHit = hitLine_Object3d(p0,p1, prjType.caliber*0.5 );
        nboxhits++;
        if(bHit) nhits++;
        return bHit;
    }
    void draw(){
        glPushMatrix();
        glTranslatef( pos.x, pos.y, pos.z );
        glScalef(R,R,R);
        //glCallList( defaultObjectHitShape );
        float g = 0.0; if(nboxhits) g=1.0;
        if(nhits){ glColor3f(1.0,g,0.0); }else{ glColor3f(0.0,g,1.0); }
        glCallList( type->ogl );
        glPopMatrix();
    }
    Target( double R_, Vec3d pos_, ObjectType* type_, int id_ ){R=R_;pos=pos_;type=type_;id=id_;  };
};

void drawBurst( Burst3d& burst, int nseg, double dt ){
    glColor3f(1.0,0.0,0.0);
    int n = burst.shots.size();
    Vec3d op;
    for(int i=0; i<n; i++){
        Draw3D::drawPointCross( burst.shots[i].pos, 0.1 );
        //if( tmpPos ) Draw3D::drawLine( tmpPos[i], burst.shots[i].pos );
        burst.shots[i].getOldPos(dt,op);
        Draw3D::drawLine( op, burst.shots[i].pos );
    }
    //if( tmpPos ){ glColor3f(1.0,1.0,0.0); Draw3D::drawPointCross( tmpPos[n-1], 0.5 ); }
    burst.shots[0].getOldPos(dt,op);
    glColor3f(1.0,1.0,0.0); Draw3D::drawPointCross( op, 0.5 );
    glColor3f(0.0,1.0,1.0); Draw3D::drawPointCross( burst.shots[0].pos, 0.5 );
    glColor3f(0.0,0.5,0.0);
    Vec3d b = burst.bbox.p+burst.bbox.hdir*burst.bbox.l;
    Draw3D::drawPointCross( burst.bbox.p, 3.0 );
    Draw3D::drawPointCross( b           , 3.0 );
    Draw3D::drawCylinderStrip_wire( nseg, burst.bbox.r, burst.bbox.r, (Vec3f)burst.bbox.p, (Vec3f)b );
}

void fireBurst( Burst3d* burst, int n, double dt, const Vec3d& pos0, const Vec3d& vel0, double vrnd, double prnd ){
    for(int i=0; i<n; i++){
        Vec3d p = pos0 - vel0*dt*i + (Vec3d){randf(-vrnd,vrnd),randf(-vrnd,vrnd),randf(-vrnd,vrnd)};
        Vec3d v = vel0 + (Vec3d){randf(-prnd,prnd),randf(-prnd,prnd),randf(-prnd,prnd)};
        burst->addShot( p, v );
    }
}

class TestAppShotHit : public AppSDL2OGL_3D { public:
    double dvel = 10.0;
    Shooter world;

    int glo_burst = 0;
    ProjectileType projType1;
    ObjectType     objType1;

    int defaultObjectShape, defaultObjectHitShape;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppShotHit( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppShotHit::TestAppShotHit( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    defaultObjectShape = glGenLists(1);
    glNewList( defaultObjectShape , GL_COMPILE );
        glEnable( GL_LIGHTING );
        glColor3f( 0.8f, 0.8f, 0.8f ); Draw3D::drawPolygons( Solids::Icosahedron_nfaces, Solids::Icosahedron_ngons, Solids::Icosahedron_faces,  Solids::Icosahedron_verts );
    glEndList();

    defaultObjectHitShape = glGenLists(1);
    glNewList( defaultObjectHitShape , GL_COMPILE );
        glDisable ( GL_LIGHTING );
        //Draw3D::drawAxis ( 3.0f );
        //glColor3f( 0.8f, 0.0f, 0.8f );
        //Draw3D::drawSphereOctLines( 16, 2.0, (Vec3f){0.0,0.0,0.0} );
        Draw3D::drawSphere_oct( 3, 1.0, (Vec3f){0.0,0.0,0.0}, true );
    glEndList();

    // https://en.wikipedia.org/wiki/5.56%C3%9745mm_NATO

    projType1.mass    = 3.56e-3;
    projType1.caliber = 5.56e-3;
    //projType1.explosiveMass = ;
    projType1.updateAux( 0.3 );

    objType1.ogl = defaultObjectHitShape;

    for(int i=0; i<500; i++){
        //Vec3d pos; pos.fromRandomCube(50.0);
        //Vec3d vel; vel.fromRandomSphereSample(); vel.mul(800.0);
        Vec3d pos; pos.fromRandomCube(150.0);   pos.z = -500.0;
        //Vec3d vel; double phi = randf(0.0,2*M_PI); vel.set(cos(phi),sin(phi),0.0); vel.mul(800.0);
        Vec3d vel = Vec3dZ*800.0;
        Burst3d* burst = new Burst3d(&projType1, i);
        fireBurst( burst, 10, 0.02, pos, vel, 1.0, 1.0 );
        world.bursts.push_back( burst );
    }

    double objR = 5.0;
    for(int i=0; i<500; i++){
        Vec3d p; p.fromRandomCube(200.0);  p.z = 100.0;
        //Object3d* o = new Object3d( objR, p, &objType1, i );
        Object3d* o = new Target( objR, p, &objType1, i );
        //o->pos(,objType1);
        world.objects.push_back(o);
    }

    cameraMoveSpeed = 1.0;
    zoom = 160.0;
    //world.perFrame = 1;
    //world.dt = 0.000001;
    world.update_world();

    world.debugFlags |= Shooter::DEBUG_FLAG_SHOTS;

    //cam.zoom   = 160;
    //cam.aspect = ASPECT_RATIO;

    printf( "SETUP DONE \n"  );
}

void TestAppShotHit::draw(){
   // printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	//if(frameCount%10==0) burst.move( 0.1, {0.0,-9.81,0.0}, airDensity, tmpPos );
	//if(frameCount%10==0)
	if(!STOP){
        world.update_world();


        double r_av=0,l_av=0,z_av,vz_av;
        int nshots=0;
        for( int i=0; i<world.bursts.size(); i++ ){
            //Burst3d& b = *world.bursts[i];
            Capsula3D& bb = world.bursts[i]->bbox;
            r_av+=bb.r;  l_av+=bb.l; z_av+=bb.p.z;
            for( Particle3d& p: world.bursts[i]->shots){ vz_av+=p.vel.z; nshots++; }
            if( world.opCount_ShotObject>5000 ){
                printf( " %i   %g %g (%g,%g,%g) (%g,%g,%g) \n", i, bb.r, bb.l, bb.p.x, bb.p.y, bb.p.z,  bb.hdir.x, bb.hdir.y, bb.hdir.z );
            }
        }
        r_av/=world.bursts.size();  l_av/=world.bursts.size(); z_av/=world.bursts.size(); vz_av/=nshots;
        //printf( " average capsual size : %g %g %g \n", r_av, l_av, z_av );
        printf( " bbox_av : %g %g %g %g ", r_av, l_av, z_av, vz_av );
        if( r_av > 100.0 ){
            STOP = true;
            cam.pos = (Vec3f)world.bursts[0]->bbox.p;
            auto& bb = world.bursts[0]->bbox;
            printf(" bbox  %g %g (%g,%g,%g) (%g,%g,%g) \n", bb.r, bb.l, bb.p.x, bb.p.y, bb.p.z,  bb.hdir.x, bb.hdir.y, bb.hdir.z);
            auto& shots = world.bursts[0]->shots;
            for(int i=0; i<shots.size(); i++){
                Particle3d p = shots[i];
                printf( " %i (%g,%g,%g) (%g,%g,%g) \n", i,  p.pos.x, p.pos.y, p.pos.z,   p.vel.x, p.vel.y, p.vel.z );
            }
        }


        if( world.opCount_ShotObject>5000 ){ exit(0); }
    }

	for( Burst3d* burst : world.bursts ){
        //burst->move( 0.1, {0.0,-9.81,0.0}, airDensity );
        drawBurst( *burst, 8, world.dt );
	}

	//glCallList( glo_burst );

    for( Object3d* o : world.objects ){
        if(o->type == &objType1 ) ((Target*)o)->draw();
    }

};

void TestAppShotHit::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
                case SDLK_c:
                    Target* o = ((Target*)world.objects[0]);
                    setToRay( (Vec3d)cam.pos, (Vec3d)cam.rot.c, o->pos);
                    o->nhits = 0;
                    o->nboxhits = 0;
                    world.update_world();
                    break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}


void TestAppShotHit::drawHUD(){
    glDisable ( GL_LIGHTING );
    glColor3f( 0.0f, 1.0f, 0.0f );
    glBegin( GL_LINES );
    float whalf = WIDTH *0.5;
    float hhalf = HEIGHT*0.5;
    glVertex3f( whalf-10,hhalf, 0 ); glVertex3f( whalf+10,hhalf, 0 );
    glVertex3f( whalf,hhalf-10, 0 ); glVertex3f( whalf,hhalf+10, 0 );
    glEnd();
}

// ===================== MAIN

TestAppShotHit * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppShotHit( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















