
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

const int nshots = 20;
Vec3d tmpPos[nshots];


void checkHit_N2( std::vector<Projectile3D*>& projectiles, std::vector<Vec3d>& objects ){
    auto it_proj = projectiles.begin();
    while( it_proj != projectiles.end() ) {
        Projectile3D * proj = *it_proj;
        bool hitted = false;
        for( Vec3d& o : objects ){

        }
        /*
        Vec3d hRay,normal;
        hRay.set_sub( proj->pos, proj->old_pos );
        double tmax = hRay.normalize();
        for( o : objects ){
            o->ray( proj->old_pos, hRay, &normal );
            if (hitted) break;
        }
        */
        if( hitted || (proj->time > projLifetime) ){
            it_proj = projectiles.erase( it_proj );
            delete proj;
        }else{
            ++it_proj;
        }
    }
}


void drawBurst( Burst3d& burst, Vec3d* tmpPos = 0 ){
    glColor3f(1.0,0.0,0.0);
    int n = burst.shots.size();
    for(int i=0; i<n; i++){
        Draw3D::drawPointCross( burst.shots[i].pos, 0.1 );
        if( tmpPos ) Draw3D::drawLine( tmpPos[i], burst.shots[i].pos );
    }
    if( tmpPos ){ glColor3f(1.0,1.0,0.0); Draw3D::drawPointCross( tmpPos[n-1], 0.5 ); }
    glColor3f(0.0,1.0,1.0); Draw3D::drawPointCross( burst.shots[0].pos, 0.5 );
    glColor3f(0.0,0.5,0.0);
    Vec3d b = burst.bbox.p+burst.bbox.hdir*burst.bbox.l;
    Draw3D::drawPointCross( burst.bbox.p, 3.0 );
    Draw3D::drawPointCross( b           , 3.0 );
    Draw3D::drawCylinderStrip_wire( 16, burst.bbox.r, burst.bbox.r, (Vec3f)burst.bbox.p, (Vec3f)b );
}

class TestAppShotHit : public AppSDL2OGL_3D {
	public:
    //MultiFight3DWorld world;
    double dvel = 10.0;
    Shooter world;

    int glo_burst = 0;
    ProjectileType projType1;
    Burst3d burst;

    //std::vector<KinematicBody*> objects;
    int nobject = 100;
    Vec3d* objects;

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
        Draw3D::drawAxis ( 3.0f );
        glColor3f( 0.8f, 0.0f, 0.8f ); Draw3D::drawSphereOctLines( 16, 2.0, (Vec3f){0.0,0.0,0.0} );
    glEndList();

    objects = new Vec3d[ nobject ];
    float Lspan = 50.0;
    for( int i=0; i<nobject; i++ ){
        objects[i].set( randf(-Lspan,Lspan), randf(-Lspan,Lspan), randf(-Lspan,Lspan) );
    }

    // https://en.wikipedia.org/wiki/5.56%C3%9745mm_NATO

    projType1.mass    = 3.56e-3;
    projType1.caliber = 5.56e-3;
    //projType1.explosiveMass = ;
    projType1.updateAux( 0.3 );


    burst.type = &projType1;
    Vec3d vel0 = {800.0,0.0,0.0};
    Vec3d pos0 = {0.0,0.0,0.0};
    double dt = 0.02;
    double vrnd = 1;
    double prnd = 1;
    for(int i=0; i<nshots; i++){
        Vec3d p = pos0 - vel0*dt*i + (Vec3d){randf(-vrnd,vrnd),randf(-vrnd,vrnd),randf(-vrnd,vrnd)};
        Vec3d v = vel0 + (Vec3d){randf(-prnd,prnd),randf(-prnd,prnd),randf(-prnd,prnd)};
        burst.addShot( p, v );
    }

    burst.move( 0.1, {0.0,-9.81,0.0}, airDensity, tmpPos );

    glo_burst = glGenLists(1);
    glNewList( glo_burst , GL_COMPILE);
    drawBurst( burst, tmpPos );
    glEndList();

    zoom = 160.0;

}

void TestAppShotHit::draw(){
   // printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	if(frameCount%10==0) burst.move( 0.1, {0.0,-9.81,0.0}, airDensity, tmpPos );

	drawBurst( burst, tmpPos );

	glCallList( glo_burst );


	/*
    //for( auto o : world.objects ) {
    for( int i=0; i<nobject; i++ ) {
        //float glMat[16];
        glPushMatrix();
        //Draw3D::toGLMat( o->lpos, o->lrot, glMat );
        //glMultMatrixf( glMat );
        Vec3d * o = objects + i;
        glTranslatef( o->x, o->y, o->z );

        //double t = raySphere( {0.0d,0.0d,0.0d}, camMat.c, 2.0, o->lpos );
        double t = raySphere( (Vec3d)cam.pos, (Vec3d)cam.rot.c, 2.0, *o );
        if( ( t>0 ) && (t < 1000.0 ) ){
            //printf( " t %f  pos (%3.3f,%3.3f,%3.3f) \n", t, o->lpos.x, o->lpos.y, o->lpos.z );
            glCallList( defaultObjectHitShape );
        }

        glCallList( defaultObjectShape );
        glPopMatrix();
    }
    */

};

void TestAppShotHit::eventHandling ( const SDL_Event& event  ){
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

/*
void TestAppShotHit::mouseHandling( ){
    int mx,my;
    SDL_GetMouseState( &mouseX, &mouseY );   mouseY=HEIGHT-mouseY;
    mouse_begin_x = (2*mouseX-WIDTH )*zoom/HEIGHT;
    mouse_begin_y = (2*mouseY-HEIGHT)*zoom/HEIGHT;
    Uint32 buttons = SDL_GetRelativeMouseState( &mx, &my);
    //if ( buttons & SDL_BUTTON(SDL_BUTTON_RIGHT)) {
        Quat4f q; q.fromTrackball( 0, 0, -mx*mouseRotSpeed, my*mouseRotSpeed );
        qCamera.qmul_T( q );
    //}
}
*/

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
















