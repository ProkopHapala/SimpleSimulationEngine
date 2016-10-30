
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
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "geom3D.h"
#include "raytrace.h"

#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

#include "Object3D.h"
#include "Warrior3D.h"
#include "Projectile3D.h"
#include "Shooter.h"


Terrain25D *  prepareTerrain()      {
    Terrain25D * terrain = new Terrain25D();
    terrain->shape = glGenLists(1);
    glNewList( terrain->shape , GL_COMPILE );
    int na=100,nb=100;
    float da=1.0,db=1.0;
    float x0=-0.5*da*na,y0=-0.5*db*nb;
    glEnable(GL_LIGHTING);
    glColor3f(0.5f,0.5f,0.5f);
    glNormal3f(0.0f,1.0f,0.0f);
    /*
    glBegin(GL_QUADS);
        glVertex3f( na*da,     0, 0 );
        glVertex3f( 0,         0, 0 );
        glVertex3f( 0,     nb*db, 0 );
        glVertex3f( na*da, nb*db, 0 );
    glEnd();
    */
    float * oldvals = new float[na*3];
    for(int ia=0; ia<na; ia++){
        glBegin(GL_TRIANGLE_STRIP);
        for(int ib=0; ib<nb; ib++){
            int i3 = 3*ib;
            Vec2d dv1,dv2;
            Vec2d p1; p1.set( (ia  )*da+x0, ib*db+y0 );
            Vec2d p2; p2.set( (ia+1)*da+x0, ib*db+y0 );
            float v1,v2;
            if( ia == 0 ){
                v1 = (float)terrain->eval( p1, dv1 );
            }else{
                v1 = oldvals[i3]; dv1.x=oldvals[i3+1]; dv1.y=oldvals[i3+2];
            }
            v2 = (float)terrain->eval( p2, dv2 );
            oldvals[i3] = v2; oldvals[i3+1] = dv2.x; oldvals[i3+2] = dv2.y;
            glNormal3f(-dv1.x,1.0,-dv1.y); glVertex3f( (float)p1.x,  v1, (float)p1.y  );
            glNormal3f(-dv2.x,1.0,-dv2.y); glVertex3f( (float)p2.x,  v2, (float)p2.y );
            printf( " %i (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f)\n", p1.x, p1.y, v1 ,  p2.x, p2.y, v2  );
        }
        glEnd();
    }

    glBegin(GL_LINES);
    for(int ia=0; ia<na; ia++){
        for(int ib=0; ib<nb; ib++){
            int i3 = 3*ib;
            Vec2d p,dv; p.set( ia*da+x0, ib*db+y0 );
            double v = (float)terrain->eval( p, dv );
            glVertex3f( (float)p.x,         v, (float)p.y );
            glVertex3f( (float)(p.x-dv.x),  v+1.0, (float)(p.y-dv.y) );
        }

    }
    glEnd();
    glEndList();
    return terrain;
}

class Tanks_single : public AppSDL2OGL_3D {
	public:
    Shooter world;
    double dvel = 10.0;

    Warrior3D *warrior1;

    int hitShape, warriorShape, objectShape;

    double camPhi,camTheta;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	virtual void keyStateHandling( const Uint8 *keys );
    virtual void mouseHandling( );
    void         camera();

	Tanks_single( int& id, int WIDTH_, int HEIGHT_ );

};

Tanks_single::Tanks_single( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    world.init_world();
    printf( "DEBUG_SHIT : %i \n", world.debug_shit );

    // ---- terrain

    //new Terrain25D();
    world.terrain = prepareTerrain();

    // ---- Objects

    int sphereShape = glGenLists(1);
    glNewList( sphereShape , GL_COMPILE );
        //glPushMatrix();

        //glDisable ( GL_LIGHTING );
        //Draw3D::drawAxis ( 3.0f );
        //glColor3f( 1.0f, 0.0f, 1.0f ); Draw3D::drawLines   ( Solids::Icosahedron_nedges, Solids::Icosahedron_edges, Solids::Icosahedron_verts                             );

        //glEnable( GL_LIGHTING );
        //glColor3f( 0.8f, 0.8f, 0.8f ); Draw3D::drawPolygons( Solids::Icosahedron_nfaces, Solids::Icosahedron_ngons, Solids::Icosahedron_faces,  Solids::Icosahedron_verts );

        glEnable( GL_LIGHTING ); glColor3f( 0.8f, 0.8f, 0.8f );   Draw3D::drawSphere_oct( 6, 1.0, {0.0,0.0,0.0} );
        //glPopMatrix();
    glEndList();

    int nobjects=10;
    double xrange = 10.0;
    for( int i=0; i<nobjects; i++){
        Object3D * o = new Object3D();
        //o->bounds.orientation.set({});
        //o->bounds.span.set(o->bounds.orientation.a.normalize(),o->bounds.orientation.a.normalize(),o->bounds.orientation.a.normalize());
        o->bounds.orientation.fromRand( {randf(0,1),randf(0,1),randf(0,1)} );
        o->bounds.span.set( randf(0.2,2.0), randf(0.2,2.0), randf(0.2,2.0) );
        //o->bounds.span.set( 1, 2.0, 0.5 );
        Mat3d m; m.set_mmul_NT( o->bounds.orientation, o->bounds.orientation );
        /*
        printf( " === %i \n", i );
        printf( " %f %f %f \n", m.ax, m.ay, m.az );
        printf( " %f %f %f \n", m.bx, m.by, m.bz );
        printf( " %f %f %f \n", m.cx, m.cy, m.cz );
        */
        //o->bounds.pos.set( randf(-xrange,xrange),randf(-xrange,xrange),randf(-xrange,xrange) );
        Vec2d p,dv; p.set( randf(-xrange,xrange),randf(-xrange,xrange) );
        double v = world.terrain->eval( p, dv );
        o->bounds.pos.set( p.x, v, p.y );

        //o->bounds.pos.set( {0.0,0.0,0.0} );
        o->shape= sphereShape;
        o->id = i;
        world.objects.push_back(o);

    }

    //camPos.set();


    warrior1 = world.makeWarrior( {0.0d,0.0d,0.0d}, {0.0d,0.0d,1.0d}, {0.0d,1.0d,0.0d}, 0 );

    warrior1->pos.set( 0.0, 0.0, -15.0 );

    zoom = 5.0;
    first_person = true;
    perspective  = true;

}

void Tanks_single::draw(){
    //delay = 200;
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    //glEnable( GL_LIGHTING );
    //glCallList( world.defaultObjectShape );
	glDisable ( GL_LIGHTING );
	Draw3D::drawAxis ( 3.0f );

    //printf( "camMat.a (%3.3f,%3.3f,%3.3f) \n", camMat.a.x, camMat.a.y, camMat.a.z );

    warrior1->gun_rot.set( camMat.c );

    world.update_world( );

    if(world.terrain) glCallList(world.terrain->shape);

    Vec3d hRay,ray0,normal;
    hRay.set(camMat.c);
    ray0.set(camPos);
    for( auto o : world.objects ) {

        glColor3f(1.0f,0.0f,0.0f); Draw3D::drawVecInPos( o->bounds.orientation.a*o->bounds.span.a*1.2, o->bounds.pos );
        glColor3f(0.0f,1.0f,0.0f); Draw3D::drawVecInPos( o->bounds.orientation.b*o->bounds.span.b*1.2, o->bounds.pos );
        glColor3f(0.0f,0.0f,1.0f); Draw3D::drawVecInPos( o->bounds.orientation.c*o->bounds.span.c*1.2, o->bounds.pos );

        float glMat[16];
        glPushMatrix();
        Draw3D::toGLMat( o->bounds.pos, o->bounds.orientation, o->bounds.span, glMat );

        glMultMatrixf( glMat );
        glCallList( o->shape );
        glPopMatrix();

        double thit;
        if ( rayPointDistance2( ray0, hRay,  o->bounds.pos, thit ) > 4.0 ) continue;
        thit = o->bounds.ray( ray0, hRay,  &normal );
        if( (thit>0)&&(thit < 0.999e+300) ){
            Vec3d hit_point; hit_point.set_add_mul( ray0, hRay, thit );
            //printf( " %i %g (%3.3f,%3.3f,%3.3f)  \n", o->id,  thit, hit_point.x, hit_point.y, hit_point.z );
            Draw3D::drawVecInPos  ( normal, hit_point );
            Draw3D::drawPointCross( hit_point, 0.1 );
        }
    }

    for( auto p : world.projectiles ) {
        Draw3D::drawLine(p->pos, p->old_pos);
        Draw3D::drawPointCross(p->pos,0.1);
    }

};





void Tanks_single::keyStateHandling( const Uint8 *keys ){

    if( warrior1 != NULL ){
        //if( keys[ SDL_SCANCODE_W ] ){ warrior1->pos.add_mul( camMat.c, +0.1 ); }
        //if( keys[ SDL_SCANCODE_S ] ){ warrior1->pos.add_mul( camMat.c, -0.1 ); }
        //if( keys[ SDL_SCANCODE_A ] ){ warrior1->pos.add_mul( camMat.a, -0.1 ); }
        //if( keys[ SDL_SCANCODE_D ] ){ warrior1->pos.add_mul( camMat.a, +0.1 ); }

        if( keys[ SDL_SCANCODE_W ] ){ warrior1->vel.add_mul( camMat.c, +0.1 ); }
        if( keys[ SDL_SCANCODE_S ] ){ warrior1->vel.add_mul( camMat.c, -0.1 ); }
        if( keys[ SDL_SCANCODE_A ] ){ warrior1->vel.add_mul( camMat.a, -0.1 ); }
        if( keys[ SDL_SCANCODE_D ] ){ warrior1->vel.add_mul( camMat.a, +0.1 ); }
        if( keys[ SDL_SCANCODE_SPACE ] ){ warrior1->vel.mul( 0.9 ); }

        //camPos.set( warrior1->pos );
        camPos.set_add( warrior1->pos, {0.0,1.0,0.0} );
    }

    //if( keys[ SDL_SCANCODE_W ] ){ camPos.add_mul( camMat.c, +0.1 ); }
	//if( keys[ SDL_SCANCODE_S ] ){ camPos.add_mul( camMat.c, -0.1 ); }
	//if( keys[ SDL_SCANCODE_A ] ){ camPos.add_mul( camMat.a, -0.1 ); }
	//if( keys[ SDL_SCANCODE_D ] ){ camPos.add_mul( camMat.a, +0.1 ); }
	//if( keys[ SDL_SCANCODE_W ] ){ camPos.z += +0.1; }
	//if( keys[ SDL_SCANCODE_S ] ){ camPos.z += -0.1; }
	//if( keys[ SDL_SCANCODE_A ] ){ camPos.x += +0.1; }
	//if( keys[ SDL_SCANCODE_D ] ){ camPos.x += -0.1; }
    if( keys[ SDL_SCANCODE_Q ] ){ qCamera.droll2(  +0.01 ); }
	if( keys[ SDL_SCANCODE_E ] ){ qCamera.droll2(  -0.01 ); }

};

void Tanks_single::mouseHandling( ){
    int mx,my;
    Uint32 buttons = SDL_GetRelativeMouseState( &mx, &my);

    camPhi   += mx*0.001;
    camTheta += my*0.001;
    camMat.a.set( sin(camPhi),0.0,cos(camPhi) );
    double ct=cos(camTheta),st=sin(camTheta);
    camMat.b.set(-camMat.a.z*st, ct, camMat.a.x*st);
    camMat.c.set(-camMat.a.z*ct,-st, camMat.a.x*ct); // up vector
    //printf("camMat:\n",);
    printf("camMat: (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f)\n",camMat.ax,camMat.ay,camMat.az,  camMat.bx,camMat.by,camMat.bz, camMat.cx,camMat.cy,camMat.cz);
    //qCamera.fromAngleAxis( camTheta, {sin(camPhi),0, cos(camPhi)} );
    //qCamera.fromAngleAxis( camTheta, {0,sin(camPhi), cos(camPhi)} );

}

void Tanks_single::camera(){
    float camMatrix[16];
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity();
    float fov = VIEW_ZOOM_DEFAULT/zoom;
    glFrustum( -ASPECT_RATIO, ASPECT_RATIO, -1, 1, 1*fov, VIEW_DEPTH*fov );
    Draw3D::toGLMatCam( {0.0d,0.0d,0.0d}, camMat, camMatrix );
    glMultMatrixf( camMatrix );
    glTranslatef ( -camPos.x, -camPos.y, -camPos.z );
    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity();
}



void Tanks_single::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_f:  warrior1->tryJump(); break;
                //case SDLK_h:  warrior1->tryJump(); break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
                case SDLK_ESCAPE:   quit(); break;
                //case SDLK_SPACE:    STOP = !STOP; printf( STOP ? " STOPED\n" : " UNSTOPED\n"); break;
                case SDLK_KP_MINUS: zoom*=VIEW_ZOOM_STEP; break;
                case SDLK_KP_PLUS:  zoom/=VIEW_ZOOM_STEP; break;
            }
            break;

        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    warrior1->trigger = true;
                    break;
            }
            break;

        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    warrior1->trigger = false;
                    break;
            }
            break;
    };
    //AppSDL2OGL::eventHandling( event );
}




void Tanks_single::drawHUD(){
    glDisable ( GL_LIGHTING );
    glColor3f( 0.0f, 1.0f, 0.0f );
    drawCrosshair( 10 );
}

// ===================== MAIN

Tanks_single * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new Tanks_single( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















