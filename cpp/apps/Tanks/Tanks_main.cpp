
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

#include "SDL_utils.h"
#include "Draw.h"
#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"
#include "testUtils.h"


#include "balistics.h"
#include "Mesh.h"
#include "Object3D.h"
#include "Warrior3D.h"
#include "Tank.h"
#include "Projectile3D.h"
#include "Shooter.h"

//int fontTex;
#include "TankHelpers.h"

class Tanks_single : public AppSDL2OGL_3D {
	public:
    Shooter world;
    double dvel = 10.0;

    //Warrior3D *warrior1;
    Tank *warrior1;

    std::vector<Object3D*> objects;

    int hitShape, warriorShape, objectShape;

    double camPhi=0,camTheta=0;

	virtual void draw   () override;
	virtual void drawHUD() override;
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  ) override;
	virtual void keyStateHandling( const Uint8 *keys )       override;
    virtual void mouseHandling   ( )                         override;
    virtual void camera          ( )                         override;

	Tanks_single( int& id, int WIDTH_, int HEIGHT_ );

};

Tanks_single::Tanks_single( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );

    world.init_world();
    world.perFrame = 3;
    printf( "DEBUG_SHIT : %i \n", world.debug_shit );

    // ---- terrain
    world.terrain = prepareTerrain();

    // ---- Objects
    int sphereShape = glGenLists(1);
    glNewList( sphereShape , GL_COMPILE );
        glEnable( GL_LIGHTING ); glColor3f( 0.8f, 0.8f, 0.8f );   Draw3D::drawSphere_oct( 6, 1.0, {0.0,0.0,0.0} );
        //glPopMatrix();
    glEndList();

    int nobjects=10;
    double xrange = 10.0;
    for( int i=0; i<nobjects; i++){
        Object3D * o = new Object3D();
        o->lrot.fromRand( {randf(0,1),randf(0,1),randf(0,1)} );
        o->grot = o->lrot;
        o->span.set( randf(0.2,2.0), randf(0.2,2.0), randf(0.2,2.0) );
        Vec2d p,dv; p.set( randf(-xrange,xrange),randf(-xrange,xrange) );
        double v = world.terrain->eval( p, dv );
        o->lpos.set( p.x, v, p.y );
        o->gpos = o->lpos;
        o->shape= sphereShape;
        o->id = i;
        objects.push_back(o);
    }

    //Tank * tank1 = new Tank();
    //warrior1 = new Warrior3D();
    warrior1 = new Tank();

    warrior1->fromFile( "data/tank1.txt" );
    warrior1->kind = 0;
    warrior1->hground = 2.5;
    warrior1->setPose( {0.0d,2.0d,0.0d}, {0.0d,0.0d,1.0d}, {0.0d,1.0d,0.0d} );
    world.registrWarrior( warrior1 );
    warrior1->makeWheels( 3, -3.0, 3.0, 3.0, -1.0, 20.0, 1.0, 0.8  );

    printf( "hull   mass : %g [kg]\n", warrior1->hull  .getArmorMass( 7890.0 ) );
    printf( "turret mass : %g [kg]\n", warrior1->turret.getArmorMass( 7890.0 ) );
    //exit(0);

    warrior1->hull  .polygonsToTriangles( true );
    warrior1->turret.polygonsToTriangles( true );

    warrior1->setInertia_box( 5.0, {3.0,3.0,6.0} );

    double maxThick = fmax( warrior1->hull.getMaxArmor(), warrior1->turret.getMaxArmor() );

    warrior1->hull.glo_armor = glGenLists(1);
    glNewList(warrior1->hull.glo_armor, GL_COMPILE);
        glDisable    ( GL_LIGHTING   );
        glShadeModel ( GL_FLAT       );
        renderArmor( warrior1->hull,   maxThick );
    glEndList();

    warrior1->turret.glo_armor = glGenLists(1);
    glNewList(warrior1->turret.glo_armor, GL_COMPILE);
        glDisable    ( GL_LIGHTING   );
        glShadeModel ( GL_FLAT       );
        renderArmor( warrior1->turret,   maxThick );
    glEndList();

    warrior1->hull.glo_captions = glGenLists(1);
    glNewList( warrior1->hull.glo_captions, GL_COMPILE);
        renderArmorCaptions( warrior1->hull,   0.15 );
    glEndList();

    warrior1->turret.glo_captions = glGenLists(1);
    glNewList(warrior1->turret.glo_captions, GL_COMPILE);
        renderArmorCaptions( warrior1->turret, 0.15 );
    glEndList();

    Tank* tank2 = new Tank();
    *tank2 = *warrior1;
    //tank2->pos.add(5.0,5.0,10.0);
    tank2->setPose( {5.0d,3.0d,10.0d}, {0.7d,0.0d,0.7d}, {0.0d,1.0d,0.0d} );

    world.registrWarrior( tank2 );
    tank2->rotateTurret( M_PI/3.0 );

    cam.zmin = 1.0;
    //cam.zoom = 100.0;
    zoom = 0.5;
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
	//glDisable ( GL_LIGHTING );
	//Draw3D::drawAxis ( 3.0f );

    //printf( "camMat.a (%3.3f,%3.3f,%3.3f) \n", camMat.a.x, camMat.a.y, camMat.a.z );

    //warrior1->gun_rot.set( camMat.c );

    glEnable( GL_LIGHTING   );
    glEnable( GL_DEPTH_TEST );
    //glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, 1 );
    if(world.terrain) glCallList(world.terrain->shape);

    //return;

    world.update_world( );

    //cam.pos = (Vec3f)warrior1->pos;
    cam.pos = (Vec3f)(warrior1->pos) +  (Vec3f){0.0,2.0,0.0} + cam.rot.c*(-4.0) ;

    Vec3d hRay,ray0,normal;
    hRay.set((Vec3d)cam.rot.c);
    ray0.set((Vec3d)cam.pos);
    for( auto o : objects ) {

        glColor3f(1.0f,0.0f,0.0f); Draw3D::drawVecInPos( o->grot.a*o->span.a*1.2, o->gpos );
        glColor3f(0.0f,1.0f,0.0f); Draw3D::drawVecInPos( o->grot.b*o->span.b*1.2, o->gpos );
        glColor3f(0.0f,0.0f,1.0f); Draw3D::drawVecInPos( o->grot.c*o->span.c*1.2, o->gpos );

        float glMat[16];
        glPushMatrix();
        Draw3D::toGLMat( o->gpos, o->grot, o->span, glMat );

        glMultMatrixf( glMat );
        glCallList( o->shape );
        glPopMatrix();

        double thit;
        if ( rayPointDistance2( ray0, hRay,  o->gpos, thit ) > 4.0 ) continue;
        thit = o->ray( ray0, hRay,  &normal );
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

    glCallList( warrior1->hull  .glo_armor    );
    glCallList( warrior1->turret.glo_armor    );

    warrior1->rotateTurretToward( (Vec3d)cam.rot.c );

    for( Warrior3D * w : world.warriors ){
        Tank * tank =  ((Tank*)w);
        Draw3D::drawShape( tank->hull.glo_armor, tank->pos, tank->rotMat );
        Mat3d grot;
        //tank->turret.globalRotT(tank->rotMat, grot);
        tank->turret.globalRot(tank->rotMat, grot);
        Draw3D::drawMatInPos( grot*60, tank->pos + (Vec3d){0.0,1.0,0.0} );
        glColor3f(1.0,1.0,1.0);
        Draw3D::drawVecInPos( tank->gun_rot, tank->pos + (Vec3d){0.0,1.0,0.0} );
        Draw3D::drawShape( tank->turret.glo_armor, tank->pos, grot  );
        drawTankWheels(tank);
    }

    //glCallList( warrior1->hull  .glo_captions ); // TODO : does not work because of Draw::billboardCamProj contains ::glGetFloatv (GL_MODELVIEW_MATRIX,  glModel);
    //glCallList( warrior1->turret.glo_captions );
    glColor3f(1.0f,0.5f,0.5f);  renderArmorCaptions( warrior1->hull,   0.15 );
    glColor3f(0.5f,0.5f,1.0f);  renderArmorCaptions( warrior1->turret, 0.15 );

    //Vec3d hRay = camMat.c;
    //Vec3d ray0 = camPos;

    hRay = (Vec3d)cam.rot.c; ray0 = (Vec3d)cam.pos;
    //hRay = warrior1->gun_rot; ray0 = warrior1->pos;

    // ray vs tank
    Tank * tank2 = (Tank*)world.warriors[1];
    int ipl; VehicleBlock* block; double effthick;
    //Vec3d normal;
    double t = tank2->ray( ray0, hRay, ipl, block, effthick, normal );
    if( ipl>=0 ){
        glColor3f(0.0f,1.0f,0.0f);
        //Draw3D::drawVecInPos( normal, ray0 + camMat.c*t );
        printf( "ray t %g \n", t );
        //Draw3D::drawPointCross( ray0 + hRay*t, 100.5 );

        Mat3d grot;
        //block->globalRotT( tank2->rotMat, grot );
        block->globalRot( tank2->rotMat, grot );

        //Draw3D::drawVecInPos( grot.dotT(block->armor[ipl].normal), ray0 + hRay*t );
        //Draw3D::drawVecInPos( grot.dotT(normal), ray0 + hRay*t );
        Draw3D::drawVecInPos( normal, ray0 + hRay*t );
        glPushMatrix();
            float glMat[16];
            Draw3D::toGLMat( tank2->pos, grot, glMat );
            glMultMatrixf( glMat );
            Draw3D::drawPolygonBorder( ipl, *block );
            Draw3D::drawPolygonNormal( ipl, *block );
        glPopMatrix();
        char str[64];

        sprintf(str,"%4.0fmm\0", effthick );
        Draw3D::drawText(str, ray0 + hRay*t, fontTex, 0.2, 0 );
        //printf( "itr %i ipl %i %g %g %g\n", itr, ipl, thick, cdot, effthick  );
    }

    glEnable     ( GL_LIGHTING   );
    glEnable     ( GL_DEPTH_TEST );
    glShadeModel ( GL_SMOOTH     );

    // Draw3D::drawAxis( 10 );

};


void Tanks_single::keyStateHandling( const Uint8 *keys ){

    if( warrior1 != NULL ){
        warrior1->power_gear[1]= 0.0; warrior1->power_gear[0]= 0.0;
        if( keys[ SDL_SCANCODE_W ] ){ warrior1->power_gear[1]= 1.0; warrior1->power_gear[0]= 1.0; }
        if( keys[ SDL_SCANCODE_S ] ){ warrior1->power_gear[1]=-1.0; warrior1->power_gear[0]=-1.0; }
        if( keys[ SDL_SCANCODE_A ] ){ warrior1->power_gear[1]= 1.0; warrior1->power_gear[0]=-1.0; }
        if( keys[ SDL_SCANCODE_D ] ){ warrior1->power_gear[1]=-1.0; warrior1->power_gear[0]= 1.0; }
        //camPos.set( warrior1->pos );
        //camPos.set_add( warrior1->pos, {0.0,2.0,0.0} );
    }
    if( keys[ SDL_SCANCODE_Q ] ){ qCamera.droll2(  +0.01 ); }
	if( keys[ SDL_SCANCODE_E ] ){ qCamera.droll2(  -0.01 ); }

};

void Tanks_single::mouseHandling( ){
    int mx,my;
    Uint32 buttons = SDL_GetRelativeMouseState( &mx, &my);

    camPhi   += mx*0.001;
    camTheta += my*0.001;
    printf( "cam.Phi,Theta : %g %g \n", camPhi, camTheta );

    cam.rot.a.set( sin(camPhi),0.0,cos(camPhi) );
    double ct=cos(camTheta),st=sin(camTheta);
    cam.rot.b.set(-cam.rot.a.z*st, ct, cam.rot.a.x*st);
    cam.rot.c.set(-cam.rot.a.z*ct,-st, cam.rot.a.x*ct); // up vector



    //printf("camMat:\n",);
    //printf("camMat: (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f)\n",camMat.ax,camMat.ay,camMat.az,  camMat.bx,camMat.by,camMat.bz, camMat.cx,camMat.cy,camMat.cz);
    //qCamera.fromAngleAxis( camTheta, {sin(camPhi),0, cos(camPhi)} );
    //qCamera.fromAngleAxis( camTheta, {0,sin(camPhi), cos(camPhi)} );

}

void Tanks_single::camera(){

    //((Quat4f)qCamera).toMatrix(cam.rot);
    cam.zoom   = zoom;
    cam.aspect = ASPECT_RATIO;
    //Cam::ortho( cam, true );
    //Cam::perspective( cam );
    if (perspective){ Cam::perspective( cam ); }
    else            { Cam::ortho( cam, true ); }

    /*
    float camMatrix[16];
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity();
    float fov = VIEW_ZOOM_DEFAULT/zoom;
    glFrustum( -ASPECT_RATIO, ASPECT_RATIO, -1, 1, 1*fov, VIEW_DEPTH*fov );
    Draw3D::toGLMatCam( {0.0d,0.0d,0.0d}, cam.rot, camMatrix );
    glMultMatrixf(  camMatrix );
    glTranslatef ( -cam.pos.x, -cam.pos.y, -cam.pos.z );
    //glTranslatef ( -10*camMat.c.x, -10*camMat.c.y, -10*camMat.c.z );
    //glTranslatef ( 0.0, -5.0, 0.0 );
    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity();
    */
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
    SDL_DisplayMode dm;
    SDL_GetDesktopDisplayMode(0, &dm);
	int junk;
	thisApp = new Tanks_single( junk , dm.w-150, dm.h-100 );
	SDL_SetWindowPosition(thisApp->window, 100, 0 );

	thisApp->loop( 1000000 );
	return 0;
}
















