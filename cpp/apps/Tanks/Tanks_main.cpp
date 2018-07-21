
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

int fontTex;

void armorColorScale( float f ){
    float mf = 1 - f;
    //float r = f*f; float g = 4*mf*f; float b = mf*mf;

    float r = 1-mf*mf; float g = f*f*f; float b = mf*mf*mf*mf;
    //float r = 1-mf*mf; float g = f*f*f; float b = f*f*f*f*f*f + 0.5*mf*mf*mf*mf;

    //float renorm = 3/(r+g+b); r*=renorm; g*=renorm; b*=renorm;
    glColor3f( r, g, b );
}

void renderArmor( const VehicleBlock& block, double maxThickness ){
    int i=0;
    for( Polygon* pl : block.polygons ){
        //printf( " pl %i npoints %i \n", i, pl->ipoints.size() );
        armorColorScale( block.armor[i].thickness/maxThickness );
        Draw3D::drawPlanarPolygon( pl->ipoints.size(), &pl->ipoints.front(), &block.points.front() );
        //Vec3d c = block.faceCog( i );
        //glColor3f(1.0f,1.0f,1.0f); Draw3D::drawVecInPos( block.armor[i].normal, c );
        i++;
    }

}

void renderArmorCaptions( const VehicleBlock& block, float sz ){
    char str[64];
    int i=0;
    for( Polygon* pl : block.polygons ){
        sprintf(str,"%i:%3.0fmm%2.1fton\0",i+1, block.armor[i].thickness, block.armor[i].mass*1e-3 );
        Vec3d c = block.faceCog( i );
        Draw3D::drawText(str, c, fontTex, sz, 0 );
        i++;
    }
}

void drawTankWheels(Tank * tank){
    Vec3d gpos;
    for(int i=0; i<tank->nwheel; i++){
        tank->getWheelPos( i, gpos );
        //printf( "%i (%g,%g,%g)\n", i, gpos.x, gpos.y, gpos.z );
        Draw3D::drawPointCross( gpos, 0.5 );
    }
}

Terrain25D *  prepareTerrain(){
//    Terrain25D * terrain = new Terrain25D();

    Terrain25D_bicubic * terrain = new Terrain25D_bicubic();
    terrain->ruler.setup( (Vec2d){10.0,10.0}*-16, (Vec2d){10.0d,10.0d} );
    terrain->allocate( {32,32} );
    terrain->makeRandom( -2.0, 2.0 );

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
            glNormal3f(-dv1.x,1.0,-dv1.y); glVertex3f( (float)p1.x,  v1, (float)p1.y );
            glNormal3f(-dv2.x,1.0,-dv2.y); glVertex3f( (float)p2.x,  v2, (float)p2.y );

            //glColor3f(v1,0.5,-v1); glVertex3f( (float)p1.x,  v1, (float)p1.y );
            //glColor3f(v2,0.5,-v2); glVertex3f( (float)p2.x,  v2, (float)p2.y );

            //printf( " %i (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f)\n", p1.x, p1.y, v1 ,  p2.x, p2.y, v2  );
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

    //Warrior3D *warrior1;
    Tank *warrior1;

    std::vector<Object3D*> objects;

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

    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );

    world.init_world();
    world.perFrame = 3;
    printf( "DEBUG_SHIT : %i \n", world.debug_shit );

    // ---- terrain

    //new Terrain25D();
    world.terrain = prepareTerrain();

    //exit(0);

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

        o->lrot.fromRand( {randf(0,1),randf(0,1),randf(0,1)} );
        o->grot = o->lrot;
        o->span.set( randf(0.2,2.0), randf(0.2,2.0), randf(0.2,2.0) );

        //o->bounds.span.set( 1, 2.0, 0.5 );
        //
        /*
        Mat3d m; m.set_mmul_NT( o->lrot, o->lrot );
        printf( " === %i \n", i );
        printf( " %f %f %f \n", m.ax, m.ay, m.az );
        printf( " %f %f %f \n", m.bx, m.by, m.bz );
        printf( " %f %f %f \n", m.cx, m.cy, m.cz );
        */
        //o->bounds.pos.set( randf(-xrange,xrange),randf(-xrange,xrange),randf(-xrange,xrange) );
        Vec2d p,dv; p.set( randf(-xrange,xrange),randf(-xrange,xrange) );
        double v = world.terrain->eval( p, dv );
        o->lpos.set( p.x, v, p.y );
        o->gpos = o->lpos;

        //o->bounds.pos.set( {0.0,0.0,0.0} );
        o->shape= sphereShape;
        o->id = i;
        //world.objects.push_back(o);
        objects.push_back(o);

    }

    //camPos.set();

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

    warrior1->initSpherical( 5.0, 10.0 );

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

    //warrior1->gun_rot.set( camMat.c );

    if(world.terrain) glCallList(world.terrain->shape);

    //return;

    world.update_world( );

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


    //Draw3D::drawShape( warrior1->pos, warrior1->qrot, warrior1->hull  .glo_armor  );
    //Draw3D::drawShape( warrior1->pos, warrior1->qrot, warrior1->turret.glo_armor  );
    for( Warrior3D * w : world.warriors ){
        Tank * tank =  ((Tank*)w);

        //Draw3D::drawShape( tank->pos, tank->qrot.get_inv(), tank->hull.glo_armor );
        //Draw3D::drawShape( tank->pos, tank->qrot, tank->hull.glo_armor );
        Draw3D::drawShape( tank->pos, tank->rotMat, tank->hull.glo_armor );

        Mat3d grot;
        //tank->turret.globalRotT(tank->rotMat, grot);
        tank->turret.globalRot(tank->rotMat, grot);

        Draw3D::drawMatInPos( grot*60, tank->pos + (Vec3d){0.0,1.0,0.0} );
        glColor3f(1.0,1.0,1.0);
        Draw3D::drawVecInPos( tank->gun_rot, tank->pos + (Vec3d){0.0,1.0,0.0} );


        Draw3D::drawShape( tank->pos, grot, tank->turret.glo_armor  );

        drawTankWheels(tank);

        /*
        Mat3d grot;
        //globalRot( const Mat3d& rot0, Mat3d& grot );
        tank->turret.globalRotT(tank->rotMat, grot);
        //Draw3D::drawShapeT( tank->pos, tank->qrot, tank->hull.glo_armor  );
        Draw3D::drawShape( tank->pos, tank->qrot.get_inv(), tank->hull.glo_armor );
        //Draw3D::drawShapeT( tank->pos, tank->qrot, tank->turret.glo_armor  );
        Draw3D::drawShape( tank->pos, grot, tank->turret.glo_armor  );
        drawTankWheels(tank);
        //Mat3d setT( const MAT& M );
        Mat3d rotMat; tank->qrot.toMatrix_T(rotMat);
        Draw3D::drawMatInPos( rotMat*60, tank->pos               );
        Draw3D::drawMatInPos( rotMat*60, tank->turret.pos );
        glColor3f(1.0f,1.0f,1.0f); Draw3D::drawVecInPos( tank->gun_rot*10.0, tank->pos );
        */

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

    warrior1->power_gear[1]= 0.0; warrior1->power_gear[0]= 0.0;

    if( warrior1 != NULL ){
        //if( keys[ SDL_SCANCODE_W ] ){ warrior1->pos.add_mul( camMat.c, +0.1 ); }
        //if( keys[ SDL_SCANCODE_S ] ){ warrior1->pos.add_mul( camMat.c, -0.1 ); }
        //if( keys[ SDL_SCANCODE_A ] ){ warrior1->pos.add_mul( camMat.a, -0.1 ); }
        //if( keys[ SDL_SCANCODE_D ] ){ warrior1->pos.add_mul( camMat.a, +0.1 ); }

        //if( keys[ SDL_SCANCODE_W ] ){ warrior1->vel.add_mul( camMat.c, +0.1 ); }
        //if( keys[ SDL_SCANCODE_S ] ){ warrior1->vel.add_mul( camMat.c, -0.1 ); }
        //if( keys[ SDL_SCANCODE_A ] ){ warrior1->vel.add_mul( camMat.a, -0.1 ); }
        //if( keys[ SDL_SCANCODE_D ] ){ warrior1->vel.add_mul( camMat.a, +0.1 ); }
        //if( keys[ SDL_SCANCODE_SPACE ] ){ warrior1->vel.mul( 0.9 ); }

        if( keys[ SDL_SCANCODE_W ] ){ warrior1->power_gear[1]= 1.0; warrior1->power_gear[0]= 1.0; }
        if( keys[ SDL_SCANCODE_S ] ){ warrior1->power_gear[1]=-1.0; warrior1->power_gear[0]=-1.0; }
        if( keys[ SDL_SCANCODE_A ] ){ warrior1->power_gear[1]= 1.0; warrior1->power_gear[0]=-1.0; }
        if( keys[ SDL_SCANCODE_D ] ){ warrior1->power_gear[1]=-1.0; warrior1->power_gear[0]= 1.0; }

        //camPos.set( warrior1->pos );
        //camPos.set_add( warrior1->pos, {0.0,2.0,0.0} );
        cam.pos = (Vec3f)(warrior1->pos) +  (Vec3f){0.0,2.0,0.0} + cam.rot.c*(-4.0) ;
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
    float camMatrix[16];
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity();
    float fov = VIEW_ZOOM_DEFAULT/zoom;
    glFrustum( -ASPECT_RATIO, ASPECT_RATIO, -1, 1, 1*fov, VIEW_DEPTH*fov );
    Draw3D::toGLMatCam( {0.0d,0.0d,0.0d}, cam.rot, camMatrix );
    glMultMatrixf( camMatrix );
    glTranslatef ( -cam.pos.x, -cam.pos.y, -cam.pos.z );
    //glTranslatef ( -10*camMat.c.x, -10*camMat.c.y, -10*camMat.c.z );
    //glTranslatef ( 0.0, -5.0, 0.0 );
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
















