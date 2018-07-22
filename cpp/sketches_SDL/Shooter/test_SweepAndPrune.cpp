
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

const int ndirs = 1;
Vec3d dirs[ndirs] = { {1.0,0.0,0.0} };

//const int nshots = 20;
//Vec3d tmpPos[nshots];

typedef std::pair<Object3d*,Object3d*> ObjectCollision;

std::vector<ObjectCollision> collisions;

class SAPTestObj : public Object3d{ public:
    int ncol=0;
    //virtual bool getShot( const Vec3d& p0, const Vec3d& p1, const ProjectileType& prjType, double dt ){
    virtual bool collideNarrow( Object3d* obj ) override {
        ncol++;
        if(obj->type==type) ((SAPTestObj*)obj)->ncol++;
        //collisions.push_back( {this,obj} );
        return true;
    };
    void draw(){
        glPushMatrix();
        glTranslatef( pos.x, pos.y, pos.z );
        glScalef(R,R,R);
        //glCallList( defaultObjectHitShape );
        //if(ncol){ glColor3f(1.0,0.0,0.0); }else{ glColor3f(0.0,0.0,1.0); }
        glCallList( type->ogl );
        glPopMatrix();
    }
    SAPTestObj( double R_, Vec3d pos_, ObjectType* type_, int id_ ){R=R_;pos=pos_;type=type_;id=id_;  };
};

class SAPTestVeh : public Vehicle3d{ public:
    int ncol=0;
    virtual bool collideNarrow( Object3d* obj ) override {
        ncol++;
        if(obj->type==type) ((SAPTestVeh*)obj)->ncol++;
        glColor3f(0.0,1.0,0.0); Draw3D::drawLine( pos, obj->pos );
        return true;
    };
    void draw(){
        glPushMatrix();
        glTranslatef( pos.x, pos.y, pos.z );
        glScalef(R,R,R);
        if(ncol){ glColor3f(1.0,0.0,0.0); }else{ glColor3f(0.0,0.0,1.0); }
        glCallList( type->ogl );
        glPopMatrix();
    }
    SAPTestVeh( double R_, Vec3d pos_, ObjectType* type_, int id_ ){R=R_;pos=pos_;type=type_;id=id_;  };
};

void drawSweep( const std::vector<SAPitem>& sweep, float y0, float dy, float WIDTH ){
    float xmin  = sweep[0]    .span.x;
    float xmax  = sweep.back().span.y;
    float sc    = WIDTH/(xmax-xmin);
    glBegin(GL_LINES);
    for( int i=0; i<sweep.size(); i++ ){
        float x0 = (sweep[i].span.x-xmin)*sc;
        float x1 = (sweep[i].span.y-xmin)*sc;
        glVertex3f(x0,y0+i*dy,0);
        glVertex3f(x1,y0+i*dy,0);
    }
    glEnd();
}

class TestAppSweepAndPrune : public AppSDL2OGL_3D { public:
    Shooter world;

    int ipicked;
    int glo_burst = 0;
    ObjectType     objType1;
    VehicleType    vehType1;

    SAPbuff sapObj,sapVeh;

    int oglSphereHitShape=0;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppSweepAndPrune( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppSweepAndPrune::TestAppSweepAndPrune( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    oglSphereHitShape = glGenLists(1);
    glNewList( oglSphereHitShape , GL_COMPILE );
        glDisable ( GL_LIGHTING );
        Draw3D::drawSphere_oct( 3, 1.0, (Vec3f){0.0,0.0,0.0}, true );
    glEndList();

    // https://en.wikipedia.org/wiki/5.56%C3%9745mm_NATO

    objType1.ogl = oglSphereHitShape;
    vehType1.ogl = oglSphereHitShape;

    double objR = 25.0;
    for(int i=0; i<25; i++){
        Vec3d p; p.fromRandomCube(200.0);  p.z = 100.0;
        Object3d* o = new SAPTestObj( objR, p, &objType1, i );
        world.objects.push_back(o);
    }

    objR = 5.0;
    for(int i=0; i<25; i++){
        Vec3d p; p.fromRandomCube(200.0);  p.z = 100.0;
        Vehicle3d* o = new SAPTestVeh( objR, p, &vehType1, i );
        world.vehicles.push_back(o);
    }

    ipicked = rand()%world.vehicles.size();

    sapObj.init( ndirs, dirs, 20 );
    sapObj.addObjects( world.objects.size(), world.objects.data() );
    sapObj.sort();

    sapVeh.init( ndirs, dirs, 20 );
    //sapVeh.addObjects( world.vehicles.size(), world.vehicles.data() );
    for( Vehicle3d* veh : world.vehicles ){ sapVeh.addObject( veh ); }
    sapVeh.sort();

    /*
    sap.sort(); printf( "DEBUG 2.3 \n");
    sap.collideSelfObjects(0);
    */

    cameraMoveSpeed = 1.0;
    zoom = 160.0;
    //world.perFrame = 1;
    //world.dt = 0.000001;
    //world.update_world();
    //world.debugFlags |= Shooter::DEBUG_FLAG_SHOTS;

    printf( "SETUP DONE \n"  );
}

void TestAppSweepAndPrune::draw(){
   // printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	for( Object3d*  o: world.objects  ){ ((SAPTestObj*)o)->ncol=0; };
	for( Vehicle3d* o: world.vehicles ){ ((SAPTestVeh*)o)->ncol=0; };

    setToRay( (Vec3d)cam.pos, (Vec3d)cam.rot.c, world.vehicles[ipicked]->pos);
    sapVeh.updateObjects();
    sapVeh.sort(); // use some fast sort-update ?
    sapVeh.collideSelfObjects(0);
    sapVeh.collideCrossObjects(0, sapObj.sweeps[0] );
    //SAPbuff::collideCrossObjects( sapVeh.sweeps[0], sapObj.sweeps[0] );
    //SAPbuff::collideCrossObjects( sapObj.sweeps[0], sapVeh.sweeps[0] );

    glColor3f(0.0,0.5,0.0);
    for( Object3d* o : world.objects ){
        if(o->type == &objType1 ) ((SAPTestObj*)o)->draw();
    }

    glColor3f(0.0,0.0,1.0);
    for( Vehicle3d* o : world.vehicles ){
        if(o->type == &vehType1 ) ((SAPTestVeh*)o)->draw();
    }
};

void TestAppSweepAndPrune::drawHUD(){
    glDisable ( GL_LIGHTING );
    glColor3f( 0.0f, 1.0f, 0.0f );

    glColor3f(0.0,0.5,0.0); drawSweep( sapObj.sweeps[0], 200, 3, WIDTH );
    glColor3f(0.0,0.0,1.0); drawSweep( sapVeh.sweeps[0], 400, 3, WIDTH );
}

void TestAppSweepAndPrune::eventHandling ( const SDL_Event& event  ){
    /*
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
            }
            break;
    };
    */
    AppSDL2OGL::eventHandling( event );
}


// ===================== MAIN

TestAppSweepAndPrune * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppSweepAndPrune( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















