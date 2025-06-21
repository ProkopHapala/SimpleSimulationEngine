
/// @file @brief  This is a creative modeling tool for generating organic, muscle-like shapes. It uses Hermite splines (`spline_hermite.h`) to define cross-sections which are then "lofted" to form a smooth 3D surface. The shape can be interactively edited by manipulating the spline control points using the `EditorGizmo.h` interface, allowing for intuitive, real-time sculpting.
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

#include "EditorGizmo.h"
#include "GUI.h"

#include "spline_hermite.h"


#include "AppSDL2OGL_3D.h"
#include "testUtils.h"


struct MusculeNode{
    Vec3d pos;
    Vec3d dir;
    Vec3d up;
    Vec3d width;

    void dir2up(){ dir.normalize(); up .makeOrtho(dir); up .normalize(); }
    void up2dir(){ up .normalize(); dir.makeOrtho(up ); dir.normalize(); }
};

class Muscule{ public:
    std::vector<MusculeNode> nodes;

    Muscule( ){};
    Muscule(int n):nodes(n){};
    Muscule(int n, Vec3d p0, Vec3d p1, double w0, double w1):nodes(n){
        initPoints( p0, p1, w0, w1 );
    }


    void dir2up(){ for(MusculeNode& o:nodes){o.dir2up();} }
    void up2dir(){ for(MusculeNode& o:nodes){o.up2dir();} }

    void initPoints( Vec3d p0, Vec3d p1, double w0, double w1=-1 ){
        int n = nodes.size();
        Vec3d  dir = p1-p0;
        double dL  = dir.normalize()/n;
        Vec3d up,side; dir.getSomeOrtho(up,side);
        if(w1<0){w1=w0;}
        double dw=(w1-w0)/n;
        for(int i=0; i<n; i++){
            nodes[i].pos = p0 + dir*(dL*i);
            nodes[i].dir = dir;
            nodes[i].up  = up;
            double w = w0+dw*i;
            nodes[i].width={w,w,1};
        }
    }

    void renderControls(){
        glBegin(GL_LINES);
        for(int i=0; i<nodes.size(); i++){
            const MusculeNode& o = nodes[i];
            Vec3d side; side.set_cross(o.dir, o.up); side.normalize();
            //printf( "node[%i] p(%g,%g,%g) dir(%g,%g,%g) up(%g,%g,%g)\n", o.pos.x,o.pos.y,o.pos.z,  o.dir.x,o.dir.y,o.dir.z,  o.up.x,o.up.y,o.up.z );
            glColor3f(0.,0,1); Draw3D::vertex(o.pos); Draw3D::vertex(o.pos + o.dir);
            glColor3f(0.,1,0); Draw3D::vertex(o.pos); Draw3D::vertex(o.pos + o.up*o.width.a);
            glColor3f(1.,0,0); Draw3D::vertex(o.pos); Draw3D::vertex(o.pos + side*o.width.b);
        }
        glEnd();
    }

    void renderSpine( int nv, bool bUp ){
        Mat3d rot;
        double inuv=1./nv;
        glBegin(GL_LINE_STRIP);
        //glBegin(GL_TRIANGLE_STRIP);
        int nnd=nodes.size();
        for(int i=1; i<nnd; i++){
            const MusculeNode& o0 = nodes[i-1];
            const MusculeNode& o1 = nodes[i  ];
            Vec3d op = o0.pos;
            Draw3D::vertex(op);
            for(int iv=1; iv<=nv; iv++){
                double t  = inuv*iv;
                double mt = 1-t;
                double c0,c1,d0,d1;
                Spline_Hermite::basis( t, c0,c1,d0,d1 );
                //rot.c = (o1.pos - o0.pos).normalized(); // dir_z azis
                Vec3d p = o1.pos*c1 + o0.pos*c0 + o1.dir*d1 + o0.dir*d0;
                rot.c = (p-op).normalized();


                Vec3d upd0=o1.up*0.5, upd1=o0.up*-0.5;
                Vec3d wd0=o1.width*0.5, wd1=o0.width*-0.5;
                if(i<(nnd-1)){ upd0.add_mul(nodes[i+1].up, 0.5); wd1.add_mul(nodes[i+1].width, 0.5); }
                if(i>1      ){ upd0.add_mul(nodes[i-2].up,-0.5); wd0.add_mul(nodes[i-2].width,-0.5); }
                rot.b = o1.up*c1 + o0.up*c0 + upd1*d1 + upd0*d0;
                Vec3d w     = o1.width*c1 + o0.width*c0 + wd1*d1 + wd0*d0;

                rot.b.makeOrtho(rot.c);
                rot.b.normalize();
                rot.a.set_cross(rot.b,rot.c);

                glColor3f(0,0,1);
                Draw3D::vertex(p);
                if(bUp){
                    glColor3f(0,1,0);
                    Draw3D::vertex(p+rot.b*w.b);
                    Draw3D::vertex(p   );
                    glColor3f(1,0,0);
                    Draw3D::vertex(p+rot.a*w.a);
                    Draw3D::vertex(p   );
                }
                op=p;
            }
        }
        glEnd();
    }


    void render( Vec2i nuv ){
        Mat3d rot;
        Vec2d dcrot; dcrot.fromAngle(2*M_PI/nuv.x);
        double inuv=1./nuv.y;
        int nnd=nodes.size();
        //glBegin(GL_POINTS);
        //glBegin(GL_TRIANGLE_STRIP);
        //glBegin(GL_LINE_LOOP);
        Vec3d ops[nuv.x+1];
        Vec3d ons[nuv.x+1];
        for(int i=1; i<nnd; i++){
            const MusculeNode& o0 = nodes[i-1];
            const MusculeNode& o1 = nodes[i  ];

            Vec3d op = o0.pos;
            for(int iv=0; iv<nuv.y; iv++){
                double t = inuv*iv;
                double c0,c1,d0,d1;
                Spline_Hermite::basis( t, c0,c1,d0,d1 );
                //rot.c = (o1.pos - o0.pos).normalized(); // dir_z azis
                Vec3d p = o1.pos*c1 + o0.pos*c0 + o1.dir*d1 + o0.dir*d0;
                Vec3d w;
                if(iv>0){
                    rot.c = (p-op).normalized();
                    Vec3d upd0=o1.up*0.5,   upd1=o0.up*-0.5;
                    Vec3d wd0=o1.width*0.5, wd1=o0.width*-0.5;
                    if(i<(nnd-1)){ upd1.add_mul(nodes[i+1].up, 0.5);  wd1.add_mul(nodes[i+1].width, 0.5); }
                    if(i>1      ){ upd0.add_mul(nodes[i-2].up,-0.5);  wd0.add_mul(nodes[i-2].width,-0.5); }
                    rot.b = o1.up*c1 + o0.up*c0 + upd1*d1 + upd0*d0;
                    //rot.b.normalize();
                    w     = o1.width*c1 + o0.width*c0 + wd1*d1 + wd0*d0;
                }else{
                    rot.c = o0.dir;
                    rot.b = o0.up;
                    w     = o0.width;
                }

                rot.b.makeOrtho(rot.c);
                rot.b.normalize();
                rot.a.set_cross(rot.b,rot.c);

                //rot.c = (o1.dir*c1 + o0.dir*c0 ).normalized(); // up vector
                //rot.a = (o1.up *c1 + o0.up *c0 ).normalized(); // up vector
                //rot.b = (o1.up *c1 + o0.up *c0 ).normalized();
                //Vec3d w = o1.width + o1.width;
                Vec2d crot{1,0};
                //glBegin(GL_LINE_LOOP);
                glBegin(GL_TRIANGLE_STRIP);
                for(int iu=0; iu<=nuv.x; iu++){

                    Vec3d pi = p + rot.a*(w.a*crot.a) + rot.b*(w.b*crot.b);
                    Vec3d ni = rot.a*crot.a + rot.b*crot.b;
                    if(iv>0){ Draw3D::normal( ons[iu] ); Draw3D::vertex(ops[iu]); }
                    Draw3D::normal( ni ); Draw3D::vertex( pi );
                    ops[iu]=pi;ons[iu]=ni;
                    crot.mul_cmplx(dcrot);
                }
                glEnd();
                op=p;
            }
        }
        //glEnd();
    }

};

Muscule mus1(3);


class TestAppMousePicking : public AppSDL2OGL_3D {
	public:
    //MultiFight3DWorld world;
    double dvel = 10.0;

    //std::vector<KinematicBody*> objects;
    int nobject = 100;
    EditorGizmo gizmo;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppMousePicking( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppMousePicking::TestAppMousePicking( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {
    gizmo.cam = &cam;

    mus1.initPoints( {-1.,-2.,0.5}, {1.,3.,0.5}, 1, 0.5 );

    mus1.nodes[0].width=(Vec3d){0.2,0.1,1.};
    mus1.nodes[1].width=(Vec3d){1.0,0.5,1.};
    mus1.nodes[2].width=(Vec3d){0.2,0.1,1.};

    for(MusculeNode& o:mus1.nodes){
        Vec3d v; v.fromRandomCube(1.0);
        //o.dir.add( v );
        o.up.add( v );
    }
    //mus1.dir2up();
    mus1.up2dir();

    /*
    int np = 30;
    gizmo.npoint = np;
    gizmo.points = new Vec3d[np];
    double sz = 5.0;
    for(int i=0; i<gizmo.npoint; i++){
        gizmo.points[i].fromRandomBox( {-sz,-sz,-sz}, {sz,sz,sz} );
    }
    */
}

void TestAppMousePicking::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    glDisable(GL_LIGHTING);
    mus1.renderControls();
    glColor3f(0,0,0);       mus1.renderSpine(5,true);

    glEnable(GL_LIGHTING);
    glColor3f(0.8,0.8,0.8); mus1.render( {8,5} );


	/*
	gizmo.draw();

	Draw3D::drawPoints( gizmo.npoint, gizmo.points ,0.1 );


	for(auto& it: gizmo.selection){
        Vec3d p = gizmo.points[it.first];
        Draw::color_of_hash( it.second + 15454 );
        Draw3D::drawPointCross( p, 0.2 );
	};
	*/

	Draw3D::drawAxis(10.0);

};




void TestAppMousePicking::eventHandling ( const SDL_Event& event  ){
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
    Vec2f pix = { 2*mouseX/float(HEIGHT) - ASPECT_RATIO,
                  2*mouseY/float(HEIGHT) - 1              };
    cam.persp = perspective;
    cam.zoom  = zoom;
    //gizmo.onEvent( pix, event );
}


void TestAppMousePicking::drawHUD(){
    /*
    glDisable ( GL_LIGHTING );
    glColor3f( 0.0f, 1.0f, 0.0f );
    glBegin( GL_LINES );
    float whalf = WIDTH *0.5;
    float hhalf = HEIGHT*0.5;
    glVertex3f( whalf-10,hhalf, 0 ); glVertex3f( whalf+10,hhalf, 0 );
    glVertex3f( whalf,hhalf-10, 0 ); glVertex3f( whalf,hhalf+10, 0 );
    glEnd();
    */
}

// ===================== MAIN

TestAppMousePicking * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppMousePicking( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
