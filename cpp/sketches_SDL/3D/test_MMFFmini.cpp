
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "testUtils.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw3D.h"
#include "SDL_utils.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"

#include "raytrace.h"
//#include "Molecule.h"
#include "MMFFmini.h"
#include "NBFF.h"
#include "MMFFBuilder.h"
#include "DynamicOpt.h"

//#include "NBSRFF.h"
//#include "IO_utils.h"

#include "AppSDL2OGL_3D.h"


// ==========================
// TestAppSoftMolDyn
// ==========================

void plotSurfPlane( Vec3d normal, double c0, Vec2d d, Vec2i n ){
    Vec3d da,db;
    normal.getSomeOrtho( da,db );
    da.mul( d.a/da.norm() );
    db.mul( d.b/db.norm() );
    //glColor3f(1.0f,0.0f,0.0f); Draw3D::drawVecInPos(normal, {0.0,0.0,0.0} );
    //glColor3f(0.0f,1.0f,0.0f); Draw3D::drawVecInPos(da*10, {0.0,0.0,0.0} );
    //glColor3f(0.0f,0.0f,1.0f); Draw3D::drawVecInPos(db*10, {0.0,0.0,0.0} );
    Draw3D::drawRectGridLines( n*2, (da*-n.a)+(db*-n.b) + normal*c0, da, db );
}

void drawAngle(int i, MMFFmini& ff){
    //Vec2i b=ff.ang2bond[i];
    float sz=0.1;
    Vec3i a=ff.ang2atom[i];
    Draw3D::drawArrow(ff.apos[a.x],ff.apos[a.y],sz);
    Draw3D::drawArrow(ff.apos[a.y],ff.apos[a.z],sz);
    //if(b.i&SIGN_MASK){ Draw3D::drawArrow(apos[a.x],apos[a.y]); }else{ Draw3D::drawArrow(apos[a.x],apos[a.y]); };
    //if(b.j&SIGN_MASK){   };
}

void drawTors(int i, MMFFmini& ff){
    float sz=0.1;
    //Vec2i b=ff.ang2bond[i];
    Quat4i t=ff.tors2atom[i];
    Draw3D::drawArrow(ff.apos[t.x],ff.apos[t.y],sz);
    Draw3D::drawArrow(ff.apos[t.y],ff.apos[t.z],sz);
    Draw3D::drawArrow(ff.apos[t.z],ff.apos[t.w],sz);
    //if(b.i&SIGN_MASK){ Draw3D::drawArrow(apos[a.x],apos[a.y]); }else{ Draw3D::drawArrow(apos[a.x],apos[a.y]); };
    //if(b.j&SIGN_MASK){   };
}

class TestAppSoftMolDyn : public AppSDL2OGL_3D {
	public:
	Molecule    mol;
	//MMFFparams  params;
    MMFFmini    ff;
    NBFF       nff;
    MMFFBuilder builder;
    DynamicOpt  opt;

    int     fontTex;
    int     ogl_sph;

    char str[256];

    Vec3d ray0;
    int ipicked  = -1, ibpicked = -1;
    int perFrame =  50;

    double drndv =  10.0;
    double drndp =  0.5;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppSoftMolDyn( int& id, int WIDTH_, int HEIGHT_ );

};



#include "NDarray.h"
using namespace NDarray;


TestAppSoftMolDyn::TestAppSoftMolDyn( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {


    /*
    printf( "----\n" );
    (1<9_r).print();
    (1<9_s).print();

    //( 4_s | 1>>9/-2_s ).print();

    printf( "----\n" );
    (-5_s).print();
    (5/-1_s).print();
    (1<5/-1_s).print();
    //(1,9,-1_s).print();
    printf( "----\n" );
    //(3>>-6_r/2|Scan(4)|Scan(1,9,-1)).print();
    (3<-6_r/2 , 4_s , 1<9/-1_s).print();
    (3<-6_r/2 , 4_s , 1<9/-1_s, 3<19/3_s).print();
    */

    VecN<double> v{1.0d,1.3d,1.5d};

    v.print();
    v.apply( [](double& f,double _)->double{ f=f*f; return 0; } );

    std::function<void(double&,double&)> f = o_max<double>;
    v.apply_(f);
    v.apply_(o_max<double>);

    //print( v.apply_( o_max ) );
    v.print();


    //A* pa = new A{ A{1}+A{3}  };
    //pa->print();


    //(A{1},A{3}).print();

    //printf( "(%i,%i)\n", t.ts[0].i, t.ts[1].i );


    //(Scan(1)|Scan(1)).print();
    //VecN a{0.5,1,0,3,5,6,65848,8487,0.1};
    //exit(0);












    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );

    double l0    = 1.5;
    double Kbond = 10.0;
    double Kang  = 3.0;
    double Ktors = 1.0;

    int tors_n = 3;

    const int natom=4,nbond=3,nang=2,ntors=1;
    Vec3d apos0[] = {
        {-2.0,0.0,0.0},  // 0
        {-1.0,2.0,0.0},  // 1
        {+1.0,2.0,0.0},  // 2
        {+2.0,0.0,1.0},  // 3
    };
    Vec2i bong2atom[] = {
        {0,1},  // 0
        {1,2},  // 1
        {2,3},  // 2
    };
    Vec2i ang2bond[] = {
        {0,1},  // 0
        {1,2},  // 1
    };
    Vec3i tors2bond[] = {
        {0,1,2},  // 0
    };
    double a0s[] ={
        2.0,
        2.0,
        2.0,
    };

    /*
    const int natom=7,nbond=7,nang=10;
    Vec3d apos0[] = {
        {-2.0,0.0,0.0},  // 0
        {-1.0,2.0,0.0},  // 1
        {+1.0,2.0,0.0},  // 2
        {+2.0,0.0,0.0},  // 3
        {+0.0,-1.0,0.0},  // 4

        {+0.0,0.0,+1.0},   // 5
        {+0.0,0.0,-1.0}   // 6
    };
    Vec2i bong2atom[] = {
        {0,1},  // 0
        {1,2},  // 1
        {2,3},  // 2
        {3,4},  // 3
        {4,0},  // 4

        {5,0},  // 5
        {6,0},  // 6
    };
    Vec2i ang2bond[] = {
        {0,1},  // 0
        {1,2},  // 1
        {2,3},  // 2
        {3,4},  // 3
        {4,0},  // 4

        {5,6},  // 5

        {0,5},  // 6
        {0,6},  // 7
        {4,5},  // 8
        {4,6}   // 9
    };
    double a0s[] ={
        2.0,
        2.0,
        2.0,
        2.0,
        1.0,

        2.0,

        2.0,
        2.0,
        2.0,
        2.0
    };
    */

    /*
    const int natom=5+2,nbond=4+3,nang=6;
    Vec3d apos0[] = {
        { 0.5, 0.5, 0.5},  // 0
        {-1.0,+1.0,+1.0},  // 1
        {+1.0,-1.0,+1.0},  // 2
        {-1.0,-1.0,-1.0},  // 3
        {+1.0,+1.0,-1.0},  // 4

        {-1.0,-1.0,-2.0},  // 5
        {+1.0,+1.0,-2.0}   // 6

        //{1.0,0.0,0.0},  // 1
        //{0.0,1.0,0.0},  // 2
        //{0.0,0.0,1.0},  // 3
        //{-1.0,-1.0,-1.0}   // 4
    };
    Vec2i bong2atom[] = {
        {0,1},  // 0
        {0,2},  // 1
        {0,3},  // 2
        {0,4},  // 3

        {5,6},  // 4
        {3,5},  // 5
        {4,6}   // 6

    };
    Vec2i ang2bond[] = {
        {0,1},
        {0,2},
        {0,3},
        {1,2},
        {1,3},
        {2,3}
    };
    double a0s[] ={
        2.0,
        2.0,
        2.0,
        2.0,
        2.0,
        2.0
    };
    */

    // ============ Build molecule

    MMFFAtom brushAtom{ };
    MMFFBond brushBond{ };

    for(int i=0;i<natom;i++){
        brushAtom.pos = apos0[i];
        builder.insertAtom(brushAtom, true);
    }
    for(int i=0;i<nbond;i++){
        brushBond.atoms=bong2atom[i];
        builder.insertBond(brushBond);
    }
    for(int i=0;i<natom;i++){
        builder.makeSPConf(i,0,0);
        //builder.makeSPConf(i);
    }



    //builder.toMMFFmini( &ff );

    //exit(0);

    ff.realloc(natom,nbond,nang,ntors);


    printf( "DEBUG 1 \n" );
    for(int i=0; i<ff.natoms; i++){
        ff.apos[i] = apos0[i];
    }
    printf( "DEBUG 2 \n" );
    for(int i=0; i<ff.nbonds; i++){
        ff.bond2atom[i] = bong2atom[i];
        ff.bond_k [i] = Kbond;
        ff.bond_l0[i] = l0;
    }
    printf( "DEBUG 3 \n" );
    for(int i=0; i<ff.nang; i++){
        ff.ang2bond[i] = ang2bond[i];
        double a0 = -a0s[i]/2.0; // NOTE: we use half-angle
        ff.ang_cs0[i] = { cos(a0), sin(a0) };
        ff.ang_k  [i] = Kang;
    }
    for(int i=0; i<ff.ntors; i++){
        ff.tors2bond[i] = tors2bond[i];
        ff.tors_k   [i] = Ktors;
        ff.tors_n   [i] = tors_n;
    }
    ff.angles_bond2atom  ();
    ff.torsions_bond2atom();
    printf( "DEBUG 4 \n" );




    nff.bindOrRealloc( ff.natoms, ff.nbonds, ff.apos, ff.aforce, 0, ff.bond2atom );
    //nff.setREQs(0,nff.n, {1.4,0.0}  )
    nff.setREQs(0       ,nff.n,{1.9080,sqrt(0.003729),0});

    opt.bindOrAlloc( 3*ff.natoms, (double*)ff.apos, 0, (double*)ff.aforce, 0 );
    //opt.setInvMass( 1.0 );
    opt.cleanVel( );

    printf( "DEBUG 5 \n" );

    ogl_sph = glGenLists(1);
    glNewList( ogl_sph, GL_COMPILE );
        //glEnable( GL_LIGHTING );
        //glColor3f( 0.8f, 0.8f, 0.8f );
        //Draw3D::drawSphere_oct(3, 0.5, {0.0,0.0,0.0} );
        Draw3D::drawSphere_oct( 2, 0.25, {0.0,0.0,0.0} );
    glEndList();

}

void TestAppSoftMolDyn::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );


	/*
	//ibpicked = world.pickBond( ray0, camMat.c , 0.5 );
    ray0 = (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y);
    Draw3D::drawPointCross( ray0, 0.1 );
    //Draw3D::drawVecInPos( camMat.c, ray0 );
    if(ipicked>=0) Draw3D::drawLine( world.apos[ipicked], ray0);
    */


    /*
    for(int i=0; i<50; i++){
        ff.apos[1]={0.0,0.0};
        ff.apos[0]={1.0,0.0};
        double a = i*2*M_PI/50;
        ff.apos[2]={cos(a),sin(a)};

        ff.eval();

        glColor3f(1.0,1.0,1.0);
        Draw3D::drawVecInPos( (ff.apos[0]-ff.apos[1])*5.0, ff.apos[1] );
        Vec2d cs = ff.ang_cs0[0];
        cs.set_mul_cmplx(cs,cs);
        Draw3D::drawVecInPos( ((Vec3d){cs.x,-cs.y,0.0})*5.0, ff.apos[1] );
    */

    //ff.apos[0].set(-2.0,0.0,0.0);
    perFrame = 1;
    //perFrame = 50;
	double F2;
	for(int itr=0; itr<perFrame; itr++){
        printf( "======= frame %i \n", frameCount );

	    //printf( "DEBUG run 1 \n" );

	    // rotate arom[0]
        //ff.apos[0] = ff.apos[1] + (ff.apos[0]-ff.apos[1]).rotate( 2*M_PI/perFrame, ff.apos[2]-ff.apos[1] );

        ff.eval();
        nff.evalLJQ_sortedMask();

        //for(int i=0; i<world.natoms; i++){ world.aforce[i].set(0.0d); }
        //printf( "DEBUG x.1 \n" );
        //world.eval_bonds(true);
        //world.eval_angles();
        //printf( "DEBUG x.2 \n" );
        //world.eval_angles();
        //printf( "DEBUG x.3 \n" );
        //world.eval_LJq_On2();

        /*
        //exit(0);
        if(ipicked>=0){
            Vec3d f = getForceSpringRay( world.apos[ipicked], (Vec3d)cam.rot.c, ray0, -1.0 );
            //printf( "f (%g,%g,%g)\n", f.x, f.y, f.z );
            world.aforce[ipicked].add( f );
        };


        for(int i=0; i<world.natoms; i++){
            world.aforce[i].add( getForceHamakerPlane( world.apos[i], {0.0,0.0,1.0}, -3.0, 0.3, 2.0 ) );
            //printf( "%g %g %g\n",  world.aforce[i].x, world.aforce[i].y, world.aforce[i].z );
        }
        */

        //exit(0);

        //for(int i=0; i<world.natoms; i++){ world.aforce[i].add({0.0,-0.01,0.0}); }
        //int ipivot = 0;
        //world.aforce[ipivot].set(0.0);
        //opt.move_LeapFrog(0.01);
        //opt.move_MDquench();

        opt.move_GD(0.01);
        //F2 = opt.move_FIRE();
        //exit(0);

    }


    glColor3f(0.0,1.0,0.0);
    //drawAngle(0, ff);
    //drawTors (0, ff);

    //glColor3f(0.6f,0.6f,0.6f); plotSurfPlane( (Vec3d){0.0,0.0,1.0}, -3.0, {3.0,3.0}, {20,20} );
    //Draw3D::drawVecInPos( (Vec3d){0.0,0.0,1.0},  (Vec3d){0.0,0.0,0.0} );

    //printf( "==== frameCount %i  |F| %g \n", frameCount, sqrt(F2) );
    //printf( "DEBUG run 2 \n" );
    // draw Bonds
    glColor3f(0.0f,0.0f,0.0f);
    /*
    for(int i=0; i<ff.nbonds; i++){
        Vec2i ib = ff.bond2atom[i];
        glColor3f(0.0f,0.0f,0.0f);
        //if(i==ibpicked) glColor3f(1.0f,0.0f,0.0f);
        Draw3D::drawLine(ff.apos[ib.x],ff.apos[ib.y]);
        sprintf(str,"%i\0",i);
        //Draw3D::drawText(str, (world.apos[ib.x]+world.apos[ib.y])*0.5, fontTex, 0.02, 0,0);
        Draw3D::drawText(str, (ff.apos[ib.x]+ff.apos[ib.y])*0.5, fontTex, 0.02, 0);
    }
    */
    //};
    // draw Atoms
    double fsc = 1.0;
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);
    for(int i=0; i<ff.natoms; i++){
        //printf( "apos[%i] (%g,%g,%g)\n", i, ff.apos[i].x,ff.apos[i].y,ff.apos[i].z );
        //glColor3f(0.0f,0.0f,0.0f); Draw3D::drawPointCross(world.apos[i],0.2);
        //printf( "aforce[%i] (%g,%g,%g) \n", i,  ff.aforce[i].x, ff.aforce[i].y, ff.aforce[i].z );
        //glColor3f(0.0f,0.0f,1.0f); Draw3D::drawVecInPos(ff.aforce[i]*fsc,ff.apos[i]);
        glColor3f(0.0f,0.0f,0.0f); Draw3D::drawPointCross(ff.apos[i],0.1);
        /*
        //glCallList( ogl_sph );
        glEnable(GL_LIGHTING);
        Mat3d mat;
        mat.setOne();
        //mat.mul();
        glColor3f(0.8f,0.8f,0.8f);
        Draw3D::drawShape(ff.apos[i],mat,ogl_sph);
        glDisable(GL_LIGHTING);
        */
    }
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);

    /*
    printf("==========\n");
    for(int i=0; i<world.natoms; i++){
        printf("iatom %i (%g,%g,%g) (%g,%g,%g) \n", i, world.apos[i].x,world.apos[i].y,world.apos[i].z, world.aforce[i].x,world.aforce[i].y,world.aforce[i].z  );
    }
    if(frameCount>=10){STOP = true;}
    */

};


void TestAppSoftMolDyn::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_p:  first_person = !first_person; break;
                //case SDLK_o:  perspective  = !perspective; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;

                //case SDLK_v: for(int i=0; i<world.natoms; i++){ ((Vec3d*)opt.vel)[i].add(randf(-drndv,drndv),randf(-drndv,drndv),randf(-drndv,drndv)); } break;
                //case SDLK_p: for(int i=0; i<world.natoms; i++){ world.apos[i].add(randf(-drndp,drndp),randf(-drndp,drndp),randf(-drndp,drndp)); } break;

                //case SDLK_LEFTBRACKET:  if(ibpicked>=0) world.bond_0[ibpicked] += 0.1; break;
                //case SDLK_RIGHTBRACKET: if(ibpicked>=0) world.bond_0[ibpicked] -= 0.1; break;

                //case SDLK_a: world.apos[1].rotate(  0.1, {0.0,0.0,1.0} ); break;
                //case SDLK_d: world.apos[1].rotate( -0.1, {0.0,0.0,1.0} ); break;

                //case SDLK_w: world.apos[1].mul( 1.1 ); break;
                //case SDLK_s: world.apos[1].mul( 0.9 ); break;

            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                //case SDL_BUTTON_LEFT:
                //    ipicked = pickParticle( world.natoms, world.apos, ray0, (Vec3d)cam.rot.c , 0.5 );
                //    break;
                //case SDL_BUTTON_RIGHT:
                //    ibpicked = world.pickBond( ray0, (Vec3d)cam.rot.c , 0.5 );
                //    printf("ibpicked %i \n", ibpicked);
                //    break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    ipicked = -1;
                    break;
                case SDL_BUTTON_RIGHT:
                    //ibpicked = -1;
                    break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}

void TestAppSoftMolDyn::drawHUD(){
    glDisable ( GL_LIGHTING );

}

// ===================== MAIN

TestAppSoftMolDyn * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppSoftMolDyn( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















