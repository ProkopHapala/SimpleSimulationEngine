
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw3D.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"

#include "Multipoles.h"
#include "grids3D.h"
#include "PBCsystem.h"
#include "MultipoleGrid.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

#define R2SAFE  1.0e-8f

constexpr int nbodies = 16;
double  charges[nbodies];
Vec3d   pos    [nbodies];


inline void addAtomicForceLJQ( const Vec3d& dp, Vec3d& f, double r0, double eps, double q ){
    //Vec3f dp; dp.set_sub( p2, p1 );
    double ir2  = 1/( dp.norm2() + R2SAFE );
    double ir   = sqrt(ir2);
    double ir2_ = ir2*r0*r0;
    double ir6  = ir2_*ir2_*ir2_;
    double fr   = ( ( 1 - ir6 )*ir6*12*eps + ir*q*-14.3996448915f )*ir2;
    f.add_mul( dp, fr );
}

/*
void eval_LJq_On2(){
    for(int i=0; i<natoms; i++){
        Vec3d ljq_i = aLJq[i];
        Vec3d pi    = apos[i];
        Vec3d f; f.set(0.0);
        for(int j=0; j<natoms; j++){
            if(i!=j){
                Vec3d& ljq_j = aLJq[j];
                double rij = ljq_i.x+ljq_j.x;
                double eij = ljq_i.y*ljq_j.y;
                double qq  = ljq_i.z*ljq_j.z;
                addAtomicForceLJQ( pi-apos[j], f, rij, -eij, qq );
            }
        }
        aforce[i].add(f);
    }
}
*/

// ======= THE CLASS

class TestAppMultipoles : public AppSDL2OGL_3D {
	public:

    MultipoleGrid grid;

    PBCsystem system1;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppMultipoles( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppMultipoles::TestAppMultipoles( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {
    double rspan = 5.0;
    for(int i=0; i<nbodies; i++){
        pos    [i].set( randf(-rspan,rspan), randf(-rspan,rspan), randf(-rspan,rspan) );
        charges[i] = randf( 0.5, 1.0 );
    }
    grid.ruler.setStep(8.0);
    grid.allocate( nbodies, {8,8,8} );
    grid.ruler.pos0.set(-32.0,-32.0,-32.0);
    //grid.setn   (8,8,8);
    printf("DEBUG 1\n");
    grid.atomsToCells( nbodies, pos, charges );
    //exit(0);

    // Test PBCsystem
    //system1.lvec.set( {1.0,0.5,0.0},{0.0,1.0,0.5},{0.5,0.0,1.0} );
    //system1.setCell( (Mat3d){ 1.0,0.5,0.0,  0.0,1.0,0.5,   0.5,0.0,1.0 } );
    //system1.setCell( (Mat3d){ 1.0,randf(-0.5,0.5),randf(-0.5,0.5),  randf(-0.5,0.5),1.0,randf(-0.5,0.5),   randf(-0.5,0.5),randf(-0.5,0.5),1.0 } );

    Mat3d lvec = (Mat3d){ 1.0,randf(-0.5,0.5),randf(-0.5,0.5),  randf(-0.5,0.5),1.0,randf(-0.5,0.5),   randf(-0.5,0.5),randf(-0.5,0.5),1.0 };
    lvec.mul(15.0);
    system1.setCell( lvec );

    system1.natoms = 50;
    system1.pos = new Vec3d[system1.natoms];
	for(int i=0; i<system1.natoms; i++){
        system1.pos[i] = system1.lvec.a*randf(0.0,1.0) + system1.lvec.b*randf(0.0,1.0) + system1.lvec.c*randf(0.0,1.0);
	}
    system1.initRuler( 6.0 );

    system1.atomsToCells();

    //printVec( system1.ruler.n );

    //exit(0);
}

void TestAppMultipoles::draw(){
    printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	/*
	for(int i=0; i<nbodies; i++){
        Draw::color_of_hash( 1545464+grid.atom2cell[i]*15464 );
        Draw3D::drawPoint( pos[i] );
        //Draw3D::drawPointCross( pos[i], 0.2 );

        Vec3i ip;
        grid.ruler.i2ixyz( grid.atom2cell[i], ip );
        Draw3D::drawLine( pos[i], grid.ruler.box2pos( ip, (Vec3d){4.0,4.0,4.0} ) );
        //Draw3D::drawLine( ,  )
	}
	*/

	//glColor3f(0.0,0.0,1.0); Draw3D::drawTriclinicBox(system1.lvec, {0.0,0.0,0.0}, {1.0,1.0,1.0} );
	//glColor3f(1.0,0.0,0.0); Draw3D::drawTriclinicBox(system1.ilvec, {0.0,0.0,0.0}, {1.0,1.0,1.0} );


	/*
	Vec3d pmin={-1.0,-2.0,-3.0};
	Vec3d pmax={2.5,1.5,0.5};
	Vec3d cmin,cmax;

	triclinicBounds( system1.ilvec, pmin, pmax, cmin, cmax );

	glColor3f(0.0,0.0,1.0); Draw3D::drawBBox( pmin, pmax );
	glColor3f(1.0,0.0,0.0); Draw3D::drawTriclinicBox(system1.ilvec, cmin, cmax );
	*/

	glColor3f(1.0,1.0,1.0);
	for(int i=0; i<system1.natoms; i++){
        Draw3D::drawPointCross( system1.pos[i] ,0.1);
	}

	/*
	glColor3f(0.0,0.0,0.0);
    for(int i=0; i<system1.nAtomPBC; i++){
        Draw3D::drawPointCross( system1.pbc_pos[i], 0.1 );
	};
	*/

	/*
	glColor3f(0.3,0.3,0.3);
	for(int ia=system1.pbc0.a; ia<system1.pbc1.a; ia++){
        for(int ib=system1.pbc0.b; ib<system1.pbc1.b; ib++){
            for(int ic=system1.pbc0.c; ic<system1.pbc1.c; ic++){
                Draw3D::drawTriclinicBox(system1.lvec, {ia,ib,ic}, {ia+1,ib+1,ic+1} );
            }
        }
    }
    */

    for(int icell=0; icell<system1.ruler.ntot; icell++){
        //int icell = system1.ruler.ixyz2i( {ia,ib,ic} );
        //printf( "\n", icell );

        if(icell != (frameCount/30)%system1.ruler.ntot ) continue;
        Draw::color_of_hash( icell + 25545 );
        int n  = system1.cellNs   [icell];
        int i0 = system1.cell2atom[icell];
        for(int i=0; i<n; i++){
            Draw3D::drawPointCross( system1.pbc_pos[i+i0] ,0.3);
        }

        Vec3i ip; system1.ruler.i2ixyz( icell, ip );
        Vec3d p = system1.ruler.box2pos( ip, {0.0,0.0,0.0} );
        double d = system1.ruler.step;
        Draw3D::drawBBox( p, p+(Vec3d){d,d,d} );
    }

	glColor3f(0.0,0.0,1.0); Draw3D::drawBBox( system1.ruler.pos0, system1.ruler.pmax );
	glColor3f(1.0,0.0,0.0); Draw3D::drawTriclinicBox(system1.lvec, system1.cmin, system1.cmax );


};


void TestAppMultipoles::eventHandling ( const SDL_Event& event  ){
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

void TestAppMultipoles::drawHUD(){
    glDisable ( GL_LIGHTING );
}

// ===================== MAIN

TestAppMultipoles * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppMultipoles( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















