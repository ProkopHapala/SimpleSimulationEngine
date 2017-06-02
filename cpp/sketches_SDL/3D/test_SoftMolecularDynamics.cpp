
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

#include "raytrace.h"
#include "Molecule.h"
#include "MMFF.h"
#include "DynamicOpt.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

// ==========================
// TestAppSoftMolDyn
// ==========================

/*
constexpr int natoms=5, nbonds=4, nang=6, ntors=0;

Vec3d  apos[natoms] = {
{ 0.0, 0.0, 0.0},
{+1.0,+1.0,+1.0},
{-1.0,-1.0,+1.0},
{+1.0,-1.0,-1.0},
{-1.0,+1.0,-1.0}
};   // atomic position

Vec2i  bond2atom[nbonds] = {{0,1},{0,2},{0,3},{0,4}};
//double bond_0   [nbonds] = {1.0,1.2,1.5,1.7};  // [A]
double bond_0   [nbonds] = {1.0,1.0,1.0,1.0};  // [A]
Vec2i  ang2bond [nang]   = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
*/

/*
constexpr int natoms=4, nbonds=3, nang=3, ntors=0;

Vec3d  apos[natoms] = {
{ 0.0, 0.0, 0.0},
{ 0.0, 0.0,+1.0},
{-1.0, 0.0, 0.0},
{+1.0, 0.0,-1.0},
};   // atomic position

Vec2i  bond2atom[nbonds] = {{0,1},{0,2},{0,3}};
double bond_0   [nbonds] = {1.0,1.0,1.0};  // [A]
Vec2i  ang2bond [nang]   = {{0,1},{1,2},{2,0}};
*/

/*
constexpr int natoms=3, nbonds=2, nang=1, ntors=0;

Vec3d  apos[natoms] = {
{ 0.0, 0.0, 0.0},
{+1.0,+1.0,0.0},
{-1.0,-1.0,0.0},
};   // atomic position

Vec2i  bond2atom[nbonds] = {{0,1},{0,2}};
double bond_0   [nbonds] = {0.8,1.2};  // [A]
Vec2i  ang2bond [nang]   = {{0,1}};
*/

// test closest point skew-lines
//Vec3d   p1,p2,dp1,dp2;
//double  t1,t2;

class TestAppSoftMolDyn : public AppSDL2OGL_3D {
	public:
	Molecule    mol;
	MMFFparams  params;
    MMFF        world;
    DynamicOpt  opt;

    int ogl_sph;
    Vec3d ray0;
    int ipicked  = -1, ibpicked = -1;
    int perFrame = 10.0;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppSoftMolDyn( int& id, int WIDTH_, int HEIGHT_ );

};


TestAppSoftMolDyn::TestAppSoftMolDyn( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    params.loadBondTypes("common_resources/BondTypes.dat");

    mol.loadMol("common_resources/propylacid.mol");
    mol.bondsOfAtoms();   mol.printAtom2Bond();
    mol.autoAngles();

    Vec3d cog = mol.getCOG_av();
    mol.addToPos( cog*-1.0d );

    world.apos      = mol.pos;
    world.bond2atom = mol.bond2atom;
    world.ang2bond  = mol.ang2bond;
    world.allocate( mol.natoms, mol.nbonds, mol.nang, 0 );
    world.ang_b2a();

    //exit(0);

    params.fillBondParams( world.nbonds, world.bond2atom, mol.bondType, mol.atomType, world.bond_0, world.bond_k );
    world.printBondParams();
    //exit(0);

    /*
    world.apos      = apos;
    world.bond2atom = bond2atom;
    world.ang2bond  = ang2bond;
    world.bond_0    = bond_0;
    world.allocate( natoms, nbonds, nang, ntors );
    world.ang_b2a();
    */

    opt.bindArrays( 3*world.natoms, (double*)world.apos, new double[3*world.natoms], (double*)world.aforce );
    opt.setInvMass( 1.0 );

    for(int i=0; i<world.nbonds; i++){
        world.bond_k[i] = 2.0;
    }

    for(int i=0; i<world.nang; i++){
        world.ang_0[i] = {1.0,0.0};
        world.ang_k[i] = 0.5;
        //Vec2i ib = world.ang2bond[i];
        //world.ang2atom [i] = (Vec3i){ world.bond2atom[ib.x].y, world.bond2atom[ib.y].y, world.bond2atom[ib.y].x };
    }

    ogl_sph = glGenLists(1);
    glNewList( ogl_sph, GL_COMPILE );
        //glEnable( GL_LIGHTING );
        //glColor3f( 0.8f, 0.8f, 0.8f );
        //Draw3D::drawSphere_oct(3, 0.5, {0.0,0.0,0.0} );
        Draw3D::drawSphere_oct(1, 0.1, {0.0,0.0,0.0} );
    glEndList();

    /*
    // test closest point skew-lines
    srand(548);
    p1  = (Vec3d){randf(-1.0,1.0),randf(-1.0,1.0),randf(-1.0,1.0)};
    p2  = (Vec3d){randf(-1.0,1.0),randf(-1.0,1.0),randf(-1.0,1.0)};
    dp1 = (Vec3d){randf(-1.0,1.0),randf(-1.0,1.0),randf(-1.0,1.0)};  dp1.normalize();
    dp2 = (Vec3d){randf(-1.0,1.0),randf(-1.0,1.0),randf(-1.0,1.0)};  dp2.normalize();
    rayLine( p1, dp1, p2, dp2, t1, t2 );
    //t2  = rayLine( p2, dp2, p1, dp1 );
    */

}

void TestAppSoftMolDyn::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	/*
	// test closest point skew-lines
	glColor3f(0.8f,0.0f,0.0f); Draw3D::drawLine( p1+(dp1*3.0), p1+(dp1*-3.0) );
	glColor3f(0.0f,0.0f,0.8f); Draw3D::drawLine( p2+(dp2*3.0), p2+(dp2*-3.0) );
	glColor3f(0.0f,0.0f,0.0f); Draw3D::drawLine( p1+(dp1*t1), p2+(dp2*t2) );
	printf( "t1 %g t2 %g \n", t1, t2);
	//Vec3d rxl; rxl.set_cross(dp1, dp2);
	//glColor3f(0.8f,0.0f,0.0f); Draw3D::drawVecInPos( rxl, p1+(dp1*t1) );
	//glColor3f(0.0f,0.0f,0.8f); Draw3D::drawVecInPos( rxl, p2+(dp2*t2) );
	return;
	*/

    ray0 = camMat.a*mouse_begin_x + camMat.b*mouse_begin_y;
    Draw3D::drawPointCross( ray0, 0.1 );
    //Draw3D::drawVecInPos( camMat.c, ray0 );
    if(ipicked>=0) Draw3D::drawLine( world.apos[ipicked], ray0);

	double F2;
	for(int itr=0; itr<perFrame; itr++){
        for(int i=0; i<world.natoms; i++){ world.aforce[i].set(0.0d); }
        world.eval_bonds();
        //world.eval_angles();
        world.eval_angcos();

        //exit(0);
        if(ipicked>=0){
            Vec3d f = getForceSpringRay( world.apos[ipicked], camMat.c, ray0, -1.0 );
            //printf( "f (%g,%g,%g)\n", f.x, f.y, f.z );
            world.aforce[ipicked].add( f );
        };

        //for(int i=0; i<world.natoms; i++){ world.aforce[i].add({0.0,-0.01,0.0}); }
        int ipivot = 0;
        world.aforce[ipivot].set(0.0);
        //opt.move_LeapFrog(0.01);
        //opt.move_MDquench();
        F2 = opt.move_FIRE();
        //exit(0);
    }

    //printf( "==== frameCount %i  |F| %g \n", frameCount, sqrt(F2) );

    for(int i=0; i<world.natoms; i++){
        //glColor3f(0.0f,0.0f,0.0f); Draw3D::drawPointCross(world.apos[i],0.2);
        glColor3f(1.0f,0.0f,0.0f); Draw3D::drawVecInPos(world.aforce[i]*30.0,world.apos[i]);

        //glCallList( ogl_sph );
        glEnable(GL_LIGHTING);
        Mat3d mat;
        mat.setOne();
        //mat.mul();
        glColor3f(0.8f,0.8f,0.8f);
        Draw3D::drawShape(world.apos[i],mat,ogl_sph);
        glDisable(GL_LIGHTING);
    }
    for(int i=0; i<world.nbonds; i++){
        Vec2i ib = world.bond2atom[i];
        glColor3f(0.0f,0.0f,0.0f);
        if(i==ibpicked) glColor3f(1.0f,0.0f,0.0f); ;
        Draw3D::drawLine(world.apos[ib.x],world.apos[ib.y]);
    }

};


void TestAppSoftMolDyn::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;


                case SDLK_a: world.apos[1].rotate(  0.1, {0.0,0.0,1.0} ); break;
                case SDLK_d: world.apos[1].rotate( -0.1, {0.0,0.0,1.0} ); break;
                case SDLK_w: world.apos[1].mul( 1.1 ); break;
                case SDLK_s: world.apos[1].mul( 0.9 ); break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    ipicked = pickParticle( world.natoms, world.apos, ray0, camMat.c , 0.5 );
                    break;
                case SDL_BUTTON_RIGHT:
                    ibpicked = world.pickBond( ray0, camMat.c , 0.5 );
                    printf("ibpicked %i \n", ibpicked);
                    break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    ipicked = -1;
                    break;
                case SDL_BUTTON_RIGHT:
                    ibpicked = -1;
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
















