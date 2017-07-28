
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw3D.h"
#include "SDL_utils.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"

#include "raytrace.h"
#include "Molecule.h"
#include "MMFF.h"
#include "DynamicOpt.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"


/*

TO DO:
 - save geom to .xyz
 - save geom to .pdb or .mol
 - add charges (read from .pdb)
 - torsion angles
 - add better model of substrate
 - include rigid body molecules (from MoleculerWorld )
 - Brute force non-bonded interactions by OpenCL
 - add some editation capabilities

 TODO Corrections:
 - vdW distances seems to be too close
 - some bonds too long
 - correct angular forcefield to repdesent kinked groups ( e.g. -OH )
*/


// ==========================
// AppMolecularEditor2
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

class AppMolecularEditor2 : public AppSDL2OGL_3D {
	public:
	Molecule    mol;
	MMFFparams  params;
    MMFF        world;
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

    double  atomSize = 0.25;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	AppMolecularEditor2( int& id, int WIDTH_, int HEIGHT_ );

};

AppMolecularEditor2::AppMolecularEditor2( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    //AtomType atyp;
    //atyp.fromString( "CA 6 4 4 1 2.00 0.09 0x11EEAA" );
    params.loadAtomTypes( "common_resources/AtomTypes.dat" );

    //for(auto kv : params.atypNames) { printf( ">>%s<< %i \n", kv.first.c_str(), kv.second ); };


    char str[1024];
    printf( "type %s \n", params.atypes[ params.atypNames.find( "C" )->second ].toString( str ) );
    printf( "type %s \n", params.atypes[ params.atypNames.find( "H" )->second ].toString( str ) );
    printf( "type %s \n", params.atypes[ params.atypNames.find( "O" )->second ].toString( str ) );
    printf( "type %s \n", params.atypes[ params.atypNames.find( "N" )->second ].toString( str ) );
    /*
    auto it = params.atypNames.find( "C" );
    if( it != params.atypNames.end() ){
        //printf( "type CA %i \n", it->second );
        printf( "type %i %s \n", it->second, params.atypes[ it->second ].toString( str ) );
    }else{
        printf("not found\n");
    }
    */

    mol.atypNames = &params.atypNames;

    //exit(0);

    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );

    params.loadBondTypes("common_resources/BondTypes.dat");

    //mol.loadMol("common_resources/propylacid.mol");
    //mol.loadMol("common_resources/precursor_OH.mol");
    mol.loadMol("common_resources/precursor_CN.mol");
    mol.bondsOfAtoms();   mol.printAtom2Bond();
    mol.autoAngles();

    Vec3d cog = mol.getCOG_av();
    mol.addToPos( cog*-1.0d );

    /*
    world.apos      = mol.pos;
    world.bond2atom = mol.bond2atom;
    world.ang2bond  = mol.ang2bond;
    world.allocate( mol.natoms, mol.nbonds, mol.nang, 0 );
    world.ang_b2a();
    //params.fillBondParams( world.nbonds, world.bond2atom, mol.bondType, mol.atomType, world.bond_0, world.bond_k );
    */

    //Vec3d pos = (Vec3d){0.0,0.0,0.0};
    Mat3d rot; rot.setOne();
    builder.insertMolecule(&mol, {0.0,0.0,0.0}, rot );
    builder.insertMolecule(&mol, {5.0,0.0,0.0}, rot );
    builder.insertMolecule(&mol, {0.0,5.0,0.0}, rot );
    builder.insertMolecule(&mol, {5.0,5.0,0.0}, rot );
    builder.toMMFF(&world, &params);

    world.ang_b2a();

    //exit(0);

    world.printBondParams();
    //exit(0);

    opt.bindArrays( 3*world.natoms, (double*)world.apos, new double[3*world.natoms], (double*)world.aforce );
    opt.setInvMass( 1.0 );
    opt.cleanVel( );

    printf( "DEBUG 2 \n" );

    for(int i=0; i<world.nbonds; i++){
        world.bond_k[i] = 2.0;
    }

    printf( "DEBUG 3 \n" );

    for(int i=0; i<world.nang; i++){
        world.ang_0[i] = {1.0,0.0};
        world.ang_k[i] = 0.5;
        //Vec2i ib = world.ang2bond[i];
        //world.ang2atom [i] = (Vec3i){ world.bond2atom[ib.x].y, world.bond2atom[ib.y].y, world.bond2atom[ib.y].x };
    }

    printf( "DEBUG 4 \n" );

    ogl_sph = glGenLists(1);
    glNewList( ogl_sph, GL_COMPILE );
        //glEnable( GL_LIGHTING );
        //glColor3f( 0.8f, 0.8f, 0.8f );
        //Draw3D::drawSphere_oct(3, 0.5, {0.0,0.0,0.0} );
        Draw3D::drawSphere_oct( 3, 1.0, {0.0,0.0,0.0} );
    glEndList();

    //printf( "bond 8 %g \n", world.bond_0[8] );
    //printf( "bond 9 %g \n", world.bond_0[9] );
    //Vec2i iat = bond2atom[8];
    //Vec2i iat = bond2atom[9];
    //exit(0);

    world.grid.n    = (Vec3i){100,100,100};
    world.grid.pos0 = (Vec3d){-5.0,-5.0,-5.0};
    world.grid.setCell( (Mat3d){ 10.0,0.0f,0.0f,  0.0,10.0f,0.0f,  0.0,0.0f,10.0f } );

    Vec3d * FF     = new Vec3d[world.grid.getNtot()];
    world.FFPauli  = new Vec3d[world.grid.getNtot()];
    world.FFLondon = new Vec3d[world.grid.getNtot()];

    world.evalGridFFexp( world.natoms, world.apos, world.aLJq, -2.0*2,  1, world.FFPauli );
    world.evalGridFFexp( world.natoms, world.apos, world.aLJq, -2.0  , -2, world.FFLondon );
    world.evalCombindGridFF( -2.0, 0.0, 2.0, 1.0, FF );
    world.grid.saveXSF( "testX.xsf", FF, 0 );
    world.grid.saveXSF( "testY.xsf", FF, 1 );
    world.grid.saveXSF( "testZ.xsf", FF, 2 );

}

void AppMolecularEditor2::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    /*
    srand(154);
    Vec3d ax  = (Vec3d){randf(-1.0,1.0),randf(-1.0,1.0),randf(-1.0,1.0)};   ax.normalize();
    Vec3d dp1 = (Vec3d){randf(-1.0,1.0),randf(-1.0,1.0),randf(-1.0,1.0)};   dp1.add_mul( ax, -ax.dot(dp1) ); dp1.normalize();
    Vec3d dp2 = (Vec3d){randf(-1.0,1.0),randf(-1.0,1.0),randf(-1.0,1.0)};   dp2.add_mul( ax, -ax.dot(dp2) ); dp2.normalize();
    printf( "%g %g %g %g \n", dp1.dot(dp1), dp2.dot(dp2), dp1.dot(ax), dp2.dot(ax) );
    Vec3d p1  = ax*randf( 0.0,1.0) + dp1*randf( -1.5,-1.0);
    Vec3d p2  = ax*randf(-1.0,0.0) + dp2*randf( -1.5,0.0);
    double t1,t2;
    double dist = rayLine( p1, dp1, p2, dp2, t1, t2 );
	// test closest point skew-lines
	//glColor3f(0.6f,0.6f,0.6f); Draw3D::drawLine( ax*1, ax*-1 );
	glColor3f(0.8f,0.0f,0.0f); Draw3D::drawLine( p1, p1+dp1 ); Draw3D::drawPointCross( p1, 0.1 );
	glColor3f(0.0f,0.0f,0.8f); Draw3D::drawLine( p2, p2+dp2 ); Draw3D::drawPointCross( p2, 0.1 );
	glColor3f(0.0f,0.0f,0.0f); Draw3D::drawLine( p1+(dp1*t1), p2+(dp2*t2) ); // Draw3D::drawLine( p1+(dp1*t1), p1 );   Draw3D::drawLine( p2, p2+(dp2*t2) );
	printf( "t1 %g t2 %g  dist %g %g \n", t1, t2, dist, (p1+(dp1*t1)-p2-(dp2*t2)).norm() );
	//Vec3d rxl; rxl.set_cross(dp1, dp2);
	//glColor3f(0.8f,0.0f,0.0f); Draw3D::drawVecInPos( rxl, p1+(dp1*t1) );
	//glColor3f(0.0f,0.0f,0.8f); Draw3D::drawVecInPos( rxl, p2+(dp2*t2) );
	return;
	*/

	//ibpicked = world.pickBond( ray0, camMat.c , 0.5 );

    ray0 = camMat.a*mouse_begin_x + camMat.b*mouse_begin_y;
    Draw3D::drawPointCross( ray0, 0.1 );
    //Draw3D::drawVecInPos( camMat.c, ray0 );
    if(ipicked>=0) Draw3D::drawLine( world.apos[ipicked], ray0);

	double F2;
	for(int itr=0; itr<perFrame; itr++){

        for(int i=0; i<world.natoms; i++){ world.aforce[i].set(0.0d); }

        //printf( "DEBUG x.1 \n" );
        world.eval_bonds(true);
        //world.eval_angles();
        //printf( "DEBUG x.2 \n" );
        world.eval_angcos();
        //printf( "DEBUG x.3 \n" );
        world.eval_LJq_On2();

        //exit(0);
        if(ipicked>=0){
            Vec3d f = getForceSpringRay( world.apos[ipicked], camMat.c, ray0, -1.0 );
            //printf( "f (%g,%g,%g)\n", f.x, f.y, f.z );
            world.aforce[ipicked].add( f );
        };


        for(int i=0; i<world.natoms; i++){
            world.aforce[i].add( getForceHamakerPlane( world.apos[i], {0.0,0.0,1.0}, -3.0, 0.3, 2.0 ) );
            //printf( "%g %g %g\n",  world.aforce[i].x, world.aforce[i].y, world.aforce[i].z );
        }

        //exit(0);

        //for(int i=0; i<world.natoms; i++){ world.aforce[i].add({0.0,-0.01,0.0}); }
        //int ipivot = 0;
        //world.aforce[ipivot].set(0.0);
        //opt.move_LeapFrog(0.01);
        //opt.move_MDquench();
        F2 = opt.move_FIRE();
        //exit(0);

    }

    glColor3f(0.6f,0.6f,0.6f); plotSurfPlane( (Vec3d){0.0,0.0,1.0}, -3.0, {3.0,3.0}, {20,20} );
    //Draw3D::drawVecInPos( (Vec3d){0.0,0.0,1.0},  (Vec3d){0.0,0.0,0.0} );

    //printf( "==== frameCount %i  |F| %g \n", frameCount, sqrt(F2) );


    for(int i=0; i<world.nbonds; i++){
        Vec2i ib = world.bond2atom[i];
        glColor3f(0.0f,0.0f,0.0f);
        if(i==ibpicked) glColor3f(1.0f,0.0f,0.0f); ;
        Draw3D::drawLine(world.apos[ib.x],world.apos[ib.y]);
        sprintf(str,"%i\0",i);
        Draw3D::drawText(str, (world.apos[ib.x]+world.apos[ib.y])*0.5, fontTex, 0.02, 0,0);
    }

    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);
    for(int i=0; i<world.natoms; i++){
        //glColor3f(0.0f,0.0f,0.0f); Draw3D::drawPointCross(world.apos[i],0.2);
        glColor3f(1.0f,0.0f,0.0f); Draw3D::drawVecInPos(world.aforce[i]*30.0,world.apos[i]);

        //glCallList( ogl_sph );
        glEnable(GL_LIGHTING);
        Mat3d mat;
        mat.setOne();
        mat.mul( atomSize*params.atypes[world.atypes[i]].RvdW );
        //glColor3f(0.8f,0.8f,0.8f);
        Draw::setRGB( params.atypes[world.atypes[i]].color );
        Draw3D::drawShape(world.apos[i],mat,ogl_sph);
        glDisable(GL_LIGHTING);
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


void AppMolecularEditor2::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_p:  first_person = !first_person; break;
                //case SDLK_o:  perspective  = !perspective; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;

                case SDLK_v: for(int i=0; i<world.natoms; i++){ ((Vec3d*)opt.vel)[i].add(randf(-drndv,drndv),randf(-drndv,drndv),randf(-drndv,drndv)); } break;
                case SDLK_p: for(int i=0; i<world.natoms; i++){ world.apos[i].add(randf(-drndp,drndp),randf(-drndp,drndp),randf(-drndp,drndp)); } break;

                case SDLK_LEFTBRACKET:  if(ibpicked>=0) world.bond_0[ibpicked] += 0.1; break;
                case SDLK_RIGHTBRACKET: if(ibpicked>=0) world.bond_0[ibpicked] -= 0.1; break;

                case SDLK_a: world.apos[1].rotate(  0.1, {0.0,0.0,1.0} ); break;
                case SDLK_d: world.apos[1].rotate( -0.1, {0.0,0.0,1.0} ); break;

                case SDLK_w: world.apos[1].mul( 1.1 ); break;
                //case SDLK_s: world.apos[1].mul( 0.9 ); break;
                case SDLK_s: printf("saving ... "); save2xyz( "out.xyz", &world, &params ); printf("... DONE "); break;
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
                    //ibpicked = -1;
                    break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}

void AppMolecularEditor2::drawHUD(){
    glDisable ( GL_LIGHTING );

}

// ===================== MAIN

AppMolecularEditor2 * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new AppMolecularEditor2( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















