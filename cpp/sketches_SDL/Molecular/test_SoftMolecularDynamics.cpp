
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

// Non-Blocking terminal input : see https://stackoverflow.com/questions/6055702/using-fgets-as-non-blocking-function-c
// Probably works only on Linux
#include <fcntl.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw3D.h"

#include "SDL_utils.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"

#include "raytrace.h"
#include "Molecule.h"
#include "MMFF.h"
#include "MMFFBuilder.h"
#include "DynamicOpt.h"

#include "Draw3D_Molecular.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

#include "repl.h"
#include "commandTree.h"

/*

ToDo (quick):
    * assign angles depending on number of bonds and atom type
    * Pin/Fix atom in place
    * use MM::Builder to automatically assign hydrogen, free-electron pairs and Pi-orbitals
    * plot real atom sizes

More complex ToDo:
    * Rotate selection around bond
    * Move atoms using Gizmo
    * Fix selection as rigid body
    * Plot energy profiles along atom movement (other atoms fixed)

*/

// ==========================
// TestAppSoftMolDyn
// ==========================

class TestAppSoftMolDyn : public AppSDL2OGL_3D {
	public:
	Molecule    mol;
	MMFFparams  params;
    MMFF        world;
    MM::Builder builder;

    DynamicOpt  opt;

    int     fontTex=0,fontTexPix=0;
    int     ogl_sph=0;

    char str[256];

    REPL::Interpreter repl;
    CommandTree commands;

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

void moveAtoms( int n, Vec3d* ps, Vec3d shift ){
    for(int i=0; i<n; i++)ps[i].add(shift);
}

TestAppSoftMolDyn::TestAppSoftMolDyn( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    REPL::init();
    repl.functions["move\n"] = [this]{ moveAtoms(world.natoms, world.apos, (Vec3d){1.0,0.0,0.0}); };


    fontTexPix = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );
    fontTex    = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );

    params.loadAtomTypes( "common_resources/AtomTypes.dat" );
    //mol.atomTypeNames = &params.atomTypeNames;
    //mol.atomTypeDict  = &params.atomTypeDict;
    params.loadBondTypes( "common_resources/BondTypes.dat");

    mol.bindParams(&params);

    mol.loadMol("common_resources/propylacid.mol");
    //mol.loadMol_old("common_resources/propylacid.mol");
    //mol.bondsOfAtoms();
    mol.autoAngles(true);   ;
    params.assignREs( mol.natoms, mol.atomType, mol.REQs );

    Vec3d cog = mol.getCOG_av();
    mol.addToPos( cog*-1.0d );

    builder.insertMolecule(&mol, {0.0,0.0,0.0}, Mat3dIdentity, false );
    builder.insertMolecule(&mol, {5.0,0.0,0.0}, Mat3dIdentity, false );
    builder.insertMolecule(&mol, {0.0,5.0,0.0}, Mat3dIdentity, false );
    builder.insertMolecule(&mol, {5.0,5.0,0.0}, Mat3dIdentity, false );
    builder.toMMFF(&world, &params );

    world.ang_b2a();

    //params.printAtomTypeDict();
    //mol.printAtom2Bond();
    //mol.printAngleInfo()
    //world.printBondParams();
    //world.printAtomInfo();

    opt.bindArrays( 3*world.natoms, (double*)world.apos, new double[3*world.natoms], (double*)world.aforce, NULL );
    opt.setInvMass( 1.0 );
    opt.cleanVel( );

    for(int i=0; i<world.nbonds; i++){
        world.bond_k[i] = 2.0;
    }

    for(int i=0; i<world.nang; i++){
        world.ang_0[i] = {1.0,0.0};
        world.ang_k[i] = 0.5;
        //Vec2i ib = world.ang2bond[i];
        //world.ang2atom [i] = (Vec3i){ world.bond2atom[ib.x].y, world.bond2atom[ib.y].y, world.bond2atom[ib.y].x };
    }

    //Draw3D::makeSphereOgl( ogl_sph, 2, 0.25 );
    Draw3D::makeSphereOgl( ogl_sph, 2, 1.0 );


    commands.root.func=[]{};
    //commands.root.leafs.insert( { "h", []{printf("Hey!\n")} }  );
    //commands.root.leafs.insert( { "o", []{printf("OK!\n")} }  );
    //commands.root.leafs.insert( { "a", []{printf("Action 1 \n")} }  );

    commands.root.addLeaf( "h", "prints Hey", []{printf("Hey!\n");} );
    commands.root.addLeaf( "o", "prints OK",  []{printf("OK!\n");}  );
    CommandNode* cnd = commands.root.addLeaf( "1", "Toolbox 1",  []{printf("Opening Toolbox 1 \n");} );
    cnd->addLeaf( "1", "Tool 1.1",  []{printf("run_tool 1.1 \n");} );
    cnd->addLeaf( "2", "Tool 1.2",  []{printf("run_tool 1.2 \n");} );
    commands.root.addLeaf("2", "Toolbox 2",  []{printf("Opening Toolbox 1 \n");} );
    cnd->addLeaf( "1", "Tool 2.1",  []{printf("run_tool 2.1 \n");} );
    cnd->addLeaf( "2", "Tool 2.2",  []{printf("run_tool 2.2 \n");} );

}


void TestAppSoftMolDyn::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    /*
	// Non-Blocking terminal input : see https://stackoverflow.com/questions/6055702/using-fgets-as-non-blocking-function-c
    int fd = fileno(stdin);
    int flags  = fcntl(fd, F_GETFL, 0);
    flags |= O_NONBLOCK;
    fcntl(fd, F_SETFL, flags);
	char inputLine[1024];
	char* s = fgets(inputLine,1024,stdin);
	if(s!=0){
        printf( "got:`%s`", s );
	}
    */
    repl.eval();

    ray0 = (Vec3d)mouseRay0(); Draw3D::drawPointCross( ray0, 0.1 ); //Draw3D::drawVecInPos( camMat.c, ray0 );
    if(ipicked>=0) Draw3D::drawLine( world.apos[ipicked], ray0);

	double F2;
	for(int itr=0; itr<perFrame; itr++){

        for(int i=0; i<world.natoms; i++){ world.aforce[i].set(0.0d); }

        //printf( "DEBUG x.1 \n" );
        world.eval_bonds(true);     // with    eval_LJq_On2
        //world.eval_bonds(false);  // without eval_LJq_On2
        //printf( "DEBUG x.2 \n" );
        //world.eval_angles();  // currently not working
        world.eval_angcos();
        //printf( "DEBUG x.3 \n" );
        world.eval_LJq_On2();

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

        //exit(0);

        //for(int i=0; i<world.natoms; i++){ world.aforce[i].add({0.0,-0.01,0.0}); }
        //int ipivot = 0;
        //world.aforce[ipivot].set(0.0);
        //opt.move_LeapFrog(0.01);
        //opt.move_MDquench();
        F2 = opt.move_FIRE();
        //exit(0);
        //printf( "==== frameCount %i  |F| %g \n", frameCount, sqrt(F2) );

    }

    glColor3f(0.6f,0.6f,0.6f); Draw3D::plotSurfPlane( (Vec3d){0.0,0.0,1.0}, -3.0, {3.0,3.0}, {20,20} );
    glColor3f(0.0f,0.0f,0.0f);
    Draw3D::drawLines ( world.nbonds, (int*)world.bond2atom, world.apos );
    Draw3D::bondLabels( world.nbonds,       world.bond2atom, world.apos, fontTex, 0.02 );
    glColor3f(1.0f,0.0f,0.0f);
    Draw3D::vecsInPoss( world.natoms, world.aforce, world.apos, 300.0              );
    Draw3D::atomsREQ  ( world.natoms, world.apos,   world.aREQ, ogl_sph, 1.0, 0.25 );

    //printf("==========\n");
    //for(int i=0; i<world.natoms; i++){
    //    printf("iatom %i (%g,%g,%g) (%g,%g,%g) \n", i, world.apos[i].x,world.apos[i].y,world.apos[i].z, world.aforce[i].x,world.aforce[i].y,world.aforce[i].z  );
    //}
    //if(frameCount>=10){STOP = true;}
};


void TestAppSoftMolDyn::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_TEXTINPUT:
            /* Add new text onto the end of our text */
            commands.eval( event.text.text );
            break;
        case SDL_KEYDOWN :
            //commands.eval( event.key. );
            //printf( "SDLK_w ", event.key. )
            switch( event.key.keysym.sym ){
                //case SDLK_p:  first_person = !first_person; break;
                //case SDLK_o:  perspective  = !perspective; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
                case SDLK_BACKSPACE: commands.goBack(); break;

                case SDLK_v: for(int i=0; i<world.natoms; i++){ ((Vec3d*)opt.vel)[i].add(randf(-drndv,drndv),randf(-drndv,drndv),randf(-drndv,drndv)); } break;
                case SDLK_p: for(int i=0; i<world.natoms; i++){ world.apos[i].add(randf(-drndp,drndp),randf(-drndp,drndp),randf(-drndp,drndp)); } break;

                case SDLK_LEFTBRACKET:  if(ibpicked>=0) world.bond_0[ibpicked] += 0.1; break;
                case SDLK_RIGHTBRACKET: if(ibpicked>=0) world.bond_0[ibpicked] -= 0.1; break;

                case SDLK_a: world.apos[1].rotate(  0.1, {0.0,0.0,1.0} ); break;
                case SDLK_d: world.apos[1].rotate( -0.1, {0.0,0.0,1.0} ); break;

                case SDLK_w: world.apos[1].mul( 1.1 ); break;
                //case SDLK_s: world.apos[1].mul( 0.9 ); break;

            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    ipicked = pickParticle( ray0, (Vec3d)cam.rot.c , 0.5, world.natoms, world.apos );
                    break;
                case SDL_BUTTON_RIGHT:
                    ibpicked = world.pickBond( ray0, (Vec3d)cam.rot.c , 0.5 );
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

void TestAppSoftMolDyn::drawHUD(){
    glDisable ( GL_LIGHTING );
    char str[1024];
    commands.curInfo( str );
    //sprintf( str, );
    //Draw3D::drawText( str, (Vec3f){10.,10.,0.0}, fontTexPix, 20, 0 );
    glTranslatef( 10 ,HEIGHT-20 ,0 );
    Draw::drawText( str, fontTexPix, 7, {100,50} );

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
















