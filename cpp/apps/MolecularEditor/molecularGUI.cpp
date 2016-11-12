
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include <SDL2/SDL_image.h>
//#include <SDL2/SDL_ttf.h>
//#include "Texture.h"

#include "Draw2D.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
//#include "DynamicOpt.h"

#include "AtomTypes.h"
#include "MoleculeType.h"
#include "MolecularWorld.h"

#include "Draw.h"
#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"
#include "SDL_utils.h"
#include "testUtils.h"

// font rendering:
//  http://www.willusher.io/sdl2%20tutorials/2013/12/18/lesson-6-true-type-fonts-with-sdl_ttf
//  http://stackoverflow.com/questions/28880562/rendering-text-with-sdl2-and-opengl


/*
int drawAtom( MoleculeType * mol, int i, int nsphere, float atomscale, uint32_t color ){
    Draw::setRGB( color );
    int nvert = Draw3D::drawSphere_oct( nsphere, atomscale*mol->typeList->vdwRs[mol->atypes[i]], mol->xyzs[i] );
    return nvert;
}

int drawBond( MoleculeType * mol, int i, int j, int nstick, float bondwidth  ){
    Vec3f ai,aj;
    convert( mol->xyzs[i], ai );
    convert( mol->xyzs[j], aj );
    int nvert = Draw3D::drawCylinderStrip( nstick, bondwidth, bondwidth, ai, aj );
    return nvert;
}
*/

int renderMoleculeCPK ( MoleculeType * mol, int nsphere, int nstick, float atomscale, float bondwidth ){
    if( mol->viewlist > 0 ) {	glDeleteLists( mol->viewlist, 1 );	}
    int nvert = 0;
    mol->viewlist = glGenLists(1);
    glNewList( mol->viewlist , GL_COMPILE );
        glShadeModel ( GL_SMOOTH );
        for( int i=0; i<mol->natoms; i++     ){
            //printf("render atom %i \n", i);
            //nvert+= drawAtom( mol, i, nsphere, atomscale, mol->typeList->colors[ mol->atypes[i] ] );
            //printf("render atom %i %i \n", mol->atypes[i] );
            uint32_t color = mol->typeList->colors[ mol->atypes[i] ];
            //printf("render atom %i type %i color %i \n", i, mol->atypes[i], color );
            Draw::setRGB( color );
            nvert += Draw3D::drawSphere_oct( nsphere, atomscale*mol->typeList->vdwRs[mol->atypes[i]], mol->xyzs[i] );
        }
        if( mol->bonds != NULL ){
            glColor3f( 0.2f, 0.2f, 0.2f );
            for( int ib=0; ib<mol->nbonds; ib+=2 ){
                //nvert+= drawBond( mol, , mol->bonds[ib+1], nstick, bondwidth );
                Vec3f ai,aj;
                convert( mol->xyzs[ mol->bonds[ib  ] ], ai );
                convert( mol->xyzs[ mol->bonds[ib+1] ], aj );
                nvert += Draw3D::drawCylinderStrip( nstick, bondwidth, bondwidth, ai, aj );
            }
        }
    glEndList();
    //printf( " nvert %i \n", nvert );
    return mol->viewlist;
}

// ============================
//   MolecularEditorApp
// ============================

class MolecularEditorApp : public AppSDL2OGL_3D {
	public:
    MolecularWorld world;

    int perFrame       = 40;
    bool converged     = false;
    double fmaxConverg = 0.00001;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );

	MolecularEditorApp( int& id, int WIDTH_, int HEIGHT_ );
};

MolecularEditorApp::MolecularEditorApp( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    //MolecularWorld( char const* filename, MoleculeType * molTypeList );
    world.fromDir( "inputs/", "atomTypes.ini", "molTypes.ini", "instances.ini" );
    world.makeFF ( );
    world.optimizer->initOpt( 0.05, 0.15 );

    for(int i=0; i<world.nMolTypes; i++){
        printf(" rendering mol %i \n", i );
        world.molTypes[i].toCOG_average();
        world.molTypes[i].findBonds( 0.6 );
        renderMoleculeCPK( &world.molTypes[i], 4, 8, 0.5, 0.2 );
    }

}

void MolecularEditorApp::draw(){

    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	//glDisable( GL_DEPTH_TEST );


	if( !converged ){
        long tick1 = getCPUticks();
        for(int iter=0; iter<perFrame; iter++){
            world.rigidOptStep( );
            double fmax = world.optimizer->getFmaxAbs( );
            printf(" opt step %i fmax %g \n", world.optimizer->stepsDone, fmax );
            if( fmax < fmaxConverg ){
                converged = true;
                printf(" converged after %i step \n", world.optimizer->stepsDone );
                break;
            }
        }
        double ticks = (getCPUticks() - tick1)/((double)perFrame);
        printf("======= %f Mticks/iter  %f ticks/interaction \n", ticks*1.0e-6, ticks/world.nInteractions );
    }
    //exit(0);


	glMatrixMode(GL_MODELVIEW);
	//glMatrixMode(GL_PROJECTION);
    for (int i=0; i<world.nmols; i++){
        if( world.instances[i]->viewlist > 0 ){
            Mat3d rotmat;
            float glMat[16];
            //rot[i].toMatrix_unitary2( rotmat );
            //rot[i].toMatrix_unitary( rotmat );
            //printf( "%i   (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f,%3.3f)\n", i,  world.pos[i].x,world.pos[i].y,world.pos[i].z,   world.rot[i].x,world.rot[i].y,world.rot[i].z,world.rot[i].w  );
            world.rot[i].toMatrix( rotmat );
            Draw3D::drawPointCross(world.pos[i],2.0);
            Draw3D::toGLMat( world.pos[i], rotmat, glMat ); // somehow not working
            //Draw3D::toGLMat( {0.0,0.0,0.0}, rotmat, glMat );
            glPushMatrix();
            glMultMatrixf( glMat );
            //glTranslatef( world.pos[i].x,world.pos[i].y,world.pos[i].z );
            //glMultTransposeMatrixf( glMat );
            //glLoadMatrixf( glMat );
            glCallList   ( world.instances[i]->viewlist );
            glPopMatrix();
        }
    };

    //exit(0);

};

void MolecularEditorApp::drawHUD(){}

void MolecularEditorApp::eventHandling ( const SDL_Event& event  ){
    //printf( "NBodyWorldApp::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_0:  formation_view_mode = 0;            printf( "view : default\n" ); break;
                //case SDLK_1:  formation_view_mode = VIEW_INJURY;  printf( "view : injury\n"  ); break;
                //case SDLK_2:  formation_view_mode = VIEW_STAMINA; printf( "view : stamina\n" ); break;
                //case SDLK_3:  formation_view_mode = VIEW_CHARGE;  printf( "view : charge\n"  ); break;
                //case SDLK_4:  formation_view_mode = VIEW_MORAL;   printf( "view : moral\n"   ); break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                break;
                case SDL_BUTTON_RIGHT:
                break;
            }
            break;
            /*
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    world.picked = NULL;
                    break;
            }
            break;
            */
    };
    AppSDL2OGL_3D::eventHandling( event );
    camStep = zoom*0.05;
}

// ===================== MAIN

MolecularEditorApp * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	thisApp = new MolecularEditorApp( junk , 800, 600 );
	thisApp->zoom = 30;
	thisApp->loop( 1000000 );
	return 0;
}
















