
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
#include "VecN.h"


#include "AppSDL2OGL_3D.h"
#include "testUtils.h"
#include "SDL_utils.h"
#include "Plot2D.h"

#include "Molecule.h"
#include "BondAdaptedMesh.h"

void drawTetrahedron(const Tetrahedron& t){
    glBegin(GL_LINES);
    Draw3D::vertex(t.a); Draw3D::vertex(t.b);
    Draw3D::vertex(t.a); Draw3D::vertex(t.c);
    Draw3D::vertex(t.a); Draw3D::vertex(t.d);
    Draw3D::vertex(t.b); Draw3D::vertex(t.c);
    Draw3D::vertex(t.c); Draw3D::vertex(t.d);
    Draw3D::vertex(t.d); Draw3D::vertex(t.b);
    glEnd();
}

// ======= THE CLASS

class TestAppRARFF: public AppSDL2OGL_3D { public:

    //RigidAtom     atom1;
    //RigidAtomType type1,type2;

    Molecule mol;
    AtomGrid agrid;

    bool bRun = false;

    int ipicked  = -1, ibpicked = -1;

    //Plot2D plot1;

    //double Emin,Emax;
    //int     npoints;
    //Vec3d*  points  =0;
    //double* Energies=0;
    //Vec3d * Forces  =0;

    int      fontTex;

    virtual void draw   ();
    virtual void drawHUD();
    //virtual void mouseHandling( );
    virtual void eventHandling   ( const SDL_Event& event  );
    //virtual void keyStateHandling( const Uint8 *keys );

    TestAppRARFF( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppRARFF::TestAppRARFF( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );

    mol.loadXYZ_bas( "data/CH4.bas" );

    double Rmax = 1.8;

    mol.findBonds_brute( Rmax );

    printf( "mol.nbonds %i \n", mol.nbonds );

    std::vector<int> neighs[mol.natoms];

    for(int i=0; i<mol.nbonds; i++){
        Vec2i& b = mol.bond2atom[i];
        neighs[b.i].push_back(b.j);
        neighs[b.j].push_back(b.i);
    }

    printf( "neighs[0].size() %i \n", neighs[0].size() );

    agrid.fromNeighbors( 0, neighs[0].size(), &neighs[0][0], mol.pos, Rmax );

    printf( "agrid.mesh npoints, nedges, npolys %i,%i,%i \n", agrid.mesh.points.size(), agrid.mesh.edges.size(), agrid.mesh.polygons.size() );

    //exit(0);

    //ff.loadFromFile_bas( "data/CH4_eFF.bas" );

    //double sz = 0.51;
    // break symmetry
    //for(int i=0; i<ff.na; i++){ ff.apos[i].add( randf(-sz,sz),randf(-sz,sz),randf(-sz,sz) );  }


    for(int i=0; i<mol.natoms; i++){
        printf( "A_pos[%i] (%g,%g,%g)\n", i, mol.pos[i].x, mol.pos[i].y, mol.pos[i].z );
    }



    //exit(0);

    //ff.move( 0.01, 0.9 );

}

void TestAppRARFF::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);


    double fsc = 1.0;
    glColor3f(0.0,0.0,0.0);
    for(int i=0; i<mol.natoms; i++){
        //printf( "apos[%i] (%g,%g,%g)\n", i, ff.apos[i].x, ff.apos[i].y, ff.apos[i].z );
        Draw3D::drawPointCross( mol.pos  [i] , 0.2 );
        //Draw3D::drawVecInPos(   ff.aforce[i]*fsc, ff.apos[i] );
        //printf( " %i %f %f %f %f  \n", i, ff.aQ[i], ff.apos[i].x,ff.apos[i].y,ff.apos[i].z );
        //printf( " %i %f %f %f %f  \n", i, ff.aQ[i], ff.aforce[i].x, ff.aforce[i].y, ff.aforce[i].z );
    }

    for(int i=0; i<mol.nbonds; i++){
        Draw3D::drawLine( mol.pos[mol.bond2atom[i].i], mol.pos[mol.bond2atom[i].j] );
    }


    glColor3f(0.0,0.0,1.0);
    for( MeshEdge& edge: agrid.mesh.edges ){
        Vec3d& a =  agrid.mesh.points[ edge.verts.a ];
        Vec3d& b =  agrid.mesh.points[ edge.verts.b ];
        //Draw3D::drawPointCross( l12.p0, 0.1 ); Draw3D::drawVecInPos( l12.hdir, l12.p0 );
        //printf( "edge %i(%g,%g,%g) %i(%g,%g,%g) \n", edge.verts.a, a.x,a.y,a.z,   edge.verts.b, b.x,b.y,b.z );
        Draw3D::drawLine( a, b );
    }

    glColor3f(0.0,1.0,0.0);
    for(int i=0; i<agrid.tetrahedrons.size(); i++){
        drawTetrahedron(agrid.tetrahedrons[i]);
    }
    //exit(0);

    for( Vec3d& p: agrid.mesh.points ){
        glEnable(GL_DEPTH_TEST);
        Draw3D::drawPointCross( p, 0.1 );
    }

    /*
    int i=0;
    for( Polygon* pl: agrid.mesh.polygons ){
        Draw::color_of_hash(i*4456464+54844); i++;
        glBegin(GL_TRIANGLE_FAN);
        for( int i : pl->ipoints ){
            Vec3d& v = agrid.mesh.points[i];
            glVertex3f( v.x,v.y,v.z  );
        }
        glEnd();
    }
    */

};


void TestAppRARFF::drawHUD(){


	glTranslatef( 100.0,100.0,0.0 );
	glScalef    ( 10.0,10.00,1.0  );
	//plot1.view();

}

/*
void TestAppRARFF::mouseHandling( ){
    int mx,my; Uint32 buttons = SDL_GetRelativeMouseState( &mx, &my);
    if ( buttons & SDL_BUTTON(SDL_BUTTON_LEFT)) {

    }
    AppSDL2OGL_3D::mouseHandling( );
};
*/

void TestAppRARFF::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
                case SDLK_SPACE: bRun = !bRun;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //ipicked = pickParticle( ff.natoms, ff.apos, ray0, (Vec3d)cam.rot.c , 0.5 );
                break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //ibpicked = pickParticle( ff.natoms, ff.apos, ray0, (Vec3d)cam.rot.c , 0.5 );
                    //printf( "dist %i %i = ", ipicked, ibpicked );
                    break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}

// ===================== MAIN

TestAppRARFF* thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppRARFF( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















