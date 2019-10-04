
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

/*
#include "Multipoles.h"
#include "PotentialFlow.h"
#include "grids3D.h"
#include "MultipoleGrid.h"
*/

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"
#include "SDL_utils.h"
#include "Plot2D.h"

//#include "MMFF.h"

#include "RARFF.h"

#define R2SAFE  1.0e-8f

// ======= THE CLASS

/*
void drawVectorArray(int n, Vec3d* ps, Vec3d* vs, double sc ){
    glBegin(GL_LINES);
    for(int i=0; i<n; i++){
        Vec3d p=ps[i];        glVertex3f(p.x,p.y,p.z);
        p.add_mul( vs[i], sc); glVertex3f(p.x,p.y,p.z);
    }
    glEnd();
}

void drawScalarArray(int n, Vec3d* ps, double* vs, double vmin, double vmax ){
    glBegin(GL_POINTS);
    double sc = 1/(vmax-vmin);
    for(int i=0; i<n; i++){
        Vec3d p=ps[i];
        double c = (vs[i]-vmin)*sc;
        glColor3f(c,c,c);
        glVertex3f(p.x,p.y,p.z);
        //printf( "i %i p(%g,%g,%g) v: %g c: %g\n", i, p.x,p.y,p.z, vs[i], c );
    }
    glEnd();
}
*/

//void drawRigidAtom( const Vec3d& pos, Vec3d* bhs ){
void drawRigidAtom( RigidAtom& atom ){
    Vec3d bhs[N_BOND_MAX];

    //rotateVectors<double>(N_BOND_MAX, atom.qrot, atom.type->bh0s, bhs );
    atom.qrot.rotateVectors(N_BOND_MAX, atom.type->bh0s, bhs, false );

    Draw3D::drawPointCross( atom.pos, 0.1 );
    //for(int i=0; i<N_BOND_MAX; i++){
    for(int i=0; i<atom.type->nbond; i++){
        //printf( "%i (%g,%g,%g)\n", i, type1.bh0s[i].x, type1.bh0s[i].y, type1.bh0s[i].z );
        //printf( "%i (%g,%g,%g)\n", i, bhs[i].x, bhs[i].y, bhs[i].z );
        Draw3D::drawVecInPos( bhs[i], atom.pos );
    }
}


template<typename Func> void numDeriv( Vec3d p, double d, Vec3d& f, Func func){
    //double e0 = Efunc(p);
    double d_=d*0.5;
    p.x+=d_; f.x = func(p); p.x-=d; f.x-=func(p); p.x+=d_;
    p.y+=d_; f.y = func(p); p.y-=d; f.y-=func(p); p.y+=d_;
    p.z+=d_; f.z = func(p); p.z-=d; f.z-=func(p); p.z+=d_;
    f.mul(1/d);
}

class TestAppRARFF: public AppSDL2OGL_3D { public:

    RigidAtom     atom1;
    RigidAtomType type1,type2;

    bool bRun = true;

    RARFF ff;

    Plot2D plot1;

    double Emin,Emax;
    int     npoints;
    Vec3d*  points  =0;
    double* Energies=0;
    Vec3d * Forces  =0;

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

    // for exp
    type1.nbond = 3;  // number bonds
    type1.rbond0 = 0.5;
    type1.acore =  4.0;
    type1.bcore = -0.7;
    type1.abond =  3.0;
    type1.bbond = -1.1;
    type1.bh0s = (Vec3d*)sp2_hs;
    type1.print();
    //exit(0);

    type2.nbond = 4;  // number bonds
    type2.rbond0 = 0.5;
    type2.acore =  4.0;
    type2.bcore = -0.7;
    type2.abond =  3.0;
    type2.bbond = -1.1;
    type2.bh0s = (Vec3d*)sp3_hs;
    type2.print();
    //exit(0);


/*
    // for overlap
    type1.nbond = 3;  // number bonds
    type1.rbond0 = 0.5;
    type1.acore =  4.0;
    type1.bcore = -0.7;
    type1.abond =  2.0;
    type1.bbond = -1.5;
*/


    //exit(0);

    /*
    ff.realloc(4);
    for(int i=0; i<ff.natom; i++){ ff.atoms[i].type=&type1; };
    ff.atoms[0].setPose( (Vec3d){0.0,0.0,0.0}, Quat4dIdentity );
    ff.atoms[1].setPose( (Vec3d){1.0,1.5,0.0}, Quat4dIdentity );
    ff.atoms[2].setPose( (Vec3d){1.0,0.0,0.0}, Quat4dIdentity );
    ff.atoms[3].setPose( (Vec3d){0.0,1.0,0.0}, Quat4dIdentity );
    */


/*
    int nang    = 6;
    ff.realloc(nang+3);
    for(int i=0; i<ff.natom; i++){ ff.atoms[i].type=&type2; };
    //double dang = 2*M_PI/(nang +0.5);
    double dang = 2*M_PI/(nang +0.0);
    for(int i=0; i<nang; i++){

        double ang = i*dang;
        double sa=sin(ang);
        double ca=cos(ang);

        double sa_=sin(ang*0.5);
        double ca_=cos(ang*0.5);
        //ff.atoms[i].type=&type1;
        ff.atoms[i].pos  = (Vec3d ){ca,sa,0};
        ff.atoms[i].qrot = (Quat4d){0.0,0.0,-sa_,ca_};
    }

    ff.atoms[6].setPose( (Vec3d){2.0,0.0,0.0}, Quat4dIdentity );  ff.atoms[6].type=&type2;     ff.atoms[7].qrot.dRot_exact(1.0, (Vec3d){0.0,0.0,1.0} );
    //ff.atoms[7].setPose( (Vec3d){0.0,2.0,0.0}, Quat4dIdentity );  ff.atoms[7].type=&type2;
    ff.atoms[7].setPose( (Vec3d){-1.0,2.0,0.0}, Quat4dIdentity );  ff.atoms[7].type=&type2;
    ff.atoms[8].setPose( (Vec3d){0.0,0.0,-2.0}, Quat4dIdentity );  ff.atoms[7].type=&type2;
*/


    srand(0);
    //srand(2);

    int nat = 12;
    ff.realloc(nat);
    for(int i=0; i<nat; i++){
        if(randf()>0.5){ ff.atoms[i].type=&type1;  }else{ ff.atoms[i].type=&type2; }
        ff.atoms[i].pos.fromRandomBox((Vec3d){-5.0,-5.0,-1.0},(Vec3d){5.0,5.0,1.0});
        ff.atoms[i].qrot.setRandomRotation();
        ff.atoms[i].cleanAux();
    }



    //ff.realloc(2);
    //ff.atoms[0].setPose( (Vec3d){1.0,0.0, -1.0}, Quat4dIdentity ); ff.atoms[0].type=&type2;
    //ff.atoms[1].setPose( (Vec3d){0.0,1.0, -1.0}, Quat4dIdentity ); ff.atoms[1].type=&type2;

    //ff.realloc(2);
    //ff.atoms[0].setPose( (Vec3d){0.2,1.0, -1.0}, Quat4dIdentity ); ff.atoms[0].type=&type1;
    //ff.atoms[1].setPose( (Vec3d){1.0,0.0, -1.0}, Quat4dFront    ); ff.atoms[1].type=&type1;


    atom1.type = &type1;
    atom1.pos  = Vec3dZero;
    atom1.qrot = Quat4dIdentity;

    int nx      = 100;
    int ny      = 100;
    npoints = nx*ny;

    _realloc(points  ,npoints);
    _realloc(Energies,npoints);
    _realloc(Forces  ,npoints);

    makeSamples({nx,ny},{-5.0,-5.0,0.2},{10.0,0.0,0.0},{0.0,10.0,0.0},points);

    RigidAtomType pairType;
    pairType.combine(type1,type1);

    printf(" >>> pairType: <<<<\n");
    pairType.print();
    //exit(0);


    Vec3d bhs[N_BOND_MAX];
    //rotateVectors<double>(N_BOND_MAX, atom1.qrot, atom1.type->bh0s, bhs );
    atom1.qrot.rotateVectors(N_BOND_MAX, atom1.type->bh0s, bhs, false );
    for(int i=0; i<npoints; i++){
        Vec3d torq;
        Vec3d dij   = points[i] - atom1.pos;
        Energies[i] = ff.pairEF( dij, pairType, bhs, Forces[i], torq );
        //printf( "%i p:(%g,%g,%g) E: %g\n", i, points[i].x, points[i].y, points[i].z, Energies[i] );
    }
    //exit(0);



    // PLOT FOCRE FIELD

    DataLine2D * line_Er = new DataLine2D(100);
    line_Er->linspan(0,10.0);
    line_Er->clr = 0xFF00FF00;

    DataLine2D * line_Fr = new DataLine2D(100);
    line_Fr->linspan(0,10.0);
    line_Fr->clr = 0xFF0000FF;

    DataLine2D * line_Frn = new DataLine2D(100);
    line_Frn->linspan(0,10.0);
    line_Frn->clr = 0xFFFF00FF;

    plot1.init();
    plot1.fontTex = fontTex;
    plot1.clrGrid = 0xFF404040;
    plot1.clrBg   = 0xFF408080;
    //plot1.lines.push_back( line1  );
    plot1.lines.push_back( line_Er  );
    plot1.lines.push_back( line_Fr  );
    plot1.lines.push_back( line_Frn );
    plot1.render();

    Vec3d p0 = (Vec3d){0.0,0.0,0.0};
    Vec3d dp = (Vec3d){1.0,0.0,0.0};   dp.normalize();

    for(int i=0; i<line_Er->n; i++){
        Vec3d f,fnum,torq;
        double x = line_Er->xs[i];
        Vec3d dij   = (p0+dp*x) - atom1.pos;
        double E = ff.pairEF( dij, pairType, bhs, f, torq );
        numDeriv( dij, 0.01, fnum, [&](Vec3d p){
            Vec3d f_,tq_;
            return ff.pairEF( p, pairType, bhs, f_, tq_ );
        } );
        //printf( "%i p(%g,%g,%g) E %g f(%g,%g,%g) fnum(%g,%g,%g) \n", i, dij.x, dij.y, dij.z, E,   f.x,f.y,f.z,  fnum.x,fnum.y,fnum.z );

        line_Er ->ys[i]  = E*10.0;
        line_Fr ->ys[i] = dp.dot(f   )*-10.0;
        line_Frn->ys[i] = dp.dot(fnum)*-10.0;

        printf( "line_Er %i x %g y %g dy %g dyn %g \n", i, line_Er->xs[i], line_Er->ys[i], line_Fr->ys[i], line_Frn->ys[i] );
    }
    //exit(0);
    plot1.render();

    // PLOT FOCRE FIELD 1D

    ff.cleanAtomForce();
    ff.interEF();
    VecN::minmax(npoints, Energies, Emin, Emax);

}

void TestAppRARFF::draw(){
    printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);

    if(bRun){
        ff.cleanAtomForce();
        ff.interEF();
        //ff.move(0.005);
        ff.moveMDdamp(0.1, 0.9);

        /*
        double cosdRot,cosdPos;
        ff.getCos( cosdRot, cosdPos );
        double damp = fmax( 0.9, fmin(cosdRot,cosdPos) );
        printf( "cosdRot %g cosdRot %g damp %g \n", cosdRot, cosdPos, damp );
        ff.moveMDdamp(0.1, damp );
        */

    }

    //Vec3d bhs[N_BOND_MAX];
    //atom1.torq = (Vec3d){0.1,0.0,0.0};
    //atom1.moveRotGD(0.8);
    //printf( "qrot (%g,%g,%g,%g)\n", atom1.qrot.x, atom1.qrot.y, atom1.qrot.z, atom1.qrot.w );

    glColor3f(1.0,1.0,1.0);
    //drawRigidAtom( atom1 );

    double fsc = 0.1;
    double tsc = 0.1;
    for(int i=0; i<ff.natom; i++){
        glColor3f(1.0,1.0,1.0); drawRigidAtom(ff.atoms[i]);
        glColor3f(1.0,0.0,0.0); Draw3D::drawVecInPos( ff.atoms[i].force*fsc, ff.atoms[i].pos  );
        glColor3f(0.0,0.0,1.0); Draw3D::drawVecInPos( ff.atoms[i].torq*tsc,  ff.atoms[i].pos  );
    };


/*
    printf("npoints %i Emin %g Emax %g \n",npoints, Emin, Emax);
    glPointSize(5);
    drawScalarArray( npoints, points, Energies, Emin, Emax );
*/
/*
    glColor3f(0.0,1.0,0.0);
    drawVectorArray( npoints, points, Forces, 0.02 );
*/

    Draw3D::drawAxis( 1.0);

    //exit(0);
};


void TestAppRARFF::drawHUD(){
/*
    glColor3f(1.0,1.0,1.0);
    txt.viewHUD( {100,220}, fontTex );

	gui.draw();

	glTranslatef( 10.0,300.0,0.0 );
	glColor3f(0.5,0.0,0.3);
    Draw::drawText( "abcdefghijklmnopqrstuvwxyz \n0123456789 \nABCDEFGHIJKLMNOPQRTSTUVWXYZ \nxvfgfgdfgdfgdfgdfgdfg", fontTex, fontSizeDef, {10,5} );
*/

	glTranslatef( 100.0,100.0,0.0 );
	glScalef    ( 10.0,10.00,1.0  );
	plot1.view();

}


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
















