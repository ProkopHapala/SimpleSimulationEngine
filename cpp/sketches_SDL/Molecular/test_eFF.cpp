
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>
#include  <functional>

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

//#include "MMFF.h"

#define R2SAFE  1.0e-8f

int i_DEBUG = 0;
#include "DynamicOpt.h"
#include "eFF.h"
#include "e2FF.h"



/*
int pickParticle( int n, Vec3d * ps, const Mat3d& cam, double R ){
    double tmin =  1e+300;
    int imin    = -1;
    for(int i=0; i<n; i++GL_LIGHTING){
        double x = cam.a.dot(ps[i]);
        doyble y = cam.b.dot(ps[i]);
        double ti = raySphere( ray0, hRay, R, ps[i] );
        if(ti<tmin){ imin=i; tmin=ti; }
    }
    return imin;
}
*/


int genFieldMap( int ogl, Vec2i ns, const Vec3d* ps, const double* Es, double vmin, double vmax ){
    //printf( "val_range: %g %g %g \n", val_range.x, val_range.y, Es[0] );
    //float clSz = 3.0;
    if(ogl) glDeleteLists(ogl,1);
    ogl = glGenLists(1);
    glNewList(ogl, GL_COMPILE);
    //glColor3f(1.0,0.0,0.0);
    //glPointSize(2.0);
    //Draw3D::drawPoints(nptot, ps, -1.0 );
    //Draw3D::drawVectorArray( nptot, ps, fs, 0.5, 0.5 );
    //Draw3D::drawScalarArray( nptot, ps, Es, 0.0, 1.0, Draw::colors_rainbow, Draw::ncolors );
    //val_range={Es[0]-clSz,Es[0]+clSz};
    //Draw3D::drawScalarGrid( {100,100}, {-5.0,-5.0,0.0}, {0.1,0.0,0.0}, {0.0,0.1,0.0}, Es, val_range.x, val_range.y, Draw::colors_RWB, Draw::ncolors );
    //Draw3D::drawScalarArray( ps, Es, val_range.x, val_range.y, Draw::colors_RWB, Draw::ncolors );
    Draw3D::drawScalarField( ns, ps, Es, vmin, vmax, Draw::colors_RWB, Draw::ncolors );
    Draw3D::drawColorScale( 20, {5.0,-5.0,0.0}, Vec3dY*10, Vec3dX*0.5, Draw::colors_RWB, Draw::ncolors );
    //exit(0);
    glEndList();
    return ogl;
}


void applyCartesianBoxForce( const Vec3d& pmin, const Vec3d& pmax,const Vec3d& k, int n, const Vec3d* ps, Vec3d* fs ){
    for(int i=0;i<n; i++){ boxForce( ps[i], fs[i], pmin, pmax, k ); }
   //for(int i=0;i<n; i++){     if(i==0){ boxForce( ps[i], fs[i], pmin, pmax, k ); printf( " atom[%i] p(%g,%g,%g) f(%g,%g,%g) \n", i, ps[i].x, ps[i].y, ps[i].z,   fs[i].x, fs[i].y, fs[i].z ); }  }
}


// ======= THE CLASS

class TestAppRARFF: public AppSDL2OGL_3D { public:

    //RigidAtom     atom1;
    //RigidAtomType type1,type2;

    bool bRun = false;

    EFF  ff;
    DynamicOpt opt;

    int perFrame = 1;

    //E2FF ff2;
    int ipicked  = -1;

    // DEBUG STUFF
    GLint ogl_fs = 0;

    Plot2D plot1;

    Vec2d Erange;
    double E0=0, Espread=1.0;
    Vec2i field_ns;
    Vec3d  * field_ps=0;
    double * field_Es=0;

    std::function<void   (const Vec3d& p, Vec3d& f)>  FFfunc;
    std::function<double (const Vec3d& p)          >  Efunc ;

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

    // ===== EVAL FF 1D Curves

    //ff.wee = 0.25;
    //ff.wae = 0.25;
    //ff.waa = 0.25;

    /*
    double w2ee = sq(ff.wee);
    double w2ae = sq(ff.wae);
    double w2aa = sq(ff.waa);
    double qaa  =  4*4;
    double qae  = -1.0;

    printf( " w2 ee,ae,aa(%g,%g,%g)  w ee,ae,aa(%g,%g,%g)  \n", w2ee, w2ae, w2aa   ,  ff.wee, ff.wae, ff.waa  );

    plot1.xsharingLines(4, 100, 0.0, 0.1);
    DataLine2D *l,*l_;

    l=plot1.lines[0]; l->clr=0xFFFF0000; l->label="Eee"; evalLine( *l, [&](double x){ Vec3d f;  return addPairEF_expQ( {x,0,0}, f, w2ee, +1.0,  ff.bEE, ff.aEE  ); } );
    l=plot1.lines[1]; l->clr=0xFFFF00FF; l->label="Eae"; evalLine( *l, [&](double x){ Vec3d f;  return addPairEF_expQ( {x,0,0}, f, w2ae, qae,   ff.bAE, ff.aAE  ); } );
    //l=plot1.lines[3]; l->clr=0xFF0080FF; l->label="Faa"; evalLine( *l, [&](double x){ Vec3d f;  return addPairEF_expQ( {x,0,0}, f, w2aa, qaa, 0,      0           ); } );

    l=plot1.lines[2]; l->clr=0xFFFF8000; l->label="Fee"; evalLine( *l, [&](double x){ Vec3d f=Vec3dZero;  addPairEF_expQ( {x,0,0}, f, w2ee, +1.0, ff.bEE, ff.aEE  ); return -f.x; } );
    l=plot1.lines[3]; l->clr=0xFFFF80FF; l->label="Fae"; evalLine( *l, [&](double x){ Vec3d f=Vec3dZero;  addPairEF_expQ( {x,0,0}, f, w2ae, qae,  ff.bAE, ff.aAE  ); return -f.x; } );
    //l=plot1.lines[2]; l->clr=0xFF0000FF; l->label="Eaa"; evalLine( *l, [&](double x){ Vec3d f=Vec3dZero;  addPairEF_expQ( {x,0,0}, f, w2aa, qaa, 0,      0           ); return f.x; } );

    plot1.update();
    plot1.autoAxes(0.5,0.5);
    printf( "axBound %g,%g %g,%g \n", plot1.axBounds.a.x, plot1.axBounds.a.y, plot1.axBounds.b.x, plot1.axBounds.b.y );
    plot1.render();
    */

        // ===== SETUP GEOM

    //ff.loadFromFile_bas( "data/CH4.bas" );
    //ff.loadFromFile_bas( "data/C2H6.bas" );
    //ff.loadFromFile_bas( "data/C_eFF.bas" );
    //ff.loadFromFile_bas( "data/CH4_eFF.bas" );
    //ff.loadFromFile_bas( "data/C2H6_e2FF.bas" );
    //ff.loadFromFile_bas( "data/C2.bas" );
    //ff.loadFromFile_bas( "data/H2.bas" );
    //ff.loadFromFile_bas( "data/C2e2.bas" );
    //ff.loadFromFile_bas( "data/H-e.bas" );
    //ff.loadFromFile_xyz( "data/C2H4_eFF_pairs.xyz" );
    //ff.loadFromFile_xyz( "data/C2H4_eFF_spin.xyz" );
    ff.loadFromFile_xyz( "data/H2_eFF_spin.xyz" );

    /*
    ff.realloc( 1, 2 );
    ff.aQ  [0] = +1.0;
    ff.apos[0] = (Vec3d){ 0.0, 0.0,0.0};
    ff.epos[0] = (Vec3d){+1.0, 0.0,0.0};
    ff.epos[1] = (Vec3d){-1.0, 0.0,0.0};
    */

    /*
    ff.realloc( 1, 4 );
    ff.aQ  [0] = +4.0;
    ff.apos[0] = (Vec3d){0.0,0.0,0.0};
    ff.epos[0] = (Vec3d){+1.0, 0.0,0.0};
    ff.epos[1] = (Vec3d){-1.0, 0.0,0.0};
    ff.epos[2] = (Vec3d){ 0.0,+1.0,0.0};
    ff.epos[3] = (Vec3d){ 0.0,-1.0,0.0};
    */

    /*
    ff.realloc( 2, 8 );
    ff.aQ  [0] = +4.0;
    ff.aQ  [1] = +4.0;
    ff.apos[0] = (Vec3d){-1.0,0.0,0.0};
    ff.apos[1] = (Vec3d){ 1.0,0.0,0.0};

    ff.epos[0] = (Vec3d){ 0.0,+0.5,0.0};
    ff.epos[1] = (Vec3d){-2.0, 0.0,0.0};
    ff.epos[2] = (Vec3d){-1.0,+1.0,0.0};
    ff.epos[3] = (Vec3d){-1.0,-1.0,0.0};

    ff.epos[4] = (Vec3d){ 0.0,-0.5,0.0};
    ff.epos[5] = (Vec3d){+2.0, 0.0,0.0};
    ff.epos[6] = (Vec3d){+1.0,+1.0,0.0};
    ff.epos[7] = (Vec3d){+1.0,-1.0,0.0};
    */

    //setGeom(ff);
    //double sz = 0.51;
    // break symmetry
    //for(int i=0; i<ff.na; i++){ ff.apos[i].add( randf(-sz,sz),randf(-sz,sz),randf(-sz,sz) );  }
    for(int i=0; i<ff.na; i++){ printf( "A_pos[%i] (%g,%g,%g)\n", i, ff.apos[i].x, ff.apos[i].y, ff.apos[i].z ); }
    for(int i=0; i<ff.ne; i++){ printf( "e_pos[%i] (%g,%g,%g)\n", i, ff.epos[i].x, ff.epos[i].y, ff.epos[i].z ); }

    ff.autoAbWs( default_aAbWs, default_eAbWs );
    //exit(0);

    opt.bindOrAlloc( 3*ff.nDOFs, (double*)ff.pDOFs, 0, (double*)ff.fDOFs, 0 );
    opt.cleanVel( );
    opt.initOpt( 0.05, 0.2 );

    // ===== EVAL FF 2D MAP



    i_DEBUG = 1;
    ff.eval();
    i_DEBUG = 0;

    Vec3d  *ps=0,*fs=0;
    double *Es=0;

    ipicked = 1;

    FFfunc = [&](const Vec3d& p, Vec3d& f)->void{
        ff.apos[ipicked] = p;
        //ff.epos[0] = p; // force on one electron
        ff.clearForce();
        ff.eval();
        //f = ff.eforce[0];
        f = ff.aforce[ipicked];
    };

    Efunc = [&](const Vec3d& p)->double{
        ff.apos[ipicked] = p; // force on one electron
        ff.clearForce();
        double E = ff.eval();
        //f = ff.eforce[0];
        return E;
    };

    //exit(0);
    //ff.move( 0.01, 0.9 );

    Vec3d pa0 = ff.apos[ipicked];
    field_ns = {100,100};
    //int nptot = sampleVecField<fffunc>( (Vec2i){60,60}, (Vec3d){-3.0,-3.0,0.0}, (Vec3d){0.1,0.0,0.0}, (Vec3d){0.0,0.1,0.0}, ps, fs );
    //int nptot = sampleVecField( fffunc, (Vec2i){60,60}, (Vec3d){-3.0,-3.0,0.0}, (Vec3d){0.1,0.0,0.0}, (Vec3d){0.0,0.1,0.0}, ps, fs );
    //int nptot = sampleVecField( fffunc, {100,100}, {-5.0,-5.0,0.0}, {0.1,0.0,0.0}, {0.0,0.1,0.0}, ps, fs );
    //int nptot =
    sampleScalarField( Efunc, field_ns, {-5.0,-5.0,+0.1}, {0.1,0.0,0.0}, {0.0,0.1,0.0}, field_ps, field_Es, Erange );
    ff.apos[ipicked] = pa0;

    E0 = field_Es[0];
    printf( "val_range: %g %g %g \n", Erange.x, Erange.y, field_Es[0] );
    Espread = 3.0;

    ogl_fs = genFieldMap( ogl_fs, field_ns, field_ps, field_Es, E0-Espread, E0+Espread );

    //setGeom(ff);
    //double sz = 0.51;
    // break symmetry
    //for(int i=0; i<ff.na; i++){ ff.apos[i].add( randf(-sz,sz),randf(-sz,sz),randf(-sz,sz) );  }
    for(int i=0; i<ff.na; i++){ printf( "A_pos[%i] (%g,%g,%g)\n", i, ff.apos[i].x, ff.apos[i].y, ff.apos[i].z ); }
    for(int i=0; i<ff.ne; i++){ printf( "e_pos[%i] (%g,%g,%g)\n", i, ff.epos[i].x, ff.epos[i].y, ff.epos[i].z ); }


    //ff.apos[0].z = 2.0;

}

void TestAppRARFF::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);

    glDisable(GL_LIGHTING);

    if(bRun){
        for(int itr=0;itr<perFrame;itr++){
            printf( " ==== frame %i i_DEBUG  %i \n", frameCount, i_DEBUG );
            double F2 = 1.0;

            ff.clearForce();
            applyCartesianBoxForce( {0.0,0.0,0.0}, {0.0,0.0,0.0}, {0,0,50.0}, ff.na, ff.apos, ff.aforce );
            //applyCartesianBoxForce( {0.0,0.0,0.0}, {0.0,0.0,0.0}, {0,0,5.0}, ff.ne, ff.epos, ff.eforce );
            //ff.evalEE();
            //ff.evalAE();
            //ff.evalAA();
            ff.eval();
            //ff.apos[0].set(.0);
            //if(bRun)ff.moveGD(0.001, 1, 1);

            F2 = opt.move_FIRE();

            printf( " |F| %g \n", sqrt(F2) );
            if(!(F2<1000000.0))perFrame=0;
        }
    }

    // --- Constrain in Z
    //double Kz = 1.0;
    //for(int i=0; i<ff.na; i++){  ff.aforce[i].z += -Kz*ff.apos[i].z;  };
    //for(int i=0; i<ff.ne; i++){  ff.eforce[i].z += -Kz*ff.epos[i].z;  };
    //printf( "na %i ne %i \n", ff.na, ff.ne );


    //ff.aforce[0].set(0.);
    //if(bRun) ff.move_GD( 0.01 );
    //if(bRun) ff.run( 1, 0.1, 0.5 );

    //Vec3d d = ff.apos[0]-ff.apos[1];

    //printf("C1-C2 %g C1-e %g C2-e %g \n", (ff.apos[0]-ff.apos[1]).norm(),
    //                                      (ff.apos[0]-ff.epos[0]).norm(),
    //                                      (ff.apos[1]-ff.epos[0]).norm() );

    glCallList(ogl_fs);
    //Draw3D::drawColorScale( 20, {0.0,0.0,0.0}, Vec3dY, Vec3dX, Draw::colors_rainbow, Draw::ncolors );
    //plot1.view();

    //printf( "apos (%g,%g,%g) \n", ff.apos[0].x, ff.apos[0].y, ff.apos[0].z );

    double Qsz = 0.05;
    double fsc = 1.0;
    glColor3f(0.0,0.0,0.0);
    for(int i=0; i<ff.na; i++){
        //printf( "apos[%i] (%g,%g,%g)\n", i, ff.apos[i].x, ff.apos[i].y, ff.apos[i].z );
        Draw3D::drawPointCross( ff.apos  [i]    , ff.aQ  [i]*Qsz );
        Draw3D::drawVecInPos(   ff.aforce[i]*fsc, ff.apos[i] );
        //printf( " %i %f %f %f %f  \n", i, ff.aQ[i], ff.apos[i].x,ff.apos[i].y,ff.apos[i].z );
        //printf( " %i %f %f %f %f  \n", i, ff.aQ[i], ff.aforce[i].x, ff.aforce[i].y, ff.aforce[i].z );
    }

    //glColor3f(1.0,1.0,1.0);
    for(int i=0; i<ff.ne; i++){
        //printf( "epos[%i] (%g,%g,%g)\n", i, ff.epos[i].x, ff.epos[i].y, ff.epos[i].z );
        if(ff.espin[i]>0){ glColor3f(0.0,0.5,1.0); }else{ glColor3f(1.0,0.5,0.0); };
        Draw3D::drawPointCross( ff.epos  [i],     Qsz );
        Draw3D::drawVecInPos(   ff.eforce[i]*fsc, ff.epos[i] );
    }

    //for(int i=0; i<ff.ne; i+=2){
    //    Draw3D::drawLine(ff.epos[i],ff.epos[i+1] );
    //}

    //exit(0);

    //ff.aforce[0].set(0.);
    //ff.aforce[1].set(0.);
    //if(bRun) ff.move_GD( 0.01 );

};


void TestAppRARFF::drawHUD(){
	//glTranslatef( 100.0, 100.0, 0.0 );
	//glScalef    (  10.0,  10.0, 1.0 );
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
    Vec3d pa0;
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;

                case SDLK_LEFTBRACKET :  Espread *= 1.2; ogl_fs=genFieldMap(ogl_fs, field_ns, field_ps, field_Es, E0-Espread, E0+Espread ); break;
                case SDLK_RIGHTBRACKET:  Espread /= 1.2; ogl_fs=genFieldMap(ogl_fs, field_ns, field_ps, field_Es, E0-Espread, E0+Espread ); break;
                case SDLK_e:
                    pa0 = ff.apos[ipicked];
                    sampleScalarField( Efunc, field_ns, {-5.0,-5.0,+0.1}, {0.1,0.0,0.0}, {0.0,0.1,0.0}, field_ps, field_Es, Erange );
                    E0 = field_Es[0];
                    ogl_fs = genFieldMap( ogl_fs, field_ns, field_ps, field_Es, E0-Espread, E0+Espread );
                    ff.apos[ipicked]= pa0;
                    break;
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
















