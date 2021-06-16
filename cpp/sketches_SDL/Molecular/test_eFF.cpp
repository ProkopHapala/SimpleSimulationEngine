
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

int DEBUG_i = 0;
int DEBUG_j = 0;

Vec3d* DEBUG_fe_ae =0;
Vec3d* DEBUG_fa_ae =0;
Vec3d* DEBUG_fe_ee =0;
Vec3d* DEBUG_fa_aa =0;


#include "DynamicOpt.h"
#include "InteractionsGauss.h"
#include "eFF.h"
//#include "e2FF.h" // old currently not working



#include "Forces.h"

#include "eFF_plots.h"


void cleanDebugForce(int ne, int na){
    for(int i=0; i<ne; i++){
        DEBUG_fe_ae[i]=Vec3dZero;
        DEBUG_fe_ee[i]=Vec3dZero;
    }
    for(int i=0; i<na; i++){
        DEBUG_fa_ae[i]=Vec3dZero;
        DEBUG_fa_aa[i]=Vec3dZero;
    }
};

void applyCartesianBoxForce( const Vec3d& pmin, const Vec3d& pmax,const Vec3d& k, int n, const Vec3d* ps, Vec3d* fs ){
    for(int i=0;i<n; i++){ boxForce( ps[i], fs[i], pmin, pmax, k ); }
   //for(int i=0;i<n; i++){     if(i==0){ boxForce( ps[i], fs[i], pmin, pmax, k ); printf( " atom[%i] p(%g,%g,%g) f(%g,%g,%g) \n", i, ps[i].x, ps[i].y, ps[i].z,   fs[i].x, fs[i].y, fs[i].z ); }  }
}

// ======= THE CLASS

class TestAppRARFF: public AppSDL2OGL_3D { public:

    bool bRun = false;

    //E2FF ff2;
    EFF  ff;
    DynamicOpt opt;

    int perFrame = 1;

    Vec2i field_ns;
    Vec2d Erange;
    double E0,Espread;
    Vec3d  * field_ps=0;
    double * field_Es=0;
    std::function<void   (const Vec3d& p, Vec3d& f)>  FFfunc;
    std::function<double (const Vec3d& p)          >  Efunc ;

    bool bMapElectron = false;
    int ipicked  = 0;

    // DEBUG STUFF
    GLint ogl_fs = 0;
    GLint oglSph = 0;

    Plot2D plot1;

    int      fontTex;

    // ---- Functions

    virtual void draw   ();
    virtual void drawHUD();
    //virtual void mouseHandling( );
    virtual void eventHandling   ( const SDL_Event& event  );
    //virtual void keyStateHandling( const Uint8 *keys );

    TestAppRARFF( int& id, int WIDTH_, int HEIGHT_ );
    void init2DMap( int n, double dx );

};

TestAppRARFF::TestAppRARFF( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );
    plot1.fontTex=fontTex;

    //checkDerivs( ff.KRSrho );   // exit(0);
    //makePlots( plot1, ff );     // exit(0);
    //makePlots2( plot1 );        // exit(0);
    //checkDerivs2();             // exit(0);

    // ===== SETUP GEOM
    //char* fname = "data/H_eFF.xyz";
    char* fname = "data/H2_eFF.xyz";
    //char* fname = "data/H2O_eFF.xyz";
    //char* fname = "data/H2_eFF_spin.xyz";
    //char* fname = "data/Ce1_eFF.xyz";
    //char* fname = "data/Ce2_eFF.xyz";
    //char* fname = "data/Ce4_eFF.xyz";
    //char* fname = "data/CH3_eFF_spin.xyz";
    //char* fname = "data/CH4_eFF_flat_spin.xyz";
    //char* fname = "data/CH4_eFF_spin.xyz";
    //char* fname = "data/C2_eFF_spin.xyz";
    //char* fname = "data/C2H4_eFF_spin.xyz";
    //char* fname = "data/C2H4_eFF_spin_.xyz";
    //char* fname = "data/C2H6_eFF_spin.xyz";
    //char* fname = "data/C2H6_eFF_spin_.xyz";
    //ff.loadFromFile_xyz( "data/C2H4_eFF_spin.xyz" );
    ff.loadFromFile_xyz( fname );

    DEBUG_fe_ae = new Vec3d[ff.ne];
    DEBUG_fa_ae = new Vec3d[ff.na];
    DEBUG_fe_ee = new Vec3d[ff.ne];
    DEBUG_fa_aa = new Vec3d[ff.na];

    //setGeom(ff);
    //double sz = 0.2;
    //for(int i=0; i<ff.na; i++){ ff.apos[i].add( randf(-sz,sz),randf(-sz,sz),randf(-sz,sz) );  }
    //for(int i=0; i<ff.ne; i++){ ff.epos[i].add( randf(-sz,sz),randf(-sz,sz),randf(-sz,sz) );  }
    //ff.autoAbWs( default_aAbWs, default_eAbWs );
    //VecN::set(ff.ne,4.0,ff.esize);

    // ==== Test Eval

    //makePlots( plot1, ff );
    //ff.loadFromFile_xyz( fname  );
    //init2DMap( 100, 0.1 );

    // ==== Bind Optimizer
    opt.bindOrAlloc( ff.nDOFs, ff.pDOFs, 0, ff.fDOFs, 0 );
    opt.cleanVel( );
    opt.initOpt( 0.05, 0.1 );
    opt.f_limit = 1000.0;

    ff.iPauliModel = 2; // dens overlap
    ff.info();

    double E = ff.eval();
    printf( "E %g | Ek %g Eee %g EeePaul %g Eaa %g Eae %g EaePaul %g \n", E, ff.Ek, ff.Eee, ff.EeePaul, ff.Eaa, ff.Eae, ff.EaePaul );
    //exit(0);

    oglSph=Draw::list(oglSph);
    Draw3D::drawSphere_oct(3,1.0d,(Vec3d){0.,0.,0.});
    glEndList();

    DEBUG
}

void TestAppRARFF::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);

    //return;

    double vminOK = 1e-6;
    double vmaxOK = 1e+3;

    //perFrame=10; // ToDo : why it does not work properly for  perFrame>1 ?
    perFrame = 1;
    double sum = 0;
    if(bRun){
        double F2 = 1.0;
        for(int itr=0;itr<perFrame;itr++){
            //printf( " ==== frame %i i_DEBUG  %i \n", frameCount, i_DEBUG );

            ff.clearForce();
            //ff.clearForce_noAlias();
            cleanDebugForce( ff.ne, ff.na);
            //VecN::sum( ff.ne*3, ff.eforce );

            //applyCartesianBoxForce( {0.0,0.0,0.0}, {0.0,0.0,0.0}, {0,0,50.0}, ff.na, ff.apos, ff.aforce );
            //applyCartesianBoxForce( {0.0,0.0,0.0}, {0.0,0.0,0.0}, {0,0,5.0},  ff.ne, ff.epos, ff.eforce );
            double E = ff.eval();
            //ff.apos[0].set(.0);
            //checkFinite( ff, vminOK, vmaxOK );

            //printf( "fa1(%g,%g,%g) fe1(%g,%g,%g)\n", fa1.x,fa1.x,fa1.x,   fe1.x,fe1.x,fe1.x );

            //VecN::set( ff.na*3, 0.0, (double*)ff.aforce );   // FIX ATOMS
            //VecN::set( ff.ne, 0.0, ff.fsize );               // FIX ELECTRON SIZE
            //if(bRun)ff.move_GD(0.001 );

            ff.move_GD( 0.001 );
            //if(bRun) ff.move_GD_noAlias( 0.0001 );

            //F2 = opt.move_FIRE();

            //checkFinite( ff, vminOK, vmaxOK );

            //printf( "frame[%i] E %g pa[0](%g,%g,%g) pe[0](%g,%g,%g) \n", frameCount, E,   ff.apos[0].x,ff.apos[0].y,ff.apos[0].z,   ff.epos[0].x,ff.epos[0].y,ff.epos[0].z );
            //printf( "frame[%i] E %g pe[0](%g,%g,%g) s %g fe[0](%g,%g,%g) fs %g \n", frameCount, E,   ff.epos[0].x,ff.epos[0].y,ff.epos[0].z,  ff.esize[0],   ff.eforce[0].x,ff.eforce[0].y,ff.eforce[0].z, ff.fsize[0] );

            printf( "frame[%i] E %g pa[0](%g,%g,%g) pe[0](%g,%g,%g) s %g \n", frameCount, E, ff.apos[0].x,ff.apos[0].y,ff.apos[0].z,   ff.epos[0].x,ff.epos[0].y,ff.epos[0].z,  ff.esize[0] );

            //printf( "E %g | Ek %g Eee %g EeePaul %g Eaa %g Eae %g EaePaul %g \n", E, ff.Ek, ff.Eee, ff.EeePaul, ff.Eaa, ff.Eae, ff.EaePaul );
            //printf( "=== %i %i frame[%i][%i] |F| %g \n", ff.na, ff.ne, frameCount, itr, sqrt(F2) );
        }
        if( F2 < 1e-6 ){
            ff.info();
            bRun=false;
        }
    }

    //printf( "e[0] r %g s %g \n", ff.epos[0].norm(), ff.esize[0] );

    //printf( "r07 r %g s %g \n", ff.epos[0].norm(), ff.esize[0] );

    // --- Constrain in Z
    //double Kz = 1.0;
    //for(int i=0; i<ff.na; i++){  ff.aforce[i].z += -Kz*ff.apos[i].z;  };
    //for(int i=0; i<ff.ne; i++){  ff.eforce[i].z += -Kz*ff.epos[i].z;  };
    //printf( "na %i ne %i \n", ff.na, ff.ne );


    //Vec3d d = ff.apos[0]-ff.apos[1];

    glCallList(ogl_fs);
    //Draw3D::drawColorScale( 20, {0.0,0.0,0.0}, Vec3dY, Vec3dX, Draw::colors_rainbow, Draw::ncolors );
    //printf( "apos (%g,%g,%g) \n", ff.apos[0].x, ff.apos[0].y, ff.apos[0].z );


    char strtmp[256];
    double Qsz = 0.05;
    double fsc = 1.0;

    for(int i=0; i<ff.na; i++){
        //printf( "apos[%i] (%g,%g,%g)\n", i, ff.apos[i].x, ff.apos[i].y, ff.apos[i].z );
        glColor3f(0.0,0.0,0.0); Draw3D::drawPointCross( ff.apos  [i]    , ff.aPars[i].x*Qsz );
        glColor3f(1.0,0.0,0.0); Draw3D::drawVecInPos  ( ff.aforce[i]*fsc, ff.apos[i] );
        //Draw3D::drawVecInPos(   ff.aforce[i]*fsc, ff.apos[i] );
        //printf( " %i %f %f %f %f  \n", i, ff.aQ[i], ff.apos[i].x,ff.apos[i].y,ff.apos[i].z );
        //printf( " %i %f %f %f %f  \n", i, ff.aQ[i], ff.aforce[i].x, ff.aforce[i].y, ff.aforce[i].z );
        //sprintf(strtmp,"%i",i);
        //Draw3D::drawText(strtmp, ff.apos[i], fontTex, 0.02, 0);

        //glColor3f(0.,0.,0.); Draw3D::drawVecInPos( ff.aforce[i],   ff.apos[i] );
        //glColor3f(0.,1.,0.); Draw3D::drawVecInPos( DEBUG_fa_ae[i], ff.apos[i] );
        //glColor3f(1.,0.,0.); Draw3D::drawVecInPos( DEBUG_fa_aa[i], ff.apos[i] );

    }

    //glColor3f(1.0,1.0,1.0);
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    //glBlendFunc(GL_DST_COLOR, GL_SRC_ALPHA);
    for(int i=0; i<ff.ne; i++){
        //printf( "epos[%i] (%g,%g,%g)\n", i, ff.epos[i].x, ff.epos[i].y, ff.epos[i].z );
        //if(ff.espin[i]>0){ glColor3f(0.0,0.5,1.0); }else{ glColor3f(1.0,0.5,0.0); };
        float alpha=0.1;
        if(ff.espin[i]>0){ glColor4f(0.0,0.0,1.0, alpha); }else{ glColor4f(1.0,0.0,0.0, alpha); };
        Draw3D::drawShape( oglSph, ff.epos[i], Mat3dIdentity*ff.esize[i],  false );
        glColor3f(1.0,0.0,0.0); Draw3D::drawVecInPos(  ff.eforce[i]*fsc ,ff.epos[i] );
        //glColor3f(1.0,0.0,0.0); Draw3D::drawPointCross( ff.epos[i], 0.1 );
        //Draw3D::drawSphere_oct(3,ff.esize[i],ff.epos[i]);

        //glColor3f(1.,1.,1.); Draw3D::drawVecInPos( ff.eforce  [i], ff.epos[i] );
        //glColor3f(1.,0.,1.); Draw3D::drawVecInPos( DEBUG_fe_ae[i], ff.epos[i] );

        //Draw3D::drawPointCross( ff.epos  [i], ff.esize[i] );
        //Draw3D::drawVecInPos(   ff.eforce[i]*fsc, ff.epos[i] );
        //Draw3D::drawSphereOctLines( 8, ff.esize[i], ff.epos[i] );
        //sprintf(strtmp,"%i",i);
        //Draw3D::drawText(strtmp, ff.epos[i], fontTex, 0.02, 0);

    }


    //for(int i=0; i<ff.ne; i+=2){
    //    Draw3D::drawLine(ff.epos[i],ff.epos[i+1] );
    //}

    //exit(0);

    //ff.aforce[0].set(0.);
    //ff.aforce[1].set(0.);
    //if(bRun) ff.move_GD( 0.01 );

    //glDisable(GL_DEPTH_TEST);
    //plot1.view();

};


void TestAppRARFF::drawHUD(){
	glTranslatef( 100.0, 250.0, 0.0 );
	glScalef    ( 100.0, 100.0, 1.0 );
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

                case SDLK_i: ff.info(); break;
                case SDLK_LEFTBRACKET :  Espread *= 1.2; ogl_fs=genFieldMap(ogl_fs, field_ns, field_ps, field_Es, E0-Espread, E0+Espread ); break;
                case SDLK_RIGHTBRACKET:  Espread /= 1.2; ogl_fs=genFieldMap(ogl_fs, field_ns, field_ps, field_Es, E0-Espread, E0+Espread ); break;
                case SDLK_e: bMapElectron=!bMapElectron; break;
                case SDLK_f:{
                    pa0 = ff.apos[ipicked];
                    sampleScalarField( Efunc, field_ns, {-5.0,-5.0,+0.1}, {0.1,0.0,0.0}, {0.0,0.1,0.0}, field_ps, field_Es, Erange );
                    E0 = field_Es[0];
                    ogl_fs = genFieldMap( ogl_fs, field_ns, field_ps, field_Es, E0-Espread, E0+Espread );
                    ff.apos[ipicked]= pa0;
                    }break;
                case SDLK_SPACE: bRun = !bRun; break;
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

void TestAppRARFF::init2DMap( int n, double dx ){

    int ipicked = 1;

    FFfunc = [&](const Vec3d& p, Vec3d& f)->void{
        if(bMapElectron){ ff.epos[ipicked] = p; }else{ ff.apos[ipicked] = p; }
         // force on one electron
        ff.clearForce();
        ff.eval();
        if(bMapElectron){ f = ff.eforce[ipicked]; }else{ f = ff.aforce[ipicked]; }
    };

    Efunc = [&](const Vec3d& p)->double{
        //ff.apos[ipicked] = p; // force on one electron
        if(bMapElectron){ ff.epos[ipicked] = p; }else{ ff.apos[ipicked] = p; }
        ff.clearForce();
        double E = ff.eval();
        return E;
    };

    field_ns = {n,n};

    Vec3d pa0;// = ff.apos[ipicked];
    if(bMapElectron){ pa0 = ff.epos[ipicked]; }else{ pa0 = ff.apos[ipicked]; }

    sampleScalarField( Efunc, field_ns, {-(dx*n)/2,-(dx*n)/2,+0.1}, {dx,0.0,0.0}, {0.0,dx,0.0}, field_ps, field_Es, Erange );
    if(bMapElectron){ ff.epos[ipicked] = pa0; }else{ ff.apos[ipicked] = pa0; }

    E0 = field_Es[0];
    printf( "val_range: %g %g %g \n", Erange.x, Erange.y, field_Es[0] );
    Espread = 3.0;

    ogl_fs = genFieldMap( ogl_fs, field_ns, field_ps, field_Es, E0-Espread, E0+Espread );

    //return ogl_fs;
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
















