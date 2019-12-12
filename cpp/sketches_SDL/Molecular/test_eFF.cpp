
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

#include "Forces.h"

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


void printFFInfo(const EFF& ff){
    printf( "=== ElectronForcefield(ne %i,na %i) \n", ff.ne, ff.na );
    double Eel   = ff.Eaa + ff.Eae + ff.Eee;
    double EPaul = ff.EeePaul + ff.EaePaul;
    double Etot  = ff.Ek + Eel + EPaul;
    printf( "Etot %g Ek %g Eel %g(ee %g, ea %g aa %g)  EPaul %g(ee %g, ae %g) \n", Etot, ff.Ek, Eel, ff.Eee,ff.Eae,ff.Eaa,   EPaul, ff.EeePaul, ff.EaePaul );
    for(int i=0; i<ff.na; i++){
        printf( "a[%i] xyzs(%g,%g,%g) fxyzs(%g,%g,%g) \n", i, ff.apos[i].x, ff.apos[i].y, ff.apos[i].z, ff.aforce[i].x, ff.aforce[i].y, ff.aforce[i].z );
    }
    for(int i=0; i<ff.ne; i++){
        printf( "e[%i] xyzs(%g,%g,%g,%g) fxyzs(%g,%g,%g,%g) \n", i, ff.epos[i].x, ff.epos[i].y, ff.epos[i].z, ff.esize[i], ff.eforce[i].x, ff.eforce[i].y, ff.eforce[i].z, ff.fsize[i] );
    }
}



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
    printf( "%i %i %li %g %g %li %i \n", ns.x,ns.y, (long)Es, vmin, vmax, (long)Draw::colors_RWB, Draw::ncolors );
    Draw3D::drawScalarField( ns, ps, Es, vmin, vmax, Draw::colors_RWB, Draw::ncolors );
    Draw3D::drawColorScale( 20, {5.0,-5.0,0.0}, Vec3dY*10, Vec3dX*0.5, Draw::colors_RWB, Draw::ncolors );
    //exit(0);
    glEndList();
    return ogl;
}


void checkDerivs( Vec3d KRSrho ){

    double qqee = 1.0;
    double r0   = 1.5;
    double sj0  = 0.7;
    double si0  = 0.58;
    double S0   = 0.0638745;

    auto func_ds = [&]( double x, double& f  )->double{
        //addKineticGauss( double s, double& fs );
        double fr,fsj;
        Vec3d fvec=Vec3dZero; f=0;
        //double e = CoulombGauss( r0, x, fr, f, qqee );           f*=x;
        //double e = getDeltaTGauss( r0*r0, x, sj0, fr, f, fsj );
        //double e = getOverlapSGauss( r0*r0, x, sj0, fr, f, fsj );
        //double e = addPauliGauss( (Vec3d){0.0,0.0,r0}, x, sj0, fvec, f, fsj, true, KRSrho );
        //double e = addDensOverlapGauss_S( (Vec3d){0.0,0.0,r0}, x, sj0, 1.0, fvec, f, fsj );
        double e = addDensOverlapGauss_P( (Vec3d){0.0,0.0,r0}, x, sj0, 1.0, fvec, f, fsj );
        return e;
    };

    auto func_dr = [&]( double x, double& f  )->double{
        //addKineticGauss( double s, double& fs );
        double fsi,fsj;
        Vec3d fvec=Vec3dZero;
        //double e = CoulombGauss  ( x, si0, f, fsi, qqee );      f*=x;
        //double e = getDeltaTGauss( x*x, si0, sj0, f, fsi, fsj );  f*=x;
        //double e = getOverlapSGauss( x*x, si0, sj0, f, fsi, fsj ); f*=x;
        //double e = PauliSGauss_anti( x, f, 0.2 );
        //double e = PauliSGauss_syn ( x, f, 0.2 );
        //double e = addPauliGauss( {0.0,0.0,x}, si0, sj0, fvec, f, fsj, true, KRSrho ); f=fvec.z;
        //double e = addDensOverlapGauss_S( (Vec3d){0.0,0.0,x}, si0, sj0, 1.0, fvec, fsi, fsj ); f=fvec.z;
        double e = addDensOverlapGauss_P( (Vec3d){0.0,0.0,x}, si0, sj0, 1.0, fvec, fsi, fsj ); f=fvec.z;
        return e;
    };

    double fE,f;

    //checkDeriv( KineticGauss, x, 0.001, fE, f );
    checkDeriv( func_ds, si0, 0.001, fE, f );
    checkDeriv( func_dr, r0 , 0.001, fE, f );
    //checkDeriv( func_dr, S0, 0.001, fE, f );


}

int makePlots( Plot2D& plot, EFF& ff ){

    int ielem = 1;
    double QQae = -1.0;
    double QQaa = +1.0;
    double QQee = QE*QE;

    double bEE     = -1.0;
    double aEE     =  2.0;
    double bEEpair = -1.0;
    double aEEpair =  0.1;

    double w2ee = 1.5;

    Vec3d  eAbw = default_eAbWs[ielem];
    Vec3d  aAbw; combineAbW( default_eAbWs[ielem] , default_eAbWs[ielem], aAbw );

    int np = 100;

    plot.xsharingLines( 3, np, 0.001, 0.05 );
    DataLine2D *l;
    double si = 1.0;
    double sj = 1.0;
    l=plot.lines[0]; l->clr=0xFFFF0000; l->label="EDens"; evalLine( *l, [&](double x){ double fsi,fsj,E; Vec3d f; E=addDensOverlapGauss_S( {0,0,x}, si, sj, 1.0, f, fsi, fsj );                     return E; } );
    l=plot.lines[1]; l->clr=0xFF0000FF; l->label="EPaul"; evalLine( *l, [&](double x){ double fsi,fsj,E; Vec3d f; E=addPauliGauss        ( {0,0,x}, si, sj,      f, fsi, fsj, true,  EFF::KRSrho ); return E; } );
    l=plot.lines[2]; l->clr=0xFF0080FF; l->label="EPaul"; evalLine( *l, [&](double x){ double fsi,fsj,E; Vec3d f; E=addPauliGauss        ( {0,0,x}, si, sj,      f, fsi, fsj, false, EFF::KRSrho ); return E; } );


    /*
    // --- derivative along pos.z
    plot.xsharingLines( 3, np, 0.001, 0.05 );
    DataLine2D *l;
    double si = 1.25;
    double sj = 0.8;

    l=plot.lines[0]; l->clr=0xFFFFFFFF; l->label="EDens"; evalLine( *l, [&](double x){ double fsi,fsj,E; Vec3d f=Vec3dZero; E=addDensOverlapGauss_S( {0,0,x}, si, sj, 1.0, f, fsi, fsj ); return E; } );
    l=plot.lines[1]; l->clr=0xFF0000FF; l->label="F1";    evalLine( *l, [&](double x){ double fsi,fsj,E; Vec3d f=Vec3dZero; E=addDensOverlapGauss_S( {0,0,x}, si, sj, 1.0, f, fsi, fsj ); return f.z; } );
    l=plot.lines[2]; l->clr=0xFF0080FF; l->label="FeeNum"; evalLine( *l, [&](double x){ double fsi,fsj,E1,E2,dx=0.001; Vec3d f=Vec3dZero;
        E1=addDensOverlapGauss_S( {0,0,x-dx}, si, sj, 1.0, f, fsi, fsj );
        E2=addDensOverlapGauss_S( {0,0,x+dx}, si, sj, 1.0, f, fsi, fsj );
        return (E2-E1)/(2*dx);
    } );
    */

    /*
    // --- derivative along si
    plot.xsharingLines( 3, np, 0.25, 0.05 );
    DataLine2D *l;
    double sj = 1.0;
    l=plot.lines[0]; l->clr=0xFFFFFFFF; l->label="EDens";  evalLine( *l, [&](double x){ double fsi,fsj,E; Vec3d f=Vec3dZero; E=addDensOverlapGauss_S( {0,0,1}, x, sj, 1.0, f, fsi, fsj ); return E; } );
    l=plot.lines[1]; l->clr=0xFF0000FF; l->label="F1";     evalLine( *l, [&](double x){ double fsi,fsj,E; Vec3d f=Vec3dZero; E=addDensOverlapGauss_S( {0,0,1}, x, sj, 1.0, f, fsi, fsj ); return fsi; } );
    l=plot.lines[2]; l->clr=0xFF0080FF; l->label="FeeNum"; evalLine( *l, [&](double x){ double fsi,fsj,E1,E2,dx=0.001; Vec3d f=Vec3dZero;
        E1=addDensOverlapGauss_S( {0,0,1}, x-dx, sj, 1.0, f, fsi, fsj );
        E2=addDensOverlapGauss_S( {0,0,1}, x+dx, sj, 1.0, f, fsi, fsj );
        return (E2-E1)/(2*dx);
    } );
    */
    //l=plot.lines[0]; l->clr=0xFF0000FF; l->label="F0";    evalLine( *l, [&](double x){ return 0; } );
    //l=plot.lines[2]; l->clr=0xFF0000FF; l->label="F2";    evalLine( *l, [&](double x){ return 0; } );

    for(int i=0;i<np;i++){
        printf( " %i %g -> %g %g %g \n", i, l->xs[i],  plot.lines[0]->ys[i], plot.lines[1]->ys[i], plot.lines[2]->ys[i] );
    }

    /*
    plot.xsharingLines( 3, np, 0.0, 0.05 );
    DataLine2D *l;
    l=plot.lines[0]; l->clr=0xFFFF0000; l->label="Eee";    evalLine( *l, [&](double x){ double fs,fr,E; E=CoulombGauss( x, 2.0, fr, fs, 1.0 ); return E;    } );
    l=plot.lines[1]; l->clr=0xFF0000FF; l->label="Fee";    evalLine( *l, [&](double x){ double fs,fr,E; E=CoulombGauss( x, 2.0, fr, fs, 1.0 ); return fr*x; } );
    l=plot.lines[2]; l->clr=0xFFFF00FF; l->label="FeeNum"; evalLine( *l, [&](double x){ double fs,fr,E1,E2,s=2.0,dx=0.001;
        E1=CoulombGauss( x-dx, s, fr, fs, 1.0 );
        E2=CoulombGauss( x+dx, s, fr, fs, 1.0 );
        return (E2-E1)/(2*dx);
    } );
    */



    //l=plot.lines[2]; l->clr=0xFFFF00FF; l->label="Eae";  evalLine( *l, [&](double x){ Vec3d f;  return addPairEF_expQ( {x,0,0}, f, eAbw.z, QQae, eAbw.y,  eAbw.x     ); } );
    //l=plot.lines[3]; l->clr=0xFF0000FF; l->label="Eaa";  evalLine( *l, [&](double x){ Vec3d f;  return addPairEF_expQ( {x,0,0}, f, aAbw.z, QQaa, aAbw.y,  aAbw.x     ); } );
    //l=plot.lines[3]; l->clr=0xFF0080FF; l->label="Faa"; evalLine( *l, [&](double x){ Vec3d f;  return addPairEF_expQ( {x,0,0}, f, w2aa, qaa, 0,      0           ); } );
    //l=plot.lines[2]; l->clr=0xFFFF8000; l->label="Fee"; evalLine( *l, [&](double x){ Vec3d f=Vec3dZero;  addPairEF_expQ( {x,0,0}, f, w2ee, +1.0, bEE, aEE  ); return -f.x; } );
    //l=plot.lines[3]; l->clr=0xFFFF80FF; l->label="Fae"; evalLine( *l, [&](double x){ Vec3d f=Vec3dZero;  addPairEF_expQ( {x,0,0}, f, w2ae, qae,  bAE, aAE  ); return -f.x; } );
    //l=plot.lines[2]; l->clr=0xFF0000FF; l->label="Eaa"; evalLine( *l, [&](double x){ Vec3d f=Vec3dZero;  addPairEF_expQ( {x,0,0}, f, w2aa, qaa, 0,      0           ); return f.x; } );


    /*

    plot.xsharingLines( 7, np, 0.0, 0.05 );
    double fc = 0.2;

    plot.lines[0]->clr = 0xFFFFFFFF; // Etot
    plot.lines[1]->clr = 0xFF0000FF; // Eaa
    plot.lines[2]->clr = 0xFFFF00FF; // Eae
    plot.lines[3]->clr = 0xFFFF0000; // Eee
    plot.lines[4]->clr = 0xFF00FF00; // Ek
    plot.lines[5]->clr = 0xFFFFFF00; // EeePaul
    plot.lines[6]->clr = 0xFF8F008F; // Eae*-0.25

    ff.apos[0]= Vec3dZero;
    ff.apos[1]= Vec3dZero;
    ff.epos[0]= Vec3dZero;
    ff.epos[1]= Vec3dZero;

    double sc=0.2;

    for(int i=0; i<np; i++){
        double x = i*0.1;

        ff.apos[0].x = -x;
        ff.apos[1].x =  x;
        ff.epos[0].x = -x;
        ff.epos[1].x =  x;

        double Etot = ff.eval();
        plot.lines[0]->ys[i] = sc*Etot;
        plot.lines[1]->ys[i] = sc*ff.Eaa;
        plot.lines[2]->ys[i] = sc*ff.Eae;
        plot.lines[6]->ys[i] = sc*ff.Eae*-0.5;
        plot.lines[3]->ys[i] = sc*ff.Eee;
        plot.lines[4]->ys[i] = sc*ff.Ek;
        plot.lines[5]->ys[i] = sc*ff.EeePaul;
        printf( "makePlots[%i] %g -> %g (%g,%g(4*%g),%g) %g,%g \n", i, x,   Etot, ff.Eaa,ff.Eae,ff.Eae*0.5,ff.Eee, ff.Ek, ff.EeePaul );
    }
    */

    plot.update();
    plot.autoAxes(0.5,0.5);
    printf( "axBound %g,%g %g,%g \n", plot.axBounds.a.x, plot.axBounds.a.y, plot.axBounds.b.x, plot.axBounds.b.y );
    plot.render();

}

int makePlots2( Plot2D& plot ){

    EFF ff;
    ff.realloc(1,1);

    int np = 60;
    plot.xsharingLines( 3, np, 0.0, 0.05 );
    double fc = 0.2;

    plot.lines[0]->clr = 0xFFFF0000; // Etot
    plot.lines[1]->clr = 0xFF00FF00; // F-pos
    plot.lines[2]->clr = 0xFF0000FF; // F-size

    ff.aQ   [0]= 4;
    ff.apos [0]= Vec3dZero;
    ff.epos [0]= Vec3dZero;
    ff.esize[0]=1.0;

    ff.autoAbWs( default_aAbWs, default_eAbWs );

    double sc=0.05;
    double Etot;
    for(int i=0; i<np; i++){
        double x = i*0.1;
        //ff.apos[0].x = -x;
        ff.epos[0].x = -x;

        ff.esize[0]=0.5; plot.lines[0]->ys[i] = sc* ff.eval();
        ff.esize[0]=1.0; plot.lines[1]->ys[i] = sc* ff.eval();
        ff.esize[0]=1.5; plot.lines[2]->ys[i] = sc* ff.eval();

        /*
        ff.clearForce();
        Etot = ff.eval();
        plot.lines[0]->ys[i] = sc*Etot;
        plot.lines[1]->ys[i] = sc*ff.eforce[0].x;
        plot.lines[2]->ys[i] = sc*ff.fsize [0];
        */
        //printf( "makePlots2[%i] %g -> E %g fe %g fize %g \n", i, x, Etot, ff.eforce[0].x, ff.fsize[0] );
        //printf( "makePlots[%i] %g -> %g (%g,%g(4*%g),%g) %g,%g \n", i, x,   Etot, ff.Eaa,ff.Eae,ff.Eae*0.5,ff.Eee, ff.Ek, ff.EeePaul );
    }


    plot.update();
    plot.autoAxes(0.5,0.5);
    printf( "axBound %g,%g %g,%g \n", plot.axBounds.a.x, plot.axBounds.a.y, plot.axBounds.b.x, plot.axBounds.b.y );
    plot.render();

    ff.dealloc();

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

    bool bMapElectron = false;
    int ipicked  = 0;

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


    void init2DMap( int n, double dx );

};

void TestAppRARFF::init2DMap( int n, double dx ){

// ===== EVAL FF 2D MAP

    //Vec3d  *ps=0,*fs=0;
    //double *Es=0;

    ipicked = 1;

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

    Vec3d pa0;// = ff.apos[ipicked];
    if(bMapElectron){ pa0 = ff.epos[ipicked]; }else{ pa0 = ff.apos[ipicked]; }
    field_ns = {n,n};
    sampleScalarField( Efunc, field_ns, {-(dx*n)/2,-(dx*n)/2,+0.1}, {dx,0.0,0.0}, {0.0,dx,0.0}, field_ps, field_Es, Erange );
    if(bMapElectron){ ff.epos[ipicked] = pa0; }else{ ff.apos[ipicked] = pa0; }

    E0 = field_Es[0];
    printf( "val_range: %g %g %g \n", Erange.x, Erange.y, field_Es[0] );
    Espread = 3.0;

    ogl_fs = genFieldMap( ogl_fs, field_ns, field_ps, field_Es, E0-Espread, E0+Espread );

}

TestAppRARFF::TestAppRARFF( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );
    plot1.fontTex=fontTex;

    checkDerivs( ff.KRSrho );   // exit(0);
    //makePlots( plot1, ff );  // return; //      exit(0);
    makePlots2( plot1 ); //return;


    // ===== SETUP GEOM
    //char* fname = "data/H_eFF.xyz";
    //char* fname = "data/H2_eFF_spin.xyz";
    //char* fname = "data/Ce1_eFF.xyz";
    //char* fname = "data/Ce2_eFF.xyz";
    //char* fname = "data/Ce4_eFF.xyz";
    //char* fname = "data/CH3_eFF_spin.xyz";
    //char* fname = "data/CH4_eFF_flat_spin.xyz";
    //char* fname = "data/CH4_eFF_spin.xyz";
    //char* fname = "data/C2_eFF_spin.xyz";
    //char* fname = "data/C2H4_eFF_spin.xyz";
    char* fname = "data/C2H4_eFF_spin_.xyz";
    //char* fname = "data/C2H6_eFF_spin.xyz";
    //ff.loadFromFile_xyz( "data/C2H4_eFF_spin.xyz" );
    ff.loadFromFile_xyz( fname );

    //setGeom(ff);
    //double sz = 0.51;
    // break symmetry
    //for(int i=0; i<ff.na; i++){ ff.apos[i].add( randf(-sz,sz),randf(-sz,sz),randf(-sz,sz) );  }
    for(int i=0; i<ff.na; i++){ printf( "A_pos[%i] (%g,%g,%g)\n", i, ff.apos[i].x, ff.apos[i].y, ff.apos[i].z ); }
    for(int i=0; i<ff.ne; i++){ printf( "e_pos[%i] (%g,%g,%g)\n", i, ff.epos[i].x, ff.epos[i].y, ff.epos[i].z ); }

    ff.autoAbWs( default_aAbWs, default_eAbWs );
    //exit(0);

    //VecN::set(ff.ne,4.0,ff.esize);

    // ==== Test Eval

    i_DEBUG = 1;
    ff.eval();
    printf( "Ek %g Eaa %g Eae %g Eee %g EeePaul %g \n", ff.Ek, ff.Eaa, ff.Eae, ff.Eee, ff.EeePaul );
    i_DEBUG = 0;
    //exit(0);

    //makePlots( plot1, ff );
    //ff.loadFromFile_xyz( fname  );

    // ==== Bind Optimizer

    //init2DMap( 100, 0.1 );

    opt.bindOrAlloc( ff.nDOFs, ff.pDOFs, 0, ff.fDOFs, 0 );
    opt.cleanVel( );
    opt.initOpt( 0.005, 0.1 );
    opt.f_limit = 1000.0;

    ff.eval();

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

    //return;

    if(bRun){
        for(int itr=0;itr<perFrame;itr++){
            printf( " ==== frame %i i_DEBUG  %i \n", frameCount, i_DEBUG );
            double F2 = 1.0;

            ff.clearForce();
            //applyCartesianBoxForce( {0.0,0.0,0.0}, {0.0,0.0,0.0}, {0,0,50.0}, ff.na, ff.apos, ff.aforce );
            //applyCartesianBoxForce( {0.0,0.0,0.0}, {0.0,0.0,0.0}, {0,0,5.0}, ff.ne, ff.epos, ff.eforce );
            //ff.evalEE();
            //ff.evalAE();
            //ff.evalAA();
            ff.eval();
            //ff.apos[0].set(.0);


            //VecN::set( ff.na*3, 0.0, (double*)ff.aforce ); // FIX ATOMS
            VecN::set( ff.ne, 0.0, ff.fsize ); // FIX ELECTRON SIZE
            //if(bRun)ff.move_GD(0.001 );

            //ff.move_GD( 0.0001 );

            F2 = opt.move_FIRE();

            printf( " |F| %g \n", sqrt(F2) );
            //if(!(F2<1000000.0))perFrame=0;
        }
    }

    //printf( "e[0] r %g s %g \n", ff.epos[0].norm(), ff.esize[0] );

    //printf( "r07 r %g s %g \n", ff.epos[0].norm(), ff.esize[0] );

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

    //printf( "apos (%g,%g,%g) \n", ff.apos[0].x, ff.apos[0].y, ff.apos[0].z );

    char strtmp[256];
    double Qsz = 0.05;
    double fsc = 1.0;
    glColor3f(0.0,0.0,0.0);
    for(int i=0; i<ff.na; i++){
        //printf( "apos[%i] (%g,%g,%g)\n", i, ff.apos[i].x, ff.apos[i].y, ff.apos[i].z );
        Draw3D::drawPointCross( ff.apos  [i]    , ff.aQ  [i]*Qsz );
        Draw3D::drawVecInPos(   ff.aforce[i]*fsc, ff.apos[i] );
        //printf( " %i %f %f %f %f  \n", i, ff.aQ[i], ff.apos[i].x,ff.apos[i].y,ff.apos[i].z );
        //printf( " %i %f %f %f %f  \n", i, ff.aQ[i], ff.aforce[i].x, ff.aforce[i].y, ff.aforce[i].z );
        sprintf(strtmp,"%i",i);
        Draw3D::drawText(strtmp, ff.apos[i], fontTex, 0.02, 0);
    }

    //glColor3f(1.0,1.0,1.0);
    for(int i=0; i<ff.ne; i++){
        //printf( "epos[%i] (%g,%g,%g)\n", i, ff.epos[i].x, ff.epos[i].y, ff.epos[i].z );
        if(ff.espin[i]>0){ glColor3f(0.0,0.5,1.0); }else{ glColor3f(1.0,0.5,0.0); };
        Draw3D::drawPointCross( ff.epos  [i], ff.esize[i] );
        Draw3D::drawVecInPos(   ff.eforce[i]*fsc, ff.epos[i] );
        Draw3D::drawSphereOctLines( 8, ff.esize[i], ff.epos[i] );
        sprintf(strtmp,"%i",i);
        Draw3D::drawText(strtmp, ff.epos[i], fontTex, 0.02, 0);
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

                case SDLK_i: printFFInfo(ff);
                case SDLK_LEFTBRACKET :  Espread *= 1.2; ogl_fs=genFieldMap(ogl_fs, field_ns, field_ps, field_Es, E0-Espread, E0+Espread ); break;
                case SDLK_RIGHTBRACKET:  Espread /= 1.2; ogl_fs=genFieldMap(ogl_fs, field_ns, field_ps, field_Es, E0-Espread, E0+Espread ); break;
                case SDLK_e: bMapElectron=!bMapElectron;
                case SDLK_f:
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
















