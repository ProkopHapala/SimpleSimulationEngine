
#ifndef eFF_plots_h
#define eFF_plots_h


#include "MMFFBuilder.h"

// ===========================================
// ============ Checking
// ===========================================

Vec3d v3sum(int n, Vec3d* fs){
    Vec3d f=Vec3dZero;
    for(int i=0; i<n; i++){ f.add(fs[i]); }
    return f;
}

Vec3d v3min(int n, Vec3d* fs){
    Vec3d f={1e+300,1e+300,1e+300};
    for(int i=0; i<n; i++){ f.setIfLower( fs[i] ); }
    return f;
}

Vec3d v3max(int n, Vec3d* fs){
    Vec3d f={-1e+300,-1e+300,-1e+300};
    for(int i=0; i<n; i++){ f.setIfGreater( fs[i] ); }
    return f;
}

bool checkScalar(const char* s, double v, double vmin, double vmax){
    bool b = false;
    if( v<vmin ){
        printf("%s %g<min(%g)\n", s, v, vmin);
        b = true;
    }
    if( v<vmax ){
        printf("%s %g<min(%g)\n", s, v, vmin);
        b = true;
    }
    return b;
}

bool checkFinite(const EFF& ff, double vmin, double vmax ){
    bool bErr = false;
    Vec3d pmin,pmax,  fmin,fmax, psum,fsum;
    pmin = v3min( ff.ne, ff.epos );  pmax = v3max( ff.ne, ff.epos );
    fmin = v3min( ff.ne, ff.eforce); fmax = v3max( ff.ne, ff.eforce );
    fsum = v3sum( ff.ne, ff.eforce); Vec3d ftot = fsum;
    psum = v3sum( ff.ne, ff.epos  ); Vec3d ptot = psum;
    //printf("pe min(%g,%g,%g) max(%g,%g,%g) | fe min(%g,%g,%g) max(%g,%g,%g) \n", pmin.x,pmin.y,pmin.z,   pmax.x,pmax.y,pmax.z,     fmin.x,fmin.y,fmin.z, fmax.x,fmax.y,fmax.z );
    //printf("fe sum(%g,%g,%g) min(%g,%g,%g) max(%g,%g,%g) \n", fsum.x,fsum.y,fsum.z,     fmin.x,fmin.y,fmin.z, fmax.x,fmax.y,fmax.z );
    printf("fe sum(%g,%g,%g) pe sum(%g,%g,%g) \n", fsum.x,fsum.y,fsum.z,  psum.x,psum.y,psum.z   );

    pmin = v3min( ff.na, ff.apos   ); pmax = v3max( ff.na, ff.apos   );
    fmin = v3min( ff.na, ff.aforce ); fmax = v3max( ff.na, ff.aforce );
    fsum = v3sum( ff.na, ff.aforce ); ptot.add(psum);
    psum = v3sum( ff.na, ff.apos   ); ftot.add(fsum);
    //printf("pa min(%g,%g,%g) max (%g,%g,%g) | fa min(%g,%g,%g) max(%g,%g,%g) \n", pmin.x,pmin.y,pmin.z,   pmax.x,pmax.y,pmax.z,    fmin.x,fmin.y,fmin.z, fmax.x,fmax.y,fmax.z );
    //printf("fa sum(%g,%g,%g) min(%g,%g,%g) max(%g,%g,%g) \n", fsum.x,fsum.y,fsum.z,     fmin.x,fmin.y,fmin.z, fmax.x,fmax.y,fmax.z );
    printf("fe sum(%g,%g,%g) pe sum(%g,%g,%g) \n", fsum.x,fsum.y,fsum.z,  psum.x,psum.y,psum.z   );
    printf( "ftot (%g,%g,%g) ptot (%g,%g,%g)\n", ftot.x, ftot.y, ftot.z,     ptot.x, ptot.y, ptot.z );

    //fmin = v3min( ff.ne, ff.epos ); pmax = v3max( ff.ne, ff.epos ); printf("pe min(%g,%g,%g) max (%g,%g,%g)\n", pmin.x,pmin.y,pmin.z,   pmax.x,pmax.y,pmax.z );
    //fmin = v3min( ff.na, ff.apos ); pmax = v3max( ff.na, ff.apos ); printf("pa min(%g,%g,%g) max (%g,%g,%g)\n", pmin.x,pmin.y,pmin.z,   pmax.x,pmax.y,pmax.z );

    //pmin = v3min( ff.ne, ff.epos ); pmax = v3max( ff.ne, ff.epos ); printf("pe min(%g,%g,%g) max (%g,%g,%g)\n", pmin.x,pmin.y,pmin.z,   pmax.x,pmax.y,pmax.z );

    //Vec3d  pe = v3sum( ff.ne, ff.epos );      bErr &= checkScalar( "fe", pe.norm(), vmin, vmax);
    //Vec3d  pa = v3sum( ff.na, ff.apos );      bErr &= checkScalar( "fa", pa.norm(), vmin, vmax);
    //double ps = VecN::sum( ff.ne, ff.esize ); bErr &= checkScalar( "fs", ps       , vmin, vmax);

    //Vec3d  fe = v3sum( ff.ne, ff.eforce );    bErr &= checkScalar( "fe", fe.norm(), vmin, vmax);
    //Vec3d  fa = v3sum( ff.na, ff.aforce );    bErr &= checkScalar( "fa", fa.norm(), vmin, vmax);
    //double fs = VecN::sum( ff.ne, ff.fsize );   bErr &= checkScalar( "fs", fs       , vmin, vmax);
    return bErr;
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

void checkDerivs2(){

    double d = 0.0001;
    double E1,E2;

    EFF ff;
    ff.realloc(2,2);
    ff.aPars [0].x= 4;
    ff.aPars [1].x= 4;
    //ff.autoAbWs( default_aAbWs, default_eAbWs );

    ff.apos [0]= {+0.7,+0.2,0.0};
    ff.apos [1]= {-0.7,-0.2,0.0};
    ff.epos [0]= {+0.1,+0.6,0.0};
    ff.epos [1]= {-0.1,-0.6,0.0};
    ff.esize[0]=0.5;
    ff.esize[1]=2.0;

    // a1_r
    ff.apos [0].x-=d; E1=ff.eval();    ff.apos[0].x+=2*d; E2=ff.eval();     ff.apos [0].x-=d; ff.clearForce(); ff.eval(); printf("f/fE  apos[0]  %g / %g | %g \n", ff.aforce[0].x, -(E2-E1)/(2*d), ff.aforce[0].y );
    //ff.apos[1].x-=d; E1=ff.eval();    ff.apos[1].x+=2*d; E2=ff.eval();    ff.apos[1].x-=d; printf("dar1 %g / %g \n",    ff.aforce[1].x, (E2-E1)/(2*d) );
    ff.epos [0].x-=d; E1=ff.eval();    ff.epos [0].x+=2*d; E2=ff.eval();    ff.epos [0].x-=d; ff.clearForce(); ff.eval(); printf("f/fE  epos[0] %g / %g | %g \n", ff.eforce[0].x, -(E2-E1)/(2*d), ff.eforce[0].y );
    ff.epos [1].x-=d; E1=ff.eval();    ff.epos [1].x+=2*d; E2=ff.eval();    ff.epos [1].x-=d; ff.clearForce(); ff.eval(); printf("f/fE  epos[1] %g / %g | %g \n", ff.eforce[1].x, -(E2-E1)/(2*d), ff.eforce[1].y );
    ff.esize[0]  -=d; E1=ff.eval();    ff.esize[0]  +=2*d; E2=ff.eval();    ff.esize[0]  -=d; ff.clearForce(); ff.eval(); printf("f/fE esize[0] %g / %g \n",      ff.fsize[0],    -(E2-E1)/(2*d) );
    ff.esize[1]  -=d; E1=ff.eval();    ff.esize[1]  +=2*d; E2=ff.eval();    ff.esize[1]  -=d; ff.clearForce(); ff.eval(); printf("f/fE esize[1] %g / %g \n",      ff.fsize[1],    -(E2-E1)/(2*d) );
    //exit(0);

}

// ===========================================
// ============ Plotting
// ===========================================

void plotAtomsPot( EFF& ff, DataLine2D *line, Vec3d p0, Vec3d dp, float sc=1.0, double s=0.0 ){
    Vec3d ps[line->n];
    for(int i=0; i<line->n; i++){  ps[i]=p0+dp*line->xs[i]; }
    //solver.orbAtPoints( io, line->n, ps, line->ys );
    ff.atomsPotAtPoints( line->n, ps, line->ys, s, 1.0 );
    for(int i=0; i<line->n; i++){  line->ys[i]*=sc; }
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



void makePlots( Plot2D& plot, EFF& ff ){

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

void makePlots2( Plot2D& plot ){

    EFF ff;
    ff.realloc(1,1);

    int np = 60;
    plot.xsharingLines( 3, np, 0.0, 0.05 );
    double fc = 0.2;

    plot.lines[0]->clr = 0xFFFF0000; // Etot
    plot.lines[1]->clr = 0xFF00FF00; // F-pos
    plot.lines[2]->clr = 0xFF0000FF; // F-size

    //ff.aQ   [0]= 4;
    ff.aPars[0].x= 4;
    ff.apos [0]= Vec3dZero;
    ff.epos [0]= Vec3dZero;
    ff.esize[0]=1.0;

    //ff.autoAbWs( default_aAbWs, default_eAbWs );

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

void gen_molecule( EFF& ff ){

    MM::Builder builder;

    // ================== Generate Atomic

    //const int natom=4,nbond=3,nang=2,ntors=1;
    const int natom=4,nbond=3;
    Vec3d apos0[] = {
        {-2.0,0.0,0.0},  // 0
        {-1.0,2.0,0.0},  // 1
        {+1.0,2.0,0.0},  // 2
        {+2.0,0.0,1.0},  // 3
    };
    Vec2i bong2atom[] = {
        {0,1},  // 0
        {1,2},  // 1
        {2,3},  // 2
    };

    // ============ Build molecule

    MM::Atom brushAtom{  6, -1,-1, Vec3dZero, MM::Atom::defaultREQ };
    MM::Bond brushBond{ -1, {-1,-1}, 1.5, 25.0 };

    builder.capBond = MM::Bond{ -1, {-1,-1}, 1.07, 15.0 };
    builder.capAtom = { 1, -1,-1, Vec3dZero, MM::Atom::defaultREQ };

    printf( "----- Atoms \n" );
    for(int i=0;i<natom;i++){
        brushAtom.pos = apos0[i];
        //builder.insertAtom(brushAtom, true);
        builder.insertAtom(brushAtom, 0);
    }
    printf( "----- Bonds \n" );
    for(int i=0;i<nbond;i++){
        brushBond.atoms=bong2atom[i];
        builder.insertBond(brushBond);
    }
    printf( "----- Confs \n" );
    for(int i=0;i<natom;i++){
        builder.makeSPConf(i,0,0);
        //builder.makeSPConf(i);
    }
    //builder.toMMFFmini( ff );
    builder.toEFF( ff, EFFparams, 0.5, 0.025 );


}

#endif
