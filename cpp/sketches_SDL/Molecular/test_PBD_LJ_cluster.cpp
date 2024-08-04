
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw3D.h"
//#include "Solids.h"
#include "raytrace.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "VecN.h"


#include "AppSDL2OGL_3D.h"
#include "testUtils.h"
#include "SDL_utils.h"
#include "Plot2D.h"

#include "Draw3D_Molecular.h"

#include "NBFF_SR.h"

#define R2SAFE  1.0e-8f

// ============= This is From LIB ReactiveFF.cpp
// /home/prokop/git/ProbeParticleModel_OclPy3/cpp/ReactiveFF.cpp


// ============= Global Variables


NBFF_SR ff;

// ============= Functions


class TestAppRARFF: public AppSDL2OGL_3D { public:

    bool bBlockAddAtom=false;
    int ipicked = -1;
    Vec3d ray0;
    Vec3d mouse_p0;
    int ogl_sph=0;

    //const char* workFileName="data/work.rff";

    Plot2D plot1;
    bool bRun = true;
    int  perFrame = 100;

    int      fontTex;

    virtual void draw   ();
    virtual void drawHUD();
    //virtual void mouseHandling( );
    virtual void eventHandling   ( const SDL_Event& event  );
    virtual void keyStateHandling( const Uint8 *keys );

    TestAppRARFF( int& id, int WIDTH_, int HEIGHT_ );

    void simulation();
    void makePotentialPlot();
    void view_BBoxes();
    void view_atoms();
    void view_bonds();
    void make_atoms( int nCluster, int nPerCluster, double xspan, double size );

};

TestAppRARFF::TestAppRARFF( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {
    fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );

    make_atoms( 4, 8, 5.0, 2.0 );

    makePotentialPlot();
    
    Draw3D::makeSphereOgl( ogl_sph, 3, 0.25 );
}

void run(int niter){
    for(int i=0; i<niter; i++){
        //ff.eval();
        //ff.solve_constr();
        //ff.move();
        ff.step( 0.5 );
    }
    
}

void TestAppRARFF::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);
    glDisable( GL_LIGHTING );

    perFrame = 10;
    //ff.bGridAccel=false;
    if(bRun){
        run(1);
        //run(3); exit(0);
        //run(perFrame);
    }

    glEnable( GL_LIGHTING );
    glColor3f( 1.0, 1.0, 1.0 );
    view_atoms();

    glDisable( GL_LIGHTING );
    glColor3f( 0.0, 0.0, 0.0 );
    view_bonds();
    //view_BBoxes();
    ray0 = (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y);
    glColor3f( 0.0, 0.0, 0.0 );
    Draw3D::drawPointCross( ray0, 0.1 );
    if(ff.ipicked>=0){
        ff.pick_ray0 = ray0;
        ff.pick_hray = (Vec3d)(cam.rot.c);
        Draw3D::drawLine( ff.apos[ff.ipicked], ray0);
    }
    Draw3D::drawAxis( 1.0);
};

void TestAppRARFF::drawHUD(){
	glTranslatef( 400.0,400.0,0.0 );
	glScalef    ( 40.0,40.0,1.0  );
	plot1.view();
}

void TestAppRARFF::keyStateHandling( const Uint8 *keys ){
    if( keys[ SDL_SCANCODE_R ] ){
        if(ipicked>=0){
            //Mat3d rot;
            //Quat4d qrot;  qrot.f=(Vec3d)cam.rot.c*0.01; qrot.normalizeW();
            //ff.qrots[ipicked].qmul(qrot);
            //ff.projectAtomBons(ipicked);
        }
    }
	AppSDL2OGL_3D::keyStateHandling( keys );
};

void TestAppRARFF::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );

    //Quat4i* caps = &capsBrush;
    //if(ipicked>=0) caps=((Quat4i*)ff.bondCaps)+ipicked;

    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;

                //case SDLK_l: ff.load(workFileName); bRun=false; break;
                //case SDLK_k: ff.save(workFileName);  break;
                //case SDLK_x: ff.saveXYZ("data/RARFF_out.xyz");  break;
                case SDLK_SPACE: bRun = !bRun;  break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;

        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:{
                    ff.ipicked = pickParticle( ray0, (Vec3d)cam.rot.c, 0.5, ff.natom, ff.apos );
                        printf( "LMB DOWN picked %i \n", ff.ipicked );
                    }break;
                case SDL_BUTTON_RIGHT:{
                    /*
                    ipicked = pickParticle( ray0, (Vec3d)cam.rot.c, 0.5, ff.natom, ff.apos, ff.ignoreAtoms );
                    printf( "remove atom %i \n", ipicked );
                    ff.ignoreAtoms[ ipicked ] = true;
                    */
                    }break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    // printf( "LMB UP picked %i bblock %i \n", ipicked, bBlockAddAtom );
                    // if( (ipicked==-1)&&(!bBlockAddAtom) ){
                    //     Vec3d mouse_p = (Vec3d)( cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y );
                    //     Vec3d dir = mouse_p-mouse_p0;
                       
                    // }
                    break;
                case SDL_BUTTON_RIGHT:
                    ff.ipicked = -1;
                    break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}

void TestAppRARFF::makePotentialPlot(){
    int np = 100;
    DataLine2D* L = 0; 
    DataLine2D* line_Er=new DataLine2D(np); L=line_Er; plot1.add(L); L->clr = 0xFF0000FF; 
    DataLine2D* line_Fr=new DataLine2D(np); L=line_Fr; plot1.add(L); L->clr = 0xFFFF0000;  // L->replace_xs( line_Er->xs );

    DataLine2D* line_Er2=new DataLine2D(np); L=line_Er2; plot1.add(L); L->clr = 0xFFFF00FF; 
    DataLine2D* line_Fr2=new DataLine2D(np); L=line_Fr2; plot1.add(L); L->clr = 0xFFFF80FF;  // L->replace_xs( line_Er->xs );

    for(int i=0;i<100;i++){
        double x = i*0.1;
        Vec3d f;  double E  = getSR( Vec3d{x,0.0,0.0}, f, ff.Rcut, ff.R0_nb, ff.E0_nb, ff.K_nb  );
        Vec3d f2; double E2 = getSR2( Vec3d{x,0.0,0.0}, f2, ff.Rcut, ff.R0_nb, ff.E0_nb, ff.K_nb, ff.Rf  );
        //printf( "[%i] r=%g R=%g fr=%g \n", i, x, E, f.x );
        line_Er->xs[i] = x; line_Er->ys[i] = E;
        line_Fr->xs[i] = x; line_Fr->ys[i] = f.x;

        line_Er2->xs[i] = x; line_Er2->ys[i] = E2;
        line_Fr2->xs[i] = x; line_Fr2->ys[i] = f2.x;

    }
    plot1.init();
    plot1.scaling.y = 10.0;
    plot1.fontTex = fontTex;
    plot1.clrGrid = 0xFF404040;
    plot1.bRenderAxes = false;
    plot1.render();
}

void TestAppRARFF::view_BBoxes(){
    //printf( "view_BBoxes() ff.pointBBs.ncell=%i \n", ff.pointBBs.ncell );
    for(int ic=0;ic<ff.pointBBs.ncell;ic++){
        int i0=ff.pointBBs.cellI0s[ic];
        int ni=ff.pointBBs.cellNs [ic];
        //printf( "ic %i io %i ni %i \n", ic, i0, ni );
        Vec3i ip;
        Draw  ::color_of_hash( 464+645*ic );
        //ff.pointBBs.i2ixyz( ic, ip );
        Draw3D::drawBBox( ff.BBs[ic].lo,  ff.BBs[ic].hi );
        //printf( "BBox[%i] p(%g,%g,%g) p(%g,%g,%g)\n", ic, ff.BBs[ic].lo.x, ff.BBs[ic].lo.y, ff.BBs[ic].lo.z, ff.BBs[ic].hi.x, ff.BBs[ic].hi.y, ff.BBs[ic].hi.z );
        
        // for(int j=i0; j<i0+ni;j++){
        //     //printf( "view_BBoxes() ic=%i j=%i i0=%i ni=%i \n", ic, j, i0, ni );
        //     int   io = ff.pointBBs.cell2obj[j];
        //     Vec3d p  = ff.apos[io];
        //     //printf( "j %i io %i p(%g,%g,%g) \n", j, io, p.x,p.y,p.z );
        //     //Draw  ::color_of_hash( 464+645*ic );
        //     Draw3D::drawPointCross( p, 0.2 );            
        // }
    }
}

void TestAppRARFF::view_atoms(){
    //glColor3f(0.0,0.0,0.0);
    double fsc = 0.1;
    double tsc = 0.1;
    //printf( "ff.natom %i \n", ff.natom );
    for(int ia=0; ia<ff.natom; ia++){
        glColor3f(1.0,1.0,1.0);
        Draw3D::drawShape( ogl_sph, ff.apos[ia], Mat3dIdentity );
        glColor3f(1.0,0.0,0.0);
        //Draw3D::drawVecInPos(  ff.aforce[ia]*100.0, ff.apos[ia] );
    };
}

void TestAppRARFF::view_bonds(){
    for(int ib=0; ib<ff.nBBs; ib++){
        int i0 = ff.pointBBs.cellI0s[ib];
        int ni = ff.pointBBs.cellNs [ib];
        int iao =  i0 + ni - 1;
        glBegin( GL_LINES );
        for(int ia=i0; ia<i0+ni; ia++){
            Draw3D::drawLine( ff.apos[ia], ff.apos[iao] );
            iao = ia;
        }
        glEnd();
    }
}

void TestAppRARFF::make_atoms( int nCluster, int nPerCluster, double xspan, double size  ){
    double vsize = 0.0;
    ff.realloc( nCluster*nPerCluster, nCluster );
    for(int ib=0; ib<nCluster; ib++){
        Vec3d c; c.fromRandomBox( Vec3d{-xspan,-xspan,-xspan} ,Vec3d{xspan,xspan,xspan});
        //const int  npi = ff.pointBBs.cellNs[ib];
        //const int  ip0 = ff.pointBBs.cellI0s[ib];
        ff.pointBBs.cellI0s[ib] = nPerCluster * ib;
        ff.pointBBs.cellNs [ib] = nPerCluster;
        Vec3d vcog=Vec3dZero;
        for(int i=0; i<nPerCluster; i++){
            int ip = ib*nPerCluster+i;
            Vec3d d; d.fromRandomBox( Vec3d{-size,-size,-size} ,Vec3d{size,size,size});
            Vec3d v; v.fromRandomBox( Vec3d{-vsize,-vsize,-vsize} ,Vec3d{vsize,vsize,vsize});
            ff.apos [ip] = c + d;
            //ff.avel [ip] = v;
            //vcog.add( v );
            ff.pointBBs.cell2obj[ff.pointBBs.cellI0s[ib]+i] =ip;
        }
        //vcog.mul( 1./nPerCluster );
        //for(int i=0; i<nPerCluster; i++){ ff.avel [ib*nPerCluster+i].sub( vcog ); }
    }
}

// ===================== MAIN

TestAppRARFF* thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppRARFF( junk , 1400, 1000 );
	thisApp->loop( 1000000 );
	return 0;
}



