
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
    void make_atoms( int natom, double xspan, double size );

};

TestAppRARFF::TestAppRARFF( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {
    fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );



    Draw3D::makeSphereOgl( ogl_sph, 3, 0.25 );

}

void run(int niter){
    for(int i=0; i<niter; i++){
        //ff.eval();
        //ff.solve_constr();
        //ff.move();
        ff.step( 0.1 );
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
        //run(perFrame);
    }else{

    }
    view_atoms();
    ray0 = (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y);
    Draw3D::drawPointCross( ray0, 0.1 );
    if(ipicked>=0) Draw3D::drawLine( ff.apos[ipicked], ray0);
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
                    //int ip = pickParticle( ray0, (Vec3d)cam.rot.c, 0.5, ff.natom, ff.apos, ff.ignoreAtoms );
                    // bBlockAddAtom=false;
                    // if( ip>=0){
                    //     if(ipicked==ip){ ipicked=-1; bBlockAddAtom=true; printf("inv\n"); }
                    //     else           { ipicked=ip;                     printf("set\n"); }
                    // }else{ ipicked=-1; }
                    // mouse_p0 = (Vec3d)( cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y );
                    // printf( "LMB DOWN picked %i/%i bblock %i \n", ipicked,ip, bBlockAddAtom );
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
                    printf( "LMB UP picked %i bblock %i \n", ipicked, bBlockAddAtom );
                    if( (ipicked==-1)&&(!bBlockAddAtom) ){
                        Vec3d mouse_p = (Vec3d)( cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y );
                        Vec3d dir = mouse_p-mouse_p0;
                       
                    }
                    break;
                case SDL_BUTTON_RIGHT:
                    ipicked = -1;
                    break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}

void TestAppRARFF::makePotentialPlot(){
    // DataLine2D * line_Er = new DataLine2D(100); line_Er->clr = 0xFF0000FF; 
    // DataLine2D * line_Fr = new DataLine2D(100); line_Fr->clr = 0xFFFF0000;   line_Fr->replace_xs( line_Er->xs );
    // DataLine2D * line_Fn = new DataLine2D(100); line_Fn->clr = 0xFF008000;   line_Fn->replace_xs( line_Er->xs );
    // testEF( ff, 0.0, 6.0, 60, line_Er->ys, line_Fr->ys,  line_Er->xs, line_Fn->ys );
    // plot1.init();
    // plot1.scaling.y = 1.0;
    // plot1.fontTex = fontTex;
    // plot1.clrGrid = 0xFF404040;
    // //plot1.clrBg   = 0xFF408080;
    // //plot1.lines.push_back( line1  );
    // plot1.lines.push_back( line_Er  );
    // plot1.lines.push_back( line_Fr  );
    // plot1.lines.push_back( line_Fn  );
    // plot1.render();
}

void TestAppRARFF::view_BBoxes(){
    for(int ic=0;ic<ff.pointBBs.ncell;ic++){
        int i0=ff.pointBBs.cellI0s[ic];
        int ni=ff.pointBBs.cellNs [ic];
        //printf( "ic %i io %i ni %i \n", ic, i0, ni );
        Vec3i ip;
        Draw  ::color_of_hash( 464+645*ic );
        //ff.pointBBs.i2ixyz( ic, ip );
        Draw3D::drawBBox( ff.BBs[ic].lo,  ff.BBs[ic].hi );
        for(int j=i0; j<i0+ni;j++){
            int io = ff.pointBBs.cell2obj[j];
            Vec3d p = ff.apos[io];
            //printf( "j %i io %i p(%g,%g,%g) \n", j, io, p.x,p.y,p.z );
            //Draw  ::color_of_hash( 464+645*ic );
            Draw3D::drawPointCross( p, 0.2 );            
        }
    }
}

void TestAppRARFF::view_atoms(){
    //glColor3f(0.0,0.0,0.0);
    double fsc = 0.1;
    double tsc = 0.1;
    //printf( "ff.natom %i \n", ff.natom );
    for(int ia=0; ia<ff.natom; ia++){
        Draw3D::drawShape( ogl_sph, ff.apos[ia], Mat3dIdentity );
    };
}

void TestAppRARFF::make_atoms( int natom, double xspan, double step  ){
    ff.realloc( natom, 16 );
    for(int i=0; i<natom; i++){
        ff.apos [i].fromRandomBox( Vec3d{-10.,-10.,-10.} ,Vec3d{10.,10.,10.});
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



