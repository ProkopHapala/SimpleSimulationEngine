
#include <cmath>
#include <cstdio>
#include "GLView.h"

#include "Draw.h"
#include "Draw2D.h"
#include "Plot2D.h"
#include "PlotScreen2D.h"
#include "testUtils.h"

#include "SDL_utils.h"

#include "gonioApprox.h"


#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

/*

// ======= Accuracy
erfx_e6(x)    MaxErr: 3.65274e-07 RMSE: 9.70469e-08 
erfx_e9(x)    MaxErr: 1.55826e-09 RMSE: 3.80278e-10 

// ======= Performance
junk;           : 0.752 ticks/call ( 7.52169e+06 1e+07 ) | 4.065e+07 
erf(x)/x        : 60.398 ticks/call ( 6.03977e+08 1e+07 ) | 3.7214e+06 
erfx_e6(x)      : 2.706 ticks/call ( 2.70644e+07 1e+07 ) | 3.7214e+06 
erfx_e9(x)      : 4.154 ticks/call ( 4.15412e+07 1e+07 ) | 3.7214e+06 

*/



// ===== Globals

Plot2D plot1;
//int fontTex;

void my_draw(){

    //printf(" my_draw ! %i %i %g \n", plot1.lines.size(), plot1.lines.back()->n, plot1.lines.back()->ys[2]  );
    glClearColor( 1.0f, 1.0f, 1.0f, 0.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glDisable( GL_DEPTH_TEST );
    //plot1.drawAxes();

    Draw2D::drawLine({-1,0},{1,0});

    plot1.drawAxes();
    plot1.view();

}

void setup(){

fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );

    Func1d myFunc = &sin;

    plot1.init();
    plot1.fontTex = fontTex;

    int nsamp = 700;
    double xmin=1.e-6;
    double dx  =0.01;
    
    DataLine2D * l_ref  = plot1.add( new DataLine2D(nsamp,xmin,dx, 0xFF000000 ) );
    DataLine2D * l0     = plot1.add( new DataLine2D(nsamp,xmin,dx, Draw::icolorScale(0.5) ) );
    DataLine2D * l1     = plot1.add( new DataLine2D(nsamp,xmin,dx, Draw::icolorScale(0.2) ) );

    DataLine2D * lerr0  = plot1.add( new DataLine2D(nsamp,xmin,dx, Draw::icolorScale(0.5) ) );
    DataLine2D * lerr1  = plot1.add( new DataLine2D(nsamp,xmin,dx, Draw::icolorScale(0.25) ) );


    for(int i=0; i<nsamp; i++){
        double x = l_ref->xs[i];
        //printf(  "%i %g (%g,%g)\n", i, phi, x, y );
        l_ref->ys[i] = erf(x)/x;


        l0->ys[i] = erfx_e6(x);
        l1->ys[i] = erfx_e9(x);

        lerr0->ys[i] = log10( fabs( l0->ys[i] - l_ref->ys[i] ) );
        lerr1->ys[i] = log10( fabs( l1->ys[i] - l_ref->ys[i] ) );
    }

    plot1.render();

    //return;

    
    const int n = 1000;
    const int m = 10000;
    double xs [n];
    for(int i=0; i<n; i++){
        xs[i]=randf(1e-6,8.0);
    }
    //VecN::arange(n,-30.0,60./n,xs);

    double c=0,s=0;
    double dn=1./n;

    TEST_ERROR_PROC_N( "erfx_e6(x)   ",{double x=xs[i]; c=erfx_e6(xs[i]); c-=erf(x)/x; STORE_ERROR(c) }, n );
    TEST_ERROR_PROC_N( "erfx_e9(x)   ",{double x=xs[i]; c=erfx_e9(xs[i]); c-=erf(x)/x; STORE_ERROR(c) }, n );


    SPEED_TEST_PROC_NM( "junk;          ", {sum+=xs[i];         }, n, m );
    SPEED_TEST_PROC_NM( "erf(x)/x       ", {double x=xs[i]; sum+=erf(x)/x; }, n, m );
    SPEED_TEST_PROC_NM( "erfx_e6(x)     ", {sum+=erfx_e6(xs[i]); }, n, m );
    SPEED_TEST_PROC_NM( "erfx_e9(x)     ", {sum+=erfx_e9(xs[i]); }, n, m );
}

int main(){

    init( 640, 480 );
    set_draw_function( my_draw );

    setup();

    run_Nframes(5000);
    
}
