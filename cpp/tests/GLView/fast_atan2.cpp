
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

Results:

ACCURACY:
atan2_a1      MaxErr: 2.90068e-05 RMSE: 1.08498e-05 
atan2_a2      MaxErr: 1.69537e-05 RMSE: 8.26705e-06 
atan2_a3      MaxErr: 1.69537e-05 RMSE: 8.26705e-06 
atan2_nvidia  MaxErr: 3.41833e-06 RMSE: 1.52343e-06 

SPEED:
junk;           : 0.923 ticks/call ( 922786 1e+06 ) | 248552 
atan2()         : 66.998 ticks/call ( 6.69982e+07 1e+06 ) | -25546.6 
atan2_a1()      : 4.034 ticks/call ( 4.03435e+06 1e+06 ) | -25546.7 
atan2_a2()      : 10.090 ticks/call ( 1.00904e+07 1e+06 ) | -25546.3 
atan2_a3()      : 6.713 ticks/call ( 6.71309e+06 1e+06 ) | -25546.3 
atan2_nvidia()  : 16.877 ticks/call ( 1.68767e+07 1e+06 ) | -25546.7 

*/




// ===== Globals

Plot2D plot1;
int fontTex;

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

    int nsamp = 1000;

    double xmin=-3*M_PI;
    double xmax= 5*M_PI;
    //Draw::icolorScale(0.1);
    
    DataLine2D * l_ref  = plot1.add( new DataLine2D(nsamp,xmin,xmax, 0xFF000000 ) );
    DataLine2D * l0     = plot1.add( new DataLine2D(nsamp,xmin,xmax, Draw::icolorScale(0.5) ) );
    DataLine2D * l1     = plot1.add( new DataLine2D(nsamp,xmin,xmax, Draw::icolorScale(0.2) ) );
    //DataLine2D * l2     = plot1.add( new DataLine2D(nsamp,xmin,xmax, Draw::icolorScale(0.4) ) );
    //DataLine2D * l3     = plot1.add( new DataLine2D(nsamp,xmin,xmax, Draw::icolorScale(0.6) ) );
    //DataLine2D * l4     = plot1.add( new DataLine2D(nsamp,xmin,xmax, Draw::icolorScale(0.8) ) );

    DataLine2D * lerr0  = plot1.add( new DataLine2D(nsamp,xmin,xmax, Draw::icolorScale(0.5) ) );
    DataLine2D * lerr1  = plot1.add( new DataLine2D(nsamp,xmin,xmax, Draw::icolorScale(0.25) ) );
    //DataLine2D * lerr2  = plot1.add( new DataLine2D(nsamp,xmin,xmax, Draw::icolorScale(0.45) ) );
    //DataLine2D * lerr3  = plot1.add( new DataLine2D(nsamp,xmin,xmax, Draw::icolorScale(0.65) ) );
    //DataLine2D * lerr4  = plot1.add( new DataLine2D(nsamp,xmin,xmax, Draw::icolorScale(0.85) ) );

    //printf( "colors %i %i %i %i \n", plot1.lines[0]->clr, plot1.lines[1]->clr, plot1.lines[2]->clr, plot1.lines[3]->clr, plot1.lines[3]->clr );

    for(int i=0; i<nsamp; i++){
        double phi = l_ref->xs[i];
        double x = cos( phi );
        double y = sin( phi );
        //printf(  "%i %g (%g,%g)\n", i, phi, x, y );
        l_ref->ys[i] = atan2(y,x);
        //l1->ys[i] = atan2_a1(y,x);
        //l2->ys[i] = atan2_a2(y,x);
        //l3->ys[i] = atan2_a3(y,x);
        //l4->ys[i] = atan2_nvidia(y,x);

        l0->ys[i] = atan2_a1(y,x);
        l1->ys[i] = atan2_t<atan2_xx_8 >(y,x);
        //l2->ys[i] = atan2_t<atan2_xx_8 >(y,x);
        //l3->ys[i] = atan2_t<atan2_xx_10>(y,x);
        //l4->ys[i] = atan2_t<atan2_xx_12>(y,x);

        lerr0->ys[i] = log10( fabs( l0->ys[i] - l_ref->ys[i] ) );
        lerr1->ys[i] = log10( fabs( l1->ys[i] - l_ref->ys[i] ) );
        //lerr2->ys[i] = log10( fabs( l2->ys[i] - l_ref->ys[i] ) );
        //lerr3->ys[i] = log10( fabs( l3->ys[i] - l_ref->ys[i] ) );
        //lerr4->ys[i] = log10( fabs( l4->ys[i] - l_ref->ys[i] ) );
    }

    plot1.render();

    //return;

    
    const int n = 1000;
    const int m = 1000;
    double xs [n];
    double ys [n];
    for(int i=0; i<n; i++){
        xs[i]=randf(-30.0,30.0);
        ys[i]=randf(-30.0,30.0);
    }
    //VecN::arange(n,-30.0,60./n,xs);

    double c=0,s=0;
    double dn=1./n;

    TEST_ERROR_PROC_N( "atan2_a1     ",{double x=xs[i]; double y=ys[i]; c=atan2_a1    (y,x); c-=atan2(y,x); STORE_ERROR(c) }, n );

    TEST_ERROR_PROC_N( "atan2_t6     ",{double x=xs[i]; double y=ys[i]; c=atan2_t<atan2_xx_6 >(y,x); c-=atan2(y,x); STORE_ERROR(c) }, n );
    TEST_ERROR_PROC_N( "atan2_t8     ",{double x=xs[i]; double y=ys[i]; c=atan2_t<atan2_xx_8 >(y,x); c-=atan2(y,x); STORE_ERROR(c) }, n );
    TEST_ERROR_PROC_N( "atan2_t10    ",{double x=xs[i]; double y=ys[i]; c=atan2_t<atan2_xx_10>(y,x); c-=atan2(y,x); STORE_ERROR(c) }, n );
    TEST_ERROR_PROC_N( "atan2_t12    ",{double x=xs[i]; double y=ys[i]; c=atan2_t<atan2_xx_12>(y,x); c-=atan2(y,x); STORE_ERROR(c) }, n );
    TEST_ERROR_PROC_N( "atan2_t14    ",{double x=xs[i]; double y=ys[i]; c=atan2_t<atan2_xx_14>(y,x); c-=atan2(y,x); STORE_ERROR(c) }, n );

    //TEST_ERROR_PROC_N( "atan2_a2     ",{double x=xs[i]; double y=ys[i]; c=atan2_a2    (y,x); c-=atan2(y,x); STORE_ERROR(c) }, n );
    //TEST_ERROR_PROC_N( "atan2_a3     ",{double x=xs[i]; double y=ys[i]; c=atan2_a3    (y,x); c-=atan2(y,x); STORE_ERROR(c) }, n );
    //TEST_ERROR_PROC_N( "atan2_nvidia ",{double x=xs[i]; double y=ys[i]; c=atan2_nvidia(y,x); c-=atan2(y,x); STORE_ERROR(c) }, n );
    //TEST_ERROR_PROC_N( "atan2_fnvidia ",{double x=xs[i]; double y=ys[i]; c=atan2_nvidia_fabs(y,x); c-=atan2(y,x); STORE_ERROR(c) }, n );

    SPEED_TEST_PROC_NM( "junk;          ", { sum+=xs[i]; sum+=ys[i]; }, n, m );

    SPEED_TEST_PROC_NM( "atan2        " , {sum+= c=atan2       (ys[i],xs[i]); }, n, m );
    SPEED_TEST_PROC_NM( "atan2_a1     " , {sum+= c=atan2_a1    (ys[i],xs[i]); }, n, m );

    SPEED_TEST_PROC_NM( "atan2_t6     " , {sum+= c=atan2_t<atan2_xx_6 >(ys[i],xs[i]); }, n, m );
    SPEED_TEST_PROC_NM( "atan2_t8     " , {sum+= c=atan2_t<atan2_xx_8 >(ys[i],xs[i]); }, n, m );
    SPEED_TEST_PROC_NM( "atan2_t10    " , {sum+= c=atan2_t<atan2_xx_10>(ys[i],xs[i]); }, n, m );
    SPEED_TEST_PROC_NM( "atan2_t12    " , {sum+= c=atan2_t<atan2_xx_12>(ys[i],xs[i]); }, n, m );
    SPEED_TEST_PROC_NM( "atan2_t14    " , {sum+= c=atan2_t<atan2_xx_14>(ys[i],xs[i]); }, n, m );

    //SPEED_TEST_PROC_NM( "atan2_a2()     " , {sum+= c=atan2_a2    (ys[i],xs[i]); }, n, m );
    //SPEED_TEST_PROC_NM( "atan2_a3()     " , {sum+= c=atan2_a3    (ys[i],xs[i]); }, n, m );
    //SPEED_TEST_PROC_NM( "atan2_nvidia() " , {sum+= c=atan2_nvidia(ys[i],xs[i]); }, n, m );
    //SPEED_TEST_PROC_NM( "atan2_fnvidia()", {sum+= c=atan2_nvidia_fabs(ys[i],xs[i]); }, n, m );

}

int main(){

    init( 640, 480 );
    set_draw_function( my_draw );

    setup();

    run_Nframes(5000);
    
}
