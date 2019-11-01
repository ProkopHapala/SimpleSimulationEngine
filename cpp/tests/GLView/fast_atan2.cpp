
#include <cmath>
#include <cstdio>
#include "GLView.h"

#include "Draw2D.h"
#include "Plot2D.h"
#include "PlotScreen2D.h"
#include "testUtils.h"

#include "SDL_utils.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

// ===== Globals

Plot2D plot1;
int fontTex;

void my_draw(){
    glClearColor( 1.0f, 1.0f, 1.0f, 0.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glDisable( GL_DEPTH_TEST );
    //plot1.drawAxes();
    plot1.view();
}

int main(){

    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );

    Func1d myFunc = &sin;

    plot1.init();
    plot1.fontTex = fontTex;

    int nsamp = 1000;

    DataLine2D * lcos_ref  = plot1.add( new DataLine2D(nsamp,-3*M_PI,5*M_PI, 0xFFFF0000) );
    DataLine2D * lcos      = plot1.add( new DataLine2D(nsamp,-3*M_PI,5*M_PI, 0xFFFFFF00) );
    DataLine2D * lcos_err  = plot1.add( new DataLine2D(nsamp,-3*M_PI,5*M_PI, 0xFFFF7F00) );

    DataLine2D * lsin_ref  = plot1.add( new DataLine2D(nsamp,-3*M_PI,5*M_PI,0xFF0000FF) );
    DataLine2D * lsin      = plot1.add( new DataLine2D(nsamp,-3*M_PI,5*M_PI,0xFF00FFFF) );
    DataLine2D * lsin_err  = plot1.add( new DataLine2D(nsamp,-3*M_PI,5*M_PI,0xFF007FFF) );

    DataLine2D * lsqrt_ref = plot1.add( new DataLine2D(nsamp,-3*M_PI,5*M_PI,0xFF00FF00) );
    DataLine2D * lsqrt     = plot1.add( new DataLine2D(nsamp,-3*M_PI,5*M_PI,0xFF007F00) );
    DataLine2D * lsqrt_err = plot1.add( new DataLine2D(nsamp,-3*M_PI,5*M_PI,0xFF407F40) );

    DataLine2D * ltan_ref  = plot1.add( new DataLine2D(nsamp,-M_PI/4,M_PI/4,0xFF00FF00) );
    DataLine2D * ltan      = plot1.add( new DataLine2D(nsamp,-M_PI/4,M_PI/4,0xFF007F00) );
    DataLine2D * ltan_err  = plot1.add( new DataLine2D(nsamp,-M_PI/4,M_PI/4,0xFF407F40) );

    for(int i=0; i<nsamp; i++){
        lcos_ref->ys[i] = cos(lcos_ref->xs[i]);
        lsin_ref->ys[i] = sin(lsin_ref->xs[i]);
        lsqrt_ref->ys[i] = 1/sqrt(lsqrt_ref->xs[i]);
        cos_sin( lcos_ref->xs[i], lcos->ys[i], lsin->ys[i] );
        //cos_sin( lcos_ref->xs[i], lcos->ys[i], lsin->ys[i], 2,1 );
        //cos_sin( lcos_ref->xs[i], lcos->ys[i], lsin->ys[i], 2,2 );
        //cos_sin( lcos_ref->xs[i], lcos->ys[i], lsin->ys[i], 2,3 );
        //cos_sin( lcos_ref->xs[i], lcos->ys[i], lsin->ys[i], 4,1 );
        //cos_sin( lcos_ref->xs[i], lcos->ys[i], lsin->ys[i], 4,2 );
        //cos_sin( lcos_ref->xs[i], lcos->ys[i], lsin->ys[i], 4,3 );
        //cos_sin( lcos_ref->xs[i], lcos->ys[i], lsin->ys[i], 6,1 );
        //cos_sin( lcos_ref->xs[i], lcos->ys[i], lsin->ys[i], 6,2 );
        //cos_sin( lcos_ref->xs[i], lcos->ys[i], lsin->ys[i], 6,3 );
        lsqrt->ys[i]    = fastInvSqrt( lsqrt_ref->xs[i] );
        lcos_err->ys[i] = log10( fabs( lcos->ys[i] - lcos_ref->ys[i] ) );
        lsin_err->ys[i] = log10( fabs( lsin->ys[i] - lsin_ref->ys[i] ) );
        lsqrt_err->ys[i] = log10( fabs( lsqrt->ys[i] - lsqrt_ref->ys[i] ) );

        double x = ltan_ref->xs[i];
        ltan->ys[i]     = x*tan_xx_12( x*x );
        ltan_ref->ys[i] = tan(ltan_ref->xs[i]);
        ltan_err->ys[i] = log10( fabs( ltan->ys[i] - ltan_ref->ys[i] ) );
        //printf( "[%i]: x %g err %g %g \n", i, lcos_ref->xs[i], lcos_err->ys[i], lsin_err->ys[i] );
    }

    plot1.render();

    //return;

    const int n = 1000;
    const int m = 1000;
    double xs[n];
    double xs_[n];
    for(int i=0; i<n; i++){
        xs [i]=randf(-30.0,30.0);
        xs_[i]=randf(-M_PI/4,M_PI/4);
    }
    //VecN::arange(n,-30.0,60./n,xs);

    double c=0,s=0;
    double dn=1./n;


    //TEST_ERROR_PROC_N( "tan() "       ,{double x=xs[i]; c=tan(x);        c-=1/sqrt(x); STORE_ERROR(c) }, n );
    TEST_ERROR_PROC_N( "tan_xx_12 () ",{double x=xs_[i]; c=x*tan_xx_12 (x*x); c-=tan(x); STORE_ERROR(c) }, n );
    TEST_ERROR_PROC_N( "tan_xx_12_() ",{double x=xs_[i]; c=x*tan_xx_12_(x*x); c-=tan(x); STORE_ERROR(c) }, n );

    TEST_ERROR_PROC_N( "fastInvSqrt() ",{double x=xs[i]; c=fastInvSqrt(x); c-=1/sqrt(x); STORE_ERROR(c) }, n );

    TEST_ERROR_PROC_N( "cos_sin()     ", {double x=xs[i]; cos_sin(x,c,s); c-=cos(x); s-=sin(x); STORE_ERROR(c); STORE_ERROR(s); }, n );

    TEST_ERROR_PROC_N( "cos_sin(2,1) ", {double x=xs[i]; cos_sin(x,c,s,2,1); c-=cos(x); s-=sin(x); STORE_ERROR(c); STORE_ERROR(s); }, n );
    TEST_ERROR_PROC_N( "cos_sin(2,2) ", {double x=xs[i]; cos_sin(x,c,s,2,2); c-=cos(x); s-=sin(x); STORE_ERROR(c); STORE_ERROR(s); }, n );
    TEST_ERROR_PROC_N( "cos_sin(2,3) ", {double x=xs[i]; cos_sin(x,c,s,2,3); c-=cos(x); s-=sin(x); STORE_ERROR(c); STORE_ERROR(s); }, n );

    TEST_ERROR_PROC_N( "cos_sin(4,1) ", {double x=xs[i]; cos_sin(x,c,s,4,1); c-=cos(x); s-=sin(x); STORE_ERROR(c); STORE_ERROR(s); }, n );
    TEST_ERROR_PROC_N( "cos_sin(4,2) ", {double x=xs[i]; cos_sin(x,c,s,4,2); c-=cos(x); s-=sin(x); STORE_ERROR(c); STORE_ERROR(s); }, n );
    TEST_ERROR_PROC_N( "cos_sin(4,3) ", {double x=xs[i]; cos_sin(x,c,s,4,3); c-=cos(x); s-=sin(x); STORE_ERROR(c); STORE_ERROR(s); }, n );

    TEST_ERROR_PROC_N( "cos_sin(6,1) ", {double x=xs[i]; cos_sin(x,c,s,6,1); c-=cos(x); s-=sin(x); STORE_ERROR(c); STORE_ERROR(s); }, n );
    TEST_ERROR_PROC_N( "cos_sin(6,2) ", {double x=xs[i]; cos_sin(x,c,s,6,2); c-=cos(x); s-=sin(x); STORE_ERROR(c); STORE_ERROR(s); }, n );
    TEST_ERROR_PROC_N( "cos_sin(6,3) ", {double x=xs[i]; cos_sin(x,c,s,6,3); c-=cos(x); s-=sin(x); STORE_ERROR(c); STORE_ERROR(s); }, n );

    SPEED_TEST_PROC_NM( "junk;        ", {double x=xs[i];s+=x;c+=x;sum=c+s;}, n, m );


    SPEED_TEST_PROC_NM( "tan() "      , {sum+= c=tan(xs_[i]);      }, n, m );
    SPEED_TEST_PROC_NM( "tan_xx_12 ()" ,{double x=xs_[i]; sum+=x*tan_xx_12 (x*x); }, n, m );
    SPEED_TEST_PROC_NM( "tan_xx_12_()" ,{double x=xs_[i]; sum+=x*tan_xx_12_(x*x); }, n, m );

    SPEED_TEST_PROC_NM( "1/sqrt()      ", {sum+=1/sqrt(xs[i]);      }, n, m );
    SPEED_TEST_PROC_NM( "fastInvSqrt() ", {sum+=fastInvSqrt(xs[i]); }, n, m );

    SPEED_TEST_PROC_NM( "cos();sin()  ", {double x=xs[i];sum+=cos(x)+sin(x);      }, n, m );
    SPEED_TEST_PROC_NM( "sincos()     ", {sincos (xs[i],&c,&s);sum+=c+s; }, n, m );
    SPEED_TEST_PROC_NM( "cos_sin( )   ", {cos_sin(xs[i],c,s);sum+=c+s; }, n, m );

    SPEED_TEST_PROC_NM( "cos_sin(2,1) ", {cos_sin(xs[i],c,s,2,1);sum+=c+s; }, n, m );
    SPEED_TEST_PROC_NM( "cos_sin(2,2) ", {cos_sin(xs[i],c,s,2,2);sum+=c+s; }, n, m );
    SPEED_TEST_PROC_NM( "cos_sin(2,3) ", {cos_sin(xs[i],c,s,2,3);sum+=c+s; }, n, m );

    SPEED_TEST_PROC_NM( "cos_sin(4,1) ", {cos_sin(xs[i],c,s,4,1);sum+=c+s; }, n, m );
    SPEED_TEST_PROC_NM( "cos_sin(4,2) ", {cos_sin(xs[i],c,s,4,2);sum+=c+s; }, n, m );
    SPEED_TEST_PROC_NM( "cos_sin(4,3) ", {cos_sin(xs[i],c,s,4,3);sum+=c+s; }, n, m );

    SPEED_TEST_PROC_NM( "cos_sin(6,1) ", {cos_sin(xs[i],c,s,6,1);sum+=c+s; }, n, m );
    SPEED_TEST_PROC_NM( "cos_sin(6,2) ", {cos_sin(xs[i],c,s,6,2);sum+=c+s; }, n, m );
    SPEED_TEST_PROC_NM( "cos_sin(6,3) ", {cos_sin(xs[i],c,s,6,3);sum+=c+s; }, n, m );


    init( 640, 480 );
    set_draw_function( my_draw );
    run_Nframes(5000);

}
