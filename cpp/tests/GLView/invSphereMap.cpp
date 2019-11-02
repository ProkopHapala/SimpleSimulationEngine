
#include <cmath>
#include <cstdio>
#include "GLView.h"


#include "fastmath.h"

#include "Draw.h"
#include "Draw2D.h"
#include "Plot2D.h"
#include "PlotScreen2D.h"
#include "testUtils.h"

#include "SDL_utils.h"

#include "gonioApprox.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

Plot2D plot1;
//int fontTex;

inline int hexface(double y, double x){
    //double t = y/x;
    int i;
    if     (x> 0.5){i=0;}
    else if(x<-0.5){i=2;} 
    else           {i=1;}
    if(y<0){ i=5-i; };
    return i;
}

inline double hex_atan_10( double x){
    double xx = x*x;
    return                     xx*(-0.002745519031537113 +xx*(-0.1122625697398381 +xx*(-0.4292923464328494 +xx*(-0.04728206847314129 ) ) ) ) 
           +x*( 0.95504016223 +xx*(-0.2935472108431604   +xx*(0.4793036343330538  +xx*(0.1959841630153277  +xx*(0.004801286948495104 ) ) ) ) );
}

inline double hex_atan_xx_22( double xx ){  
    /*
    double x4 = xx*xx;
    return 0.954926185539 
        +x4*(0.1896742126887878 
        +x4*(0.0869966813196642 
        +x4*(0.02262655250936047 
        +x4*(0.001587753332013546 
        +x4*(1.149172690378952e-05 ) ) ) ) ) 
        +xx*(-0.3181862047974173 
        +x4*(-0.1298472375382908 
        +x4*(-0.05024162161949579 
        +x4*(-0.007345327319912838 
        +x4*(-0.0002026406963035743 ) ) ) ) );
    */
    return    0.954926185539 
        +xx*(-0.3181862047974173 
        +xx*( 0.1896742126887878 
        +xx*(-0.1298472375382908 
        +xx*( 0.0869966813196642 
        +xx*(-0.05024162161949579 
        +xx*( 0.02262655250936047 
        +xx*(-0.007345327319912838 
        +xx*( 0.001587753332013546 
        +xx*(-0.0002026406963035743 
        +xx*( 1.149172690378952e-05 )))))))))); 
}

inline double hex_atan_xx_10( double xx ){  
    return  -0.954928347643 +xx*(0.3182205591492981 +xx*(-0.1892837387809918 +xx*(0.1229241138904131 +xx*(-0.05607123646283267 ) ) ) ) ;
}

inline double hexcoord(double y, double x, int& i){
    //double t = y/x;
    //int i;
    double d;
    const bool by = y<0;
    if     (x> 0.5){
        i = 0;
        double a = y/x;
        d = a*hex_atan_xx_22( a*a );
        d = by?d+1:d;
    }else if(x<-0.5){
        i = 2;
        double a = y/x;
        d = a*hex_atan_xx_22( a*a );
        d = by?d:d+1;
    }else{ // |x|<0.5 y!=0
        i = 1;
        double a = x/y;
        d = a*hex_atan_xx_10( a*a )+0.5;
    }
    if(by){ i = 5-i; };
    return d;
}


/*

A = ( 0.5, 0.8 )
B = ( 1.0, 0.0 )

A*a + b*B = (x,y)

x = 0.5*a + 1.0*b
y = 0.8*a + 0.0*b

x = 0.5*a + b
y = 0.8*a

a = y/0.8

c = -0.5/0.8

x + c*y = b


*/


inline double hexcoord_lin(double y, double x, int& i){
    //double t = y/x;
    //int i;
    double d;
    //double cos60 = 0.86602540378;
    const bool by = y<0;
    double abs_y  = fabs(y);
    if     (x> 0.5){
        i = 0;
        //double a =     abs_y*(  1/0.86602540378);
        //double b = x + abs_y*(-0.5/0.86602540378);
        //d = a/(a+b);
        d = y/( abs_y*0.5 + x*0.86602540378 );
        d = by?d+1:d;
        //d=0;
    }else if(x<-0.5){
        i = 2;
        d = y/( abs_y*-0.5 + x*0.86602540378 );
        d = by?d:d+1;
        //d=0;
    }else{ 
        i = 1;
        d = (0.5-0.86602540378*x/y);
    }
    if(y<0){ i = 5-i; };
    return d;
}

inline double hex_atan_lin_xx_10( double xx ){ 
    return 1.1026577831 +xx*(-0.4900688562766967 +xx*(0.3919957429499588 +xx*(-0.3721192362110416 +xx*(0.3733397771719061 +xx*(-0.337901579282555 +xx*(0.1871287576348613 ) ) ) ) ) );
};

inline double hexcoord_10(double y, double x, int& i){
    double a;
    const bool by = y<0;
    double abs_y  = fabs(y);
    if     (x> 0.5){
        i = 0;
        a = y/( abs_y*0.5 + x*0.86602540378 );
        a = by?a+0.5:a-0.5;
    }else if(x<-0.5){
        i = 2;
        a = y/( abs_y*-0.5 + x*0.86602540378 );
        a = by?a-0.5:a+0.5;
    }else{
        i = 1;
        a = (-0.86602540378*x/y);
    }
    double d = a*hex_atan_lin_xx_10( a*a ) + 0.5;
    if(y<0){ i = 5-i; };
    return d;
}

inline void hexSphere( Vec3d p, Vec2d& d, Vec2i& ip ){
    double r2xy = p.x*p.x + p.z*p.z;
    d.y = hexcoord( r2xy, p.z, ip.y );
    d.x = hexcoord( p.x,  p.y, ip.x );
}

void my_draw(){

    //printf(" my_draw ! %i %i %g \n", plot1.lines.size(), plot1.lines.back()->n, plot1.lines.back()->ys[2]  );
    glClearColor( 1.0f, 1.0f, 1.0f, 0.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glDisable( GL_DEPTH_TEST );
    
    //plot1.drawAxes();
    plot1.view();

    //Draw2D::drawLine({-1,0},{1,0});

    /*
    int n = 100;
    double dx = 2*M_PI/n;

    for(int i=0; i<100; i++){
        double phi = dx*i;
        double x = sin(phi);
        double y = cos(phi);
        //Draw::setRGB(  );
        Draw::colorScale( (int)(hexface(y,x)/6.0) );
        glVertex3f( x, y, 0.0 );
    }
    */

}

void setup(){

    //fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );

    plot1.init();
    //plot1.fontTex = fontTex;

    int nsamp = 1000;
    double xmin=-2*M_PI;
    double xmax= 2*M_PI;
    //Draw::icolorScale(0.1);

    DataLine2D * lref = plot1.add( new DataLine2D(nsamp,xmin,xmax, 0xFF000000 ) );
    DataLine2D * l1    = plot1.add( new DataLine2D(nsamp,xmin,xmax, 0xFFFF0000 ) );
    DataLine2D * l2    = plot1.add( new DataLine2D(nsamp,xmin,xmax, 0xFF0000FF ) );

    DataLine2D * le1   = plot1.add( new DataLine2D(nsamp,xmin,xmax, 0xFFFF8000 ) );
    DataLine2D * le2   = plot1.add( new DataLine2D(nsamp,xmin,xmax, 0xFF0080FF ) );

    int ii;
    for(int i=0; i<nsamp; i++){
        double phi = lref->xs[i];
        double x = cos( phi );
        double y = sin( phi );
        //printf(  "%i %g (%g,%g)\n", i, phi, x, y );
        //lref->ys[i] = (int)( fastFract(phi/(2*M_PI))*6.0 );
        lref->ys[i] = fastFract(3*phi/M_PI );
        //lref->ys[i] = fastFract( (3.0/M_PI)*atan2(y,x) );
        //l   ->ys[i] = hexface(y,x);
        l1   ->ys[i] = hexcoord ( y, x, ii );
        //l   ->ys[i] = hexcoord_lin( y, x, ii );
        l2   ->ys[i] = hexcoord_10( y, x, ii );

        le1->ys[i] = log10( fabs( l1->ys[i] - lref->ys[i] ) );
        le2->ys[i] = log10( fabs( l2->ys[i] - lref->ys[i] ) );
    }

    plot1.render();


    const int n = 1000;
    const int m = 1000;
    double xs [n];
    double ys [n];
    for(int i=0; i<n; i++){
        //xs[i]=randf(-30.0,30.0);
        //ys[i]=randf(-30.0,30.0);
        double phi = randf(-30,30);
        xs[i]=cos(phi);
        ys[i]=sin(phi);
    }
    //VecN::arange(n,-30.0,60./n,xs);

    double c=0,s=0;
    double dn=1./n;

    TEST_ERROR_PROC_N( "hexcoord     ",{double x=xs[i]; double y=ys[i]; c=hexcoord   (y,x,ii); c-=fastFract( (3.0/M_PI)*atan2(y,x) ); STORE_ERROR(c) }, n );
    TEST_ERROR_PROC_N( "hexcoord_10  ",{double x=xs[i]; double y=ys[i]; c=hexcoord_10(y,x,ii); c-=fastFract( (3.0/M_PI)*atan2(y,x) ); STORE_ERROR(c) }, n );

    SPEED_TEST_PROC_NM( "junk;        ", { sum+=xs[i]; sum+=ys[i]; }, n, m );
    SPEED_TEST_PROC_NM( "atan2        " , {sum+= atan2   (ys[i],xs[i]); }, n, m );
    SPEED_TEST_PROC_NM( "hexface      " , {sum+= hexface (ys[i],xs[i]); }, n, m );
    SPEED_TEST_PROC_NM( "hexcoord     " , {sum+= hexcoord(ys[i],xs[i], ii); }, n, m );
    SPEED_TEST_PROC_NM( "hexcoord_lin " , {sum+= hexcoord_lin(ys[i],xs[i], ii); }, n, m );
    SPEED_TEST_PROC_NM( "hexcoord_10 " , {sum+= hexcoord_10(ys[i],xs[i], ii); }, n, m );


}

int main(){
    init( 800, 600 );
    set_draw_function( my_draw );
    setup();
    run_Nframes(5000);
}
