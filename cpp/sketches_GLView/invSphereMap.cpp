
#include <cmath>
#include <cstdio>
#include "GLView.h"


#include "fastmath.h"

#include "Draw.h"
#include "Draw2D.h"
#include "Draw3D.h"
#include "Plot2D.h"
#include "PlotScreen2D.h"
#include "testUtils.h"

#include "SDL_utils.h"

#include "gonioApprox.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>


int ogl=0;
Plot2D plot1;
//int fontTex;



// Fibonachi sampling Sphere
// http://extremelearning.com.au/evenly-distributing-points-on-a-sphere/
// https://arxiv.org/pdf/0912.4540.pdf
// https://blog.wolfram.com/2011/07/28/how-i-made-wine-glasses-from-sunflowers/

// Fibonachi Spiral Disk - Texture map-  https://www.shadertoy.com/view/ldd3Wn
// O = texture(iChannel0, vec2(5,8)*( exp2(-length(I=(I+I-R)/R.y)) + atan(I,I.yx)/6.28) );


// Inverse Finbonachi on sphere
// https://www.shadertoy.com/view/lllXz4
// https://www.shadertoy.com/view/XlfyWl
// Spherical Fibonacci Mapping
// http://lgdv.cs.fau.de/publications/publication/Pub.2015.tech.IMMD.IMMD9.spheri/
// Authors: Benjamin Keinert, Matthias Innmann, Michael Sänger, Marc Stamminger
// (code copied from: https://www.shadertoy.com/view/4t2XWK)
/*
float sf2id(vec3 p, float n) 
{
    float phi = min(atan(p.y, p.x), PI), cosTheta = p.z;
    
    float k  = max(2.0, floor( log(n * PI * sqrt(5.0) * (1.0 - cosTheta*cosTheta))/ log(PHI*PHI)));
    float Fk = pow(PHI, k)/sqrt(5.0);
    
    vec2 F = vec2( round(Fk), round(Fk * PHI) );

    vec2 ka = -2.0*F/n;
    vec2 kb = 2.0*PI*madfrac(F+1.0, PHI-1.0) - 2.0*PI*(PHI-1.0);    
    mat2 iB = mat2( ka.y, -ka.x, -kb.y, kb.x ) / (ka.y*kb.x - ka.x*kb.y);

    vec2 c = floor( iB * vec2(phi, cosTheta - (1.0-1.0/n)));
    float d = 8.0;
    float j = 0.0;
    for( int s=0; s<4; s++ ) 
    {
        vec2 uv = vec2( float(s-2*(s/2)), float(s/2) );
        
        float cosTheta = dot(ka, uv + c) + (1.0-1.0/n);
        
        cosTheta = clamp(cosTheta, -1.0, 1.0)*2.0 - cosTheta;
        float i = floor(n*0.5 - cosTheta*n*0.5);
        float phi = 2.0*PI*madfrac(i, PHI-1.0);
        cosTheta = 1.0 - (2.0*i + 1.)/n;
        float sinTheta = sqrt(1.0 - cosTheta*cosTheta);
        
        vec3 q = vec3( cos(phi)*sinTheta, sin(phi)*sinTheta, cosTheta);
        float squaredDistance = dot(q-p, q-p);
        if (squaredDistance < d) 
        {
            d = squaredDistance;
            j = i;
        }
    }
    return j;
}

vec3 id2sf( float i, float n) 
{
    float phi = 2.0*PI*madfrac(i,PHI);
    float zi = 1.0 - (2.0*i+1.)/n;
    float sinTheta = sqrt( 1.0 - zi*zi);
    return vec3( cos(phi)*sinTheta, sin(phi)*sinTheta, zi);
}
*/



// ==================== Functions

void SphereSamplesFibonachi( int n, Vec3d* ps ){
    constexpr const double golden_angle = M_PI * ( 3.0 - sqrt(5.0) );
    //theta  = golden_angle * np.arange(n)
    double dz = 2./n;
    double z  = -1+dz*0.5;
    for(int i=0; i<n; i++){
        //z      = np.linspace(1.0 - 1.0/n, 1.0/n - 1.0, n)
        double r   = sqrt( 1. - z*z );
        double phi = golden_angle * i;
        ps[i]      = {r*cos(phi),r*sin(phi),z};
        z+=dz;
    }
}

inline int hexface(double y, double x){
    int i=0;
    const bool by = y<0;
    if(by){ y=-y; };
    if     (x<-0.57735026919*y){i=2;}
    else if(x<+0.57735026919*y){i=1;} 
    if(by){ i=5-i; };
    return i;
}

inline double hexcoord_lin(double y, double x, int& i){
    double d, abs_y;
    x*=0.86602540378;
    const bool by = y<0;
    if(by){ abs_y=y*-0.5; }else{ abs_y=y*0.5; }
    double l1 = x - abs_y;
    double l2 = x + abs_y;
    if      (l2<0){
        i=2;
        d = y/l1;
        d = by?d:d+1;
    }else if(l1<0){
        i=1;
        d = (0.5-x/y);
    }else{
        i=0;
        d = y/l2;
        d = by?d+1:d;
    }
    if(by){ i=5-i; };
    return d;
}

inline double hex_atan_xx_lin_10( double xx ){ 
    return 1.1026577831 +xx*(-0.4900688562766967 +xx*(0.3919957429499588 +xx*(-0.3721192362110416 +xx*(0.3733397771719061 +xx*(-0.337901579282555 +xx*(0.1871287576348613 ) ) ) ) ) );
}

inline double hexcoord_10(double y, double x, int& i){
    double a, abs_y;
    x*=0.86602540378;
    const bool by = y<0;
    if(by){ abs_y=y*-0.5; }else{ abs_y=y*0.5; }
    double l1 = x - abs_y;
    double l2 = x + abs_y;
    if      (l2<0){
        i=2;
        a = y/l1;
        a = by?a-0.5:a+0.5;
    }else if(l1<0){
        i=1;
        a = -x/y;
    }else{
        i=0;
        a = y/l2;
        a = by?a+0.5:a-0.5;
    }
    double d = a*hex_atan_xx_lin_10( a*a ) + 0.5;
    if(by){ i=5-i; };
    return d;
}

inline void hexSphere_lin( const Vec3d& p, Vec2d& d, Vec2i& ip ){
    double r2xy = p.x*p.x + p.y*p.y;
    d.y = hexcoord_lin( r2xy, p.z, ip.y );
    d.x = hexcoord_lin( p.y,  p.x, ip.x );
}

inline void hexSphere_10( const Vec3d& p, Vec2d& d, Vec2i& ip ){
    double r2xy = p.x*p.x + p.y*p.y;
    d.y = hexcoord_10( r2xy, p.z, ip.y );
    d.x = hexcoord_10( p.y,  p.x, ip.x );
}

inline int octFace( const Vec3d& p ){
    return (p.x>0)|((p.y>0)<<1)|((p.z>0)<<2);
}

inline int octCoord(const Vec3d& p, Vec3d& d ){
    int    i=0;
    //double r=0;
    if(p.x>0){ i|=1; d.x=p.x; }else{ d.x=-p.x; };
    if(p.y>0){ i|=2; d.y=p.y; }else{ d.y=-p.y; };
    if(p.z>0){ i|=4; d.z=p.z; }else{ d.z=-p.z; };
    double invr=1/(d.x+d.y+d.z);
    d.x*=invr; d.y*=invr; d.z*=invr;
    return i;
}

// ==================== MAIN

void my_draw(){

    //printf(" my_draw ! %i %i %g \n", plot1.lines.size(), plot1.lines.back()->n, plot1.lines.back()->ys[2]  );
    glClearColor( 1.0f, 1.0f, 1.0f, 0.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    //glDisable( GL_DEPTH_TEST );
    
    //plot1.view(false);

    glCallList(ogl);

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
    //DataLine2D * l     = plot1.add( new DataLine2D(nsamp,xmin,xmax, 0xFFFF0000 ) );
    DataLine2D * l1    = plot1.add( new DataLine2D(nsamp,xmin,xmax, 0xFFFF0000 ) );
    //DataLine2D * l2    = plot1.add( new DataLine2D(nsamp,xmin,xmax, 0xFF0000FF ) );

    //DataLine2D * le1   = plot1.add( new DataLine2D(nsamp,xmin,xmax, 0xFFFF8000 ) );
    //DataLine2D * le2   = plot1.add( new DataLine2D(nsamp,xmin,xmax, 0xFF0080FF ) );

    int ii;
    for(int i=0; i<nsamp; i++){
        double phi = lref->xs[i];
        double r = randf(0.1,5.);
        double x = cos( phi )*r;
        double y = sin( phi )*r;
        //printf(  "%i %g (%g,%g)\n", i, phi, x, y );
        //lref->ys[i] = (int)( fastFract(phi/(2*M_PI))*6.0 );
        lref->ys[i] = fastFract(3*phi/M_PI );
        //lref->ys[i] = fastFract( (3.0/M_PI)*atan2(y,x) );
        //l   ->ys[i] = hexface(y,x);
        //l1   ->ys[i] = hexcoord ( y, x, ii );
        l1   ->ys[i] = hexcoord_lin( y, x, ii );
        //l2   ->ys[i] = hexcoord_10( y, x, ii );

        //le1->ys[i] = log10( fabs( l1->ys[i] - lref->ys[i] ) );
        //le2->ys[i] = log10( fabs( l2->ys[i] - lref->ys[i] ) );
    }

    plot1.render();


    const int n = 1000;
    const int m = 10000;
    double xs [n];
    double ys [n];
    for(int i=0; i<n; i++){
        xs[i]=randf(-30.0,30.0);
        ys[i]=randf(-30.0,30.0);
        //double phi = randf(-30,30);
        //xs[i]=cos(phi);
        //ys[i]=sin(phi);
    }
    //VecN::arange(n,-30.0,60./n,xs);

    double c=0,s=0;
    double dn=1./n;

    TEST_ERROR_PROC_N( "hexcoord_lin ",{double x=xs[i]; double y=ys[i]; c=hexcoord_lin(y,x,ii); c-=fastFract( (3.0/M_PI)*atan2(y,x) ); STORE_ERROR(c) }, n );
    TEST_ERROR_PROC_N( "hexcoord_10  ",{double x=xs[i]; double y=ys[i]; c=hexcoord_10 (y,x,ii); c-=fastFract( (3.0/M_PI)*atan2(y,x) ); STORE_ERROR(c) }, n );

    SPEED_TEST_PROC_NM( "junk;        ",  {sum+=xs[i]; sum+=ys[i]; }, n, m );
    SPEED_TEST_PROC_NM( "atan2        " , {sum+= atan2   (ys[i],xs[i]); }, n, m );
    SPEED_TEST_PROC_NM( "hexface      " , {sum+= hexface (ys[i],xs[i]); }, n, m );
    SPEED_TEST_PROC_NM( "hexcoord_lin " , {sum+= hexcoord_lin(ys[i],xs[i], ii); }, n, m );
    SPEED_TEST_PROC_NM( "hexcoord_10  " , {sum+= hexcoord_10(ys[i],xs[i], ii); }, n, m );

    std::vector<Vec3d> ps{100000};
    SphereSamplesFibonachi(ps.size(),&ps[0]);

    ogl = glGenLists(1);
    glNewList(ogl,GL_COMPILE);
    glColor3f(0,0,0);
    glBegin( GL_POINTS );
    for(const Vec3d& p: ps ){
        //Vec2d d;
        Vec3d d;
        Vec2i ip; 
        //hexSphere_lin( p, d, ip );
        //hexSphere_10( p, d, ip );
        
        int i = octCoord( p, d );
        glColor3f (d.x,d.y,0.5);

        //Draw::color_of_hash( 1646454 + ip.x*6 + ip.y );

        //printf("%g,%g,%g \n", p.x,p.y,p.z );
        glVertex3f(p.x,p.y,p.z);
        Draw3D::vertex(p);
    }
    glEnd();
    //Draw3D::drawLine((Vec3d){0.0,0.0,0.0},{1,1,1});
    //Draw3D::drawPoints(ps.size(),&ps[0],-1);
    glEndList();

}

int main(){
    init( 800, 600 );
    set_draw_function( my_draw );
    setup();
    run_Nframes(5000);
}
