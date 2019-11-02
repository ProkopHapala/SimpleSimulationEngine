
#ifndef  Draw3D_Surf_h
#define  Draw3D_Surf_h

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "Vec2.h"
#include "Draw.h"

#include <math.h>
#include <cstdlib>
#include <stdio.h>
#include <stdint.h>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

#include "Mesh.h"
#include "Draw3D.h"

//#include <SDL2/SDL_opengl.h>

namespace Draw3D{

template<void colorFunc(double f)>
void drawTriaglePatch( Vec2i i0, Vec2i n, int NX, double * height, double vmin, double vmax ){
    Vec2f a,b,p;
    a.set( 1.0d, 0.0d           ); //a.mul(scale);
    b.set( 0.5d, 0.86602540378d ); //b.mul(scale);
    //glDisable(GL_SMOOTH);
    //int ii = 0;
    double renorm=1.0d/(vmax-vmin);
    for (int iy=0; iy<n.y-1; iy++){
        glBegin( GL_TRIANGLE_STRIP );
        int ii = (i0.y+iy)*NX + i0.x;
        for (int ix=0; ix<n.x; ix++){
            p.set( ix*a.x+iy*b.x, ix*a.y+iy*b.y );
            colorFunc( (height[ii     ]-vmin)*renorm ); glVertex3f( p.x    , p.y    , 0 );
            colorFunc( (height[ii + NX]-vmin)*renorm ); glVertex3f( p.x+b.x, p.y+b.y, 0 );
            ii++;
        }
        glEnd();
    }
}

template<typename UVfunc>
Vec3f getUVFuncNormal( Vec2f uv, float eps,  UVfunc func ){
    Vec2f o;
    Vec3f nor,da,db;
    /*
    o=uv; o.a+=eps; if(o.a>1)o.a=1; da.set(func(o));
    o=uv; o.a-=eps; if(o.a<0)o.a=0; da.sub(func(o));
    o=uv; o.b+=eps; if(o.b>1)o.b=1; db.set(func(o));
    o=uv; o.b-=eps; if(o.b<0)o.b=0; db.sub(func(o));
    */
    o=uv; o.a+=eps; da.set(func(o));
    o=uv; o.a-=eps; da.sub(func(o));
    o=uv; o.b+=eps; db.set(func(o));
    o=uv; o.b-=eps; db.sub(func(o));
    nor.set_cross(db,da); nor.normalize();
    return nor;
}

template<typename UVfunc>
void drawSmoothUVFunc( Vec2i n, Vec2f UVmin, Vec2f UVmax, float voff, UVfunc func ){

    float eps = 0.001;
    Vec2f duv = UVmax-UVmin; duv.mul( {1.0f/n.a,1.0f/n.b} );
    //int i0=mesh.vpos.size()/3;
    printf( "n (%i,%i) duv (%f,%f) \n", n.a, n.b, duv.a, duv.b  );
    Vec2f uv;
    for(int ia=0;ia<=n.a;ia++){
        glBegin(GL_TRIANGLE_STRIP);
        for(int ib=0;ib<=n.b;ib++){
            //int i=mesh.vpos.size();
            uv.a = UVmin.a+duv.a*ia;
            uv.b = UVmin.b+duv.b*ib + voff*duv.b*ia;
            Vec3f p,nr;

            //printf( " %i %i: uv (%f,%f) p (%f,%f,%f)\n", ia,ib, uv.a, uv.b,  p.x,p.y,p.z  );
            p  = func(uv);
            nr = getUVFuncNormal(uv,eps,func);
            glNormal3f(nr.x,nr.y,nr.z);
            glVertex3f(p.x,p.y,p.z);

            uv.a += duv.a;
            uv.b += voff*duv.b;

            p  = func(uv);
            nr = getUVFuncNormal(uv,eps,func);
            glNormal3f(nr.x,nr.y,nr.z);
            glVertex3f(p.x,p.y,p.z);

        }
        glEnd();
    }

}

template<typename UVfunc>
void drawWireUVFunc( Vec2i n, Vec2f UVmin, Vec2f UVmax, float voff, UVfunc func ){

    float eps = 0.001;
    Vec2f duv = UVmax-UVmin; duv.mul( {1.0f/n.a,1.0f/n.b} );
    //int i0=mesh.vpos.size()/3;
    Vec2f uv;
    for(int ia=0;ia<=n.a;ia++){
        // strips
        glBegin(GL_LINE_LOOP);
        for(int ib=0;ib<=n.b;ib++){
            uv.a = UVmin.a+duv.a*ia;
            uv.b = UVmin.b+duv.b*ib +  voff*duv.b*ia;
            Vec3f p;
            p  = func(uv);
            glVertex3f(p.x,p.y,p.z);

        }
        glEnd();

        //
        if(ia<n.a){
            glBegin(GL_LINE_LOOP);
            for(int ib=0;ib<=n.b;ib++){
                uv.a = UVmin.a+duv.a*ia;
                uv.b = UVmin.b+duv.b*ib + voff*duv.b*ia;
                Vec3f p;

                p  = func(uv);
                glVertex3f(p.x,p.y,p.z);

                uv.a += duv.a;
                uv.b += voff*duv.b;

                p  = func(uv);
                glVertex3f(p.x,p.y,p.z);

            }
            glEnd();
        }
    }
}


template<typename UVfunc>
void drawExtrudedWireUVFunc( Vec2i n, float thick, Vec2f UVmin, Vec2f UVmax, float voff, UVfunc func ){

    float eps = 0.001;
    Vec2f duv = UVmax-UVmin; duv.mul( {1.0f/n.a,1.0f/n.b} );
    //int i0=mesh.vpos.size()/3;
    Vec2f uv;
    Vec3f p,nr,nrl;
    for(int ia=0;ia<=n.a;ia++){
        // strips
        glBegin(GL_TRIANGLE_STRIP);
        for(int ib=0;ib<=n.b;ib++){
            uv.a = UVmin.a+duv.a*ia;
            uv.b = UVmin.b+duv.b*ib + voff*duv.b*ia;
            p  = func(uv);
            nr = getUVFuncNormal(uv,eps,func);
            nrl.set_cross(func(uv+(Vec2f){0.0,duv.b*0.5})-p, nr); nrl.normalize();
            glNormal3f(nrl.x,nrl.y,nrl.z);
            glVertex3f(p.x,p.y,p.z);
            p.add_mul(nr,thick);
            glVertex3f(p.x,p.y,p.z);
        }
        glEnd();
        if(ia<n.a){
            glBegin(GL_TRIANGLE_STRIP);
            for(int ib=0;ib<=n.b;ib++){
                uv.a = UVmin.a+duv.a*ia;
                uv.b = UVmin.b+duv.b*ib + voff*duv.b*ia;

                p  = func(uv);
                nr = getUVFuncNormal(uv,eps,func);
                nrl.set_cross(func(uv+(Vec2f){0.0,duv.b*0.5})-p, nr); nrl.normalize();
                glNormal3f(nrl.x,nrl.y,nrl.z);
                glVertex3f(p.x,p.y,p.z);
                p.add_mul(nr,thick);
                glVertex3f(p.x,p.y,p.z);

                uv.a += duv.a;
                uv.b += voff*duv.b;

                p  = func(uv);
                nr = getUVFuncNormal(uv,eps,func);
                nrl.set_cross(func(uv+(Vec2f){0.0,duv.b*0.5})-p, nr); nrl.normalize();
                glNormal3f(nrl.x,nrl.y,nrl.z);
                glVertex3f(p.x,p.y,p.z);
                p.add_mul(nr,thick);
                glVertex3f(p.x,p.y,p.z);

            }
            glEnd();
        }
    }
}

inline Vec3f HarmonicTubeUVfunc( Vec2f p, float R1, float R2, float L, float freq, float amp ){
    Vec2f csb; csb.fromAngle(p.b);
    float R = ((1-p.a)*R1 + p.a*R2)*(1.0+amp*cos(freq*p.a));
    //float R = (1-p.a)*R1 + p.a*R2;
    return (Vec3f){ csb.a*R, csb.b*R, L*p.a };
}

inline Vec3f ConeUVfunc( Vec2f p, float R1, float R2, float L ){
    Vec2f csb; csb.fromAngle(p.b);
    float R = (1-p.a)*R1 + p.a*R2;
    return (Vec3f){csb.a*R,csb.b*R,L*p.a };
}

inline Vec3f SphereUVfunc( Vec2f p, float R ){
    Vec2f csa; csa.fromAngle(p.a);
    Vec2f csb; csb.fromAngle(p.b);
    return (Vec3f){csa.a*csb.a*R,csa.a*csb.b*R,csa.b*R };
}

inline Vec3f TorusUVfunc( Vec2f p, float r, float R ){
    Vec2f csa; csa.fromAngle(p.a);
    Vec2f csb; csb.fromAngle(p.b);
    return (Vec3f){csb.a*(R+r*csa.a),csb.b*(R+r*csa.a),r*csa.b};
}

inline Vec3f TeardropUVfunc( Vec2f p, float R1, float R2, float L ){
    Vec2f csa; csa.fromAngle(p.a*M_PI-M_PI*0.5);
    Vec2f csb; csb.fromAngle(p.b);
    float f =  0.5-csa.b*0.5;
    float R = (1-f)*R1 + f*R2;
    return (Vec3f){csa.a*csb.a*R,csa.a*csb.b*R,csa.b*R-L*f };
}

inline Vec3f ParabolaUVfunc( Vec2f p, float K ){
    Vec2f csb; csb.fromAngle(p.b);
    float r = p.a;
    float l = p.a*p.a*K;
    return (Vec3f){csb.a*r,csb.b*r,l };
}

inline Vec3f HyperbolaRUVfunc( Vec2f p, float R, float K ){
    Vec2f csb; csb.fromAngle(p.b);
    float r = p.a;
    float l = K*sqrt( p.a*p.a + R*R); // - K*R;
    return (Vec3f){csb.a*r,csb.b*r,l };
}

inline Vec3f HyperbolaLUVfunc( Vec2f p, float R, float K ){
    Vec2f csb; csb.fromAngle(p.b);
    float l = p.a;
    float r = K*sqrt( p.a*p.a + R*R);
    return (Vec3f){csb.a*r,csb.b*r,l };
}

inline void drawUV_HarmonicTube( Vec2i n, Vec2f UVmin, Vec2f UVmax, float R1, float R2, float L, float voff, float freq, float amp, bool wire ){
    auto uvfunc = [&](Vec2f uv){return HarmonicTubeUVfunc(uv,R1,R2,L,freq,amp);};
    if(wire){ drawWireUVFunc( n, UVmin, UVmax, voff, uvfunc ); }
    else    { drawSmoothUVFunc( n, UVmin, UVmax,voff, uvfunc ); }
}

inline void drawUV_Cone( Vec2i n, Vec2f UVmin, Vec2f UVmax, float R1, float R2, float L, float voff, bool wire ){
    auto uvfunc = [&](Vec2f uv){return ConeUVfunc(uv,R1,R2,L);};
    if(wire){ drawWireUVFunc( n, UVmin, UVmax, voff, uvfunc ); }
    else    { drawSmoothUVFunc( n, UVmin, UVmax,voff, uvfunc ); }
}

inline void drawUV_Sphere( Vec2i n, Vec2f UVmin, Vec2f UVmax, float R, float voff, bool wire ){
    auto uvfunc = [&](Vec2f uv){return SphereUVfunc(uv,R);};
    if(wire){ drawWireUVFunc  ( n, UVmin, UVmax, voff, uvfunc ); }
    else    { drawSmoothUVFunc( n, UVmin, UVmax,voff, uvfunc ); }
}

inline void drawUV_Torus( Vec2i n, Vec2f UVmin, Vec2f UVmax, float r, float R, float voff, bool wire ){
    auto uvfunc = [&](Vec2f uv){return TorusUVfunc(uv,r,R);};
    if(wire){ drawWireUVFunc( n, UVmin, UVmax,voff, uvfunc ); }
    else    { drawSmoothUVFunc( n, UVmin, UVmax,voff, uvfunc ); }
}

inline void drawUV_Teardrop( Vec2i n, Vec2f UVmin, Vec2f UVmax, float R1, float R2, float L, float voff, bool wire ){
    auto uvfunc = [&](Vec2f uv){return TeardropUVfunc(uv,R1,R2,L);};
    if(wire){ drawWireUVFunc( n, UVmin, UVmax,voff, uvfunc ); }
    else    { drawSmoothUVFunc( n, UVmin, UVmax,voff, uvfunc ); }
}

inline void drawUV_Parabola( Vec2i n, Vec2f UVmin, Vec2f UVmax, float R, float L, float voff, bool wire ){
    float K = L/(R*R);
    UVmin.a*=R; UVmax.a*=R;
    //printf( "drawUV_Parabola: R %f L %f K %f \n", R, L, K  );
    auto uvfunc = [&](Vec2f uv){return ParabolaUVfunc(uv,K);};
    if(wire){ drawWireUVFunc( n, UVmin, UVmax,voff, uvfunc ); }
    else    { drawSmoothUVFunc( n, UVmin, UVmax,voff, uvfunc ); }
}

inline void drawUV_Parabola_ExtrudedWire( Vec2i n,Vec2f UVmin, Vec2f UVmax, float R, float L, float voff, double thick ){
    float K = L/(R*R);
    UVmin.a*=R; UVmax.a*=R;
    auto uvfunc = [&](Vec2f uv){return ParabolaUVfunc(uv,K);};
    drawExtrudedWireUVFunc( n,thick, UVmin, UVmax,voff, uvfunc );
}

inline void drawUV_Hyperbola( Vec2i n, Vec2f UVmin, Vec2f UVmax, float r, float R, float L, float voff, bool wire ){
    //printf( "drawUV_Hyperbola: r %f R %f L %f \n", r, R, L  );
    if(r>0){
        float K = R/L;
        UVmin.a*=L; UVmax.a*=L;
        auto uvfunc = [&](Vec2f uv){return HyperbolaLUVfunc(uv,r,K);};
        if(wire){ drawWireUVFunc( n, UVmin, UVmax,voff, uvfunc ); }
        else    { drawSmoothUVFunc( n, UVmin, UVmax,voff, uvfunc ); }
    }else{
        r=-r;
        float K = L/R;
        UVmin.a*=R; UVmax.a*=R;
        auto uvfunc = [&](Vec2f uv){return HyperbolaRUVfunc(uv,r,K);};
        if(wire){ drawWireUVFunc( n, UVmin, UVmax,voff, uvfunc ); }
        else    { drawSmoothUVFunc( n, UVmin, UVmax,voff, uvfunc ); }
    }
}

}; // namespace Draw3D


#endif

