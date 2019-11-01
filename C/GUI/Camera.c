
#ifndef  Camera_c
#define  Camera_c

#include "Camera.h"

#include "Vec3math.h"
#include "Mat3math.h"

/*
#define  TYPE  float  
#define  TPREF f3  
#define  VEC   float3 
#include "Vec3math_bare.h"

#define  TYPE float  
#define  TPREF f3x3  
#define  VEC  float3 
#define  MAT  float3x3 
#include "Mat3math_bare.h"
*/

Camera Camera_default(){
    Camera cam;
    cam.pos  = (float3){0.0f,0.0f,-50.0f};
    cam.rot  = f3x3_Eye;
    cam.zoom  = 10.0f;
    cam.aspect= 1.0;
    cam.zmin  = 10.0;
    cam.zmax  = 10000.0;
    return cam;
}

/*
    inline void lookAt( Vec3f p, float R ){ pos = p + rot.c*-R; }
    inline void lookAt( Vec3d p, float R ){ Vec3f p_; convert(p,p_); lookAt(p_,R); }

    inline float getTgX()const{ return 1.0/(zoom*aspect); }
    inline float getTgY()const{ return 1.0/(zoom);            }

    inline void word2screenOrtho( const Vec3f& pWord, Vec3f& pScreen ) const {
        Vec3f p; p.set_sub(pWord,pos);
        rot.dot_to( p, p );
        pScreen.x = p.x/(2*zoom*aspect);
        pScreen.y = p.y/(2*zoom);
        pScreen.z = (p.z-zmin)/(zmax-zmin);
    }

    inline Vec2f word2pixOrtho( const Vec3f& pWord, const Vec2f& resolution ) const {
        Vec3f p; p.set_sub(pWord,pos);
        rot.dot_to( p, p );
        return Vec2f{ resolution.x*(0.5f+p.x/(2*zoom*aspect)),
                      resolution.y*(0.5f+p.y/(2*zoom)) };
    }

    inline void word2screenPersp( const Vec3f& pWord, Vec3f& pScreen ) const {
        Vec3f p; p.set_sub(pWord,pos);
        rot.dot_to( p, p );
        float resc = zmin/(2*p.z*zoom);
        pScreen.x = p.x*resc/aspect;
        pScreen.y = p.y*resc;
        //pScreen.z = p.z/zmin;        // cz
        
        //(2*zmin)/w      0       0              0            0
        //0          (2*zmin)/h   0              0            0
        //0               0       (zmax+zmin)/(zmax-zmin)    (2*zmin*zmax)/(zmax-zmin)
        //0               0       0               0           -1
        //------
        //x_  =  ((2*zmin)/w)  * x
        //y_  =  ((2*zmin)/h ) * y
        //z_  =   (zmax+zmin)/(zmax-zmin)
    }

    inline Vec2f word2pixPersp( const Vec3f& pWord, const Vec2f& resolution ) const {
        Vec3f p; p.set_sub(pWord,pos);
        rot.dot_to( p, p );
        float resc = zmin/(2*p.z*zoom);
        return (Vec2f){
            resolution.x*( 0.5f + p.x*resc/aspect ),
            resolution.y*( 0.5f + p.y*resc        ) };
    }

    inline void pix2rayOrtho( const Vec2f& pix, Vec3f& ro ) const {
        float resc = 1/zoom;
        ro = rot.a*(pix.a*resc) + rot.b*(pix.b*resc);
    }

    inline Vec3f pix2rayPersp( const Vec2f& pix, Vec3f& rd ) const {
        float resc = 1/zoom;
        rd = rot.a*(pix.a*resc) + rot.b*(pix.b*resc);
    }

    inline bool pointInFrustrum( Vec3f p ) const {
        p.sub(pos);
        Vec3f c;
        rot.dot_to( p, c );
        float tgx = c.x*zoom*aspect;
        float tgy = c.y*zoom;
        float cz  = c.z*zmin;
        return (tgx>-cz)&&(tgx<cz) && (tgy>-cz)&&(tgy<cz) && (c.z>zmin)&&(c.z<zmax);
    }

    inline bool sphereInFrustrum( Vec3f p, float R ) const {
        p.sub(pos);
        Vec3f c;
        rot.dot_to( p, c );
        float my = c.z*zmin/zoom;
        float mx = my/aspect + R;  my+=R;
        return (c.x>-mx)&&(c.x<mx) && (c.y>-my)&&(c.y<my) && ((c.z+R)>zmin)&&((c.z-R)<zmax);
    }

*/

#endif