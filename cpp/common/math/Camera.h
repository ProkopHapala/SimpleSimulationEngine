
#ifndef  Camera_h
#define  Camera_h

#include <math.h>
#include <cstdlib>
#include <stdint.h>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

template<typename T>
class CameraT{ public:
    Vec3T<T>  pos    = (Vec3T<T>){0.0f,0.0f,-50.0f};
    Mat3T<T>  rot    = (Mat3T<T>)Mat3dIdentity;
    T zoom    = 10.0;
    T  aspect = 1.0;
    T  zmin   = 10.0;
    T  zmax   = 10000.0;
    bool   persp  = true;

    inline void lookAt( Vec3T<T> p, T R ){ pos = p + rot.c*-R; }
    //inline void lookAt( Vec3T<T> p, T R ){ Vec3T<T> p_; convert(p,p_); lookAt(p_,R); }

    inline T getTgX()const{ return 1.0/(zoom*aspect); }
    inline T getTgY()const{ return 1.0/(zoom);            }

    inline void word2screenOrtho( const Vec3T<T>& pWord, Vec3T<T>& pScreen ) const {
        Vec3T<T> p; p.set_sub(pWord,pos);
        rot.dot_to( p, p );
        pScreen.x = p.x/(2*zoom*aspect);
        pScreen.y = p.y/(2*zoom);
        pScreen.z = (p.z-zmin)/(zmax-zmin);
    }

    inline Vec2T<T> word2pixOrtho( const Vec3T<T>& pWord, const Vec2f& resolution ) const {
        Vec3T<T> p; p.set_sub(pWord,pos);
        rot.dot_to( p, p );
        return Vec2f{ resolution.x*(0.5+p.x/(2*zoom*aspect)),
                      resolution.y*(0.5+p.y/(2*zoom)) };
    }

    inline void word2screenPersp( const Vec3T<T>& pWord, Vec3T<T>& pScreen ) const {
        Vec3T<T> p; p.set_sub(pWord,pos);
        rot.dot_to( p, p );
        T  resc = zmin/(2*p.z*zoom);
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
    }

    inline Vec2T<T> word2pixPersp( const Vec3T<T>& pWord, const Vec2f& resolution ) const {
        Vec3T<T> p; p.set_sub(pWord,pos);
        rot.dot_to( p, p );
        T  resc = zmin/(2*p.z*zoom);
        return (Vec2f){
            resolution.x*( 0.5 + p.x*resc/aspect ),
            resolution.y*( 0.5 + p.y*resc        ) };
    }

    inline void pix2rayOrtho( const Vec2f& pix, Vec3T<T>& ro ) const {
        //T  resc = 1/zoom;
        T  resc = zoom;
        ro = rot.a*(pix.a*resc) + rot.b*(pix.b*resc);
    }

    inline void pix2rayPersp( const Vec2f& pix, Vec3T<T>& rd ) const {
        T  resc = 1/zoom;
        rd = rot.a*(pix.a*resc) + rot.b*(pix.b*resc);
    }

    //inline Vec3T<T> pix2ray( const Vec2f& pix, Vec3T<T>& rd, Vec3T<T>& ro ){
    inline void pix2ray( const Vec2f& pix, Vec3T<T>& rd, Vec3T<T>& ro ){
        if(persp){
            ro = pos;
            pix2rayPersp( pix, rd );
        }else{
            rd = rot.c;
            pix2rayOrtho( pix, ro );
        }
    }

    inline bool pointInFrustrum( Vec3T<T> p ) const {
        p.sub(pos);
        Vec3T<T> c;
        rot.dot_to( p, c );
        T  tgx = c.x*zoom*aspect;
        T  tgy = c.y*zoom;
        T  cz  = c.z*zmin;
        return (tgx>-cz)&&(tgx<cz) && (tgy>-cz)&&(tgy<cz) && (c.z>zmin)&&(c.z<zmax);
    }

    inline bool sphereInFrustrum( Vec3T<T> p, T  R ) const {
        p.sub(pos);
        Vec3T<T> c;
        rot.dot_to( p, c );
        T  my = c.z*zmin/zoom;
        T  mx = my/aspect + R;  my+=R;
        return (c.x>-mx)&&(c.x<mx) && (c.y>-my)&&(c.y<my) && ((c.z+R)>zmin)&&((c.z-R)<zmax);
    }

};

using Camera   = CameraT<float>;
using Camera_d = CameraT<double>;

#endif

