
#ifndef AOrotations_h
#define AOrotations_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "Mat4.h"
#include "spline_hermite.h"

struct AOPairType{
    //double * Is = 0; // integrals
    double nI       = 0;
    double stepI    = 1;
    double invStepI = 1;
    double * Iss=0;
    double * Izz=0;
    double * Iyy=0;
    double * Isz=0;
    double * Izs=0;

    double dot(const Quat4d& c1, const Quat4d& c2, const Vec3d& hdir, double r ){
        //Vec3d up,side;
        // non-zero integrals are ss, sz, zs,  zz, yy, xx,
        Mat3d rot;
        rot.c=hdir;
        rot.c.getSomeOrtho( rot.b, rot.a );

        Quat4d cb1,cb2;
        rot.dot_to(c1.p,cb1.p);   cb1.s = c1.s;
        rot.dot_to(c2.p,cb2.p);   cb2.s = c2.s;

        double sr = r * invStepI;
        int    is = (int)sr;
        double ds = sr-is;

        double fss = Spline_Hermite::val( ds, Iss+is );
        double fsz = Spline_Hermite::val( ds, Isz+is );
        double fzs = Spline_Hermite::val( ds, Izs+is );
        double fzz = Spline_Hermite::val( ds, Izz+is );
        double fyy = Spline_Hermite::val( ds, Iyy+is );

        return fss*cb1.s  *cb2.s   + fsz*cb1.s  *cb2.p.z  + fzs*cb1.p.z*cb2.s +
               fzz*cb1.p.z*cb2.p.z + fyy*cb1.p.y*cb2.p.y  + fyy*cb1.p.x*cb2.p.x;

    }

    void outer(const Vec3d& hdir, double r, Mat4d& M ){
        //Vec3d up,side;
        // non-zero integrals are ss, sz, zs,  zz, yy, xx,
        //Mat3d rot;
        //rot.c=hdir;
        //rot.c.getSomeOrtho( rot.b, rot.a );



        //Quat4d cb1,cb2;
        //rot.dot_to(c1.p,cb1.p);   cb1.s = c1.s;
        //rot.dot_to(c2.p,cb2.p);   cb2.s = c2.s;

        double sr = r * invStepI;
        int    is = (int)sr;
        double ds = sr-is;

        double fss = Spline_Hermite::val( ds, Iss+is );
        double fsz = Spline_Hermite::val( ds, Isz+is );
        double fzs = Spline_Hermite::val( ds, Izs+is );
        double frr = Spline_Hermite::val( ds, Izz+is );
        double fpp = Spline_Hermite::val( ds, Iyy+is );

        M.s.s  = fss;

        M.s.x  = fsz * hdir.x;
        M.s.y  = fsz * hdir.y;
        M.s.z  = fsz * hdir.z;

        M.px.s = fzs * hdir.x;
        M.py.s = fzs * hdir.y;
        M.pz.s = fzs * hdir.z;

        Vec3d ha,hb;
        hdir.getSomeOrtho(ha,hb);

        Mat3Sd Uc,Uab;
        //Ua.from_outer(ha);
        //Ub.from_outer(hb);
        Uc .from_outer(hdir);
        Uab.from_outer(ha);
        Uab.add_outer (hb);

        M.px.x = Uc.xx * frr  +  Uab.xx * fpp;
        M.px.y = Uc.xy * frr  +  Uab.xy * fpp;
        M.px.z = Uc.xz * frr  +  Uab.xz * fpp;

        M.py.x = Uc.yx * frr  +  Uab.yx * fpp;
        M.py.y = Uc.yy * frr  +  Uab.yy * fpp;
        M.py.z = Uc.yz * frr  +  Uab.yz * fpp;

        M.pz.x = Uc.zx * frr  +  Uab.zx * fpp;
        M.pz.y = Uc.zy * frr  +  Uab.zy * fpp;
        M.pz.z = Uc.zz * frr  +  Uab.zz * fpp;

    }

};

#endif



