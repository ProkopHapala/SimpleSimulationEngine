#ifndef Terrain25D_h
#define Terrain25D_h

#include "fastmath.h"
#include "Vec2.h"

class Terrain25D {
	public:

	Vec2d pos;
	Vec2d a;
	Vec2d b;
    double zmin = -1.0;
    double zmax =  1.0;

	int    shape;

    inline void boundingPlanes( double zRay, double z0, double tmax, double& tstart, double& tend,  double zmin, double zmax ){
        double tzmin  = (zmin-z0)/zRay;
        double tzmax  = (zmax-z0)/zRay;
        if(tzmin<tzmax){
            tstart=fmin(0,tzmin); tend=fmax(tmax,tzmax);
        }else{
            tstart=fmin(0,tzmax); tend=fmax(tmax,tzmin);
        };
    }
    inline void boundingPlanes( double zRay, double z0, double tmax, double& tstart, double& tend ){
        boundingPlanes( zRay, z0, tmax, tstart, tend, zmin, zmax );
    }


    virtual double eval(const Vec2d& pos, Vec2d& deriv ){
        constexpr double scx = 0.1;
        constexpr double scy = 0.2;
        double cax = cos(scx*pos.x);
        double sax = sin(scx*pos.x);
        double cay = cos(scy*pos.y);
        double say = sin(scy*pos.y);
        deriv.set( scx*cax + say, sax + scy*cay );
        return sax + say;
    };

    virtual double ray( const Vec3d& hRay, const Vec3d& ray0, double tmax, Vec3d& normal ){
        constexpr double scx = 0.1;
        constexpr double scy = 0.2;

        // bounding planes
        //if( hRay.z == 0 ) && (ray0.z)
        double tstart,tend;
        boundingPlanes( hRay.z, ray0.z, tmax, tstart, tend );
        /*
        double tzmin  = (zmin-ray0.z)/hRay.z;
        double tzmax  = (zmax-ray0.z)/hRay.z;
        double tstart,tend;
        if(tzmin<tzmax){ tstart=fmin(0,tzmin); tend=fmax(tmax,tzmax);  }else{ tstart=fmin(0,tzmax); tend=fmax(tmax,tzmin);  };
        */

        double tspan = (tend-tstart);
        Vec3d p; p.set_lincomb(1,ray0, tspan,hRay);
        double dt = fabs( (hRay.z+(scx+scy)*0.72) );
        Vec3d dRay; dRay.set_mul( hRay, dt );
        int n = tspan/dt;
        Vec2d deriv;
        double old_val = eval( {p.x, p.y}, deriv );
        for(int i=0; i<n; i++){
            p.add(dRay);
            double val = eval( {p.x, p.y}, deriv );
            if( val*old_val < 0 ){
                normal.set( deriv.x, deriv.y, -1 );
                return tstart+i*dt;
            }
        }
        return 1e+300;
    };

};

#endif
