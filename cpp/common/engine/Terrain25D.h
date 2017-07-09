#ifndef Terrain25D_h
#define Terrain25D_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"

#include "Ruler2DFast.h"

class Terrain25D { public:
    // --- variables
	//Vec2d pos;
	//Vec2d a;
	//Vec2d b;
    double zmin = -1.0;
    double zmax =  1.0;
    double drag        = -0.5d;
    double restitution = -0.8d;
	int    shape;

	// --- function declarations
	virtual double eval( const Vec2d& pos, Vec2d& deriv );
	virtual double ray ( const Vec3d& hRay, const Vec3d& ray0, double tmax, Vec3d& normal );

	// --- inline function

    inline void boundingPlanes( double zRay, double z0, double tmax, double& tstart, double& tend,  double zmin, double zmax ){
        double tzmin  = (zmin-z0)/zRay;
        double tzmax  = (zmax-z0)/zRay;
        if(tzmin<tzmax){ tstart=fmin(0,tzmin); tend=fmax(tmax,tzmax); }
        else           { tstart=fmin(0,tzmax); tend=fmax(tmax,tzmin); };
    }

    inline void boundingPlanes( double zRay, double z0, double tmax, double& tstart, double& tend ){
        boundingPlanes( zRay, z0, tmax, tstart, tend, zmin, zmax );
    }
};

class Terrain25D_bicubic : public Terrain25D { public:
    // --- variables
    Ruler2DFast ruler;
    double * heights;
    // --- function declarations
	virtual double eval( const Vec2d& pos, Vec2d& deriv );
	//virtual double ray ( const Vec3d& hRay, const Vec3d& ray0, double tmax, Vec3d& normal );

	void allocate( Vec2i n ){ ruler.setN( n ); if(heights) delete [] heights; heights = new double[ruler.ntot]; };
	void makeRandom( double zmin_, double zmax_ ){
        zmin=zmin_; zmax=zmax_;
        for(int i=0; i<ruler.ntot; i++ ){ heights[i] = randf(zmin, zmax); }
	}

};

#endif
