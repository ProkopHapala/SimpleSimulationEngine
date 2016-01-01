#ifndef AeroSurf2D_h
#define AeroSurf2D_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "fastmath.h"
#include "Vec2.h"

class RigidBody2D;

class KinematicBody2D {
	public:
	double phi;
	Vec2d  rot;
	Vec2d  pos;
	inline void setAngle( double phi_ ){  phi = phi_; rot.fromAngle( phi_ ); }
};

class AeroSurf2D : public KinematicBody2D {
	public:
	static constexpr double SAFETY_v = 1e-6;
	double area;  // [m^2]
	double CD0    = 0.02;  
	double dCD    = 0.9;  
	double dCDS   = 0.9;  
	double dCL    = 6.28;
	double dCLS   = 2.82743338823;
	double sStall = 0.16;
	double wStall = 0.08;
	
	int      nPolar;
	double * cDs;
	double * cLs;

	bool loadPolar    ( char const* filename );
	void plot_polar   ( double x0, double y0, double fscale, double phi0 );
	virtual void draw ( RigidBody2D& platform );
	double fromString ( char * s );
	char * toString   ( );
	void assertAeroForce( RigidBody2D& platform, const Vec2d& vel, double density );

	inline void polarModel( double ca, double sa, double& CD, double& CL ){
		//double abs_sa = fabs( sa ); 
		double abs_sa = (sa>0)?sa:-sa;
		double wS     = trashold_cub( abs_sa, sStall, sStall+wStall );
		double mS     = 1 - wS; 
		CD     = CD0 + ( mS*dCD*abs_sa + wS*dCDS        ) * abs_sa;
		//CL     =       ( mS*dCL        + wS*dCLS*abs_ca ) * sa;
		if( ca <0 ){ ca=-ca; sa=-sa; };
		CL     =       ( mS*dCL        + wS*dCLS*ca ) * sa;
	}

};

#endif
