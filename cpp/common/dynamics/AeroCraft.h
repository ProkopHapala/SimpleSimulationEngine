
#ifndef AeroCraft_h
#define AeroCraft_h

#include "fastmath.h"
#include "Body.h"
#include "AeroSurf.h"

//#include <SDL2/SDL_opengl.h>
//#include <Draw3D.h>

class Propeler{
    public:
    Vec3d lpos;
    Vec3d dir;
    double area;
    double power;
    double efficiency;
    double CD;
    double vstatic; // velocity of stream for properler stationary with respect to air

    inline double getThrust( double v0 ){
        // derived from those equations
        // dm = S * v0
        // F = dm * Dv
        // P = F * v0 + 0.5*dm*(Dv**2) = S*(v0**2)*Dv + 0.5*S*v0*(Dv**2)
        //printf( "v0 %g \n", v0 );
        double density = 1.22;
        double dm = density*area*(v0+vstatic);
        double a  = 0.5*dm;
        double b  = dm*v0;
        double c  = -power;
        double Dv1,Dv2;
        quadratic_roots( a, b, c,  Dv1, Dv2 );
        return dm * Dv2 * efficiency - dm*v0*CD;
    }

    inline void fromString( const char * str ){
        sscanf ( str, " %lf %lf %lf    %lf %lf %lf    %lf %lf %lf %lf",
                    &lpos.x, &lpos.y, &lpos.z,
                    &dir.x, &dir.y, &dir.z,
                    &area, &power, &efficiency, &CD
                );
        dir.normalize();
        double density = 1.22;
        vstatic = pow(power/(density*area),0.333333);
        printf ( "%lf %lf %lf %lf\n", area, power, efficiency, CD );
    }

};

class AeroCraft : public RigidBody { public:

	int nPanels = 0, nPropelers= 0;
	AeroSurface * panels      =NULL;
	Propeler    * propelers   =NULL;

	AeroSurface * leftAirelon =NULL;
	AeroSurface * rightAirelon=NULL;
	AeroSurface * elevator    =NULL;
	AeroSurface * rudder      =NULL;

	Vec3d totalThrust;

	// ==== function declarations

	//virtual void render();
	int fromFile( const char * fname );
	void applyAeroForces( const Vec3d& vwind );


	double getTotalPower() const;

	// ==== inline functions

};


#endif
