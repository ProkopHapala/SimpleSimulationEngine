#ifndef Fluid2D_h
#define Fluid2D_h

#include <vector>
#include "Noise.h"
#include "arrayAlgs.h"
#include "GridIndex2D.h"
//#include "Grid2D.h"
#include "Grid2DAlgs.h"

//#include "Fluid2D.cpp"

/*
References:
http://www.dgp.toronto.edu/people/stam/reality/index.html
http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/GDC03.pdf
http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/ns.pdf
/home/prokop/Dropbox/MyDevSW/Processing/AlgorithmicArtGeneration/PlasticOrogennyParticles
*/

inline void acum( int n, double* a, double* b, double dt){ for(int i=0;i<n;i++){  a[i] += b[i] * dt; } };
inline void setarr(int n, double* a, double val){ for(int i=0;i<n;i++){  a[i]=val; } }
//void swap_ptr(void*& a, void*& b){ void*& tmp=a; a=b; b=tmp; };

class Fluid2D : public Grid2DAlg { public:

    //int ndiffuse  = 20;
    //int npressure = 20;
    int ndiffuse  = 5;
    int npressure = 5;
    double	visc = 0.0001*4;
	double	diff = 3.0;
	double N_INVERSE = (double)1.0/(double)ntot;
    double dx,dy;
    double dt;

	// QUESTION :   would it be easier to use double * vx,vy instead of Vec2d* v ????
	//Vec2d  *vel=0,*vel_= 0;
	double *vx=0,*vy=0,*vx_=0,*vy_=0;
	double *div,*p;   // not really neede, can be relplaxed by vy_,vx_
	double *dens=0,*dens_=0;
    double *source=0;

    // ===== Function Declarations

    void allocate( Vec2i ns_ );
    void corners          ( double* u);
    void boundary_zero    ( double* u);
    void boundary_reflect ( double* u);
    void boundary_absorb  ( double* u);
    void boundary_periodic( double* u);
    void set_bnd     ( int b, double* u);
    void advect      ( double* u, double* u0, double* vx, double* vy, double dt);
    void diverg      ( double* div, double* p, double* vx, double* vy );
    void pressureBlur( double* div, double* p );
    void accelerate  ( double* p, double* vx, double* vy, double dt );
    void diffuse     ( double* u, double* u0, double diff, double dt );
    void fluidStep   ( double dt);

    // ===== inline Functions

    double interpBilinear( Vec2d p, double* u ){
        double x=_clamp(p.x, 0.5, n.x-1.5);  // boundary condition
        double y=_clamp(p.y, 0.5, n.y-1.5);
        int ix0 = (int)x;     //int ix1 = ix+1;
        int iy0 = (int)y;     //int iy1 = iy+1;
        int i0  = ip2i( {ix0,iy0} );
        int i1  = i0 + n.x;
        double fx = x-ix0;  double mx = 1-fx;
        double fy = y-iy0;  double my = 1-fy;
        return my * (mx*u[i0] + fx*u[i0+1])
             + fy * (mx*u[i1] + fx*u[i1+1]);
    }

};

#endif
