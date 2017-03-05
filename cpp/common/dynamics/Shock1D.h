
#ifndef Shock1D_h
#define Shock1D_h

#include "fastmath.h"


/*
 Simulation of shock waves in layered system ( 1 infinite planar 2 cylindrical 3 spherical symmetry )
  - simulation of implosion in atomic bombs and 
  - shaped charges
  - explosions and ablation acceleration - including nuclear
*/




/*
EOS update
 - we should compute update pressure in two situations
   - (1) adiabatic update - volume of cell changed without exchange of heat (or particles)
   - (2) isochoric update - heat was created in cell (e.g. by reaction, or radiation heating/cooling) without change of volume
- but in general case particle number and gamma (heat capacities) may be modified due to ionization of atoms etc.
 => we need to compute internal energy and redistribute it over the new particles
   pV = nRT = NkT = U/alpha  
   - isochoric pressure depends only on total internal energy, if the same ammount of energy is distributed over larger number of particles pressure is unchanged
   - but some part of energy (latent heat) is consumet for that phase transition (e.g. ionization) 
*/

#define SYMMETRY   3

class ShockMaterial{
    public:
    double dens0       = 1.0;
    double p0          = 0.0;
    double BulkModulus = 1e+9;
    double alpha       = 3.0/2.0; // number of DOFs divided by 2;   alpha = 1/(gamma-1)
    double gamma       = 5.0/3.0; // heat capacity ratio  // https://en.wikipedia.org/wiki/Heat_capacity_ratio1
    
    // Stiffened equation of state // http://www.sklogwiki.org/SklogWiki/index.php/Stiffened_equation_of_state
    //inline double p_addiabatic( double density ){   
    //    return pow( density/ShockMaterial->dens0, gamma ); 
    //};
    
};

class ShockVolume{
    public:
    ShockMaterial * material;
    
    double mass;
    
    // internal energy
    double V;
    
    // ---- updates using pressure
    double p;
    
    inline double isochoricUpdate( double Q ){
        // pV=nRT; U=apha*nRT; U = alpha * pV 
        // U_ = (U + Q);   p_ = U_/(V*apha)
        double Va = V * material->alpha;
        p = (p*Va + Q)/Va;
        return p; 
    } 
    
    inline void addiabaticUpdate( double V_ ){
        // double W = p*V*( pow(V_/V, 1-gamma ) - 1 )
        p = p*pow( V_/V, material->gamma );
        V = V_;
    }
    
    // https://en.wikipedia.org/wiki/Equation_of_state
    inline double getPressure( double V_ ){
        addiabaticUpdate( V_ );
        double pressure = p + material->p0;   // cohessive pressure
        //pressure += degeneracyPressure( V_ );
        //-- adiabatic   https://en.wikipedia.org/wiki/Adiabatic_process
        //pressure += Cadiabat*pow(density,gamma);
        //-- solid      https://en.wikipedia.org/wiki/Bulk_modulus
        //pressure +=   BulkModulus*(density-dens0); 
        return pressure;
    }
    
    /*
     // ---- updates using internal energy
    double U; 
    inline isochoricUpdate( double Q ){
        double U_ = U + Q;
    }
    inline double addiabaticUpdate( double V_ ){
        double p = U/( V * material->alpha );    // pV = nRT = NkT = U/alpha;  U = alpha * nRT 
        K  = pow( V, gamma );                    // 
        double W = p;                            // http://hyperphysics.phy-astr.gsu.edu/hbase/thermo/adiab.html
        U += W; 
        V = V;
        return p_;
    }
    */
    
    inline double update( double dt ){
        // there can be e.g. chemical or nucler reactions within cell    
    }
    
    inline void init( double V_ ){
        double V    = V;
        double mass = material->dens0 * V;
    }
};

class ShockSystem1D{
	public:
	// params
	//constexpr symmetry = 1; // 0 - plane;  1 - cylinder;  3  - sphere 
	
	double outer_pressure = 0.0;
	
	int       nlayers;
	//double  * mass;
	//double    * imass;
	//int     * materials;   
	ShockVolume * cells; 
	
	//constexpr nMatMax;
	//int       nMat;    
	//ShockMaterial  matData[nMatMax];  // contains parameters for equations of state
	
	// state varibles
	double  * bounds     = NULL;  // n (n+1 but one is in 0)
	double  * bforce     = NULL;  // ----,,-------
	//double  * pos      = NULL;  // or boundary ?
	double  * velocity   = NULL;  // of COG or boundary ?

	// axuliary - should be stored or just computed on-the-fly ? maybe just for debugging?  
	//double  * density    = NULL;  // or volume ?
	//double  * E          = NULL;

	//double  * Temperature;
	//double  * volume;
	
	// ==== function declaration
	
	void evalForce();
	void move     (double dt);
	void update   (double dt);
	
	void allocate( int nlayers_ );
	void init();

	// ==== inline function implementation
	
    inline void boundProperties( double r, double& volume, double& area ){
        switch(SYMMETRY){
            case( 3 ): {          //  sphere
                double r2 = r*r;
                volume = 4.18879020479*r2*r;  
                area   = 12.5663706144*r2;   
                break;}
            case( 2 ): {          //  cylinder
                volume = 3.14159265359*r*r;  
                area   = 6.28318530718*r;
                break;  }   
            default:  {        //  1 plane
                area   = 1;
                volume = r;
            } 
        } 
    }
    
    inline double get_dR( double r1, double r2, double& dr1, double& dr2 ){
        switch(SYMMETRY){
            case 3: {          //  sphere
                // Sphere dcog/dr1:   3*r1**3/(r1**3 - r2**3) - 9*r1**2*(r1**4 - r2**4)/(4*(r1**3 - r2**3)**2)
                // Sphere dcog/dr2:  -3*r2**3/(r1**3 - r2**3) + 9*r2**2*(r1**4 - r2**4)/(4*(r1**3 - r2**3)**2)
                double rr2 = r2*r2; double rrr2=rr2*r2;
                double rr1 = r1*r1; double rrr1=rr1*r1;
                double im3 = 1/(rrr1 - rrr2);
                double m4  = rr1*rr1 - rr2*rr2;
                dr1   = ( 3*rrr1 -  2.25*rr1*m4*im3)*im3;  
                dr2   = (-3*rrr2 +  2.25*rr2*m4*im3)*im3;
                }
            case 2: {          //  cylinder
                //Cylinder dcog/dr1:  2*r1**2/(r1**2 - r2**2) - 4*r1*(r1**3 - r2**3)/(3*(r1**2 - r2**2)**2)
                //Cylinder dcog/dr2:  -2*r2**2/(r1**2 - r2**2) + 4*r2*(r1**3 - r2**3)/(3*(r1**2 - r2**2)**2
                double rr2 = r2*r2;
                double rr1 = r1*r1;
                double im2 = 1/(rr1 - rr2);
                double m3  = rr1*r1 - rr2*r2;
                dr1   = ( 3*rr1 -  1.33333333333*r1*m3*im2)*im2;  
                dr2   = (-3*rr2 +  1.33333333333*r2*m3*im2)*im2;
                }   
            default: {         //  1 plane
                dr1=1; dr2=1;
                }
        } 
    }
    
};

#endif

