
#ifndef CLCFSF_h
#define CLCFSF_h

/*

####################################################################
#### Compact Linear Combination of Floating Spherical Functions ####
####################################################################

## Idea :
 - Bonding orbitals (BO) are highly localized occupied Molecular orbitas (MOs)
    - This high localization makes evaluation of many interaction very efficient (linear scaling and highly peraelized), almost like a forcefield
    - Using only spherical functions makes formulation and evaluation of many integrals very simple. Although nothing prevent us from extending to non-spherical function in future
    - It is strongly inspired by eFF (electron force-field) and original method of "Floating Gaussian Orbitals" ... one can say it is generalization of these methods
 - We do not Diagonalize any matrix.
    - Diagonalization of Mamiltonian in quantum solutions is in fact just a way of minimizing the variational functional ( which is obvious from Car-Parrignelo method and from using iterative solver such as Conjugate-Gradient for solving sparse hamiltonian )
    - In our case the optimization problem is rather non-linear ( Energy is not linear with respect to movement of basis functions ), therefore we use non-linear optimization algorithm such as FIRE not only for ions but also for electron degrees of freedom
    - We do multiple steps of linear optimization (changing expansion coefs for fixed geometry) per one step of non-linear one (moving atoms and floating functions). although non-linear optimization is posible, linear one is actually faster.
 - How to ensure orthogonality of occupied orbitals?
    - normaly orthogonality is ensured by diagonalization of Hamiltonian, but we want to avoid this cost.
    - We use "Pauli repulsion" or "overlpa froce" instead
 - Since main computational cost in many DFT codes is evaluation of integrals (overlap, kinetic, electrostatic, exchange correlation) we try to minimize number of these integrals
   - we go opposite direction than "frozen/contracted gaussian orbitals". In FGO gaussians of several widths are combined with fixed linear combination coefficient under one object assigned to single variational parameter in the Hamiltonian matrix. This way size of Hamiltonian is limited.
   - But in our case we do not solve Hamoltonian matrix (O(N^3) cost) instead we just move it by a linear solver (O(N) cost). Therefore, in our case the calculation of the integrals is the time limiting step.
 - Electrostatic potential is calculated by multipole expansion of orbital density with distance. Thanks to high localization we can do this expansion very fast (Like fast multipole method)
    - Electron denisty of each orbital (which is linear combination of basis functions) is re-projected on few axuliary basis functions
        rho(r)=phi(r) = sum_i{ c_i chi_i(r)} sum_j{ c_j chi_j(r)}
 - The main problem is evaluation of exchange-correlation. Which is non-linear and cannot be easily experesed as combination of orbitals.
    - This can be done on a real-space grid, but that is very costly perhaps
    - We can use some approximation from Fireball ... e.g. McWeda approximation of Pavel Jelinek
    - We can oalso sort of fit some empirical approximation using machine learning etc.


    IDEA 1) - multiple function per one center. Each can have different cutoff.
       - Variation of coeficients is much faster than variation of position ( re-evaluation of integrals is not required )


*/

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Forces.h"
//#include "GridFF.h"

#include "spline_hermite.h"


struct PairInteractionTable{
    double* Eoverlap=0; // <i|j>
    double* Ekinetic=0; // <i|Lapalce|j>
    // From interpolation of spline (higher order) we can get also derivative i.e. forces
    //double* Foverlap=0; // d<i|j>/dRij
    //double* Fkinetic=0; // <i|Lapalce|j>/dRij
};


class CLCFSF{

    double Rcut    =6.0;  // cutoff beyond which two basis functions chi has no overlap
    double RcutOrb =9.0;  // cutoff beoud which orbital (i.e. localized cluster of basis functions) has no overlap

    double Rcut2   =Rcut*Rcut;

    int natoms=0  // number of atoms (nuclei, ions)
    int perOrb=0; // number of spherical functions per orbital
    int nOrb  =0; // number of single-electron orbitals in system

    Vec3d*  epos=0;   // position of spherical function for expansion of orbitals
    double* ecoefs=0; // expansion coefficient of expansion of given orbital

    // density approx
    Vec3d*  edip =0;   // Axuliary array to store dipoles
    double* erho =0;   // temporary array to store density projection on pair overlap functions

    // ToDo: implement different basis function types
    //int* etypes = 0; // currently not used

    double* eEs=0;
    double* epfs=0;   //   force acting on position of orbitals
    double* ecfs=0;   //   force acting on combination coefficnet of orbitals

    Vec3d*  apos =0;  // positioon of atoms
    Vec3d*  aQs  =0;  // charge of atom
    int*    atype=0;  // type of atom (in particular its pseudo-potential)



    //PairInteractionTable** interactions;
    //PairInteractionTable** eInteractions;

    // Interaction tables   [itype][jtype][isample]
    //double*** ISs = 0;   // < chi_i(r) |    1    | chi_i(r) >  Overlap interactin between two  spherical bais functions
    //double*** IKs = 0;   // < chi_i(r) | Laplace | chi_i(r) >  Kinetic energy interactin between two spherical functions

    double* Wfs = 0;  // |Chi(r)>   radial wave function table
    double* ISs = 0;  // < chi_i(r) |    1    | chi_i(r) >  Overlap interactin between two  spherical bais functions
    double* IKs = 0;  // < chi_i(r) | Laplace | chi_i(r) >  Kinetic energy interactin between two spherical functions

    double*** IrhoV = 0; // < rho(r) V(r) >  interactin between spherical_eletron_density_function and spherical electrostatic_potential_function

    inline int getOrbOffset(int iorb){ return irob*nOrb; }

    /*
    inline double interpolate( Vec3d Rij, double* Is, Vec3d& fij ){
        // ToDo: Rcut may be read from
        double r2 = Rij.norm2();
        if(r2>Rcut2){ fij = Vec3dZero; return 0; }
        double r = sqrt(r2);
        double val,dval;
        Spline_Hermite::valAndDeriv( r, Is, val, dval );
        fij.set_mul( Rij, dval/r );
        return val;
    }
    */

    inline double evalShotRange( Vec3d Rij, double* Is, int i, int j ){

        // ToDo: Rcut may be read from
        double r2 = Rij.norm2();
        if(r2>Rcut2){ fij = Vec3dZero; return 0; }

        double r = sqrt(r2);
        double S,dS,  K,dK;
        int    i   = (int)r;
        T      x   =  r - i;
        Spline_Hermite::valdval( x, S, dS, ISs[i], ISs[i+1], ISs[i+2], ISs[i+3] ); // Overlap
        Spline_Hermite::valdval( x, K, dK, IKs[i], IKs[i+1], IKs[i+2], IKs[i+3] ); // Kinetic energy
        //Pauli(); // Pauli is like S^2 ? ... that is not linear ... should we rather use density overlap ?
        double dP = dS; // ToDo : More sophiticated model for Pauli Repulsion

        double ci = ecoefs[i];
        double cj = ecoefs[j];
        double cij = ci*cj;
        Vec3d fij; fij.set_mul( Rij, cij*(dK + dP)/r );

        epfs[i].add(fij);
        epfs[j].sub(fij);
        return val*cij;
    }

    void evalShortRange(){ // evaluate Energy components given by direct wave-function overlap ( below cutoff Rcut )
        for(int io=0; io<nOrb; io++){
            int i0 = getOrbOffset(io);
            for(int jo=0; jo<nOrb; jo++){
                int j0 = getOrbOffset(jo);
                if( !OrbitalOverps( i, j ) ) continue; // optimization, not to calculate needless interactions
                for(int j=j0; j<j0+perOrb; j++){ // ToDo : we may allow different number of orbitals per electron later (?)
                    Vec3d eposj = epos[j];
                    for(int i=i0; i<i0+perOrb; j++){
                        Vec3d fij;
                        //interpolate( epos[j] - epos[i], Is[ityp][jtyp], f ); // ToDo: different basis function types
                        evalShotRange( epos[j] - epos[i], Is[ityp][jtyp], f );

                    }
                }
            }
        }
    }

    double orb2dens(int io, double* coefs, double* rhos, Vec3d& dips){ // project orbital on axuliary density functions
       int i0=getOrbOffset(io);
        double Q=0;
        for(int i=i0; i<i0+perOrb; i++){
            double ci = coefs[i];
            double qii = ci*ci;
            *rhoAux=qii;
            Q += ci*ci;
            for(int j=i0; j<i; j++){

                double r2 = Rij.norm2();
                if(r2>Rcut2) continue;

                double S,dS;
                Spline_Hermite::valderiv( sqrt(r2), S, dS, ISs ); // Overlap

                double qij = S * ci*coefs[j];
                *rhos += qij;
                rhos++;
                // ToDo : Store qij to list of axuliary functions
                Q +=  qij; //
            }
        }
        double renorm = sqrt(1./Q);
        for(int i=i0; i<i0+perOrb; i++){ coefs[j] *=renorm;  };
        return Q;
    }

    double projectOrbs(){   // project density of all orbitals onto axuliary charge representation ( charges, dipoles and axuliary functions )
        int nqOrb = perOrb*(perOrb+1)/2;
        int irho0=0;
        for(int io=0; io<nOrb; io++){
            int i0 = getOrbOffset(jo);
            eQs[io] = orb2dens(io, coefs+i0, erhos+irho0, edips[io] );
            irho0+=nqOrb;
        }
    }

    void evalElectrostatics( ){
        for(int io=0; io<nOrb; io++){
            int i0 = getOrbOffset(io);
            for(int jo=0; jo<nOrb; jo++){
                int j0 = getOrbOffset(jo);
                if( OrbitalOverps( i, j ) ){
                    // Calculate detailed short-range electrostatics

                }else{
                    // Calculate approximate long-range electrostatics (using charges and dipoles)

                }
            }
        }
    }

    void evalExchangeCorrelation(){
        // not sure how exactly to do that now - no using real space grid since it is slow
    }

    double eval(){
        projectOrbs();
        evalElectrostatics();
    }


};

#endif
