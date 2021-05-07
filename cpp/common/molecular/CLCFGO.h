
#ifndef CLCFGO_h
#define CLCFGO_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"

#include "spline_hermite.h"
#include "GaussianBasis.h"
#include "Forces.h"

#include "Grid.h"
//#include "AOIntegrals.h"

#include "InteractionsGauss.h"
#include "GaussianBasis.h"

/// \defgroup Molecular  Molecular

/*!
\ingroup Molecular
 @{

##################################################################
#### Compact Linear Combination of Floating Gaussian Orbitals ####
##################################################################

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
     \f$   rho(r)=phi(r) = sum_i{ c_i chi_i(r)} sum_j{ c_j chi_j(r)}  \f$
 - The main problem is evaluation of exchange-correlation. Which is non-linear and cannot be easily experesed as combination of orbitals.
    - This can be done on a real-space grid, but that is very costly perhaps
    - We can use some approximation from Fireball ... e.g. McWeda approximation of Pavel Jelinek
    - We can oalso sort of fit some empirical approximation using machine learning etc.


    IDEA 1) - multiple function per one center. Each can have different cutoff.
       - Variation of coeficients is much faster than variation of position ( re-evaluation of integrals is not required )


 ### Make alternative with Gaussian Orbitals (instead of numerica basis)
    - It makes it possible to express auxuliary density functions exactly

   \f$     rho_ij(r)                                                                \f$
   \f$        = Gauss(r-R_i) * Gauss(r-R_j)                                         \f$
   \f$        = \exp( -w_i*( (x-x_i)^2 + (y-y_i)^2 + (y-y_i)^2 ) ) * \exp( -w_j*( (x-x_j)^2 + (y-y_j)^2 + (y-y_i)^2 ) ) \f$
   \f$        = \exp( -w_i*( (x-x_i)^2 + (y-y_i)^2 + (y-y_i)^2 ) -wj*( (x-x_j)^2 + (y-y_j)^2 + (y-y_j)^2 )  )         \f$
   \f$        x : wi*(x-x_i)^2 + wj*(x-x_j)^2   =    (w_i+w_j)x^2   - 2*(w_i*x_i+w_j*x_j)x - (w_i*x_i^2+w_j*x_j^2)        \f$
    - http://www.tina-vision.net/docs/memos/2003-003.pdf
*/
class CLCFGO{ public:

// ToDo : Later properly
constexpr static const Quat4d default_AtomParams[] = {
//  Q   sQ   sP   cP
{ 0.,  1.0, 1.0, 0.0 }, // 0
{ 1.,  0.1, 0.1, 0.0 }, // 1 H
{ 0.,  1.0, 1.0, 2.0 }, // 2 He
{ 1.,  0.1, 0.1, 2.0 }, // 3 Li
{ 2.,  0.1, 0.1, 2.0 }, // 4 Be
{ 3.,  0.1, 0.1, 2.0 }, // 5 B
{ 4.,  0.1, 0.1, 2.0 }, // 6 C
{ 5.,  0.1, 0.1, 2.0 }, // 7 N
{ 6.,  0.1, 0.1, 2.0 }, // 8 O
{ 7.,  0.1, 0.1, 2.0 }, // 9 F
};

    bool bNormalize     = true;
    bool bEvalKinetic   = true;
    bool bEvalCoulomb   = true;
    bool bEvalPauli     = true;
    bool bEvalExchange  = true;
    bool bEvalAECoulomb = true;
    bool bEvalAEPauli   = true;
    bool bEvalAE        = true;
    bool bEvalAA        = true;
    int  iPauliModel    = 1;

    bool bRescaleKinetic = true;

    double Ek=0, Eee=0,EeePaul=0,EeeExch=0, Eae=0,EaePaul=0, Eaa=0; ///< different kinds of energy

    bool bOptAtom = true;
    bool bOptEPos = true;
    bool bOptSize = true;
    bool bOptCoef = true;

    double KPauliOverlap = 50.0; // ToDo : This is just "bulgarian constant" for now
    double KPauliKin     = 50.0; // ToDo : Not sure if we should use this - perhaps this model of pauli energy should work "ab inition"
    constexpr static const Vec3d KRSrho = { 1.125, 0.9, 0.2 }; ///< eFF universal parameters

    double Rcut    =6.0;  ///< cutoff beyond which two basis functions chi has no overlap
    double RcutOrb =9.0;  ///< cutoff beoud which orbital (i.e. localized cluster of basis functions) has no overlap
    double Rcut2     =Rcut*Rcut;
    double RcutOrb2  =RcutOrb*RcutOrb;

    int natypes =0;

    int natom =0; ///< number of atoms (nuclei, ions)
    int perOrb=0; //!< Brief number of spherical functions per orbital
    int perOrb2=0;
    int nOrb  =0; //!< Brief number of single-electron orbitals in system
    // this is evaluated automaticaly
    int nBas  =0; ///< number of basis functions
    int nqOrb =0; ///< number of charges (axuliary density elements) per orbital
    int nQtot =0; ///< total number of charge elements

    // atoms (ions)
    Vec3d*  apos   =0;  ///< [A] positioon of atoms
    Vec3d*  aforce =0;  ///< [eV/A] force on atoms
    double* aQs    =0;  ///< [e] charge of atomic core ( Q_nucleus - Q_core_electrons )
    double* aPcoef =0;  ///< [eV] coeficient of pauli repulsion of core electrons
    double* aPsize =0;  ///< [A] size of core in Pauli interaction
    double* aQsize =0;  ///< [A] size of core in coulombic interaction
    int*    atype  =0;  ///< type of atom (in particular IKinetic pseudo-potential)
    //Vec3d  * aAbWs =0; ///< atomic   parameters (amplitude, decay, width)
    //Vec3d  * eAbWs =0; ///< electron parameters (amplitude, decay, width)

    // orbitals
    Vec3d*  opos  =0;   ///< [A] store positions for the whole orbital
    Vec3d*  odip  =0;   ///< [eA] Axuliary array to store dipoles for the whole orbital
    double* oEs   =0;   ///< [eV] orbital energies
    double* oQs   =0;   ///< [e] total charge in orbital before renormalization (just DEBUG?)
    int*    onq   =0;   ///< number of axuliary density functions per orbital
    int*    ospin =0;

    // --- Wave-function components for each orbital
    Vec3d*  epos  =0; ///< [A] position of spherical function for expansion of orbitals
    double* esize =0; ///< [A] spread of gassian basisfunction
    double* ecoef =0; ///< [e^.5] expansion coefficient of expansion of given orbital
    // --- Forces acting on wave-functions components
    Vec3d*  efpos  =0; ///< [eV/A]  force acting on position of orbitals   d_E/d_epos[i]
    double* efsize =0; ///< [eV/A]  force acting on combination coefficnet of orbitals  d_E/d_esize[i]
    double* efcoef =0; ///< [e^.5V] force acting on size of gaussians d_E/d_ecoef[i]

    // --- Auxuliary electron density expansion basis functions
    Vec3d * rhoP  =0; ///< [A] position of axuliary electron density function
    double* rhoQ  =0; ///< [e] charge in axuliary electron density function
    double* rhoS  =0; ///< [A] spread of axuliary electron density function
    // --- Forces acting on auxuliary density basis functions
    Vec3d * rhofP =0; ///< [eV/A] force on position of axuliary electron density function d_E/d_rhoP[i]
    double* rhofQ =0; ///< [V] force on charge in axuliary electron density function d_E/d_rhoQ[i]
    double* rhofS =0; ///< [eV/A] force on spread of axuliary electron density function d_E/d_rhoS[i]
    //double* rhoEQ =0; /// coulomb energy

    // ======= Functions

    void realloc( int natom_, int nOrb_, int perOrb_, int natypes_ ){
        // ToDo : We may automatize this somehow ????
        // atoms
        if( natom != natom_ ){
            natom = natom_;
            _realloc( apos  ,natom );
            _realloc( aforce,natom );
            _realloc( aQs   ,natom );
            //_realloc( aAbWs ,natom);
            //_realloc( eAbWs ,natom);
            _realloc( aPcoef ,natom);
            _realloc( aPsize ,natom);
            _realloc( aQsize ,natom);
            _realloc( atype ,natom ); // not used  now
        }
        if( (nOrb != nOrb_)||(perOrb != perOrb_) ){
            nOrb    = nOrb_;
            perOrb  = perOrb_;
            nBas    = nOrb * perOrb;
            perOrb2 = perOrb*perOrb;
            nqOrb   = perOrb*(perOrb+1)/2;
            nQtot   = nOrb*nqOrb;

            // orbitals
            _realloc( ospin,nOrb);
            _realloc( opos, nOrb);
            _realloc( odip, nOrb);
            _realloc( oEs , nOrb);
            _realloc( oQs , nOrb);
            _realloc( onq , nOrb);

            // --- Wave-function components for each orbital
            _realloc( epos , nBas );
            _realloc( esize, nBas );
            _realloc( ecoef, nBas );
            // --- Forces acting on wave-functions components
            _realloc( efpos ,   nBas );
            _realloc( efsize,   nBas );
            _realloc( efcoef,   nBas );

            // --- Auxuliary electron density expansion basis functions
            _realloc( rhoP, nQtot );
            _realloc( rhoQ, nQtot );
            _realloc( rhoS, nQtot );
            // --- Forces acting on auxuliary density basis functions
            _realloc( rhofP, nQtot );
            _realloc( rhofQ, nQtot );
            _realloc( rhofS, nQtot );
            //_realloc( rhoEQ, nQtot );
        }
    }

    void setDefaultValues( ){
        // ToDo : We may automatize this somehow ????
        // atoms
        for(int i=0; i<natom;  i++){
            apos  [i]=Vec3dZero;
            aforce[i]=Vec3dZero;
            aQs   [i]=1.;
            aQsize[i]=1.0;
            aPsize[i]=1.0;
            aPcoef[i]=0.0;
            atype [i]=0;
        }
        for(int i=0; i<nOrb;  i++){
            opos[i]=Vec3dZero;
            odip[i]=Vec3dZero;
            oEs [i]=0;
            oQs [i]=0;
            onq [i]=0;
            ospin [i]=1;
        }
        for(int i=0; i<nBas;  i++){
            // --- Wave-function components for each orbital
            epos [i]=Vec3dZero;
            esize[i]=1.;
            ecoef[i]=1.;
            // --- Forces acting on wave-functions components
            efpos [i]=Vec3dZero;
            efsize[i]=1.;
            efcoef[i]=1.;
        }
        for(int i=0; i<nQtot;  i++){
            // --- Auxuliary electron density expansion basis functions
            rhoP[i]=Vec3dZero;
            rhoQ[i]=0;
            rhoS[i]=1.;
            // --- Forces acting on auxuliary density basis functions
            rhofP[i]=Vec3dZero;
            rhofQ[i]=0;
            rhofS[i]=0;
            //rhoEQ[i]=0;
        }
    }

    void turnAllEvalSwitches(bool b){
        bNormalize     = b;
        bEvalKinetic   = b;
        bEvalCoulomb   = b;
        bEvalPauli     = b;
        bEvalExchange  = b;
        bEvalAECoulomb = b;
        bEvalAEPauli   = b;
        bEvalAE        = b;
        bEvalAA        = b;
    }
    void turnAllOptSwitches(bool b){
        bOptAtom = b;
        bOptEPos = b;
        bOptSize = b;
        bOptCoef = b;
    }
    void turnAllSwitches(bool b){ turnAllEvalSwitches(b); turnAllOptSwitches(b); };

    inline int getOrbOffset(int iorb)const{ return iorb*perOrb;  }
    inline int getRhoOffset(int iorb)const{ return iorb*nqOrb; }

    inline void setRcut( double Rcut_ ){
        Rcut   = Rcut_;
    }

    inline bool checkOrbitalCutoff(int i, int j)const{
        double r2 = (opos[i]-opos[j]).norm2();
        return r2<RcutOrb2;
    }

    inline void setElectron( int io, int j, Vec3d p, double s, double c ){
        int i = io*perOrb + j; epos[i]=p; esize[i]=s; ecoef[i]=c;
    }

    inline void setAtom( int i, Vec3d p, double q, double sQ, double sP, double cP ){
        apos[i]=p; aQs[i]=q; aQsize[i]=sQ; aPsize[i]=sP; aPcoef[i]=cP;
    }

    inline void setAtom( int i, int itype, Vec3d p ){
        const Quat4d& par = default_AtomParams[itype];
        apos[i]=p; aQs[i]=par.x; aQsize[i]=par.y; aPsize[i]=par.z; aPcoef[i]=par.w;
    }

    void clearAuxDens(){
        // We do not need to clear this, it is set fully in projection
        //for(int i=0; i<nQtot; i++){ rhoP[i] = 0; };
        //for(int i=0; i<nQtot; i++){ rhoQ[i] = 0; };
        //for(int i=0; i<nQtot; i++){ rhoS[i] = 0; };
        // We need only clear forces which are assembled
        for(int i=0; i<nQtot; i++){ rhofP[i] = Vec3dZero; };
        for(int i=0; i<nQtot; i++){ rhofQ[i] = 0; };
        for(int i=0; i<nQtot; i++){ rhofS[i] = 0; };
        //for(int i=0; i<nQtot; i++){ rhoEQ[i] = 0; };
    }

    inline void cleanForces(){
        for(int i=0; i<natom; i++){ aforce[i] = Vec3dZero; }
        for(int i=0; i<nBas;  i++){ efpos [i] = Vec3dZero; }
        for(int i=0; i<nBas;  i++){ efsize[i] = 0;         }
        for(int i=0; i<nBas;  i++){ efcoef[i] = 0;         }
        //clearAuxDens();
    }

    void reportOrbitals()const{
        for(int io=0; io<nOrb; io++){ for(int ib=0; ib<perOrb; ib++){
            int i=io*perOrb+ib;
            printf( "orb[%i,%i=%i] p(%g,%g,%g) s %g c %g \n", io, ib, i, epos[i].x,epos[i].y,epos[i].z, esize[i], ecoef[i] );
        } }
    }

    void reportCharges(int io)const{
        int i0    = getOrbOffset(io);
        int irho0 = getRhoOffset(io);
        int ij=0;
        for(int i=i0; i<i0+perOrb; i++){
            printf( "Q[%i|%i,%i|%i] %g|%g\n", io,i,i,ij, rhoQ[irho0+ij], ecoef[i] );
            ij++;
            for(int j=i0; j<i; j++){
                printf( "Q[%i|%i,%i|%i] %g|%g,%g\n", io,i,j,ij, rhoQ[irho0+ij], ecoef[i], ecoef[j] );
                ij++;
            }
        }
    }
    void reportCharges(){ for(int io=0; io<nOrb; io++){ reportCharges(io); } }




// ===========================================================================================================================
// ==================== Project Orbitals ( Normalize, Eval Kinetic Energy, Project to Auxiliary Density Basis )
// ===========================================================================================================================

    //double projectOrb(int io, Vec3d* Ps, double* Qs, double* Ss, Vec3d& dip, bool bNormalize ){ // project orbital on axuliary density functions
    double projectOrb(int io, Vec3d& dip ){

        bool bMakeRho = bEvalCoulomb || bEvalAECoulomb || bNormalize;   // in order to normalize we must calculate total charge in orbital
        //printf( "projectOrb() bMakeRho %i  | bEvalCoulomb %i bEvalAECoulomb %i bNormalize %i \n",  bMakeRho , bEvalCoulomb , bEvalAECoulomb , bNormalize );

        int i0    = getOrbOffset(io);
        int irho0 = getRhoOffset(io);
        Vec3d*   Ps =rhoP+irho0;
        double*  Qs =rhoQ+irho0;
        double*  Ss =rhoS+irho0;
        double Q=0;
        //double DT=0; // kinetic energy change by orthognalization
        double Ek=0; // kinetic energy
        dip         = Vec3dZero;
        Vec3d qcog  = Vec3dZero;
        //Vec3d oqcog = opos[io]; // Do We need this ?
        int ii=0;

        for(int i=i0; i<i0+perOrb; i++){

            Vec3d  pi  = epos [i];
            double ci  = ecoef[i];
            double si  = esize[i];
            double qii = ci*ci; // overlap = 1

            if(bMakeRho){
                Qs[ii]  = qii;
                Ps[ii]  = pi;
                Ss[ii]  = si*M_SQRT1_2;
                Q      += qii;
                qcog.add_mul( pi, qii );
                //printf( "projectOrb[%i|%i]:bMakeRho %i  qii %g  i0 %i \n", io,i, bMakeRho, qii, i0 );
            }

            double fr,fsi,fsj;
            if(bEvalKinetic){
                //double fEki;
                //Ek += qii*addKineticGauss( si*M_SQRT2, fEki );
                //Ek += qii*Gauss::kinetic( si );
                Ek += qii*Gauss::kinetic_s(  0.0, si, si,   fr, fsi, fsj );
                efsize[i]+= -2*fsi*qii;
            }

            ii++;

            for(int j=i0; j<i; j++){
                Vec3d pj  = epos[j];
                Vec3d Rij = pj-pi;
                double r2 = Rij.norm2();
                if(r2>Rcut2) continue;

                double cj  = ecoef[j];
                double sj  = esize[j];
                double cij = ci*cj*2;

                // --- Evaluate Normalization, Kinetic & Pauli Energy

                //const double resz = M_SQRT2; // TODO : PROBLEM !!!!!!!   getOverlapSGauss and getDeltaTGauss are made for density gaussians not wave-function gaussians => we need to rescale "sigma" (size)
                //double dSr, dSsi, dSsj;
                //double Sij  = getOverlapSGauss( r2, si*resz, sj*resz, dSr, dSsi, dSsj );

                //double dTr, dTsi, dTsj;
                //double DTij = getDeltaTGauss  ( r2, si*resz, sj*resz, dTr, dTsi, dTsj ); // This is not normal kinetic energy this is change of kinetic energy due to orthogonalization
                 //DT += DTij*cij;


                // TODO WARRNING : We calculate this BEFORE NORMALIZATION !!! we should do it after !!!!

                if(bMakeRho){
                    //printf( "projectOrb[%i|%i,%i]:bMakeRho %i \n", io,i,j, bMakeRho );
                    // --- Project on auxuliary density functions
                    Vec3d  pij;
                    double sij;
                    //double Cij = Gauss::product3D_s( si, pi, sj, pj, sij, pij );
                    double Sij = Gauss::product3D_s_new( si, pi, sj, pj, sij, pij );
                    double qij = Sij*cij;
                    Qs[ii] = qij;
                    Ps[ii] = pij;
                    Ss[ii] = sij;

                    Q += qij;
                    qcog.add_mul( pij, qij );
                }

                if(bEvalKinetic){
                    //double Ekij = Gauss::kinetic(  r2, si, sj ) * 2; // TODO : <i|Lapalace|j> between the two gaussians
                    double Kij = Gauss:: kinetic_s(  r2, si, sj,   fr, fsi, fsj );   // fr*=2; fsi*=2, fsj*=2;
                    Ek += Kij*cij;
                    Vec3d fij = Rij*(fr*cij);
                    efpos [i].add( fij ); efpos[j].sub( fij );
                    efsize[i]-= fsi*cij ; efsize[j]-= fsj*cij;
                    efcoef[i]-= Kij*cj ;  efcoef[j]-= Kij*ci;
                    //if(DEBUG_iter==DEBUG_log_iter){ printf(" Kij %g cij %g Ekij \n", Kij, cij, Kij*cij ); }
                }

                ii++;
                // ToDo : Store qij to list of axuliary functions
            }
        }
        onq[io] = ii;
        oQs[io] = Q;
        opos[io].set_mul( qcog, 1./Q );
        if(bNormalize){
            double renorm  = sqrt(1./Q);
            double renorm2 = renorm*renorm;
            //printf( "project orb[$i]: Q %g renorm %g renorm2 %g \n", io, Q, renorm, renorm2 );
            for(int i=i0; i<i0+perOrb; i++){ ecoef[i] *=renorm;  };
            for(int i= 0; i<ii       ; i++){ Qs   [i] *=renorm2; };
        }
        // ToDo: Renormalize also  rhos?
        return Ek;
    }

    double projectOrbs(bool bNormalize){   // project density of all orbitals onto axuliary charge representation ( charges, dipoles and axuliary functions )
        int nqOrb = perOrb*(perOrb+1)/2;
        int i0 =0;
        int ii0=0;
        Ek = 0;
        for(int io=0; io<nOrb; io++){
            //int i0  = getOrbOffset(jo);
            //oQs[io] =
            Ek += projectOrb( io, odip[io] );
            //projectOrb(  io, rhoP+ii0, rhoQ+ii0, rhoS+ii0, odip[io], true );
            //projectOrb(io, ecoefs+i0, erho+irho0, erhoP+irho0, odip[io] );
            ii0+=nqOrb;
            i0 +=perOrb;
        }
        return Ek;
    }


// ===========================================================================================================================
// ==================== Eval Coulombin interaction  -  Using Auxuliary density basis
// ===========================================================================================================================

    void toRho( int i, int j, int ij ){
        /// NOTE : this function is mostly for debugging of  projectOrb()
        Vec3d  pi  = epos [i];
        double ci  = ecoef[i];
        double si  = esize[i];
        Vec3d  pj  = epos [j];
        double cj  = ecoef[j];
        double sj  = esize[j];

        Vec3d Rij = pj-pi;
        double r2 = Rij.norm2();

        //double dSr, dSsi, dSsj;
        //const double resz = M_SQRT2; // TODO : PROBLEM !!!!!!!   getOverlapSGauss and getDeltaTGauss are made for density gaussians not wave-function gaussians => we need to rescale "sigma" (size)
        //double Sij_  = getOverlapSGauss( r2, si*resz, sj*resz, dSr, dSsi, dSsj );

        // ToDo : we should not need   getOverlapSGauss,    Gauss::product3D_s  should calculate Sij

        Vec3d  pij;
        double sij;
        //double Cij = Gauss::product3D_s( si, pi, sj, pj, sij, pij );
        double Sij = Gauss::product3D_s_new( si, pi, sj, pj, sij, pij );
        double cij = ci *cj;
        //printf(  "ci*cj %g ci %g cj %g \n", cij, ci, cj );
        double qij = Sij*cij*2; // TODO CHECK: should there by realy coefficeint 2.0 ?  .... test by grid !
        //double qij = Sij;
        //double qij = Sij*cij*4.85;
        //double qij = Sij_*cij*4;
        //double qij = Cij*cij*2;
        rhoQ[ij] = qij;
        rhoP[ij] = pij;
        rhoS[ij] = sij;
    }

    //Vec3d fromRho( int i, int j, int ij ){
    void fromRhoDiag( int i, int ij ){
        /// This function function cares about diagonal density terms rho_ii = wi * wi
        //Vec3d  pi  = epos [i];
        //double ci  = ecoef[i];
        //double si  = esize[i];
        Vec3d  Fpi  = rhofP[ij];
        double Fsi  = rhofS[ij];
        double dEdQ = rhofQ[ij];
        double Q    = rhoQ[ij];
        //double Eqi = rhoEQ[ij];
        efpos [i].add( Fpi );
        efsize[i] += Fsi*Q * M_SQRT1_2;
        if(DEBUG_iter==DEBUG_log_iter){
            //printf( "fromRho[%i]ii[%i] Fpi(%g,%g,%g) Fqi %g | Qi %g \n", i, ij, Fpi.x,Fpi.y,Fpi.z, rhoEQ[ij],    rhofQ[ij] );
            //printf( "fromRho[%i]ii[%i] E %g qij %g Fs %g \n", i, ij, dEdQ*Q, Q, Fsi*Q );
        }
    }

    //Vec3d fromRho( int i, int j, int ij ){
    void fromRho( int i, int j, int ij ){
    //void fromRho( int i, int j, int ij, double& aij, double& dCsi, double& dCsj, Vec3d& dCdp ){
        /// NOTE : this function is mostly for debugging of  assembleOrbForces()

        Vec3d  pi  = epos [i];
        double ci  = ecoef[i];
        double si  = esize[i];

        Vec3d  pj  = epos [j];
        double cj  = ecoef[j];
        double sj  = esize[j];

        //printf(  ":ci*cj %g ci %g cj %g \n", ci*cj, ci, cj );

        Vec3d Rij = pj-pi;
        double r2 = Rij.norm2();

        double Sij;
        Vec3d  p;
        double s;
        double dssi,dssj;
        Vec3d  dxsi,dxsj;
        double dxxi,dxxj;
        double dSr;
        double dSsi,dSsj;


        // NOTE: we must compute it twice
        // 1) in toRho() to obtain charges and positions
        // 2) here to obtain derivatives dXxi,dXxj,dCr
        // We cannot connect it because CoublombElement needs do (1)before and (2)after itself
        Sij = Gauss::product3D_s_deriv(
            si,   pi,
            sj,   pj,
            s ,   p ,
            dssi, dssj,
            dxsi, dxsj,
            dxxi, dxxj,
            dSsi, dSsj, dSr
        );

        Vec3d  Fp   = rhofP[ij];  //   fij = Rij*(-fr*qij);   ... from CoulombElement
        double Fs   = rhofS[ij];
        double dEdQ = rhofQ[ij];  // TODO: Eqi is probably redudent ( Fqi is the same )  !!!

        double cij = 2*ci*cj;

        // TODO : we should make sure there is a recoil ( forces point to opposite direction )
        double fsi = Fp.dot( dxsi ) + Fs*dssi + dEdQ*dSsi*cij;
        double fsj = Fp.dot( dxsj ) + Fs*dssj + dEdQ*dSsj*cij;
        Vec3d  fxi = Fp*dxxi;
        Vec3d  fxj = Fp*dxxj;

        Vec3d  dSdp = Rij*(dSr*cij);
        DEBUG_dQdp = dSdp;
        Vec3d Fq   = dSdp*dEdQ;
        Vec3d Fxi  = fxi + Fq;
        Vec3d Fxj  = fxj + Fq;

        efpos [i].add( Fxi );
        efpos [j].add( Fxj );
        efsize[i] += fsi;
        efsize[j] += fsj;

        if(DEBUG_iter==DEBUG_log_iter){
            //printf( "fromRho[%i,%i] eqj %g E %g Fs %g dSsi %g dCsi %g cij %g \n", i,j, dEdQ, dEdQ*rhoQ[ij], Fs, dssi, dSsi, cij );
            //printf( "fromRho[%i,%i] E %g qij %g Fp %g fp*dxxi %g Fq %g \n", i,j, dEdQ*rhoQ[ij], rhoQ[ij], Fxi.x, (Fp*dxxi).x, Fq.x  );
            //printf( "fromRho[%i,%i] dS %g  dSr %g cij %g dEdQ %g Fq.x %g F %g F[i] %g F[j] %g \n", i,j, dSdp.x, dSr*Rij.x, ci*cj, dEdQ, Fq.x, Fxi.x, efpos[i].x, efpos[j].x );
        }
    }

    void assembleOrbForces_fromRho(int io ){
        // NOTE : This is less efficient but more modular version of assembleOrbForces which use fromRho function per each element
        //        It is used for debuging since it is more modular, but after everything works original assembleOrbForces should be prefered
        //   line_Fana->ys[i]  = 0.5*solver.efpos[0].x + E_*line_dQi_ana->ys[i];
        //   Fana              = efpos                + EK*dQdp; ..... TODO
        int i0    = getOrbOffset(io);
        int irho0 = getRhoOffset(io);
        int ii    = irho0;
        for(int i=i0; i<i0+perOrb; i++){
            fromRhoDiag( i, ii );  // Why is this zero ????
            ii++;
            for(int j=i0; j<i; j++){
                fromRho( i, j, ii );
                ii++;
            }
        }
    }

    double CoublombElement( int i, int j ){

        Vec3d  pi = rhoP[i];
        double qi = rhoQ[i];
        double si = rhoS[i];

        Vec3d  pj = rhoP[j];
        double qj = rhoQ[j];
        double sj = rhoS[j];

        Vec3d Rij = pj-pi;
        double r2 = Rij.norm2();
        Vec3d fp; double fs;
        double E = Gauss::Coulomb( Rij, r2, si, sj, qi*qj, fp, fs );
        /*
        double r    = sqrt(r2 + R2SAFE);
        double s2   = si*si + sj*sj;
        double s    = sqrt(s2);
        double fr,fs;
        //double Eqq  = CoulombGauss( r, s*2, fr, fs, qij );
        //double E  = Gauss::Coulomb( r, s*2, fr, fs );   // Q :  Should there be the constant s*2 ????
        double E  = Gauss::Coulomb( r, s, fr, fs );       // WARRNING  :  removed the contant s*2 to s  ... is it correct ?
        //printf( "CoublombElement[%i,%i] q(%g,%g) E %g fs %g fr %g s %g r %g \n", i,j, qi,qj, E, fs, fr, s, r );
        Vec3d fp = Rij*(fr*qij);
        fs *= qij;
        */
        rhofP[i].add(fp);   rhofP[j].sub(fp);
        rhofS[i] -= fs*si;  rhofS[j] -= fs*sj; // Q: ??? Should not this be switched (i<->j)  rhofS[i] -= fs*sj instead of rhofS[i] -= fs*si ???
        rhofQ[i] += E*qj;   rhofQ[j] += E*qi;  // ToDo : need to be made more stable ... different (qi,qj)

        return E;
    }


    double CoulombOrbPair( int io, int jo ){
        /// Calculate Elecrostatic Energy&Force for projected density of orbitas @io and @jo
        int i0 = getRhoOffset(io);
        int j0 = getRhoOffset(jo);
        int nio = onq[io];
        int njo = onq[jo];
        double Ecoul=0;
        //printf( "CoulombOrbPair[%i,%i] (%i:%i) (%i:%i) \n", io, jo, i0,i0+nio, j0,j0+njo );
        for(int i=i0; i<i0+nio; i++){
        //int i=2;{ // DEBUG
            Vec3d  pi = rhoP[i];
            double qi = rhoQ[i];
            double si = rhoS[i];
            for(int j=j0; j<j0+njo; j++){
            //int j=j0+2;{ // DEBUG
                Vec3d  pj = rhoP[j];
                Vec3d Rij = pj-pi;
                double r2 = Rij.norm2();
                //if(r2>Rcut2) continue;    // ToDo : do this check
                double qj = rhoQ[j];
                double sj = rhoS[j];
                Vec3d fp; double fs,qij = qi*qj;
                double E = Gauss::Coulomb( Rij, r2, si, sj, qij, fp, fs );
                /*
                double r    = sqrt(r2 + R2SAFE);
                double s2   = si*si + sj*sj;
                double s    = sqrt(s2);
                double qij = qi*qj;
                double fr,fs;
                // see    InteractionsGauss.h :: addCoulombGauss( const Vec3d& dR, double si, double sj, Vec3d& f, double& fsi, double& fsj, double qq ){
                //double Eqq  = CoulombGauss( r, s*2, fr, fs, qij );
                //fs*=4;
                double E  = Gauss::Coulomb( r, s, fr, fs ); // NOTE : remove s*2 ... hope it is fine ?
                // --- Derivatives (Forces)
                Vec3d fp = Rij*(fr * qij);
                fs       *=           qij;
                //fs      *= M_SQRT1_2 * qij;
                */
                rhofP[i].add(fp);    rhofP[j].sub(fp);
                rhofS[i] -= fs*si;   rhofS[j] -= fs*sj;
                rhofQ[i] += E*qj;    rhofQ[j] += E*qi; // ToDo : need to be made more stable ... different (qi,qj)
                Ecoul    += E*qij;

                //printf( "CoulombOrbPair[%i,%i|%i,%i] E %g E,q(%g,%g) r %g s(%g,%g) \n",io,jo,i,j, Ecoul, E, qij, sqrt(r2), si, sj );

                if(DEBUG_iter=DEBUG_log_iter){
                    //printf( "CoublombElement[%i,%i] q(%g,%g) E %g fs %g fr %g s %g r %g \n", i,j, qi,qj, E, fs, fr, s, r );
                    //printf( "CoublombElement[%i,%i] q(%g,%g) E %g fs %g fr %g s %g r %g \n", i,j, qi,qj, E*qij, fp.x, fs*si, s, r );
                    //printf( "CoublombElement[%i,%i] q(%g,%g) E %g fs %g fr %g s %g r %g \n", i,j, qi,qj, E, fs, fr, s, r );
                    //printf( "CoulombOrbPair[%i,%i][%i,%i] e %g E %g s %g(%g,%g) q %g(%g,%g) r %g fr %g \n", io,jo, i,j,  E, E*qi*qj, s,si,sj, qij,qi,qj, r, fr );
                    //printf( "CoulombOrbPair[%i,%i] E %g r %g \n", i-i0,j-j0,E*qij,r,  );
                }
            }
        }

        //printf( " Ecoul[%i,%i] %g \n", io, jo, Ecoul );
        //printf( "CoulombOrbPair Eorb %g \n", Ecoul );
        return Ecoul;
    }


/*

### Coulomb Force Derivatives

        evalElectrostatICoulomb();
        for(int i=0; i<nOrb; i++) assembleOrbForces(i);

   Problem is derivative of Q
     qi = Sab(ra-rb) * Ca * Cb
     qj = Scd(rc-rd) * Cc * Cd                   ... where Ca,Cb,Cc,Cd are expansion coeficients in given basis and S is overlap integral for given pair of basis functions
   EQij(r) = qi(ra,rb) * qj(rc,rd) * Kij(ri,rj)  ... where qi,qj are charges in some overlap cloud and Kij is Coulomb Matrix kernel between the two clouds
   Force calculated by derivatives as:
     FQ_xa = dEQ/dxa =  (qi*qj) * (dKij/ri)/(dri/dxa) + ( Kij*qj ) * (dqi/dxa)

*/

    double evalElectrostatICoulomb( ){
        /// Calculate detailed short-range electrostatICoulomb
        double r2safe = 0.1;
        //double E = 0;
        Eee = 0;
        for(int io=0; io<nOrb; io++){
            int i0 = getRhoOffset(io);
            Vec3d opi  = opos[io];
            Vec3d dipi = odip[io];
            for(int jo=io+1; jo<nOrb; jo++){
                Vec3d dop = opos[jo] - opi;
                double r2 = dop.norm2();
                //printf( "evalElectrostatICoulomb[%i,%i]  %g <? %g \n", io, jo, r2,RcutOrb2 );

                double dEee = CoulombOrbPair( io, jo ); Eee+=dEee;
                printf( "evalElectrostatICoulomb[%i,%i] Eee %g r2(%g)<?R2cutOrb2(%g) \n", io, jo, dEee, r2, RcutOrb2 );

                /*
                // ==== Long Range Electrostatics Approximation
                if( r2<RcutOrb2 ){
                    int j0 = getRhoOffset(jo);
                    // ToDo : nrho,nV may vary
                    E+= CoulombOrbPair( io, jo );
                }else{
                    double ir2 = 1/(r2+r2safe);
                    double ir  = sqrt(ir2);
                    // Calculate approximate long-range electrostatICoulomb (using charges and dipoles)
                    const Vec3d& dipj = odip[jo];
                    E += ir +                                   // charge-charge
                      + dop.dot(dipi) + dop.dot(dipj) *ir*ir2  // diple-charge
                      + dipi.dot(dipj);                        // dipole-dipole
                }
                */
            }
        }
        return Eee;
    }


// ===========================================================================================================================
// ==================== Eval Cross-Terms ( Contact of two different orbitals, such as Pauli and Exchange interaction )
// ===========================================================================================================================

    double pairOverlapDervis( int io, int jo, Gauss::Blob* Bs, Gauss::PairDeriv* dBs ){
        // ToDo : This may be integrated within   pairOverlapAndKineticDervis() using only switch
        int i0=io*perOrb;
        int j0=jo*perOrb;
        int ij=0;
        // --- Project   rho_ij(r)  to temp axuliary basis
        //if(DEBUG_iter==DEBUG_log_iter) reportOrbitals();
        double S=0;
        for(int i=i0; i<i0+perOrb; i++){
            Vec3d  pi  = epos [i];
            double ci  = ecoef[i];
            double si  = esize[i];
            for(int j=j0; j<j0+perOrb; j++){
                Vec3d pj  = epos[j];
                Vec3d Rij = pj-pi;
                double r2 = Rij.norm2();
                //if(r2>Rcut2) continue;
                double cj  = ecoef[j];
                double sj  = esize[j];
                double cij = ci*cj;
                //pairs[ij].get(si,pi, sj,pj);
                Gauss::product3DDeriv(si,pi,sj,pj, Bs[ij], dBs[ij]);
                Bs[ij].charge *= cij;
                S += Bs[ij].charge;
                //if(DEBUG_iter==DEBUG_log_iter) printf( "Exchange:Project[%i,%i] ss(%g,%g) cij %g Sij %g Qij %g r %g x(%g,%g) \n", i, j, si, sj, cij, pairs[ij].C, cij*pairs[ij].C, sqrt(r2), pi.x, pj.x );
                //printf( "pairOverlapDervis[%i,%i]%i\n", i, j, ij );
                ij++;
            }
        }
        return S;
    }

    double pairKineticDervis( int io, int jo, Quat4d* TDs ){
        // TODO : this can be merged with  pairOverlapDervis() for better performance
        // ToDo : This is not required anymore since we have    pairOverlapAndKineticDervis()
        int i0=io*perOrb;
        int j0=jo*perOrb;
        int ij=0;
        // --- Project   rho_ij(r)  to temp axuliary basis
        //if(DEBUG_iter==DEBUG_log_iter) reportOrbitals();
        double T=0;
        for(int i=i0; i<i0+perOrb; i++){
            Vec3d  pi  = epos [i];
            double ci  = ecoef[i];
            double si  = esize[i];
            for(int j=j0; j<j0+perOrb; j++){
                Vec3d pj  = epos[j];
                Vec3d Rij = pj-pi;
                double r2 = Rij.norm2();
                //if(r2>Rcut2) continue;
                double cj  = ecoef[j];
                double sj  = esize[j];
                double cij = ci*cj;
                //pairs[ij].get(si,pi, sj,pj);
                //Gauss::product3DDeriv(si,pi,sj,pj, Bs[ij], dBs[ij]);
                Quat4d& TD = TDs[ij];
                //r2 *= KRSrho.x*KRSrho.x;
                //si *= KRSrho.y;
                //sj *= KRSrho.y;
                double Tij = Gauss:: kinetic_s(  r2, si, sj,  TD.z, TD.x, TD.y );
                TD.e  = Tij;
                T    += Tij;
                //Bs[ij].charge *= cij;
                //S += Bs[ij].charge;
                //if(DEBUG_iter==DEBUG_log_iter) printf( "Exchange:Project[%i,%i] ss(%g,%g) cij %g Sij %g Qij %g r %g x(%g,%g) \n", i, j, si, sj, cij, pairs[ij].C, cij*pairs[ij].C, sqrt(r2), pi.x, pj.x );
                //printf( "pairOverlapDervis[%i,%i]%i\n", i, j, ij );
                ij++;
            }
        }
        return T;
    }

    double pairOverlapAndKineticDervis( int io, int jo, Gauss::Blob* Bs, Gauss::PairDeriv* dBs, Quat4d* TDs, double& S ){
        // ToDo: If we put here a switch we can use it in place of    pairOverlapAndKineticDervis()
        int i0=io*perOrb;
        int j0=jo*perOrb;
        int ij=0;
        // --- Project   rho_ij(r)  to temp axuliary basis
        //if(DEBUG_iter==DEBUG_log_iter) reportOrbitals();
        double T=0;  S=0;
        for(int i=i0; i<i0+perOrb; i++){
            Vec3d  pi  = epos [i];
            double ci  = ecoef[i];
            double si  = esize[i];
            for(int j=j0; j<j0+perOrb; j++){
                Vec3d pj  = epos[j];
                Vec3d Rij = pj-pi;
                double r2 = Rij.norm2();
                //if(r2>Rcut2) continue;
                double cj  = ecoef[j];
                double sj  = esize[j];
                double cij = ci*cj;
                //pairs[ij].get(si,pi, sj,pj);
                // ---- Overlap
                Gauss::product3DDeriv(si,pi,sj,pj, Bs[ij], dBs[ij]);
                Bs[ij].charge *= cij;
                S += Bs[ij].charge;
                // ---- Kinetic
                Quat4d& TD = TDs[ij];
                if(bRescaleKinetic){
                    r2 *= KRSrho.x*KRSrho.x;
                    si *= KRSrho.y*M_SQRT1_2;
                    sj *= KRSrho.y*M_SQRT1_2;
                }
                double Tij = Gauss::kinetic_s(  r2, si, sj,  TD.z, TD.x, TD.y );
                TD.e  = Tij;
                T    += Tij;
                //if(DEBUG_iter==DEBUG_log_iter) printf( "Exchange:Project[%i,%i] ss(%g,%g) cij %g Sij %g Qij %g r %g x(%g,%g) \n", i, j, si, sj, cij, pairs[ij].C, cij*pairs[ij].C, sqrt(r2), pi.x, pj.x );
                //printf( "pairOverlapDervis[%i,%i]%i\n", i, j, ij );
                ij++;
            }
        }
        return T;
    }

    double evalCoulombPair( int ni, int nj, Gauss::Blob* Bis, Gauss::Blob* Bjs, Gauss::Blob* dBis, Gauss::Blob* dBjs ){
        // ToDo : can we integrate this with    pairOverlapAndKineticDervis()    ???????
        //        perhaps yes but the cost of duplicating    Bs and dBs
        double E = 0;
        for(int i=0; i<ni; i++){
            const Gauss::Blob&  Bi =  Bis[i];
                  Gauss::Blob& dBi = dBis[i];
            for(int j=0; j<nj; j++){   // TODO : perhaps we should make special version for Bis==Bjs  with  j<=i
            //int j=1; {
                const Gauss::Blob&  Bj =  Bjs[j];
                      Gauss::Blob& dBj = dBjs[j];
                Vec3d Rij = Bi.pos - Bj.pos;
                double r2 = Rij.norm2();
                double r  = sqrt(r2 + R2SAFE );
                double s2 = Bi.size*Bi.size + Bj.size*Bj.size;
                double s  = sqrt(s2);
                double qij = Bi.charge*Bj.charge * 2 ;
                //printf( "//Exchange:Coulomb[%i,%i] r %g s %g \n", i, j, r, s  );
                double fr,fs;
                double e  = Gauss::Coulomb( r, s, fr, fs ); // NOTE : remove s*2 ... hope it is fine ?
                //printf( "Exchange:Coulomb[%i,%i] qij %g E %g fr %g fs %g \n", i, j, r, s, qij, E*qij, fs, fs );
                // --- Derivatives (Forces)
                Vec3d fp = Rij*(fr * qij);
                fs       *=           qij;
                //if(DEBUG_iter==DEBUG_log_iter) printf( "ExchangeOrbPair:Coulomb[%i,%i] qij %g e %g fx %g fs %g | r %g s %g fr %g fs %g  \n", i,j, qij, e, fp.x, fs,    r, s, fr, fs );
                dBi.pos.add(fp);              dBj.pos.sub(fp);
                dBi.size   -= fs*Bi.size;     dBj.size   -= fs*Bj.size;
                dBi.charge += e*Bj.charge*2;  dBj.charge += e*Bi.charge*2; // ToDo : need to be made more stable ... different (qi,qj)
                E      += e*qij;
            }
        }
        return E;
    }

    double forceFromOverlaps( int io, int jo, Gauss::PairDeriv* dBs, Gauss::Blob* fBs ){
        int i0=io*perOrb;
        int j0=jo*perOrb;
        int ij=0;
        // --- Project   rho_ij(r)  to temp axuliary basis
        //if(DEBUG_iter==DEBUG_log_iter) reportOrbitals();
        double S=0;
        ij=0;
        //printf( "forceFromOverlaps perOrb %i nqOrb %i  (%i,%i)\n", perOrb, nqOrb, i0,j0 );
        for(int i=i0; i<i0+perOrb; i++){
            //printf( "forceFromOverlaps[%i]\n", i );
            Vec3d  pi  = epos [i];
            double ci  = ecoef[i];
            for(int j=j0; j<j0+perOrb; j++){
                //printf( "forceFromOverlaps[%i,%i]%i nqOrb %g\n", i, j, ij, nqOrb );
                Vec3d pj  = epos[j];
                Vec3d Rij = pj-pi;
                double r2 = Rij.norm2();
                //if(r2>Rcut2) continue;
                double cj  = ecoef[j];
                double cij = ci*cj;
                //if(DEBUG_iter==DEBUG_log_iter) printf( "Exchange:fromRho[%i,%i] x(%g,%g) ", i, j, pi.x, pj.x );
                //pair.backForce( Rij, cij, fQs[ij], fPs[ij], fSs[ij],  efpos[i], efpos[j], efsize[i], efsize[j] );
                Gauss::productBackForce( fBs[ij], dBs[ij], Rij, cij, efpos[i], efpos[j], efsize[i], efsize[j] );
                ij++;
            }
        }
        //printf("forceFromOverlaps() DONE!");
        return S;
    }

    // ToDo : This may be perhaps integrated into  ***  forceFromOverlaps() ***
    void applyPairForce( int io, int jo, Gauss::Blob* Bs, Gauss::PairDeriv* dBs, double Amp ){
        //printf( "applyPairForce [%i,%i] Amp %g \n", io, jo, Amp );
        int i0=io*perOrb;
        int j0=jo*perOrb;
        int ij=0;
        // --- Project   rho_ij(r)  to temp axuliary basis
        //if(DEBUG_iter==DEBUG_log_iter) reportOrbitals();
        //double E=0;
        for(int i=i0; i<i0+perOrb; i++){
            Vec3d  pi  = epos [i];
            double ci  = ecoef[i];
            double si  = esize[i];
            for(int j=j0; j<j0+perOrb; j++){
                Vec3d pj  = epos[j];
                Vec3d Rij = pj-pi;
                double r2 = Rij.norm2();
                //if(r2>Rcut2) continue;
                double cj  = ecoef[j];
                double sj  = esize[j];
                double cij = ci*cj*Amp;
                //pairs[ij].get(si,pi, sj,pj);
                const Gauss::PairDeriv& dB = dBs[ij];
                double eij=Amp*Bs[ij].charge;
                //E   += eij;
                Vec3d Fp = Rij*(dB.dCr*cij);
                efpos [i].add( Fp ); efpos[j].sub( Fp );
                efsize[i]-= dB.dCsi*cij; efsize[j]-= dB.dCsj*cij;
                //efcoef[i]-= eij*cj ;  efcoef[j]-= eij*ci;
                efcoef[i]-= eij/ci ;  efcoef[j]-= eij/cj;
                //if(DEBUG_iter==DEBUG_log_iter) printf( "applyPairForce[%i,%i]%i  dCr %g dCsi %g dCsj %g cij*Amp %g \n", i, j, ij, dB.dCr, dB.dCsi, dB.dCsj, cij );
                //if(DEBUG_iter==DEBUG_log_iter) printf( "applyPairForce[%i,%i]%i Qij %g eij %g dCr %g Fp.x %g cij %g Amp %g \n", i, j, ij, Bs[ij].charge, eij, dB.dCr, Fp.x, cij, Amp );
                //printf( "applyPairForce[%i,%i]%i\n", i, j, ij );
                ij++;
            }
        }
        //return E;
    }

    // ToDo : This may be perhaps integrated into  ***  forceFromOverlaps() ***
    double pauliCrossKinetic( int io, int jo, Gauss::Blob* Bs, Gauss::PairDeriv* dBs, Quat4d* TDs, double S, double T, Vec3d KRSrho, bool anti ){
        // This is Pauli potential derived from Kinetic energy;
        // See Eq.3a,b in article:
        // Su, J. T., & Goddard, W. A. (2009), The Journal of Chemical Physics, 131(24), 244501.
        // The dynamics of highly excited electronic systems: Applications of the electron force field.
        // See  https://doi.org/10.1063/1.3272671    or    http://aip.scitation.org/doi/10.1063/1.3272671
        //
        // Same spin (Eq.3a)    Ep = T * ( (1-a)*S/(1+S) + S/(1-S) )
        // anti spin (Eq.3b)    Ep = T * (    a *S/(1+S)           )   // Negligible ?
        //
        //See : /home/prokop/git/SimpleSimulationEngine/cpp/common/molecular/InteractionsGauss.h
        // addPauliGauss( const Vec3d& dR, double si, double sj, Vec3d& f, double& fsi, double& fsj, bool anti, const Vec3d& KRSrho )

        anti = true; // DEBUG WARRNING !!!!!!!!

        T *= KPauliKin;

        double eS,fS;
        if(anti){ eS = PauliSGauss_anti( S, fS, KRSrho.z ); }
        else    { eS = PauliSGauss_syn ( S, fS, KRSrho.z ); }
        double TfS = T*fS;

        int i0=io*perOrb;
        int j0=jo*perOrb;
        int ij=0;
        // --- Project   rho_ij(r)  to temp axuliary basis
        //if(DEBUG_iter==DEBUG_log_iter) reportOrbitals();
        double E=0;
        for(int i=i0; i<i0+perOrb; i++){
            Vec3d  pi  = epos [i];
            double ci  = ecoef[i];
            double si  = esize[i];
            for(int j=j0; j<j0+perOrb; j++){
                Vec3d pj  = epos[j];
                Vec3d Rij = pj-pi;
                double r2 = Rij.norm2();
                //if(r2>Rcut2) continue;
                double cj  = ecoef[j];
                double sj  = esize[j];
                //double cij = ci*cj*Amp;
                double cij = ci*cj;
                //pairs[ij].get(si,pi, sj,pj);
                const Gauss::PairDeriv& dB = dBs[ij];
                const Quat4d&           TD = TDs[ij];

                // ToDo : this is in    addPauliGauss()    @     /home/prokop/git/SimpleSimulationEngine/cpp/common/molecular/InteractionsGauss.h

                //double dTr, dTsi, dTsj;
                //double T = Gauss:: kinetic_s(  r2, si, sj,   dTr, dTsi, dTsj );

                double dTsi = TD.x* KPauliKin;
                double dTsj = TD.y* KPauliKin;
                double dTr  = TD.z* KPauliKin;

                //double fsi =          ( dTsi*eS + TfS*dB.dCsi )*KRSrho.y;
                //double fsj =          ( dTsj*eS + TfS*dB.dCsj )*KRSrho.y;
                //Vec3d  fp  =  Rij * ( ( dTr *eS - TfS*dB.dCr  )*KRSrho.x*KRSrho.x ); // second *KRSrho.x because dR is not multiplied

                double fsi =        (    dTsi*eS + TfS*dB.dCsi   );
                double fsj =        (    dTsj*eS + TfS*dB.dCsj   );
                Vec3d  fp  =  Rij * ( -( dTr *eS - TfS*dB.dCr  ) ); // second *KRSrho.x because dR is not multiplied

                efpos [i].add( fp ); efpos[j].sub( fp );
                efsize[i]+= fsi*cij; efsize[j]+= fsj*cij;
                efcoef[i]+= T*cj   ; efcoef[j]+= T*ci;
                E += T*cij; // TODO - this should be total kinetic energy for whole orbital not for single basis !!!!
                //if(DEBUG_iter==DEBUG_log_iter) printf( "pauliCrossKinetic[%i,%i]%i Qij %g eij %g dCr %g Fp.x %g cij %g \n", i, j, ij, cij*S, T*cij, dB.dCr, fp.x, cij );
                ij++;
            }
        }
        printf( "E %g T %g eS %g S %g \n", T*eS, T, eS, S );
        return T*eS;
        //return E;
    }

    double pauliKineticChnageVB( int io, int jo ){
        // This is Pauli potential derived from change of Kinetic Energy due to orbital orthogonalization;
        // Detla_T = Tsym + Tanti - T11 - T12
        // where T11 = <psi_1|T|psi_2>, Tsym=<psi_sym||psi_sym>, Tanti=<psi_anti||psi_anti>
        // psi_sym = (psi_1+psi_2)/|psi_1+psi_1| = (psi_1+psi_2)/2(1+S12)
        // psi_sym = (psi_1+psi_2)/|psi_1+psi_1| = (psi_1+psi_2)/2(1-S12)
        // From this it can be derived
        // Detla_T =  ( T11 + T22 - T12/S12 )/( S12^2/(1-S12^2) )
        printf( "EeePaul[%i,%i] ", io, jo );
        double T = 0;
        double S = 0;
        int i0=io*perOrb;
        int j0=jo*perOrb;
        int ij=0;
        for(int i=i0; i<i0+perOrb; i++){
            Vec3d  pi  = epos [i];
            double ci  = ecoef[i];
            double si  = esize[i];
            for(int j=j0; j<j0+perOrb; j++){
                Vec3d pj  = epos[j];
                Vec3d Rij = pj-pi;
                double r2 = Rij.norm2();
                //if(r2>Rcut2) continue;
                double cj  = ecoef[j];
                double sj  = esize[j];
                //double cij = ci*cj*Amp;
                double cij = ci*cj;
                //pairs[ij].get(si,pi, sj,pj);
                //const Gauss::PairDeriv& dB = dBs[ij];
                //const Quat4d&           TD = TDs[ij];
                double  dTr, dTsi, dTsj, dSr, dSsi, dSsj;
                { //                Overlap Integral
                    double dS = Gauss::overlap_s_deriv(  si,pi,  sj,pj,  dSr, dSsi, dSsj );
                    dS*=cij; S+=dS;  // integrate S12 = <psi_1|psi_2>
                }
                if(iPauliModel>0){ // Kinetic Energy integral
                    double dT = Gauss::kinetic_s      (  r2, si, sj,     dTr, dTsi, dTsj );
                    dT*=cij; T+=dT;  // integrate T12 = <psi_1|T|psi_2>
                }
                ij++;
            }
        }
        double E;
        if(iPauliModel==2){
            double T11 = oEs[io];
            double T22 = oEs[jo];
            double S2 = S*S;
            E = (T11 + T22 - 2*T/S)*(S2/(1-S2));
            // ToDo : Here we will need derivatives of kinetic energy according to forces
            //   This can be done if we calculate kinetic forces first and
            E *= const_K_eVA;
        }else{
            E = S*S*KPauliOverlap;
        }
        printf( "E %g \n", E );
        return E;
        //return E;
    }

    double evalCrossOrb(int io, int jo){
        double dEpaul=0, dEexch=0;
        const double R2Safe = sq(0.1);
        // we first project on axuliary density
        Gauss::PairDeriv dBs[perOrb2];
        Gauss::Blob       Bs[perOrb2];
        Gauss::Blob      fBs[perOrb2];
        // transform to auxuliatery overlap-density basis and calculate corresponding derivative transform (Jacobian)
        double T=0,S=0;
        Quat4d TDs[perOrb2];
        for(int i=0; i<perOrb2; i++){ fBs[i].setZero(); }
        bool anti = ( ospin[io] != ospin[jo] );
        if( bEvalPauli && (iPauliModel==1) ){
            printf( "EeePaul[%i,%i] ", io, jo );
            T  = pairOverlapAndKineticDervis( io, jo, Bs, dBs, TDs, S );
            dEpaul = pauliCrossKinetic          ( io, jo, Bs, dBs, TDs, S, T, KRSrho,  anti  );
            //printf( "EeePaul[%i,%i] %g ", io, jo, dEpaul );
        }else{
            if( bEvalPauli && (!anti) ){
                Gauss::iDEBUG=1;
                printf( "EeePaul[%i,%i] ", io, jo );
            }
            S = pairOverlapDervis( io, jo, Bs, dBs );
            Gauss::iDEBUG=0;
            if( bEvalPauli && (!anti) ){
                dEpaul = S*S*KPauliOverlap;
                //printf( "EeePaul[%i,%i]= %g \n", io, jo, dEpaul );
                applyPairForce   ( io, jo, Bs, dBs, S*2*KPauliOverlap );
            }
        }
        // evaluate coulombic terms in axuilary density basis ( i.e. between charge blobs )
        if( bEvalExchange && (!anti) ) dEexch = evalCoulombPair( perOrb2, perOrb2, Bs, Bs, fBs, fBs );
        //printf( "nqOrb %i \n", perOrb2 );
        //for(int i=0; i<nqOrb; i++){ printf( "fBs[%i] fq %g fs %g fp(%g,%g,%g) \n", i, fBs[i].charge, fBs[i].size,    fBs[i].pos.x,fBs[i].pos.x,fBs[i].pos.x ); }
        // transform forces back to wave-function basis, using previously calculated derivatives
        //printf( "EeePaul[%i,%i]= %g \n", io, jo, dEpaul );
        EeePaul+=dEpaul;
        EeeExch+=dEexch;
        forceFromOverlaps( io, jo, dBs, fBs );
        return dEexch + dEpaul;
    }

    double evalCrossTerms(){ // evaluate Energy components given by direct wave-function overlap ( below cutoff Rcut )
        double E = 0;
        for(int io=0; io<nOrb; io++){ for(int jo=0; jo<io; jo++){
        //for(int io=0; io<nOrb; io++){ for(int jo=io+1; jo<nOrb; jo++){
                //if( !checkOrbitalCutoff(io,jo) ) continue; // optimization, not to calculate needless interactions
                E += evalCrossOrb( io, jo);
            }
        }
        //printf( "evalExchange() E=%g \n", E  );
        return E;
    }

    double evalPauli(){ // evaluate Energy components given by direct wave-function overlap ( below cutoff Rcut )
        double E = 0;
        for(int io=0; io<nOrb; io++){
            for(int jo=0; jo<io; jo++){
                E += pauliKineticChnageVB(io,jo);
            }
        }
        return E;
    }

    void evalExchangeCorrelation(){
        // not sure how exactly to do that now - no using real space grid since it is slow
    }


// ==================================================================
// ==================== Interactions with atomic cores
// ==================================================================


double evalAEoverlap( int jo, Gauss::PairDeriv*  dBs, double* Ss, double si, const Vec3d&  pi ){  // Interaction of atomic core with electron wavefunction (Pauli)
    ///const Vec3d  pi   = apos [ia];
    //const Vec3d  abwi = eAbWs[ia]; // ToDo: We actually need only the width
    //double si = eAbWs[ia].z;
    //double si = aPsize[ia];
    //printf( "evalAEoverlap() \n" );
    double S=0;
    int j0=jo*perOrb;
    for(int jj=0; jj<perOrb; jj++){
        int j=j0+jj;
        Vec3d pj  = epos[j];
        Vec3d Rij = pj-pi;
        double r2 = Rij.norm2();
        //if(r2>Rcut2) continue;
        double cj  = ecoef[j];
        double sj  = esize[j];
        Gauss::Blob Bjunk;
        double Si = Gauss::product3DDeriv( si,pi, sj,pj, Bjunk, dBs[jj] );
        Ss[jj] = Si;
        S += Si*cj;
    }
    return S;
}

void applyPairForceAE( int ia, int jo, Gauss::PairDeriv* dBs,  double* Ss, double Amp ){
    int j0=jo*perOrb;
    const Vec3d  pi   = apos [ia];
    //const Vec3d  abwi = eAbWs[ia]; // ToDo: We actually need only the width
    //printf( "applyPairForceAE() \n" );
    //for(int j=j0; j<j0+perOrb; j++){
    for(int jj=0; jj<perOrb; jj++){
        int j=j0+jj;
        Vec3d pj  = epos[j];
        Vec3d Rij = pj-pi;
        double r2 = Rij.norm2();
        //if(r2>Rcut2) continue;
        double cj  = ecoef[j];
        double sj  = esize[j];
        double cij = cj*Amp;
        const Gauss::PairDeriv& dB = dBs   [jj];
        double eij                 = Amp*Ss[jj];
        Vec3d fp = Rij*( dB.dCr*cij);
        aforce[ia].add( fp );
        efpos [j] .sub( fp );
        efsize[j] += dB.dCsj*cij;
        efcoef[j] += eij/cj;
        //if(DEBUG_iter==DEBUG_log_iter) printf( "applyPairForceAE[%i,%i] qij %g e %g fp.x %g | r %g s %g \n", ia,j, cij, eij, fp.x,   sqrt(r2), sj  );
        //ij++;
    }
    //return E;
}

double evalArho( int ia, int jo ){ // Interaction of atomic core with electron density  (Coulomb)
    //printf( " evalArho DEBUG \n");
    const Vec3d  pi = apos  [ia];
    const double qi = aQs   [ia];
    const double si = aQsize[ia];
    //const double si = eAbWs[ia].z;
    int j0  = getRhoOffset(jo);
    int njo = onq[jo];
    double E=0;
    for(int j=j0; j<j0+njo; j++){
    //int j=j0+2;{ // DEBUG
        Vec3d  pj = rhoP[j];
        Vec3d Rij = pj-pi;
        double r2 = Rij.norm2();
        //if(r2>Rcut2) continue;
        double qj = rhoQ[j];
        double sj = rhoS[j];
        Vec3d fp; double fs,qij=qi*qj;
        double e = Gauss::Coulomb( Rij, r2, si, sj, qij, fp, fs );
        /*
        double r    = sqrt(r2 + R2SAFE);
        double s2   = si*si + sj*sj;
        double s    = sqrt(s2);
        double qij = -qi*qj;
        double fr,fs;
        double e  = Gauss::Coulomb( r, s, fr, fs ); // NOTE : remove s*2 ... hope it is fine ?
        // --- Derivatives (Forces)
        Vec3d fp = Rij*( fr * qij );
        fs       *=           qij;
        //fs      *= M_SQRT1_2 * qij;
        */
        aforce[ia].add(fp);
        rhofP [j] .sub(fp);
        rhofS [j] += fs*sj;
        //rhofQ [j] += e*qi*0.0; // Not sure why this is not used here
        E         += e*qij;
        printf( "evalArho[%i,%i] E %g e,q(%g,%g) r %g s(%g,%g) \n", ia,j, e*qij,e,qij, sqrt(r2), si, sj );
        //if(DEBUG_iter==DEBUG_log_iter) printf( "evalArho[%i,%i] qij %g e %g fp.x %g fr %g fs %g | r %g s %g \n", ia,j, qij, e, fp.x,   fr, fs,    s, r );
    }
    return E;
}

double evalAE(){
    Gauss::PairDeriv dBs[perOrb];
    double            Ss[perOrb];
    Eae=0,EaePaul=0;
    for(int ia=0; ia<natom; ia++){
        //const Vec3d  pi   = apos [ia];
        //const double qqi  = aQs  [ia];
        //const Vec3d  abwi = eAbWs[ia];
        for(int jo=0; jo<nOrb; jo++){
            if(bEvalAEPauli){// --- Pauli Overlap
                double S   = evalAEoverlap( jo, dBs, Ss, aPsize[ia], apos[ia] );
                double K   = aPcoef[ia];
                EaePaul   += K*S*S;
                double Amp = K*S*2;
                applyPairForceAE ( ia, jo, dBs, Ss, Amp );
            }
            // ---- Coulomb
            if(bEvalAECoulomb){
                Eae      += evalArho( ia, jo );
            }
        }
    }
    return Eae + EaePaul;
}

/// evaluate Atom-Atom forces
double evalAA(){
    Eaa=0;
    for(int i=0; i<natom; i++){
        const Vec3d  pi   = apos[i];
        const double qi   = aQs[i];
        for(int j=0; j<i; j++){
            Vec3d f = Vec3dZero; // HERE WAS THE ERROR (missing initialization !!!! )
            Vec3d  abw;
            const Vec3d dR  = apos[j] - pi;
            Eaa +=  addAtomicForceQ( dR, f, qi*aQs[j] );
            aforce[j].sub(f);
            aforce[i].add(f);
        }
    }
    return Eaa;
}

// ==================================================================
// ==================== Eval Total Energy And Forces
// ==================================================================

    double eval(){
        double E=0;
        cleanForces();
        if( bEvalCoulomb || (bEvalAECoulomb && bEvalAE) ) clearAuxDens();
        //projectOrbs( true );
        double Ek = projectOrbs( false ); // here we calculate kinetic energy of each orbital and project them to auxuliary charge density basis
        if( bEvalKinetic  ) E+=Ek;
        //reportCharges();
        if( bEvalPauli    ) E += evalPauli();
        //if( bEvalExchange ) E += evalExchange();
        if( bEvalCoulomb  ) E += evalElectrostatICoulomb(); // repulsion between aux-density clouds => should not distinguish density terms here
        if( bEvalAE       ) E += evalAE();
        if( bEvalAA       ) E += evalAA();
        for(int io=0; io<nOrb; io++){
            assembleOrbForces_fromRho(io);
        }
        //E += evalCrossTerms(); // ToDo : This should replace evalPauli and evalExchange()
        //printf( "eval() E=%g \n", E  );
        return E;
    }

    double moveGD( double dt){
        double F2a  = 0;
        double F2ep = 0;
        double F2es = 0;
        double F2ec = 0;
        if(bOptAtom)for(int i=0; i<natom;i++){ apos [i].add_mul( aforce[i], dt );    F2a +=aforce[i].norm2(); }
        if(bOptEPos)for(int i=0; i<nBas; i++){ epos [i].add_mul( efpos [i], dt );    F2ep+=efpos [i].norm2(); }
        if(bOptSize)for(int i=0; i<nBas; i++){ esize[i] +=       efsize[i]* dt  ;    F2es+=sq(efsize[i]);     }
        if(bOptCoef)for(int i=0; i<nBas; i++){ ecoef[i] +=       efcoef[i]* dt  ;    F2ec+=sq(efcoef[i]);     }
        // ToDo : We should out project directions which breaks normalization (!!!) but we can do it later - it is mostly importaint for dynamics, gradient desncent should be fine
        return F2a + F2ep + F2es + F2ec;
    }

// ==================================================================
// ==================== On Grid Utils
// ==================================================================

    /// ToDo : Coulomb integral C{i,j}  on grid can be calculated like this:
    //  for i in MOs :
    //      wf_j = orb2grid( i, grid )  // potencial or density
    //      for j in MOs:
    //          Iij = 0
    //          for kbas in basis_decomposition[j]:
    //              for pos in grid:
    //                  Iij += kbas(pos)^2 * wf_j[pos]^2

    double atomsPotAtPoint( const Vec3d& pos, double s, double Q )const{
        //Gauss::PairDeriv dBs[perOrb];
        //double            Ss[perOrb];
        double E =0;
        Vec3d fp; double fs;
        for(int ia=0; ia<natom; ia++){
            if(bEvalAEPauli){
                Vec3d pij; double sij;
                double S = Gauss::product3D_s_new( aPsize[ia], apos[ia], s, pos, sij, pij );
                //printf( "atomsPotAtPoint[%i] Pauli S %g c %g pij.x %g sij %g aPsize %g \n", ia, S, aPcoef[ia], pij.x, sij, aPsize[ia] );
                E += aPcoef[ia]*S*S;
            }
            if(bEvalAECoulomb){
                Vec3d Rij = pos - apos[ia];
                // Coulomb( const Vec3d& Rij, double r2, double si, double sj, double qij, Vec3d& fp, double& fs ){
                E += Gauss::Coulomb( Rij, Rij.norm2(), aQsize[ia], s, 1, fp, fs )*aQs[ia]*Q;
            }
        }
        return E;
    }

    double* atomsPotAtPoints( int n, Vec3d* ps, double* out=0, double s=0.0, double Q=1.0 )const{
        if(out==0){ out = new double[n]; };
        for(int i=0; i<n; i++){
            out[i] = atomsPotAtPoint( ps[i], s, Q );
            //printf( "atomsPotAtPoints[%i/%i] %g @(%g,%g,%g) Q[0] %g \n", i, n, out[i], ps[i].x,ps[i].y,ps[i].z, aQs[0] );
        }
        return out;
    }

    double orbAtPoint( int io, const Vec3d& pos )const{
        int i0     = getOrbOffset( io );
        double wfsum=0;
        for(int i=0; i<perOrb; i++){
            int i_=i0+i;
            Vec3d dR  = pos - epos[i_];
            double r2 = dR.norm2();
            wfsum += Gauss::bas3D_r2( r2, esize[i_] ) * ecoef[i_];
        }
        return wfsum;
    }

    double* orbAtPoints( int io, int n, Vec3d* ps, double* out=0 )const{
        if(out==0){ out = new double[n]; };
        for(int i=0; i<n; i++){
            out[i] = orbAtPoint( io, ps[i] );
            //printf( "orbAtPoints[%i/%i] %g @(%g,%g,%g) \n", i, n, out[i], ps[i].x,ps[i].y,ps[i].z );
        }
        return out;
    }

    double orb2grid( int io, const GridShape& gridShape, double* buff )const{
        int i0     = getOrbOffset( io );
        Vec3d*  Ps = epos +i0;
        double* Cs = ecoef+i0;
        double* Ss = esize+i0;
        //printf( "DEBUG orb2grid io %i i0 %i perOrb %i \n", io, i0, perOrb );
        //for(int i=0; i<perOrb; i++){ printf( "wf[%i] C(%e) P(%g,%g,%g) s %g 1/|f^2| %g \n", i, Cs[i], Ps[i].x,Ps[i].y,Ps[i].z, Ss[i], Gauss::sqnorm3Ds( Ss[i] ) ); }
        //printf( "DEBUG norm3Ds(1) %g %g \n", Gauss::norm3Ds( 1.0 ), pow( M_PI*2,-3./2) );
        return evalOnGrid( gridShape, [&](int ig, const Vec3d& pos, double& res){
            double wfsum = 0.0;
            for(int i=0; i<perOrb; i++){
                Vec3d dR  = pos - Ps[i];
                double r2 = dR.norm2();
                double si = Ss[i];
                wfsum += Gauss::bas3D_r2( r2, si ) * Cs[i];
                //printf( "ig [%i] dR(%g,%g,%g) wf %g  r2 si C %g %g %g \n", ig, dR.x,dR.y,dR.z, wfsum, r2, si, Cs[i] );
                // ToDo : Fast Gaussian ?
            }
            res += wfsum*wfsum;
            //printf( "ig [%i] wf %g \n", ig, wfsum );
            buff[ig] = wfsum;
        });
    }

    void orb2xsf( const GridShape& grid, int iorb, const char* fname )const{
        int ng = grid.n.totprod();
        double  dV = grid.voxelVolume();
        double* buff = new double[ ng ];
        double Q = orb2grid( 0, grid, buff );
        //printf( "orb2xsf Q  %g \n", Q );
        grid.saveXSF( fname, buff, -1 );
    }

    double rho2grid( int io, const GridShape& gridShape, double* buff )const{
        int i0 = getRhoOffset(io);
        int ni = onq[io];
        //printf( "DEBUG rho2grid i0 %i ni %i \n", i0, ni  );
        //printf( "DEBUG DensOverlapOrbPair i0,j0  %i,%i ni,nj %i,%i \n", i0,j0, ni, nj );
        double renorm = Gauss::norm3Ds(1);
        Vec3d*  Ps = rhoP + i0;
        double* Cs = rhoQ + i0;
        double* Ss = rhoS + i0;
        //printf( "DEBUG orb2grid io %i i0 %i perOrb %i \n", io, i0, perOrb );
        //for(int i=0; i<ni; i++){ printf( "rho[%i] C(%e) P(%g,%g,%g) s %g 1/|f^2| %g \n", i, Cs[i], Ps[i].x,Ps[i].y,Ps[i].z, Ss[i], Gauss::sqnorm3Ds( Ss[i] ) ); }
        //printf( "DEBUG norm3Ds(1) %g %g \n", Gauss::norm3Ds( 1.0 ), pow( M_PI*2,-3./2) );
        return evalOnGrid( gridShape, [&](int ig, const Vec3d& pos, double res){
            double rho_sum = 0.0;
            for(int i=0; i<ni; i++){
                Vec3d dR  = pos - Ps[i];
                double r2 = dR.norm2();
                double si = Ss[i];
                double ci = Cs[i];
                //if(i==2) ci*=1.28;
                rho_sum += Gauss::rho3D_r2( r2, si ) * ci;
                // ToDo : Fast Gaussian ?
            }
            buff[ig] = rho_sum;
        });
    }

    double hartree2grid( int io, const GridShape& gridShape, double* buff )const{
        int i0 = getRhoOffset(io);
        int ni = onq[io];
        //printf( "DEBUG rho2grid i0 %i ni %i \n", i0, ni  );
        //printf( "DEBUG DensOverlapOrbPair i0,j0  %i,%i ni,nj %i,%i \n", i0,j0, ni, nj );
        double renorm = Gauss::norm3Ds(1);
        Vec3d*  Ps = rhoP + i0;
        double* Cs = rhoQ + i0;
        double* Ss = rhoS + i0;
        //printf( "DEBUG orb2grid io %i i0 %i perOrb %i \n", io, i0, perOrb );
        //for(int i=0; i<ni; i++){ printf( "rho[%i] C(%e) P(%g,%g,%g) s %g 1/|f^2| %g \n", i, Cs[i], Ps[i].x,Ps[i].y,Ps[i].z, Ss[i], Gauss::sqnorm3Ds( Ss[i] ) ); }
        //printf( "DEBUG norm3Ds(1) %g %g \n", Gauss::norm3Ds( 1.0 ), pow( M_PI*2,-3./2) );
        return evalOnGrid( gridShape, [&](int ig, const Vec3d& pos, double res){
            double v_sum = 0.0;
            for(int i=0; i<ni; i++){
                Vec3d dR  = pos - Ps[i];
                double r2 = dR.norm2();
                double si = Ss[i];
                double ci = Cs[i];
                //if(i==2) ci*=1.28;
                //rho_sum += Gauss::elpot( sqrt(r2)/si ) * ci;
                v_sum += erfx_e6( sqrt(r2)/(si*M_SQRT2) ) * (ci * const_El_eVA);
                // ToDo : Fast Gaussian ?
            }
            buff[ig] = v_sum;
        });
    }



bool loadFromFile( char const* filename, bool bCheck ){
    //printf(" filename: >>%s<< \n", filename );
    FILE * pFile;
    pFile = fopen (filename,"r");
    if( pFile == NULL ){
        printf("ERROR in CLCFGO::loadFromFile(%s) : No such file !!! \n", filename );
        return -1;
    }
    int ntot;
    const int nbuff = 1024;
    char buff[nbuff]; char* line;
    //fscanf (pFile, " %i \n", &ntot );
    int natom_=0, nOrb_=0, perOrb_=0; bool bClosedShell=0;
    line=fgets(buff,nbuff,pFile);
    sscanf (line, "%i %i %i\n", &natom_, &nOrb_, &perOrb_, &bClosedShell );
    //printf("na %i ne %i perORb %i \n", natom, nOrb, perOrb_);
    //printf("na %i ne %i perORb %i \n", natom_, nOrb_, perOrb_ );
    if(bClosedShell) nOrb_*=2;
    realloc( natom_, nOrb_, perOrb_, 1 );
    double Qasum = 0.0;
    for(int i=0; i<natom; i++){
        double x,y,z;
        double Q,sQ,sP,cP;
        fgets( buff, nbuff, pFile); //printf( "fgets: >%s<\n", buf );
        int nw = sscanf (buff, "%lf %lf %lf %lf %lf %lf %lf", &x, &y, &z,    &Q, &sQ, &sP, &cP );
        //printf( "atom[%i] p(%g,%g,%g) Q %g sQ %g sP %g cP %g \n", i, x, y, z,    Q, sQ, sP, cP );
        apos  [i]=(Vec3d){x,y,z};
        aQs   [i]=Q;
        aQsize[i]=sQ;
        aPsize[i]=sP;
        aPcoef[i]=cP;
        Qasum += Q;
    }
    int nBasRead = nBas;
    if( bClosedShell ) nBasRead/=2;
    for(int i=0; i<nBasRead; i++){
        double x,y,z;
        double s,c;
        int spin;
        fgets( buff, nbuff, pFile); // printf( "fgets: >%s<\n", buf );
        int nw = sscanf (buff, "%lf %lf %lf %lf %lf %i", &x, &y, &z,  &s, &c, &spin );
        epos [i]=(Vec3d){x,y,z};
        esize[i]=s;
        ecoef[i]=c;
        int io=i/perOrb;
        if( !bClosedShell ){ if(nw>5)ospin[io]=spin; }else{ ospin[io]=1; };
        //printf( "ebasis[%i,%i|%i] p(%g,%g,%g) s %g c %g spin %i | nw %i io %i \n", i/perOrb, i%perOrb,i, x, y, z,  s, c, spin,  nw, io  );
    }
    if( bClosedShell ){
        for(int i=0; i<nBasRead; i++){
            int j = i+nBasRead;
            epos [j]=epos[i];
            esize[j]=esize[i];
            ecoef[j]=ecoef[i];
            ospin[j/perOrb]=-1;
        }
    }
    //printf( "Qtot = %g (%g - 2*%i) \n",  Qasum - nOrb, Qasum, nOrb );
    fclose (pFile);
    return 0;
}

void printSetup(){
    printf("===CLCFGO Setup :\n");
    printf("iPauliModel %i \n", iPauliModel );
    printf("bRescaleKinetic %i \n", bRescaleKinetic );
    printf("KPauliKin %g KPauliOverlap %g \n", KPauliKin, KPauliOverlap );
    printf("KRSrho (%g,%g,%g) \n", KRSrho.x, KRSrho.y, KRSrho.z );
    printf("Rcut %g RcutOrb %g \n",  Rcut,  RcutOrb );
    printf("bOptEPos %i \n",bOptAtom );
    printf("bOptEPos %i \n",bOptEPos );
    printf("bOptSize %i \n",bOptSize );
    printf("bOptEPos %i \n",bOptCoef );
    printf("bNormalize     %i \n",bNormalize     );
    printf("bEvalKinetic   %i \n",bEvalKinetic   );
    printf("bEvalCoulomb   %i \n",bEvalCoulomb   );
    printf("bEvalPauli     %i \n",bEvalPauli     );
    printf("bEvalExchange  %i \n",bEvalExchange  );
    printf("bEvalAE        %i \n",bEvalAE        );
    printf("bEvalAECoulomb %i \n",bEvalAECoulomb );
    printf("bEvalAEPauli   %i \n",bEvalAEPauli   );
    printf("bEvalAA        %i \n",bEvalAEPauli   );
}

void printAtoms(){
    printf( "===CLCFGO::printAtoms()\n");
    for(int i=0; i<natom; i++){
        printf( "eFF::atom[%i] p(%g,%g,%g)[A] Q(%g[e],%g[A])  Pauli(%g[eV],%g[A]) \n", i, apos[i].x,apos[i].y,apos[i].z,  aQs[i], aQsize[i], aPcoef[i], aPsize[i] );
    }
}

void printElectrons(){
    printf( "===CLCFGO::printElectrons()\n");
    for(int io=0; io<nOrb; io++){
        printf( ">> orb[%i] spin %i \n", io, ospin[io] );
        for(int j=0; j<perOrb; j++){
            int ie=io*perOrb+j;
            printf( "e[%i,%i|%i] p(%g,%g,%g)[A] size %g coef %g \n", io,j,ie, epos[ie].x,epos[ie].y,epos[ie].z, esize[ie], ecoef[ie] );
        }
    }
}



// ===========================================================================================================================
// ==================== Old Versions of Functions - Should be removed
// ===========================================================================================================================

// ===========================================================================================================================
// ==================== Old Version of overlap and Pauli Evaluation - ToDo : Remove This
// ===========================================================================================================================

/*

    double DensOverlapOrbPair( int io, int jo ){
        int i0 = getRhoOffset(io);
        int j0 = getRhoOffset(jo);
        int ni = onq[io];
        int nj = onq[jo];
        //printf( "DEBUG DensOverlapOrbPair i0,j0  %i,%i ni,nj %i,%i \n", i0,j0, ni, nj );
        double Srho = 0;

        double renorm = Gauss::norm3Ds(1);
        for(int i=i0; i<i0+ni; i++){
            Vec3d  pi = rhoP[i];
            double qi = rhoQ[i];
            double si = rhoS[i];
            for(int j=j0; j<j0+nj; j++){
                Vec3d  pj = rhoP[j];
                Vec3d Rij = pj-pi;
                double r2 = Rij.norm2();

                if(r2>Rcut2) continue;

                double qj = rhoQ[j];
                double sj = rhoS[j];

                double dSr,dSsi,dSsj;
                //double Sij = getOverlapSGauss( r2, si*M_SQRT2, sj*M_SQRT2, dSr, dSsi, dSsj )*M_SQRT1_2;
                double Sij = getOverlapSGauss( r2, si*M_SQRT2, sj*M_SQRT2, dSr, dSsi, dSsj )*renorm;
                //double Sij = getOverlapSGauss( r2, si*2, sj*2, dSr, dSsi, dSsj );
                //double Sij = getOverlapSGauss( r2, si, sj, dSr, dSsi, dSsj );

                //printf( "i,j %i %i  r %g si,j(%g,%g) qi*qj*S(%g,%g|%g):%g \n", i,j, sqrt(r2),  si,sj, qi,qj,Sij, Sij*qi*qj );

                Srho += Sij*qi*qj;

            }
        }
        return Srho;
    }

    //inline double evalShotRange( Vec3d Rij, double* Is, int i, int j ){
    inline double evalOverlap( int io, int jo ){
        /// Pauli Energy between two orbitals depend on IKinetic overlap  ~S^2
        ///   See addPauliGauss() in common/molecular/InteractionsGauss.h
        int i0 = getOrbOffset(io);
        int j0 = getOrbOffset(jo);
         /// ToDo : we may allow different number of orbitals per electron later (?)
        double Ssum = 0;
        //printf("============ i0 %i j0 %i perOrb %i \n", i0, j0, perOrb);
        for(int i=i0; i<i0+perOrb; i++){
            Vec3d  pi = epos[i];
            double si = esize[i];
            double ci = ecoef[i];

            //Vec3d& fpi = efpos [i];
            //double fsi = efsize[i];
            //double fci = efcoef[i];
            for(int j=j0; j<j0+perOrb; j++){
                /// ToDo: Rcut may be read from
                Vec3d Rij = epos[j]-pi;
                double r2 = Rij.norm2();
                if(r2>Rcut2)continue; //{ fij = Vec3dZero; return 0; }
                double sj = esize[j];
                double cj = ecoef[j];
                double dSr, dSsi, dSsj;
                //double Sij = getOverlapSGauss( r2, si, sj, dSr, dSsi, dSsj );
                //double Sij = getOverlapSGauss( r2, si*2, sj*2, dSr, dSsi, dSsj );
                double Sij = getOverlapSGauss( r2, si*M_SQRT2, sj*M_SQRT2, dSr, dSsi, dSsj );   dSr=-dSr; //dSsi=dSsi; dSsj=dSsj;
                /// NOTE : si,sj scaled by sqrt(2) because they are made for density widths not wave-functions width
                //printf( "[%i,%i]x[%i,%i] r,S %g %g ci,j  %g %g -> %g \n", io,i, jo,j, sqrt(r2), Sij, ci, cj, Sij*ci*cj );
                double cij = ci*cj;
                Ssum += Sij*cij;

                dSsi*=M_SQRT2;
                dSsj*=M_SQRT2;

                // --- Derivatives ( i.e. Forces )
                Vec3d fij = Rij*(dSr*cij);
                efpos [i].add( fij     ); efpos[j].sub ( fij );
                efsize[i]+= dSsi*cij ; efsize[j]+= dSsj*cij ;
                efcoef[i]+= Sij*cj   ; efcoef[j]+= Sij*ci  ;
            }
        }
        return Ssum;
    }

    double orbNorm(int io )const{
        int i0=getOrbOffset(io);
        double Ssum=0;
        int ii=0;
        for(int i=i0; i<i0+perOrb; i++){
            Vec3d  pi  = epos [i];
            double si  = esize[i];
            double ci  = ecoef[i];
            for(int j=i0; j<i; j++){
                Vec3d Rij = epos[j]-pi;
                double r2 = Rij.norm2();
                if(r2>Rcut2)continue; //{ fij = Vec3dZero; return 0; }
                double sj = esize[j];
                double cj = ecoef[j];
                double dSr, dSsi, dSsj;
                double Sij = getOverlapSGauss( r2, si*M_SQRT2, sj*M_SQRT2, dSr, dSsi, dSsj );
                // NOTE : si,sj scaled by sqrt(2) because they are made for density widths not wave-functions width
                Ssum += Sij*ci*cj;
            }
        }
        return Ssum;
    }

    //inline double evalShotRange( Vec3d Rij, double* Is, int i, int j ){
    inline double evalPauli( int io, int jo ){
        /// Pauli Energy between two orbitals depend on IKinetic overlap  ~S^2
        ///  See addPauliGauss() in common/molecular/InteractionsGauss.h
        int i0 = getOrbOffset(io);
        int j0 = getOrbOffset(jo);
         // ToDo : we may allow different number of orbitals per electron later (?)
         double E = 0;
         for(int i=i0; i<i0+perOrb; i++){
            Vec3d  pi = epos[i];
            double si = esize[i];
            for(int j=j0; j<j0+perOrb; j++){
                // ToDo: Rcut may be read from
                Vec3d Rij = epos[j]-pi;
                double r2 = Rij.norm2();
                if(r2>Rcut2)continue; //{ fij = Vec3dZero; return 0; }
                double sj = esize[j];

                double dSr, dSsi, dSsj;

                double Sij = getOverlapSGauss( r2, si, sj, dSr, dSsi, dSsj );

                double ci = ecoef[i];
                double cj = ecoef[j];
                double cij = ci*cj;

                // TODO: We need model for Pauli !!!!
                double eP = 0;
                double fP = 0; //
                Vec3d fij;

                //fij.set_mul( Rij, cij*fP/r ); // ToDo:  evaluate fP properly

                // ToDo:  ecfs should be set as well
                //ecfs[i] +=dcP;
                //ecfs[j] +=dcP;
                //epfs[i].add(fij);
                //epfs[j].sub(fij);
                //return eP*cij;
                E += eP*cij;
            }
        }
        return E;
    }


    inline double CrossOverlap( int io, int jo ){
        // NOTE : Later we need to evaluate some function F(S(x)) ( e.g. polynominal of cross-overlap  )
        //          => we need chain derivs   dF(S(x))/dx = (dF/dS)*(dS/dx)
        //        Most probably it will be integral of overlap density Sum{rho} = Sum{Sij^2}
        //         in such case derivative    dF(S(xi))/dx =  Stot * (dSij/dxi)
        int i0 = getOrbOffset(io);
        int j0 = getOrbOffset(jo);
         double E = 0;
         double SOij = 0.0;
         double CAmp = 1.0;
         for(int i=i0; i<i0+perOrb; i++){
            Vec3d  pi = epos [i];
            double si = esize[i];
            double ci = ecoef[i];
            for(int j=j0; j<j0+perOrb; j++){
                Vec3d pj  = epos[j];
                Vec3d Rij = pj-pi;
                double r2 = Rij.norm2();
                if(r2>Rcut2)continue; //{ fij = Vec3dZero; return 0; }
                double sj = esize[j];
                double cj = ecoef[j];
                double cij = ci*cj;

                double Sij;
                Vec3d  p;
                double s;
                double dssi,dssj;
                Vec3d  dxsi,dxsj;
                double dxxi,dxxj;
                double dSr;
                double dSsi,dSsj;

                Sij = Gauss::product3D_s_deriv(
                    si,   pi,
                    sj,   pj,
                    s ,   p ,
                    dssi, dssj,
                    dxsi, dxsj,
                    dxxi, dxxj,
                    dSsi, dSsj, dSr
                );

                double eij=CAmp*Sij;
                E   += eij*cij;
                Vec3d Fp = Rij*(dSr*cij);
                efpos [i].add( Fp ); efpos[j].sub( Fp );
                efsize[i]+= dSsi*cij; efsize[j]+= dSsj*cij;
                efcoef[i]+= eij*cj ;  efcoef[j]+= eij*ci;

            }
        }
        return E;
    }

    double evalPauli(){ // evaluate Energy components given by direct wave-function overlap ( below cutoff Rcut )
        double E = 0;
        for(int io=0; io<nOrb; io++){
            for(int jo=0; jo<nOrb; jo++){
                //if( !checkOrbitalCutoff(io,jo) ) continue; // optimization, not to calculate needless interactions
                //E += evalPauli( io, jo );
                E += CrossOverlap( io, jo ); // WARRNING : This is just DEBUG, later we need to use some function of overlap
            }
        }
        return E;
    }

    double ExchangeOrbPair(int io, int jo){

        // ToDo : Remove This   --   This was replaced by other functions like
        //    pairOverlapDervis()
        //    pairKineticDervis()

        // It may be actually easier to evaluate Hartree-Fock like exchange, because it is linear
        // see:   https://en.wikipedia.org/wiki/Exchange_interaction#Exchange_of_spatial_coordinates
        // https://chemistry.stackexchange.com/questions/61176/what-is-the-exchange-interaction
        // => Exchange integral should be rather easy to compute, just by calculating overlap between the two functions   K = <rho_ij(r1)| 1/|r1-r2| |rho_ij(r2)>
        // scheme:
        //  1) for pair of orbital, construct auxiliary rho_ij (by expansion on auxuliary functions)
        //  2) callculate self-convolution of this auxuliary function

        //if( ~(0==io) ) return 0;
        //if(DEBUG_iter==DEBUG_log_iter)printf( " -- ExchangeOrbPair[%i,%i] \n", io, jo );

        const double R2Safe = sq(0.1);
        // we first project on axuliary density

        double cijs[nqOrb];

        Vec3d  fPs[nqOrb];
        double fSs[nqOrb];
        double fQs[nqOrb];

        Gauss::PairDerivs pairs[nqOrb];

        int i0=io*perOrb;
        int j0=jo*perOrb;
        int ij=0;
        // --- Project   rho_ij(r)  to temp axuliary basis
        //if(DEBUG_iter==DEBUG_log_iter) reportOrbitals();
        for(int i=i0; i<i0+perOrb; i++){
            Vec3d  pi  = epos [i];
            double ci  = ecoef[i];
            double si  = esize[i];
            for(int j=j0; j<j0+perOrb; j++){
                Vec3d pj  = epos[j];
                Vec3d Rij = pj-pi;
                double r2 = Rij.norm2();
                //if(r2>Rcut2) continue;
                double cj  = ecoef[j];
                double sj  = esize[j];
                double cij = ci*cj;
                cijs [ij] = cij;
                pairs[ij].get(si,pi, sj,pj);
                //if(DEBUG_iter==DEBUG_log_iter) printf( "Exchange:Project[%i,%i] ss(%g,%g) cij %g Sij %g Qij %g r %g x(%g,%g) \n", i, j, si, sj, cij, pairs[ij].C, cij*pairs[ij].C, sqrt(r2), pi.x, pj.x );
                ij++;
            }
        }
        // --- self-coulomb  of   rho_ij(r)
        double E = 0;
        for(int i=0; i<nqOrb; i++){
            fPs[i] = Vec3dZero;
            fSs[i] = 0;
            fQs[i] = 0;
        }
        for(int i=0; i<nqOrb; i++){
        //int i=0;{
            //Vec3d  pi = posij[i];
            //double qi = rhoij[i];
            //double si = szij [i];
            const Gauss::PairDerivs& A = pairs[i];
            double ci = cijs[i];
            double qi = ci*A.C;
            for(int j=0; j<=i; j++){
            //int j=1; {
                const Gauss::PairDerivs& B = pairs[j];
                Vec3d Rij = B.p - A.p;
                double r2 = Rij.norm2();
                double cj = cijs[j];
                double qj = cj*B.C;

                double r  = sqrt(r2 + R2SAFE );
                double s2 = A.s*A.s + B.s*B.s;
                double s  = sqrt(s2);
                double fr,fs;

                double qij = qi*qj*2;

                //printf( "//Exchange:Coulomb[%i,%i] r %g s %g \n", i, j, r, s  );
                double e  = Gauss::Coulomb( r, s, fr, fs ); // NOTE : remove s*2 ... hope it is fine ?
                //printf( "Exchange:Coulomb[%i,%i] qij %g E %g fr %g fs %g \n", i, j, r, s, qij, E*qij, fs, fs );

                // --- Derivatives (Forces)
                Vec3d fp = Rij*(-fr * qij);
                fs       *=           qij;

                //if(DEBUG_iter==DEBUG_log_iter) printf( "ExchangeOrbPair:Coulomb[%i,%i] qij %g e %g fx %g fs %g | r %g s %g fr %g fs %g  \n", i,j, qij, e, fp.x, fs,    r, s, fr, fs );

                fPs[i].add(fp);    fPs[j].sub(fp);
                fSs[i] -= fs*A.s;  fSs[j] -= fs*B.s;
                fQs[i] += e*qj*2;  fQs[j] += e*qi*2; // ToDo : need to be made more stable ... different (qi,qj)
                E      += e*qij;
            }
        }

        //if(DEBUG_iter==DEBUG_log_iter){for(int i=0; i<nqOrb; i++){printf( "[%i] e %g fx %g | x %g q %g \n", i,  fQs[i], fPs[i].x,     pairs[i].p.x, pairs[i].C*cijs[i] );}}

        //printf( " Exchange:Coulomb DONE \n" );
        //printf( " Exchange:Coulomb DONE \n" );
        // --- From Rho
        ij=0;
        for(int i=i0; i<i0+perOrb; i++){
            Vec3d  pi  = epos [i];
            double ci  = ecoef[i];
            double si  = esize[i];
            for(int j=j0; j<j0+perOrb; j++){
                Vec3d pj  = epos[j];
                Vec3d Rij = pj-pi;
                double r2 = Rij.norm2();
                //if(r2>Rcut2) continue;
                double cj  = ecoef[j];
                double sj  = esize[j];
                double cij = ci*cj;
                //cijs[ij] = cij;
                Gauss::PairDerivs& pair = pairs[ij];
                //if(DEBUG_iter==DEBUG_log_iter) printf( "Exchange:fromRho[%i,%i] x(%g,%g) ", i, j, pi.x, pj.x );
                pair.backForce( Rij, cij, fQs[ij], fPs[ij], fSs[ij],  efpos[i], efpos[j], efsize[i], efsize[j] );
                ij++;
            }
        }
        //printf( " Exchange:FromRho DONE \n" );
        //if(DEBUG_iter==DEBUG_log_iter)printf( "ExchangeOrbPair() E=%g \n", E  );
        return E;
    }

    double evalExchange(){ // evaluate Energy components given by direct wave-function overlap ( below cutoff Rcut )
        double E = 0;
        for(int io=0; io<nOrb; io++){ for(int jo=0; jo<io; jo++){
        //for(int io=0; io<nOrb; io++){ for(int jo=io+1; jo<nOrb; jo++){
                //if( !checkOrbitalCutoff(io,jo) ) continue; // optimization, not to calculate needless interactions
                E += ExchangeOrbPair( io, jo);
            }
        }
        //printf( "evalExchange() E=%g \n", E  );
        return E;
    }


*/








};

/// @}

#endif
