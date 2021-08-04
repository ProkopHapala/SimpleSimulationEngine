
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
{ 0.,  0.1, 0.1, 2.0 }, // 2 He
{ 1.,  0.1, 0.1, 2.0 }, // 3 Li
{ 2.,  0.1, 0.1, 2.0 }, // 4 Be
{ 3.,  0.1, 0.1, 2.0 }, // 5 B
{ 4.,  0.1, 0.1, 2.0 }, // 6 C
{ 5.,  0.1, 0.1, 2.0 }, // 7 N
{ 6.,  0.1, 0.1, 2.0 }, // 8 O
{ 7.,  0.1, 0.1, 2.0 }, // 9 F
};

constexpr static const Vec3d KRSrho = { 1.125, 0.9, 0.2 }; ///< eFF universal parameters

    //double DEBUG_rescale_Q =0;

    bool bDealoc=false;

    bool bNormalize     = true;
    bool bNormForce     = false; // remove/constrain normal force
    bool bEvalKinetic   = true;
    bool bEvalCoulomb   = true;
    bool bEvalPauli     = true;
    bool bEvalExchange  = true;
    bool bEvalAECoulomb = true;
    bool bEvalAEPauli   = true;
    bool bEvalAE        = true;
    bool bEvalAA        = true;

    bool bRescaleKinetic = true;

    bool bOptAtom = true;
    bool bOptEPos = true;
    bool bOptSize = true;
    bool bOptCoef = true;

    int  iPauliModel    = 2;

    double KPauliOverlap = 50.0; // ToDo : This is just "bulgarian constant" for now
    double KPauliKin     = 50.0; // ToDo : Not sure if we should use this - perhaps this model of pauli energy should work "ab inition"

    double Rcut    =6.0;  ///< cutoff beyond which two basis functions chi has no overlap
    double RcutOrb =9.0;  ///< cutoff beoud which orbital (i.e. localized cluster of basis functions) has no overlap
    double Rcut2     =Rcut*Rcut;
    double RcutOrb2  =RcutOrb*RcutOrb;

    double minSize=0.1,maxSize=3.0;

    double Ek=0,Eee=0,EeePaul=0,EeeExch=0,Eae=0,EaePaul=0,Eaa=0, Etot=0; ///< different kinds of energy

    int natypes=0;
    int natom  =0; ///< number of atoms (nuclei, ions)
    int nOrb   =0; //!< Brief number of single-electron orbitals in system
    int nBas   =0; ///< number of basis functions
    int ndofs=0;
    int perOrb =0; //!< Brief number of spherical functions per orbital
    int perOrb2=0;
    // this is evaluated automaticaly
    int nqOrb  =0; ///< number of charges (axuliary density elements) per orbital
    int nQtot  =0; ///< total number of charge elements

    double* dofs=0;
    double* fdofs=0;

    // atoms (ions)
    Vec3d*  apos   =0;  ///< [A] positioon of atoms
    Vec3d*  aforce =0;  ///< [eV/A] force on atoms
    /*
    double* aQs    =0;  ///< [e] charge of atomic core ( Q_nucleus - Q_core_electrons )
    double* aPcoef =0;  ///< [eV] coeficient of pauli repulsion of core electrons
    double* aPsize =0;  ///< [A] size of core in Pauli interaction
    double* aQsize =0;  ///< [A] size of core in coulombic interaction
    */
    Quat4d * aPars = 0; ///< electron params { x=Q,y=sQ,z=sP,w=cP }
    int*    atype  =0;  ///< type of atom (in particular IKinetic pseudo-potential)
    //Vec3d  * aAbWs =0; ///< atomic   parameters (amplitude, decay, width)
    //Vec3d  * eAbWs =0; ///< electron parameters (amplitude, decay, width)






    // orbitals
    Vec3d*  opos  =0;   ///< [A] store positions for the whole orbital
    Vec3d*  odip  =0;   ///< [eA] Axuliary array to store dipoles for the whole orbital
    double* oEs   =0;   ///< [eV] orbital energies
    double* oQresc=0;   ///< probably not needed byt makes debugging of outProjectNormalForces() easier
    double* oQs   =0;   ///< [e] total charge in orbital before renormalization (just DEBUG?)
    int*    onq   =0;   ///< number of axuliary density functions per orbital, should be equal to nQorb
    int*    ospin =0;
    int*    ofix  =0;
    // ToDo : all this blobs ( pos, size, coef should be probably grouped together in Blob{pos,size,coef} )
    // --- Wave-function components for each orbital
    Vec3d*  epos  =0; ///< [A] position of spherical function for expansion of orbitals
    double* esize =0; ///< [A] spread of gassian basisfunction
    double* ecoef =0; ///< [e^.5] expansion coefficient of expansion of given orbital
    // --- Forces acting on wave-functions components
    Vec3d*  efpos  =0; ///< [eV/A]  force acting on position of orbitals   d_E/d_epos[i]
    double* efsize =0; ///< [eV/A]  force acting on combination coefficnet of orbitals  d_E/d_esize[i]
    double* efcoef =0; ///< [e^.5V] force acting on size of gaussians d_E/d_ecoef[i]
    // --- normal forces ( to constrain normality of wavefunction )
    Vec3d*  enfpos  =0;
    double* enfsize =0;
    double* enfcoef =0;
    // --- Kinetic Energy Forces
    Gauss::Blob* fTs = 0;
    //Vec3d*  eTfpos  =0;
    //double* eTfsize =0;
    //double* eTfcoef =0;

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
        bDealoc=true;
        // ToDo : We may automatize this somehow ????
        // atoms
        nBas    = nOrb_ * perOrb_;
        ndofs = natom_*3 + nBas*(3+1+1);
        _realloc(  dofs  ,ndofs );
        _realloc( fdofs  ,ndofs );
        printf( " nAtom %i nOrb %i perOrb %i nBas %i nDOFs %i \n", natom_, nOrb_, perOrb_, nBas, ndofs );

        apos  =(Vec3d*) dofs;
        aforce=(Vec3d*)fdofs;
        epos  =(Vec3d*) (  dofs + (natom_*3)          );
        esize =         (  dofs + (natom_*3) + 3*nBas );
        ecoef =         (  dofs + (natom_*3) + 4*nBas );
        efpos =(Vec3d*) ( fdofs + (natom_*3)          );
        efsize=         ( fdofs + (natom_*3) + 3*nBas );
        efcoef=         ( fdofs + (natom_*3) + 4*nBas );


        printf( "CLCFGO::realloc *X %li %li \n", dofs+(natom_*3), dofs+((natom_*3)+(5*nBas)-1) );

        if( natom != natom_ ){
            natom = natom_;
            //_realloc( apos  ,natom );
            //_realloc( aforce,natom );
            //_realloc( aQs   ,natom );
            //_realloc( aAbWs ,natom);
            //_realloc( eAbWs ,natom);
            //_realloc( aPcoef ,natom);
            //_realloc( aPsize ,natom);
            //_realloc( aQsize ,natom);
            _realloc( aPars, natom );
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
            _realloc( oQresc , nOrb); ///< probably not needed byt makes debugging of outProjectNormalForces() easier
            _realloc( oQs , nOrb);
            _realloc( onq , nOrb);
            _realloc( ofix, nOrb);
            // --- Wave-function components for each orbital
            //_realloc( epos , nBas );
            //_realloc( esize, nBas );
            //_realloc( ecoef, nBas );
            // --- Forces acting on wave-functions components
            //_realloc( efpos ,   nBas );
            //_realloc( efsize,   nBas );
            //_realloc( efcoef,   nBas );
            // --- normal forces ( to constrain normality of wavefunction )
            _realloc( enfpos ,   nBas );
            _realloc( enfsize,   nBas );
            _realloc( enfcoef,   nBas );
            // --- Kinetic Energy Forces
            _realloc( fTs ,   nBas );
            //_realloc( eTfpos ,   nBas );
            //_realloc( eTfsize,   nBas );
            //_realloc( eTfcoef,   nBas );

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

    void dealloc(){
        // atoms
        delete [] apos  ;
        delete [] aforce;
        //delete [] aQs   ;
        //delete [] aPcoef;
        //delete [] aPsize;
        //delete [] aQsize;
        delete [] aPars;
        delete [] atype ; // not used  now
        // orbitals
        delete [] ospin;
        delete [] opos;
        delete [] odip;
        delete [] oEs ;
        delete [] oQresc; // ToDo: probably not needed, but makes debugging easier
        delete [] oQs ;
        delete [] onq ;
        delete [] ofix;
        // --- Wave-function components for each orbital
        delete [] epos ;
        delete [] esize;
        delete [] ecoef;
        // --- Forces acting on wave-functions components
        delete [] efpos ;
        delete [] efsize;
        delete [] efcoef;
        // --- normal forces ( to constrain normality of wavefunction )
        delete [] enfpos ;
        delete [] enfsize;
        delete [] enfcoef;
        // --- Kinetic Energy Forces
        delete [] fTs;
        //delete [] eTfpos ;
        //delete [] eTfsize;
        //delete [] eTfcoef;
        // --- Auxuliary electron density expansion basis functions
        delete [] rhoP;
        delete [] rhoQ;
        delete [] rhoS;
        // --- Forces acting on auxuliary density basis functions
        delete [] rhofP;
        delete [] rhofQ;
        delete [] rhofS;
        //_realloc( rhoEQ, nQtot );
    }

    ~CLCFGO(){ if(bDealoc)dealloc(); }

    void setDefaultValues( ){
        // ToDo : We may automatize this somehow ????
        // atoms
        for(int i=0; i<natom;  i++){
            apos  [i]=Vec3dZero;
            aforce[i]=Vec3dZero;
            //aQs   [i]=1.;
            //aQsize[i]=1.0;
            //aPsize[i]=1.0;
            //aPcoef[i]=0.0;
            aPars [i] = { 1.0, 0.0, 0.1, 0.0 };   /// { x=Q,y=sQ,z=sP,w=cP }
            atype [i]=0;
        }
        for(int i=0; i<nOrb;  i++){
            opos[i]=Vec3dZero;
            odip[i]=Vec3dZero;
            oEs [i]=0;
            oQresc[i] = 1; // ToDo: probably not needed, but makes debugging easier
            oQs [i]=0;
            onq [i]=nqOrb;
            ospin[i]=1;
            ofix [i]=0;
        }
        for(int i=0; i<nBas;  i++){
            // --- Wave-function components for each orbital
            epos [i]=Vec3dZero;
            esize[i]=1;
            ecoef[i]=1;
            // --- Forces acting on wave-functions components
            efpos [i]=Vec3dZero;
            efsize[i]=0;
            efcoef[i]=0;
            // --- normal forces ( to constrain normality of wavefunction )
            enfpos [i]=Vec3dZero;
            enfsize[i]=0;
            enfcoef[i]=0;
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
        bNormForce     = b;
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
        apos[i]=p; aPars[i].set(q,sQ,sP,cP);
        //aQs[i]=q; aQsize[i]=sQ; aPsize[i]=sP; aPcoef[i]=cP;
    }

    inline void setAtom( int i, int itype, Vec3d p ){
        const Quat4d& par = default_AtomParams[itype];
        apos[i]=p; aPars[i] = par;
        //aQs[i]=par.x; aQsize[i]=par.y; aPsize[i]=par.z; aPcoef[i]=par.w;
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


    void cleanEForces(){
        for(int i=0; i<nBas;  i++){ efpos [i] = Vec3dZero; }
        for(int i=0; i<nBas;  i++){ efsize[i] = 0;         }
        for(int i=0; i<nBas;  i++){ efcoef[i] = 0;         }
    }

    inline void cleanForces(){
        for(int i=0; i<nOrb; i++){ oEs[i]=0; };
        for(int i=0; i<natom; i++){ aforce[i] = Vec3dZero; }
        cleanEForces();
        if(bNormForce){
            for(int i=0; i<nBas;  i++){ enfpos [i] = Vec3dZero; }
            for(int i=0; i<nBas;  i++){ enfsize[i] = 0;         }
            for(int i=0; i<nBas;  i++){ enfcoef[i] = 0;         }
        }
        if(bEvalPauli){
            for(int i=0; i<nBas;  i++){ fTs[i].p = Vec3dZero; }
            for(int i=0; i<nBas;  i++){ fTs[i].s = 0;         }
            for(int i=0; i<nBas;  i++){ fTs[i].c = 0;         }
        }
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
// ==================== Normalize Orb - This is to make sure we do nothing wrong by normalization after evaluation of kinetic energy
// ===========================================================================================================================

    void clampSizes(){
        for(int i=0; i<nBas; i++){
            //printf( " esize[%i] %g -> ", i, esize[i] );
            esize[i] = _clamp(esize[i], minSize, maxSize);
            //printf( " %g \n", esize[i] );
        }
    }

    void swapBas( int i, int j){
        _swap( epos [i], epos [j] );
        _swap( esize[i], esize[j] );
        _swap( ecoef[i], ecoef[j] );
    }

    void checkOrder(){
        for(int io=0; io<nOrb; io++){
            int i0=getOrbOffset(io);
            for(int ii=1; ii<perOrb; ii++){
                int i=i0+ii;
                double c0 = ecoef[i-1];
                double c1 = ecoef[i  ];
                if( (c0*c0)<(c1*c1) ){ swapBas(i,i-1); }
            }
        }
    }

    void normalizeOrbSigns(){
        for(int io=0; io<nOrb; io++){
            int i0=getOrbOffset(io);
            double c0 = ecoef[i0];
            if(c0<0){
                for(int ii=0; ii<perOrb; ii++){
                    double& ec = ecoef[i0+ii];
                    ec = -ec;
                }
            }
        }
    }

    double normalizeOrb(int io ){
        int i0    = getOrbOffset(io);
        double Q=0;
        for(int i=i0; i<i0+perOrb; i++){
            Vec3d  pi  = epos [i];
            double ci  = ecoef[i];
            double si  = esize[i];
            Q      += ci*ci; // overlap = 1
            for(int j=i0; j<i; j++){
                Vec3d pj  = epos[j];
                Vec3d Rij = pj-pi;
                double r2 = Rij.norm2();
                if(r2>Rcut2) continue;
                double cj  = ecoef[j];
                double sj  = esize[j];
                //double sij; Vec3d pij;
                //double Sij = Gauss::product3D_s_new( si, pi, sj, pj, sij, pij ); // ToDo : use more efficient/cheaper method here ( without sij,pij etc. )

                _Gauss_sij_aux( si, sj );
                _Gauss_overlap( r2, si, sj );
                double qij = S*ci*cj*2;

                Q += S*(ci*cj)*2;
            }
        }
        double renorm  = sqrt(1./Q);
        //printf( "normalizeOrb(%i) Q %g renorm %g \n", io, Q, renorm );
        //DEBUG_rescale_Q = renorm;
        oQresc[io] = renorm; // ToDo: probably not needed, but makes debugging easier
        double renorm2 = renorm*renorm;
        for(int i=i0; i<i0+perOrb; i++){ ecoef[i] *=renorm;  };
        return Ek;
    }

    void normalizeOrbs(){
        for(int io=0; io<nOrb; io++) normalizeOrb(io);
    }

// ===========================================================================================================================
// ==================== Project Orbitals ( Normalize, Eval Kinetic Energy, Project to Auxiliary Density Basis )
// ===========================================================================================================================

    //double projectOrb(int io, Vec3d* Ps, double* Qs, double* Ss, Vec3d& dip, bool bNormalize ){ // project orbital on axuliary density functions
    double projectOrb(int io, Vec3d& dip ){

        int i0      = getOrbOffset(io);
        int irho0   = getRhoOffset(io);
        Vec3d*   Ps = rhoP+irho0;
        double*  Qs = rhoQ+irho0;
        double*  Ss = rhoS+irho0;

        double   Q  = 0;
        double   Ek = 0; // kinetic energy
        dip         = Vec3dZero;
        Vec3d qcog  = Vec3dZero;
        int ii=0;

        bool bMakeRho = bEvalCoulomb || bEvalAECoulomb || bNormalize;   // in order to normalize we must calculate total charge in orbital

        for(int i=i0; i<i0+perOrb; i++){

            Vec3d  pi  = epos [i];
            double ci  = ecoef[i];
            double si  = esize[i];
            double qii = ci*ci; // overlap = 1
            Q      += qii;
            //if(ii==0) printf( "\n qii(%g|%g)", qii, ci );

            if( bMakeRho ){
                Qs[ii]  = qii;
                Ps[ii]  = pi;
                Ss[ii]  = si*M_SQRT1_2; //  si is multiplied by sqrt(1/2) when put to rho_ii
                //Qs[ii]  = 0; // !!!! DEBUG !!!! - remove diagonal part of density ( ignore fromRhoDiag )
                //if( fabs(Qs[ii])>1e-16 )printf( "orb[%i] rho[%i|%i] q %g s %g p(%g,%g,%g) \n", io, ii,i, Qs[ii], Ss[ii], Ps[ii].x,Ps[ii].y,Ps[ii].z );
                qcog.add_mul( pi, qii );
            }

            double fsi;
            double Tii = Gauss::kinetic_r0_derivs( si, fsi );
            Ek       += qii*   Tii;
            // Copy to Kinetic - we will use them for Pauli
            //eTfsize[i]+= qii*-2*fsi;
            //eTfcoef[i]+= Tii*-2*ci ; // ToDo : This needs to be checked
            //fTs[i].set( Vec3dZero, qii*-2*fsi, Tii*-2*ci  );
            fTs[i].s += qii*-2*fsi;
            fTs[i].c += Tii*-2*ci ;
            if( bEvalKinetic ){
                efsize[i]+= qii*-2*fsi;
                efcoef[i]+= Tii*-2*ci ; // ToDo : This needs to be checked
            }
            if( bNormForce ){
                //  Sii      = ci^2 <Gi|Gj> ( where <Gi|Gj>=1 because they are normalized )
                // =>  d_Sii/d_si = 0    and   d_Sii/d_ci = 2*ci
                enfcoef[i]  += -2*ci;
                //enfcoef[i]  +=  2*ci;
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

                _Gauss_sij_aux( si, sj );
                _Gauss_overlap( r2, si, sj );
                double qij = S*cij;
                Q     += qij;

                if( bNormForce ){
                    Vec3d fij = Rij*(dS_dr*cij);
                    enfpos [i].add( fij )  ; enfpos [j].sub( fij )   ;
                    enfsize[i]-= dS_dsi*cij; enfsize[j]-= dS_dsj*cij;
                    enfcoef[i]-= S*cj*2    ; enfcoef[j]-= S*ci*2   ;
                    //enfpos [i].sub( fij )  ; enfpos [j].add( fij )   ;
                    //enfsize[i]+= dS_dsi*cij; enfsize[j]+= dS_dsj*cij;
                    //enfcoef[i]+= S*cj*2    ; enfcoef[j]+= S*ci*2   ;
                }

                if( bMakeRho ){
                    _Gauss_product(pi,pj,si,sj)
                    Qs[ii] = qij;
                    Ps[ii] = pos_ij;
                    Ss[ii] = size_ij;
                    //Qs[ii] = 0; // !!!!! DEBUG !!!!! - only diagonal
                    //if( fabs(Qs[ii])>1e-16 )printf( "orb[%i] rho[%i|%i,%i] q %g s %g p(%g,%g,%g) \n", io, ii,i,j, Qs[ii], Ss[ii], Ps[ii].x,Ps[ii].y,Ps[ii].z );
                    qcog.add_mul( pos_ij, qij );
                }

                _Gauss_tau    ( r2, si, sj );
                _Gauss_kinetic( r2, si, sj );
                Ek += T*cij;
                Vec3d fij = Rij*(dT_dr*cij);
                fTs[i].add( fij   , -dT_dsi*cij, T*cj*-2 );
                fTs[j].add( fij*-1, -dT_dsj*cij, T*ci*-2 );
                if( bEvalKinetic ){
                    efpos [i].add( fij )  ; efpos [j].sub( fij )  ;
                    efsize[i]-= dT_dsi*cij; efsize[j]-= dT_dsj*cij;
                    efcoef[i]-= T*cj*2    ; efcoef[j]-= T*ci*2    ;
                    // Copy to Kinetic - we will use them for Pauli
                    //eTfpos [i].add( fij )  ; eTfpos [j].sub( fij )  ;
                    //eTfsize[i]-= dT_dsi*cij; eTfsize[j]-= dT_dsj*cij;
                    //eTfcoef[i]-= T*cj*2    ; eTfcoef[j]-= T*ci*2    ;

                }
                ii++;

            }
        }
        if( bNormalize && ( fabs(Q-1)>1e-8 ) ){  printf( "ERROR in CLCFGO::projectOrb(): psi_%i is not normalized |psi|^2 = %g \n", io, Q ); exit(0); }
        //printf( "projectOrb[%i] Q %g Ek %g \n", io, Q, Ek );
        onq [io] = ii;
        opos[io] = qcog;
        oQs [io] = Q;
        if(bEvalKinetic)oEs [io] = Ek;
        if( !bEvalKinetic ) Ek=0;
        return Ek;
    }

// ===========================================================================================================================
// ==================== Project Orbitals ( Normalize, Eval Kinetic Energy, Project to Auxiliary Density Basis )
// ===========================================================================================================================

    //double projectOrb(int io, Vec3d* Ps, double* Qs, double* Ss, Vec3d& dip, bool bNormalize ){ // project orbital on axuliary density functions
    double projectOrb_norm(int io, Vec3d& dip ){

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

            //printf(  " epos[%i,%i] (%g,%g,%g)\n", io, i, pi.x, pi.y, pi.z );

            if(bMakeRho){
                Qs[ii]  = qii;
                Ps[ii]  = pi;
                Ss[ii]  = si*M_SQRT1_2; //  si is multiplied by sqrt(1/2) when put to rho_ii
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
                    //double qij = Sij*cij;
                    double qij = Sij*cij*2; // because qij=qji   (a+b)*(a+b) = a^2 + b^2 + 2*ab
                    //qij*=M_SQRT2*0.96; // This does not help for all distances - at close it produce to high, at distance too low => Sij must decay slower
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
        //printf( "orb[%i] onq %i nQorb %i \n", io, onq[io], nqOrb );
        //printf( "orb[%i] onq %i nQorb %i \n", io, onq[io], nqOrb );
        oQs[io] = Q;
        opos[io].set_mul( qcog, 1./Q );
        if(bNormalize){
            double renorm  = sqrt(1./Q);
            double renorm2 = renorm*renorm;
            //if(perOrb>1)printf( "io[%i] Q %g  s(%g,%g) r %g \n", io, Q, ecoef[0], ecoef[1], (epos[1]-epos[0]).norm() );
            if(bEvalKinetic){ // ToDo: we should  renormalize also coefs

                Ek *= renorm2;
                for(int i=i0; i<i0+perOrb; i++){
                    efpos [i].mul(renorm2);
                    efsize[i]*=renorm2;
                    //efsize[i]*=renorm;
                    efcoef[i]*=renorm;
                }
            }
            //printf( "project orb[$i]: Q %g renorm %g renorm2 %g \n", io, Q, renorm, renorm2 );
            for(int i=i0; i<i0+perOrb; i++){ ecoef[i] *=renorm;  };
            for(int i= 0; i<ii       ; i++){ Qs   [i] *=renorm2; };
        }
        // ToDo: Renormalize also  rhos?
        return Ek;
    }

    double projectOrbs( ){   // project density of all orbitals onto axuliary charge representation ( charges, dipoles and axuliary functions )
        //printf( "   bNormalize %i \n", bNormalize );
        int nqOrb = perOrb*(perOrb+1)/2;
        int i0 =0;
        int ii0=0;
        Ek = 0;
        for(int io=0; io<nOrb; io++){
            //int i0  = getOrbOffset(jo);
            //Ek += projectOrb( io, odip[io] );
            //projectOrb(  io, rhoP+ii0, rhoQ+ii0, rhoS+ii0, odip[io], true );
            //projectOrb(io, ecoefs+i0, erho+irho0, erhoP+irho0, odip[io] );

            if(bNormalize) normalizeOrb( io );
            Ek          += projectOrb  ( io, odip[io] );

            ii0+=nqOrb;
            i0 +=perOrb;
        }
        //printf( "C++ ca,cb %g %g \n", ecoef[0], ecoef[1] );
        return Ek;
    }

    double outProjectNormalForces( int io ){
        //return 0;
        /// out project component of force which break normalization
        // We calculate derivative of normalized energy
        // Ei_    = Ei/Qi
        // dEi/dx = d(Ei/Qi)/dx = (1/Qi)(dEi/dx) - (Ei/Qi^2)(dQi/dx)
        // Assuming functions are already normalized Qi=1
        // dEi/dx = (dEi/dx) - Ei*(dQi/dx)
        int i0   = getOrbOffset(io);
        double c = -oEs[io];
        // ---- out-project f -= ds*c
        for(int i=i0; i<i0+perOrb; i++){
            //if(i==0){ printf( "F: %g -= %g * %g \n", efpos[i].x, enfpos[i].x, c   ); }
            //if(i==0){ printf( "C++ F: %g -= %g * %g \n", efcoef[i], enfcoef[i], c   ); }
            efpos [i] .add_mul( enfpos [i], c );
            efsize[i] +=        enfsize[i]* c  ;
            efcoef[i] +=        enfcoef[i]* c  ;
            //efcoef[i] *= DEBUG_rescale_Q;
            efcoef[i] *= oQresc[io]; // probably not needed but makes debugging easier
        }
        //printf( "c %g \n", c );
        return c;
    }

    double outProjectNormalForces( ){
        double csum=0;
        for(int io=0; io<nOrb; io++){ csum += outProjectNormalForces( io ); }
        return csum;
    }

// ========================================================================================================
// ==================== Pauli Repulsion Models
// ========================================================================================================

    void forceOrb( int io, double K, const Gauss::Blob* Bs ){
        int i0=io*perOrb;
        for(int i=0; i<perOrb; i++){
            int j=i0+i;
            //printf( "forceOrb[%i,%i]  K %g fs %g  \n", io,i, K, efsize[i] );
            Bs[i].applyForceScaled( K, efpos[j], efsize[j], efcoef[j] );
        }
    }

    // ToDo : probably we don't have to do this pairwise, we can sum it and use only Blobs
    void forceOrbPair( int io, int jo, double K, const Gauss::PairInt* Is ){
        int i0=io*perOrb;
        int j0=jo*perOrb;
        int ij=0;
        for(int i=i0; i<i0+perOrb; i++){
            double ci  = ecoef[i];
            for(int j=j0; j<j0+perOrb; j++){
                double cj  = ecoef[j];
                const Gauss::PairInt& I = Is[ij];
                I.applyForceScaled( K, efpos[i], efpos[j], efsize[i], efsize[j] );
                //double KS = -K;
                efcoef[i]+=-K*cj*0.64;
                efcoef[j]+=-K*ci*0.64;
                ij++;
            }
        }
    }

    double pauliAtomOrb( int ia, int io ){
        double Ssum = 0;
        int i0=io*perOrb;
        Gauss::Blob DiS[perOrb];
        for(int i=0; i<perOrb;i++){ DiS[i].setZero(); }
        Vec3d  pj    = apos  [ia];
        //Quat4d& apar = aPars[ia];
        //double cj  = aPars[ia].w;
        double cj=1;
        double sj  = aPars[ia].z;
        Vec3d  fpj = Vec3dZero;
        for(int i_=0; i_<perOrb; i_++){
            int i=i0+i_;
            Vec3d  pi = epos [i];
            double ci = ecoef[i];
            double si = esize[i];
            Vec3d Rij = pj-pi;
            double r2 = Rij.norm2();
            //if(r2>Rcut2) continue;
            //double cij = ci*cj;
            double cij = ci*cj;
            _Gauss_sij_aux(si,sj)
            _Gauss_overlap(r2,si,sj)
            double Scij = S*cij;
            Ssum += Scij;
            Vec3d fp = Rij* dS_dr*cij;
            DiS[i_].add( fp, -dS_dsi*cij, -S*cj );
            fpj    .sub( fp );
        }
        double Amp = KPauliOverlap*aPars[ia].w;
        double E = Ssum*Ssum*Amp;
        double f =    2*Ssum*Amp;
        forceOrb( io, f, DiS );
        aforce[ia].add_mul(fpj, f );
        //oEs[io]+=E*0.5; // why 0.5 ?
        oEs[io]+=E;
        return E;
    }

    double pauliOrbPair( int io, int jo ){
        // This is Pauli potential derived from change of Kinetic Energy due to orbital orthogonalization;
        // Detla_T = Tsym + Tanti - T11 - T12
        // where T11 = <psi_1|T|psi_2>, Tsym=<psi_sym||psi_sym>, Tanti=<psi_anti||psi_anti>
        // psi_sym = (psi_1+psi_2)/|psi_1+psi_1| = (psi_1+psi_2)/2(1+S12)
        // psi_sym = (psi_1+psi_2)/|psi_1+psi_1| = (psi_1+psi_2)/2(1-S12)
        // From this it can be derived
        // Detla_T =  ( T11 + T22 - T12/S12 )/( S12^2/(1-S12^2) )
        //printf( "pauliKineticChnageVB[%i,%i] ", io, jo );
        double Tsum = 0;
        double Ssum = 0;
        int i0=io*perOrb;
        int j0=jo*perOrb;
        int ij=0;

        Gauss::Blob DiS[perOrb];
        Gauss::Blob DjS[perOrb];
        Gauss::Blob DiT[perOrb];
        Gauss::Blob DjT[perOrb];
        for(int i=0; i<perOrb;i++){
            DiS[i].setZero();
            DjS[i].setZero();
            DiT[i].setZero();
            DjT[i].setZero();
        }

        //for(int i=i0; i<i0+perOrb; i++){
        for(int i_=0; i_<perOrb; i_++){
            int i=i0+i_;
            Vec3d  pi  = epos [i];
            double ci  = ecoef[i];
            double si  = esize[i];
            //for(int j=j0; j<j0+perOrb; j++){
            for(int j_=0; j_<perOrb; j_++){
                int j=j0+j_;
                Vec3d pj  = epos[j];
                Vec3d Rij = pj-pi;
                double r2 = Rij.norm2();
                //if(r2>Rcut2) continue;
                double cj  = ecoef[j];
                double sj  = esize[j];

                double cij = ci*cj;
                _Gauss_sij_aux(si,sj)
                _Gauss_overlap(r2,si,sj)
                double Scij = S*cij;
                Ssum += Scij;
                DiS[i_].add( Rij* dS_dr*cij, -dS_dsi*cij, -S*cj );
                DjS[j_].add( Rij*-dS_dr*cij, -dS_dsj*cij, -S*ci );
                //printf( "pauliOrbPair[%i,%i|%i,%i] r %g pi(%g,%g,%g) pj(%g,%g,%g) \n", io,jo, i,j, sqrt(r2),   pi.x,pi.y,pi.z,  pj.x,pj.y,pj.z );
                //printf( "pauliOrbPair[%i,%i|%i,%i] r %g si %g sj %g dSr %g dSsi %g dSsj %g S %g Scij %g Ssum %g \n", io,jo, i,j, sqrt(r2), si, sj, dS_dr, dS_dsi, dS_dsj, S, Scij, Ssum );
                // ToDo : Kinetic and Overlap share much of calculations => make sense to calculate them together in one function
                if(iPauliModel>0){ // Kinetic Energy integral
                    _Gauss_tau    ( r2,si,sj )
                    _Gauss_kinetic( r2,si,sj )
                    double Tcij = T*cij;
                    Tsum -= Tcij;
                    DiT[i_].add( Rij*-dT_dr*cij, dT_dsi*cij, T*cj );
                    DjT[j_].add( Rij* dT_dr*cij, dT_dsj*cij, T*ci );
                }
                ij++;
            }
        }
        double E;
        //printf("iPauliModel %i \n", iPauliModel);
        if(iPauliModel==2){ // Orthogonalization Kinetic energy Valence Bond KE:  Ep = ( Sij^2/(1-Sij^2) )* ( Tii + Tjj - 2*Tij/Sij )
            const double S2SAFE = 0.0001;
            double T11 = oEs[io]; // why 0.5 ?
            double T22 = oEs[jo];
            double S2  = Ssum*Ssum;   // Ssum<1
            double D   = 1/(1-S2  + S2SAFE);
            double D2  = D*D;
            double  fS = Ssum*D;
            double fS2 = Ssum*fS;
            E          = (T11 + T22)*fS2 -2*Tsum*fS;
            double kS  = ( (T11 + T22)*2*Ssum  -2*Tsum*(1+S2) )*D2;
            double kT  = -2*fS;
            forceOrb( io, kS, DiS );
            forceOrb( jo, kS, DjS );
            forceOrb( io, kT, DiT );
            forceOrb( jo, kT, DjT );
            forceOrb( io,fS2, fTs+i0 );  // d_T( Tii*S^2/(1+S^2) )
            forceOrb( jo,fS2, fTs+j0 );  // d_T( Tjj*S^2/(1+S^2) )
            //printf( "pauliOrbPair[%i,%i] E %g S %g T %g D %g | %g %g %g \n", io, jo, E, Ssum, Tsum, D,  (T11 + T22)*fS2,   2*Tsum*fS, fS );
        }else if(iPauliModel==3){ // Juat for debugging
            E = Tsum; // Just Cross-Kinetic T12 = <psi1|T|psi2>
            forceOrb( io, 1, DiT );
            forceOrb( jo, 1, DjT );
        }else if(iPauliModel==4){
            E = Ssum; // Just Cross-Overlap S12 = <psi1|psi2>
            forceOrb( io, 1, DiS );
            forceOrb( jo, 1, DjS );
        }else if(iPauliModel==5){   //   Ep = ( Sij/(1-Sij^2) )* Tij   =  ( Sij^2/(1-Sij^2) )* ( Tij/Sij )
            double S2 = Ssum*Ssum;
            double D  = 1/(1-S2);
            double  fS  = Ssum*D;
            double dfST = Tsum*(1+S2)*D*D;
            E         = fS*Tsum;
            forceOrb( io, dfST, DiS );
            forceOrb( jo, dfST, DjS );
            forceOrb( io, fS,   DiT );
            forceOrb( jo, fS,   DjT );
        }else if(iPauliModel==6){   //   Ep = Sij*Tij
            E = Ssum*Tsum;
            forceOrb( io, Tsum, DiS );
            forceOrb( jo, Tsum, DjS );
            forceOrb( io, Ssum, DiT );
            forceOrb( jo, Ssum, DjT );
        }else{ // Pauli Model 0 :   E = K*S^2
            E = Ssum*Ssum*KPauliOverlap;
            double f = 2*Ssum*KPauliOverlap;
            //printf( "pauliOrbPair(%i,%i) E %g f %g Ssum %g KPauli %g iPauliModel %i\n", io, jo, E, f, Ssum, KPauliOverlap, iPauliModel );
            forceOrb( io, f, DiS );
            forceOrb( jo, f, DjS );
        }
        //if(E>1e+4)printf( "pauliOrbPair[%i,%i] E %g S %g T %g \n", io, jo, E, Ssum, Tsum );
        //printf( "pauliOrbPair[%i,%i] E %g S %g T %g \n", io, jo, E, Ssum, Tsum );
        //return E;
        //return Ssum;
        //oEs[io]+=E*0.5; // why 0.5 ?
        //oEs[jo]+=E*0.5;
        oEs[io]+=E; // why 0.5 ?
        oEs[jo]+=E;
        return E;
    }

    double evalPauli(){ // evaluate Energy components given by direct wave-function overlap ( below cutoff Rcut )
        //printf( "======== evalPauli() \n" );
        EeePaul = 0;
        for(int io=0; io<nOrb; io++){
            for(int jo=0; jo<io; jo++){
                //printf( "evalPauli()[%i,%i] spin %i %i \n", io, jo, ospin[io], ospin[jo] );
                if( ospin[io]==ospin[jo] ) EeePaul += pauliOrbPair(io,jo);
            }
        }
        return EeePaul;
    }

// ===========================================================================================================================
// ==================== Coulomb interaction  -  Using Auxuliary density basis
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

        //Vec3d  pij;
        //double sij;
        //double Cij = Gauss::product3D_s( si, pi, sj, pj, sij, pij );
        //double Sij = Gauss::product3D_s_new( si, pi, sj, pj, sij, pij );

        _Gauss_sij_aux(si,sj)
        _Gauss_overlap(r2,si,sj)
        _Gauss_product(pi,pj,si,sj)

        double cij = ci *cj;
        //printf(  "ci*cj %g ci %g cj %g \n", cij, ci, cj );
        double qij = S*cij*2; // TODO CHECK: should there by realy coefficeint 2.0 ?  .... test by grid !
        //double qij = Sij;
        //double qij = Sij*cij*4.85;
        //double qij = Sij_*cij*4;
        //double qij = Cij*cij*2;
        rhoQ[ij] = qij;
        rhoP[ij] = pos_ij;
        rhoS[ij] = size_ij;
    }

    //Vec3d fromRho( int i, int j, int ij ){
    void fromRhoDiag( int i, int ij ){
        // Q = ci*ci*Sii
        // dE(pa,sa,qa)/dsi = (dE/dpa)*(dpa/dsi) + (dE/dsa)*(dsa/dsi) + (dE/dqa)*(dqa/dsi)
        //  (dqa/dsi)=0         because normalized ???? really even if we have normalization of ?
        //  (dpa/dsi)=0         because in center
        //  (dsa/dsi)=sqrt(0.5) because exp(-r2/2si^2)*exp(-r2/2si^2) = exp(-r2/si^2)
        // dE(pa,sa,qa)/dpi = (dE/dpa)*(dpa/dpi) + (dE/dsa)*(dsa/dpi) + (dE/dqa)*(dqa/dpi)
        //  (dqa/dpi)=0         because in center
        //  (dpa/dpi)=1         because in center
        //  (dsa/dpi)=0         because in center
        /// This function function cares about diagonal density terms rho_ii = wi * wi
        Vec3d  Fpi  = rhofP[ij];
        double Fsi  = rhofS[ij];
        double dEdQ = rhofQ[ij];
        double   ci = ecoef[i];
        efpos [i].add( Fpi );
        efsize[i] += Fsi*M_SQRT1_2;
        efcoef[i] += dEdQ*ci*2;
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
        double dEdQ = rhofQ[ij];

        double cij = 2*ci*cj;

        // BACKUP
        //double fsi = Fp.dot( dxsi ) + Fs*dssi + dEdQ*dSsi*cij;
        //double fsj = Fp.dot( dxsj ) + Fs*dssj + dEdQ*dSsj*cij;

        // Derivation:
        // dE/dpi = (dE/dQij)*(dQij/dpi) + (dE/dpij)*(dpij/dpi)  + (dE/dsij)*(dsij/dpi) // dsij/dpi = 0
        // dE/dpi = (dE/dQij)*(dQij/dpi) + (dE/dpij)*(dpij/dpi)
        // TODO : we should make sure there is a recoil ( forces point to opposite direction )

        Vec3d dSdp = Rij*(dSr*cij);
        //DEBUG_dQdp = dSdp;
        Vec3d Fq   = dSdp*dEdQ;

        efpos [i].add( Fp*dxxi - Fq );
        efpos [j].add( Fp*dxxj + Fq );
        efsize[i] += Fp.dot( dxsi ) + Fs*dssi + dEdQ*dSsi*cij;
        efsize[j] += Fp.dot( dxsj ) + Fs*dssj + dEdQ*dSsj*cij;

        efcoef[i] += dEdQ*cj*Sij*2;
        efcoef[j] += dEdQ*ci*Sij*2;

        // ToDo : What about efcoef[i], efcoef[j] ?

        if(DEBUG_iter==DEBUG_log_iter){
            //if(j==0)printf( "fromRho x %g  F %g =(Fpi %g * dXxi %g)+(dSdp %g * dEdQ %g)\n", pj.x, Fxj.x,  Fp.x, dxxj, dSdp.x, dEdQ );
            //printf( "fromRho x %g  dSdp %g =( dSr %g * cij %g )\n", pi.x, dSdp.x, dSr*Rij.x, cij );
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
                double fr,fs, qij = qi*qj; // NOTE : there should not be factor of 2 (2*qi*qj) because all pairs qi,qj are evaluated independently
                //printf( "CoulombOrbPair[%i|%i,%i|%i] q(%g,%g) s(%g,%g) \n", io,i, jo,j, qi,qj, si,sj );
                //double E = Gauss::Coulomb( Rij, r2, si, sj, qij, fp, fs );

                if( fabs(qij)<1e-16) continue;
                double e = Gauss::Coulomb( sqrt(r2), sqrt(si*si+sj*sj), fr, fs );
                Vec3d fp = Rij*(fr*qij);
                fs*=qij;

                rhofP[i].add(fp);    rhofP[j].sub(fp);
                //rhofS[i] += fs*si*qi;   rhofS[j] += fs*sj*qj;
                rhofS[i] += fs*si;   rhofS[j] += fs*sj;
                //rhofS[i] += fs;   rhofS[j] += fs;
                rhofQ[i] += -e*qj;   rhofQ[j] += -e*qi;
                Ecoul    += e*qij;
                //if( fabs(qij)>1e-16 )printf( "CoulombOrbPair[%i,%i|%i,%i] q %g r %g E %g s(%g,%g) \n",io,jo,i,j, qij, sqrt(r2), e*qij, si, sj );
                //if( fabs(qij)>1e-16 )printf( "CoulombOrbPair[%i,%i|%i,%i] r %g qij(%g,%g) E %g s(%g,%g) \n",io,jo,i,j, sqrt(r2), qi,qj,  E, si, sj );
                //if( fabs(qij)>1e-16 )printf( "CoulombOrbPair[%i,%i|%i,%i] r %g qij(%g,%g) E %g fp %g \n",io,jo,i,j, sqrt(r2), qi,qj,  E, fp.x );
                //if( fabs(qij)>1e-16 )printf( "CoulombOrbPair[%i,%i|%i,%i] r %g | dEdq %g = ( e %g * qj %g ) E %g \n",io,jo,i,j, sqrt(r2), e*qj*-2, e, qj, e*qij );
                //if( fabs(qij)>1e-16 )printf( "CoulombOrbPair[%i,%i|%i,%i] r %g E %g fp %g fr %g qij %g Eij.x %g \n", io,jo,i,j,sqrt(r2),e*qij,fp.x, fr, qij, Rij.x );
                if(DEBUG_iter=DEBUG_log_iter){
                    //printf( "CoulombOrbPair[%i,%i] E %g r %g \n", i-i0,j-j0,E*qij,r,  );

                }
            }
        }
        //printf( " Ecoul[%i,%i] %g \n", io, jo, Ecoul );
        // ToDo : perhaps there should be factor 2 ( i-j and j-i interaction ?)
        //oEs[io]+=Ecoul*0.5; // why 0.5 ?
        //oEs[jo]+=Ecoul*0.5;
        oEs[io]+=Ecoul;
        oEs[jo]+=Ecoul;
        return Ecoul;
    }

    double evalElectrostatICoulomb( ){
        /// Calculate detailed short-range electrostatICoulomb
        //printf("DEBUG evalElectrostatICoulomb \n");
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
                //printf("dEee[%i,%i] %g \n", io,jo, Eee );
                //printf( "evalElectrostatICoulomb[%i,%i] Eee %g r2(%g)<?R2cutOrb2(%g) \n", io, jo, dEee, r2, RcutOrb2 );

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
        //printf("Eee %g \n", Eee );
        return Eee;
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
                Bs[ij].c *= cij;
                S += Bs[ij].c;
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
                //Bs[ij].c *= cij;
                //S += Bs[ij].c;
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
                Bs[ij].c *= cij;
                S += Bs[ij].c;
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
        //        perhaps yes but tr Gauss:Coulomb()he cost of duplicating    Bs and dBs
        double E = 0;
        for(int i=0; i<ni; i++){
            const Gauss::Blob&  Bi =  Bis[i];
                  Gauss::Blob& dBi = dBis[i];
            for(int j=0; j<nj; j++){   // TODO : perhaps we should make special version for Bis==Bjs  with  j<=i
            //int j=1; {
                const Gauss::Blob&  Bj =  Bjs[j];
                      Gauss::Blob& dBj = dBjs[j];
                Vec3d Rij = Bi.p - Bj.p;
                double r2 = Rij.norm2();
                double r  = sqrt(r2 + R2SAFE );
                double s2 = Bi.s*Bi.s + Bj.s*Bj.s;
                double s  = sqrt(s2);
                double qij = Bi.c*Bj.c * 2 ;
                //printf( "//Exchange:Coulomb[%i,%i] r %g s %g \n", i, j, r, s  );
                double fr,fs;
                double e  = Gauss::Coulomb( r, s, fr, fs ); // NOTE : remove s*2 ... hope it is fine ?
                //printf( "Exchange:Coulomb[%i,%i] qij %g E %g fr %g fs %g \n", i, j, r, s, qij, E*qij, fs, fs );
                // --- Derivatives (Forces)
                Vec3d fp  = Rij*(fr * qij);
                fs       *=           qij;
                //if(DEBUG_iter==DEBUG_log_iter) printf( "ExchangeOrbPair:Coulomb[%i,%i] qij %g e %g fx %g fs %g | r %g s %g fr %g fs %g  \n", i,j, qij, e, fp.x, fs,    r, s, fr, fs );
                dBi.p.add(fp);      dBj.p.sub(fp);
                dBi.s -= fs*Bi.s;   dBj.s -= fs*Bj.s;
                dBi.c += e *Bj.c*2; dBj.c += e *Bi.c*2; // ToDo : need to be made more stable ... different (qi,qj)
                E     += e *qij;
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
                double eij=Amp*Bs[ij].c;
                //E   += eij;
                Vec3d Fp = Rij*(dB.dCr*cij);
                efpos [i].add( Fp ); efpos[j].sub( Fp );
                efsize[i]-= dB.dCsi*cij; efsize[j]-= dB.dCsj*cij;
                //efcoef[i]-= eij*cj ;  efcoef[j]-= eij*ci;
                efcoef[i]-= eij/ci ;  efcoef[j]-= eij/cj;
                //if(DEBUG_iter==DEBUG_log_iter) printf( "applyPairForce[%i,%i]%i  dCr %g dCsi %g dCsj %g cij*Amp %g \n", i, j, ij, dB.dCr, dB.dCsi, dB.dCsj, cij );
                //if(DEBUG_iter==DEBUG_log_iter) printf( "applyPairForce[%i,%i]%i Qij %g eij %g dCr %g Fp.x %g cij %g Amp %g \n", i, j, ij, Bs[ij].c, eij, dB.dCr, Fp.x, cij, Amp );
                //printf( "applyPairForce[%i,%i]%i\n", i, j, ij );
                ij++;
            }
        }
        //return E;
    }

// ===========================================================================================================================
// ==================== Eval Cross-Terms ( Contact of two different orbitals, such as Pauli and Exchange interaction )
// ===========================================================================================================================

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
            //printf( "EeePaul[%i,%i] ", io, jo );
            T      = pairOverlapAndKineticDervis( io, jo, Bs, dBs, TDs, S );
            dEpaul = pauliCrossKinetic          ( io, jo, Bs, dBs, TDs, S, T, KRSrho,  anti  );
            //printf( "EeePaul[%i,%i] %g ", io, jo, dEpaul );
        }else{
            if( bEvalPauli && (!anti) ){
                Gauss::iDEBUG=1;
                //printf( "EeePaul[%i,%i] ", io, jo );
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
        //for(int i=0; i<nqOrb; i++){ printf( "fBs[%i] fq %g fs %g fp(%g,%g,%g) \n", i, fBs[i].c, fBs[i].s,    fBs[i].p.x,fBs[i].p.x,fBs[i].p.x ); }
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
    //const double qi = aQs   [ia];
    //const double si = aQsize[ia];
    const double qi = aPars[ia].x;
    const double si = aPars[ia].y;
    //printf( "si %g \n", si );
    //const double si = eAbWs[ia].z;
    int j0  = getRhoOffset(jo);
    int njo = onq[jo];
    double E=0;
    //printf( "evalArho s: " );
    for(int j=j0; j<j0+njo; j++){
    //int j=j0+2;{ // DEBUG
        //if(j!=0) continue;
        Vec3d  pj = rhoP[j];
        Vec3d Rij = pj-pi;
        double r2 = Rij.norm2();
        //if(r2>Rcut2) continue;
        double qj = rhoQ[j];
        double sj = rhoS[j];

        double fr,fs, qij = qi*qj; // NOTE : there should not be factor of 2 (2*qi*qj) because all pairs qi,qj are evaluated independently
        if( fabs(qij)<1e-16) continue;
        double e = Gauss::Coulomb( sqrt(r2), sqrt(si*si+sj*sj), fr, fs );
        //printf( "e %g \n", e   );
        Vec3d fp = Rij*(-fr*qij);
        fs*=qij;

        aforce[ia].sub(fp);
        rhofP[j]  .add(fp);
        rhofS[j] += fs*sj;
        rhofQ[j] += -e*qi;
        E        +=  e*qij;
        //printf( "C++ r %g E %g e %g qij %g qi %g qj %g \n", sqrt(r2), e*qij, e, qij, qi, qj );
        //printf( "evalArho[%i,%i] E %g e,q(%g,%g) r %g s(%g,%g) \n", ia,j, e*qij,e,qij, sqrt(r2), si, sj );
        //if(DEBUG_iter==DEBUG_log_iter) printf( "evalArho[%i,%i] qij %g e %g fp.x %g fr %g fs %g | r %g s %g \n", ia,j, qij, e, fp.x,   fr, fs,    s, r );
    }

    //printf( "\n" );
    //oEs[jo]+=E*0.5; // why 0.5 ?
    oEs[jo]+=E;
    //printf( "evalArho E %g \n", E );
    return E;
}

double evalAE(){
    Gauss::PairDeriv dBs[perOrb];
    double            Ss[perOrb];
    Eae=0,EaePaul=0;
    for(int ia=0; ia<natom; ia++){
        for(int jo=0; jo<nOrb; jo++){
            if(bEvalAEPauli){// --- Pauli Overlap
                EaePaul += pauliAtomOrb( ia, jo );
                //printf( "EaePaul %g bEvalAEPauli %i \n", EaePaul, bEvalAEPauli );
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
        //const double qi   = aQs[i];
        const double qi   = aPars[i].x;
        for(int j=0; j<i; j++){
            Vec3d f = Vec3dZero; // HERE WAS THE ERROR (missing initialization !!!! )
            //Vec3d  abw;
            //double qj = aQs[j];
            double qj = aPars[j].x;
            const Vec3d dR  = apos[j] - pi;
            Eaa +=  addAtomicForceQ( dR, f, qi*qj );
            aforce[j].sub(f);
            aforce[i].add(f);
            //aforce[j].add(f);
            //aforce[i].sub(f);
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
        double Ek = projectOrbs( ); // here we calculate kinetic energy of each orbital and project them to auxuliary charge density basis
        if( bEvalKinetic  ) E+=Ek;
        //printf( "E1 %g \n", E );
        //reportCharges();
        if( bEvalPauli    ) E += evalPauli();
        //printf( "E2 %g \n", E );
        //if( bEvalExchange ) E += evalExchange();
        if( bEvalCoulomb  ) E += evalElectrostatICoulomb(); // repulsion between aux-density clouds => should not distinguish density terms here
        //printf( "E3 %g \n", E );
        if( bEvalAE       ) E += evalAE();
        //printf( "E4 %g \n", E );
        if( bEvalAA       ) E += evalAA();
        //printf( "E5 %g \n", E );
        for(int io=0; io<nOrb; io++){
            assembleOrbForces_fromRho(io);
        }
        //printf( "E6 %g \n", E );
        //E += evalCrossTerms(); // ToDo : This should replace evalPauli and evalExchange()
        if( bNormForce && bNormalize){ outProjectNormalForces(); }
        //printf( "eval() E=%g Ek %g Eee %g EeePaul %g EeeExch %g Eae %g EaePaul %g Eaa %g \n", E, Ek,Eee,EeePaul,EeeExch,Eae,EaePaul,Eaa   );
        //printf( "E7 %g \n", E );
        Etot=E;
        return E;
    }


    double orthogonalizeStep( int niter = 1 ){
        /// Exact Iterative Symetric Orthogonalization
        /// to make absolutely sure the orbitals are orthogonal, we do several step of orthogonalization per each step of
        /// according to https://scicomp.stackexchange.com/a/37683/4696
        ///  phi_i <= phi_i - 0.5 * Sum_j phi_j  <phi_i|phi_j>
        ///  therefore we compute force f = Sum_j phi_j  <phi_i|phi_j>
        /// and then perform gradient descent   phi_i <= phi_i + f*0.5     (dt=0.5)
        //printf( " 0 epos (%g,%g,%g) efpos (%g,%g,%g) \n", epos[0].x,epos[0].y,epos[0].z, efpos[0].x,efpos[0].y,efpos[0].z );
        double K  = KPauliOverlap;
        int iMode = iPauliModel;
        KPauliOverlap = 0.5;
        iPauliModel = 0;
        double F2 = 0;
        for(int i=0; i<niter; i++){
            //cleanForces();
            clampSizes();
            normalizeOrbs();
            cleanEForces();
            //if( bEvalCoulomb || (bEvalAECoulomb && bEvalAE) ) clearAuxDens();
            //double Ek = projectOrbs( );
            //printElectrons();
            evalPauli(); // ToDo : perhaps we should output maximum ramianing overlap ?
            //printElectronForces();
            //if( bNormForce && bNormalize){ outProjectNormalForces(); }
            F2 += moveElectrons( 0.5 );
        }
        //printf( " 3 epos (%g,%g,%g) efpos (%g,%g,%g) \n", epos[0].x,epos[0].y,epos[0].z, efpos[0].x,efpos[0].y,efpos[0].z );
        KPauliOverlap = K;
        iPauliModel   = iMode;
        return F2;
    }

    void forceInfo(){
        double F2a  = 0;
        double F2ep = 0;
        double F2es = 0;
        double F2ec = 0;
        Vec3d  fcog = Vec3dZero;
        if(bOptAtom)for(int i=0; i<natom;i++){ F2a +=aforce[i].norm2(); fcog.add(aforce[i]); }
        if(bOptEPos)for(int i=0; i<nBas; i++){ F2ep+=efpos [i].norm2(); fcog.add(efpos [i]); }
        if(bOptSize)for(int i=0; i<nBas; i++){ F2es+=sq(efsize[i]);     }
        if(bOptCoef)for(int i=0; i<nBas; i++){ F2ec+=sq(efcoef[i]);     }
        // ToDo : We should out project directions which breaks normalization (!!!) but we can do it later - it is mostly importaint for dynamics, gradient desncent should be fine
        //printf( " Force a %g e:p,s,c( %g, %g, %g ) fcog(%g,%g,%g) \n", sqrt(F2a), sqrt(F2ep), sqrt(F2es), sqrt(F2ec), fcog.x,fcog.y,fcog.z );
        //return F2a + F2ep + F2es + F2ec;
    }

    double moveElectrons( double dt ){
        double F2ep = 0;
        double F2es = 0;
        double F2ec = 0; //DEBUG
        //printf( " 1 epos (%g,%g,%g) efpos (%g,%g,%g) \n", epos[0].x,epos[0].y,epos[0].z, efpos[0].x,efpos[0].y,efpos[0].z );
        for(int io=0; io<nOrb; io++){
            int i0=getOrbOffset(io);
            //printf("ofix[%i] %i \n", io, ofix[io] );
            if(ofix[io]>0) continue;
            for(int i=i0; i<i0+perOrb; i++){
                if(bOptEPos){ epos [i].add_mul( efpos [i], dt );    F2ep+=efpos [i].norm2(); }
                if(bOptSize){ esize[i] +=       efsize[i]* dt  ;    F2es+=sq(efsize[i]);     }
                if(bOptCoef){ ecoef[i] +=       efcoef[i]* dt  ;    F2ec+=sq(efcoef[i]);     }
            }
        }
        //printf( " 2 epos (%g,%g,%g) efpos (%g,%g,%g) \n", epos[0].x,epos[0].y,epos[0].z, efpos[0].x,efpos[0].y,efpos[0].z );
        //printf( " F2ep %g F2es %g F2ec %g \n", F2ep, F2es, F2ec );
        return F2ep + F2es + F2ec;
    }

    double moveGD( double dt){
        double F2a  = 0;
        if(bOptAtom)for(int i=0; i<natom;i++){ apos [i].add_mul( aforce[i], dt );    F2a +=aforce[i].norm2(); } //DEBUG
        double F2e = moveElectrons( dt );
        //if(bOptEPos)for(int i=0; i<nBas; i++){ epos [i].add_mul( efpos [i], dt );    F2ep+=efpos [i].norm2(); } //DEBUG
        //if(bOptSize)for(int i=0; i<nBas; i++){ esize[i] +=       efsize[i]* dt  ;    F2es+=sq(efsize[i]);     } //DEBUG
        //if(bOptCoef)for(int i=0; i<nBas; i++){ ecoef[i] +=       efcoef[i]* dt  ;    F2ec+=sq(efcoef[i]);     } //DEBUG
        // ToDo : We should out project directions which breaks normalization (!!!) but we can do it later - it is mostly importaint for dynamics, gradient desncent should be fine
        return F2a + F2e;
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
                //Vec3d pij; double sij;
                double r2 = (apos[ia] - pos).norm2();
                double si  = aPars[ia].z;
                //double cij = aPars[ia].w;
                _Gauss_sij_aux(si,s);
                _Gauss_overlap(r2,si,s);
                double Ssum = S;
                E          += Ssum*Ssum*KPauliOverlap*aPars[ia].w;
                //double f    =    2*Ssum*KPauliOverlap;
                //double S = Gauss::product3D_s_new( aPars[ia].z, apos[ia], s, pos, sij, pij );
                //printf( "atomsPotAtPoint[%i] Pauli S %g c %g pij.x %g sij %g aPsize %g \n", ia, S, aPcoef[ia], pij.x, sij, aPsize[ia] );
                //E += aPars[ia].w*S*S;
            }
            if(bEvalAECoulomb){
                Vec3d Rij = pos - apos[ia];
                // Coulomb( const Vec3d& Rij, double r2, double si, double sj, double qij, Vec3d& fp, double& fs ){
                E += Gauss::Coulomb( Rij, Rij.norm2(), aPars[ia].y, s, 1, fp, fs )*aPars[ia].x*Q;
            }
        }
        return E;
    }


    /*
    double atomsPotAtPoint( const Vec3d& pos, double s, double Q )const{
        //Gauss::PairDeriv dBs[perOrb];
        //double            Ss[perOrb];
        double E =0;
        Vec3d fp; double fs;
        for(int ia=0; ia<natom; ia++){
            if(bEvalAEPauli){
                Vec3d pij; double sij;
                double S = Gauss::product3D_s_new( aPars[ia].z, apos[ia], s, pos, sij, pij );
                //printf( "atomsPotAtPoint[%i] Pauli S %g c %g pij.x %g sij %g aPsize %g \n", ia, S, aPcoef[ia], pij.x, sij, aPsize[ia] );
                E += aPars[ia].w*S*S;
            }
            if(bEvalAECoulomb){
                Vec3d Rij = pos - apos[ia];
                // Coulomb( const Vec3d& Rij, double r2, double si, double sj, double qij, Vec3d& fp, double& fs ){
                E += Gauss::Coulomb( Rij, Rij.norm2(), aPars[ia].y, s, 1, fp, fs )*aPars[ia].x*Q;
            }
        }
        return E;
    }
    */
    double* atomsPotAtPoints( int n, Vec3d* ps, double* out=0, double s=0.0, double Q=1.0 )const{
        if(out==0){ out = new double[n]; };
        for(int i=0; i<n; i++){
            out[i] = atomsPotAtPoint( ps[i], s, Q );
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
        }
        return out;
    }

    double rhoAtPoint( int io, const Vec3d& pos )const{
        //int i0 = getOrbOffset( io );
        int i0   = getRhoOffset(io);
        int ni   = onq[io];
        double rho_sum=0;
        for(int i=0; i<ni; i++){
            int i_=i0+i;
            Vec3d dR  = pos - rhoP[i_];
            double r2 = dR.norm2();
            rho_sum += Gauss::rho3D_r2( r2, rhoS[i_] ) * rhoQ[i_];
        }
        return rho_sum;
    }
    double* rhoAtPoints( int io, int n, Vec3d* ps, double* out=0 )const{
        if(out==0){ out = new double[n]; };
        for(int i=0; i<n; i++){
            out[i] = rhoAtPoint( io, ps[i] );
        }
        return out;
    }

    double hartreeAtPoint( int io, const Vec3d& pos )const{
        //int i0   = getOrbOffset( io );
        int i0     = getRhoOffset(io);
        int perOrb = onq[io];
        double v_sum=0;
        for(int i=0; i<perOrb; i++){
            int i_=i0+i;
            //Vec3d dR  = pos - epos[i_];
            Vec3d dR  = pos - rhoP[i_];
            double r2 = dR.norm2() + 1e-16;
            double r = sqrt(r2);
            //v_sum +=  (ecoef[i_] * const_El_eVA) * erf_6_plus( r/( esize[i_]        ) )/r; // This is using spread/size of wavefunction blob
            //v_sum   +=  (rhoQ [i_] * const_El_eVA) * erf_6_plus( r/( rhoS [i_]*M_SQRT2) )/r; // This is using spread/size of density      blob
            double fr,fs;
            v_sum += rhoQ [i_] * Gauss::Coulomb( r, rhoS[i_], fr, fs );
            //v_sum += rhoQ [i_] * Gauss::Coulomb( r, rhoS[i_]*M_SQRT2, fr, fs );  // this factor M_SQRT2 should not be necessary if rhoS[i_] is width of charge density blobs (not wave functions)
        }
        return v_sum;
    }
    double* hartreeAtPoints( int io, int n, Vec3d* ps, double* out=0 )const{
        if(out==0){ out = new double[n]; };
        for(int i=0; i<n; i++){
            out[i] = hartreeAtPoint( io, ps[i] );
        }
        return out;
    }

    double orb2grid( int io, const GridShape& gridShape, double* buff )const{
        return evalOnGrid( gridShape, [&](int ig, const Vec3d& pos, double& res){ buff[ig] = orbAtPoint( io, pos ); });
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
        return evalOnGrid( gridShape, [&](int ig, const Vec3d& pos, double res){ buff[ig] = rhoAtPoint( io, pos ); });
    }

    double hartree2grid( int io, const GridShape& gridShape, double* buff )const{
        return evalOnGrid( gridShape, [&](int ig, const Vec3d& pos, double res){ buff[ig] = hartreeAtPoint( io, pos ); });
    }


bool loadFromFile( char const* filename ){
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
    sscanf (line, "%i %i %i %b\n", &natom_, &nOrb_, &perOrb_, &bClosedShell );
    //printf( ">>%s<< | na %i no %i perOrb %i  bClosedShell %i \n", line, natom_, nOrb_, perOrb_, bClosedShell );
    //printf("na %i ne %i perORb %i \n", natom, nOrb, perOrb_);
    //printf("na %i ne %i perORb %i \n", natom_, nOrb_, perOrb_ );
    if(bClosedShell) nOrb_*=2;
    realloc( natom_, nOrb_, perOrb_, 1 );
    //printf("na %i ne %i perOrb %i | na_ %i ne_ %i perOrb_ %i \n", natom, nOrb, perOrb,  natom_, nOrb_, perOrb_);
    setDefaultValues( );
    double Qasum = 0.0;
    for(int i=0; i<natom; i++){
        double x,y,z;
        double Q,sQ,sP,cP;
        fgets( buff, nbuff, pFile); //printf( "fgets: >%s<\n", buf );
        int nw = sscanf (buff, "%lf %lf %lf %lf %lf %lf %lf", &x, &y, &z,    &Q, &sQ, &sP, &cP );
        //printf( "atom[%i] p(%g,%g,%g) Q %g sQ %g sP %g cP %g \n", i, x, y, z,    Q, sQ, sP, cP );
        apos  [i]=(Vec3d){x,y,z};
        //aQs   [i]=Q;
        //aQsize[i]=sQ;
        //aPsize[i]=sP;
        //aPcoef[i]=cP;
        aPars[i].set(Q,sQ,sP,cP);
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
        //printf( "bas[%i|%i,%i] p(%g,%g,%g) s %g c %g spin %i \n", i, i/perOrb, i%perOrb, x, y, z,    s, c, spin );
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

bool saveToFile( char const* filename ){
    //printf(" filename: >>%s<< \n", filename );
    FILE * pFile;
    pFile = fopen (filename,"w");
    if( pFile == NULL ){
        printf("ERROR in CLCFGO::saveToFile(%s) : No such file !!! \n", filename );
        return -1;
    }
    fprintf( pFile, "%i %i %i %i\n", natom, nOrb, perOrb, 0 );
    for(int i=0; i<natom; i++){
        fprintf( pFile, "%3.6f %3.6lf %3.6lf   %3.6lf %3.6f %3.6f %3.6f \n", apos[i].x, apos[i].y, apos[i].z,    aPars[i].x, aPars[i].y, aPars[i].z, aPars[i].w );
    }
    int i=0;
    for(int io=0; io<nOrb;   io++){
        for(int ii=0; ii<perOrb; ii++){
            fprintf( pFile, "%3.6f %3.6f %3.6f  %2.6f  %2.6f  %i", epos[i].x, epos[i].y, epos[i].z,  esize[i], ecoef[i], ospin[io] );
            if(ii==0){ fprintf( pFile, " // Orb #%i", io ); }
            fprintf( pFile, "\n" );
            i++;
        }
    }
    fclose (pFile);
    return 0;
}


void printSetup(){
    printf("===CLCFGO Setup :\n");
    printf("const_K_eVA %g const_El_eVA %g \n", const_K_eVA, const_El_eVA );
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

void printEnergies(){
    printf( "Etot %g | Ek %g Eee,p(%g,%g) Eae,p(%g,%g) Eaa %g \n", Etot, Ek, Eee,EeePaul, Eae,EaePaul, Eaa );
}

void printAtoms(){
    printf( "===CLCFGO::printAtoms()\n");
    for(int i=0; i<natom; i++){
        //printf( "eFF::atom[%i] p(%g,%g,%g)[A] Q(%g[e],%g[A])  Pauli(%g[eV],%g[A]) \n", i, apos[i].x,apos[i].y,apos[i].z,  aQs[i], aQsize[i], aPcoef[i], aPsize[i] );
        printf( "eFF::atom[%i] p(%g,%g,%g)[A] Q(%g[e],%g[A])  Pauli(%g[eV],%g[A]) \n", i, apos[i].x,apos[i].y,apos[i].z,  aPars[i].x, aPars[i].y, aPars[i].w, aPars[i].z );
    }
}

void printElectrons(){
    printf( "===CLCFGO::printElectrons()\n");
    for(int io=0; io<nOrb; io++){
        //printf( ">> orb[%i] spin %i \n", io, ospin[io] );
        for(int j=0; j<perOrb; j++){
            int ie=io*perOrb+j;
            printf( "e[%i,%i|%i] p(%g,%g,%g)[A] size %g coef %g \n", io,j,ie, epos[ie].x,epos[ie].y,epos[ie].z, esize[ie], ecoef[ie] );
        }
    }
}

void printElectronForces(){
    printf( "===CLCFGO::printElectronForces()\n");
    for(int io=0; io<nOrb; io++){
        for(int j=0; j<perOrb; j++){
            int ie=io*perOrb+j;
            printf( "F e[%i,%i|%i] p(%g,%g,%g)[A] size %g coef %g \n", io,j,ie, efpos[ie].x,efpos[ie].y,efpos[ie].z, efsize[ie], efcoef[ie] );
        }
    }
}

int bas2str(char* str, int ie){
    return sprintf( str, "c %5.2f s %5.2f p(%5.2f,%5.2f,%5.2f)\n",  ecoef[ie], esize[ie], epos[ie].x,epos[ie].y,epos[ie].z );
}

int Eterms2str(char* str){
    // Ek=0,Eee=0,EeePaul=0,EeeExch=0,Eae=0,EaePaul=0,Eaa=0, Etot=0;
    return sprintf( str, "Etot %3.3f Ek %3.3f Eee,P(%3.3f,%3.3f) Eae,P(%3.3f,%3.3f) Eaa %g )\n", Etot, Ek, Eee, EeePaul, Eae, EaePaul, Eaa );
}

char* orb2str(char* str0, int io){
    char* str=str0;
    int i0 = getOrbOffset(io);
    for(int ii=0; ii<perOrb; ii++){
        int i=i0+ii;
        str += bas2str(str, i);
    }
    return str;
}

char* orbs2str(char* str0){
    char* str=str0;
    for(int i=0; i<nOrb; i++){
        str+=sprintf(str,"\norb_%i E %g \n", i, oEs[i] );
        str=orb2str(str, i);
    }
    return str;
}

};

/// @}

#endif
