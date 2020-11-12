
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

    bool bEvalKinetic  = true;
    bool bEvalCoulomb  = true;
    bool bEvalPauli    = true;
    bool bEvalExchange = true;

    double Rcut    =6.0;  ///< cutoff beyond which two basis functions chi has no overlap
    double RcutOrb =9.0;  ///< cutoff beoud which orbital (i.e. localized cluster of basis functions) has no overlap
    double Rcut2     =Rcut*Rcut;
    double RcutOrb2  =RcutOrb*RcutOrb;


    int natypes =0;

    int natom =0; ///< number of atoms (nuclei, ions)
    int perOrb=0; //!< Brief number of spherical functions per orbital
    int nOrb  =0; //!< Brief number of single-electron orbitals in system
    // this is evaluated automaticaly
    int nBas  =0; ///< number of basis functions
    int nqOrb =0; ///< number of charges (axuliary density elements) per orbital
    int nQtot =0; ///< total number of charge elements

    // atoms (ions)
    Vec3d*  apos   =0;  ///< positioon of atoms
    Vec3d*  aforce =0;  ///< positioon of atoms
    double* aQs    =0;  ///< charge of atom
    int*    atype  =0;  ///< type of atom (in particular IKinetic pseudo-potential)

    // orbitals
    Vec3d*  opos =0;   ///< store positions for the whole orbital
    Vec3d*  odip =0;   ///< Axuliary array to store dipoles for the whole orbital
    double* oEs  =0;   ///< orbital energies
    double* oQs  =0;   ///< total charge in orbital before renormalization (just DEBUG?)
    int*    onq  =0;   ///< number of axuliary density functions per orbital

    // --- Wave-function components for each orbital
    Vec3d*  epos  =0; ///< position of spherical function for expansion of orbitals
    double* esize =0;
    double* ecoef =0;  ///< expansion coefficient of expansion of given orbital
    // --- Forces acting on wave-functions components
    Vec3d*  efpos  =0; ///<   force acting on position of orbitals
    double* efsize =0; ///<   force acting on combination coefficnet of orbitals
    double* efcoef =0; ///<  force acting on size of gaussians

    // --- Auxuliary electron density expansion basis functions
    Vec3d * rhoP  =0; ///< position of density axuliary functio
    double* rhoQ  =0; ///< temporary array to store density projection on pair overlap functions
    double* rhoS  =0;
    // --- Forces acting on auxuliary density basis functions
    Vec3d * rhofP =0; ///< position of density axuliary functio
    double* rhofQ =0; ///< temporary array to store density projection on pair overlap functions
    double* rhofS =0;
    //double* rhoEQ =0; /// coulomb energy

    // ======= Functions

    void realloc( int natom_, int nOrb_, int perOrb_, int natypes_ ){
        // atoms
        if( natom != natom_ ){
            natom = natom_;
            _realloc( apos  ,natom );
            _realloc( aforce,natom );
            _realloc( aQs   ,natom );
            _realloc( atype ,natom ); // not used  now
        }
        if( (nOrb != nOrb_)||(perOrb != perOrb_) ){
            nOrb   = nOrb_;
            perOrb = perOrb_;
            nBas   = nOrb * perOrb;
            nqOrb  = perOrb*(perOrb+1)/2;
            nQtot  = nOrb*nqOrb;

            // orbitals
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
        // atoms
        for(int i=0; i<natom;  i++){
            apos  [i]=Vec3dZero; 
            aforce[i]=Vec3dZero;
            aQs   [i]=1.; 
            atype [i]=0;
        }
        for(int i=0; i<nOrb;  i++){
            opos[i]=Vec3dZero;
            odip[i]=Vec3dZero;
            oEs [i]=0;
            oQs [i]=0;
            onq [i]=0;
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

    inline int getOrbOffset(int iorb)const{ return iorb*nOrb;  }
    inline int getRhoOffset(int iorb)const{ return iorb*nqOrb; }

    inline void setRcut( double Rcut_ ){
        Rcut   = Rcut_;
    }

    inline bool checkOrbitalCutoff(int i, int j)const{
        double r2 = (opos[i]-opos[j]).norm2();
        return r2<RcutOrb2;
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

    double orbNorm(int io ){
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
                Vec3d Fp = Rij*(-dSr*cij);
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

    void reportCharges(int io){
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


    //double projectOrb(int io, Vec3d* Ps, double* Qs, double* Ss, Vec3d& dip, bool bNormalize ){ // project orbital on axuliary density functions
    double projectOrb(int io, Vec3d& dip, bool bNormalize ){
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
        Vec3d oqcog = opos[io];
        int ii=0;

        for(int i=i0; i<i0+perOrb; i++){

            Vec3d  pi  = epos [i];
            double ci  = ecoef[i];
            double si  = esize[i];
            double qii = ci*ci; // overlap = 1
            
            if(bEvalCoulomb){
                qcog.add_mul( pi, qii );
                /// ToDo: MUST USE PRODUCT OF GAUSSIANS !!!!   gaussProduct3D( double wi, const Vec3d& pi, double wj, const Vec3d& pj,  double& wij, Vec3d& pij ){
                Q      += qii;
                Qs[ii]  = qii;
                Ps[ii]  = pi;
                Ss[ii]  = si*M_SQRT1_2;
                //Qs[ii]  = 0; // DEBUG
            }

            double fr,fsi,fsj;
            if(bEvalKinetic){
                //double fEki;
                //Ek += qii*addKineticGauss( si*M_SQRT2, fEki );
                //Ek += qii*Gauss::kinetic( si );
                Ek += qii*Gauss:: kinetic_s(  0.0, si, si,   fr, fsi, fsj );
                efsize[i]+= 2*fsi*qii;
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

                if(bEvalCoulomb){
                    // --- Project on auxuliary density functions
                    Vec3d  pij;
                    double sij;
                    //double Cij = Gauss::product3D_s( si, pi, sj, pj, sij, pij );
                    double Sij = Gauss::product3D_s_new( si, pi, sj, pj, sij, pij );
                    double qij = Sij*cij;
                    qcog.add_mul( pij, qij );
                    Q         += qij;

                    Qs[ii] = qij;
                    Ps[ii] = pij;   // center of axuliary overlap density function in the middle between the two wavefunctions
                    Ss[ii] = sij;
                }
               
                if(bEvalKinetic){
                    //double Ekij = Gauss::kinetic(  r2, si, sj ) * 2; // TODO : <i|Lapalace|j> between the two gaussians
                    double Kij = Gauss:: kinetic_s(  r2, si, sj,   fr, fsi, fsj );   // fr*=2; fsi*=2, fsj*=2;
                    Ek += Kij*cij;
                    Vec3d fij = Rij*(fr*cij);
                    efpos [i].add( fij ); efpos[j].sub( fij );
                    efsize[i]+= fsi*cij ; efsize[j]+= fsj*cij;
                    efcoef[i]+= Kij*cj ;  efcoef[j]+= Kij*ci;
                    if(DEBUG_iter==DEBUG_log_iter){
                        printf(" Kij %g cij %g Ekij \n", Kij, cij, Kij*cij );
                    }
                }

                ii++;
                // ToDo : Store qij to list of axuliary functions
            }
        }
        onq[io] = ii;
        oQs[io] = Q;
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
        double Ek = 0;
        for(int io=0; io<nOrb; io++){
            //int i0  = getOrbOffset(jo);
            //oQs[io] =
            Ek += projectOrb( io, odip[io], bNormalize );
            //projectOrb(  io, rhoP+ii0, rhoQ+ii0, rhoS+ii0, odip[io], true );
            //projectOrb(io, ecoefs+i0, erho+irho0, erhoP+irho0, odip[io] );
            ii0+=nqOrb;
            i0 +=perOrb;
        }
        return Ek;
    }

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
        efsize[i] += Fsi*Q;
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

        double qij = qi*qj;

        Vec3d Rij = pj-pi;
        double r2 = Rij.norm2();

        double r    = sqrt(r2 + R2SAFE);
        double s2   = si*si + sj*sj;
        double s    = sqrt(s2);

        double fr,fs;
        //double Eqq  = CoulombGauss( r, s*2, fr, fs, qij );

        //double E  = Gauss::Coulomb( r, s*2, fr, fs );   // Q :  Should there be the constant s*2 ????
        double E  = Gauss::Coulomb( r, s, fr, fs );       // WARRNING  :  removed the contant s*2 to s  ... is it correct ?

        //printf( "CoublombElement[%i,%i] q(%g,%g) E %g fs %g fr %g s %g r %g \n", i,j, qi,qj, E, fs, fr, s, r );

        Vec3d fp = Rij*(-fr*qij);
        fs *= qij * M_SQRT1_2;
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
                //if(r2>Rcut2) continue;

                double qj = rhoQ[j];
                double sj = rhoS[j];

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
                Vec3d fp = Rij*(-fr * qij);
                fs       *=           qij;
                rhofP[i].add(fp);    rhofP[j].sub(fp);
                rhofS[i] -= fs*si;   rhofS[j] -= fs*sj;
                rhofQ[i] += E*qj;    rhofQ[j] += E*qi; // ToDo : need to be made more stable ... different (qi,qj)
                Ecoul    += E*qij;

                if(DEBUG_iter=DEBUG_log_iter){
                    //printf( "CoublombElement[%i,%i] q(%g,%g) E %g fs %g fr %g s %g r %g \n", i,j, qi,qj, E, fs, fr, s, r );
                    //printf( "CoublombElement[%i,%i] q(%g,%g) E %g fs %g fr %g s %g r %g \n", i,j, qi,qj, E*qij, fp.x, fs*si, s, r );
                    //printf( "CoublombElement[%i,%i] q(%g,%g) E %g fs %g fr %g s %g r %g \n", i,j, qi,qj, E, fs, fr, s, r );
                    //printf( "CoulombOrbPair[%i,%i][%i,%i] e %g E %g s %g(%g,%g) q %g(%g,%g) r %g fr %g \n", io,jo, i,j,  E, E*qi*qj, s,si,sj, qij,qi,qj, r, fr );
                    //printf( "CoulombOrbPair[%i,%i] E %g r %g \n", i-i0,j-j0,E*qij,r,  );
                }
            }
        }
        
        //printf( " Ecoul %g \n", Ecoul );
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
        double E = 0;
        for(int io=0; io<nOrb; io++){
            int i0 = getRhoOffset(io);
            Vec3d opi  = opos[io];
            Vec3d dipi = odip[io];
            for(int jo=io+1; jo<nOrb; jo++){
                Vec3d dop = opos[jo] - opi;
                double r2 = dop.norm2();
                //printf( "evalElectrostatICoulomb[%i,%i]  %g <? %g \n", io, jo, r2,RcutOrb2 );

                E += CoulombOrbPair( io, jo );

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
        for(int i=0; i<nOrb; i++){ 
            assembleOrbForces_fromRho(i); 
        }
        return E;
    }

    double ExchangeOrbPair(int io, int jo){
        // It may be actually easier to evaluate Hartree-Fock like exchange, because it is linear
        // see:   https://en.wikipedia.org/wiki/Exchange_interaction#Exchange_of_spatial_coordinates
        // https://chemistry.stackexchange.com/questions/61176/what-is-the-exchange-interaction
        // => Exchange integral should be rather easy to compute, just by calculating overlap between the two functions   K = <rho_ij(r1)| 1/|r1-r2| |rho_ij(r2)>
        // scheme:
        //  1) for pair of orbital, construct auxiliary rho_ij (by expansion on auxuliary functions)
        //  2) callculate self-convolution of this auxuliary function

        //if( ~(0==io) ) return 0;
        if(DEBUG_iter==DEBUG_log_iter)printf( " -- ExchangeOrbPair[%i,%i] \n", io, jo );

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

                //if(DEBUG_iter==DEBUG_log_iter) printf( "Exchange:Project[%i,%i] cij %g Sij %g Qij %g r %g x(%g,%g) \n", i, j, cij, pairs[ij].C, cij*pairs[ij].C, sqrt(r2), pi.x, pj.x );
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
        //for(int i=0; i<nqOrb; i++){
        int i=0;{
            //Vec3d  pi = posij[i];
            //double qi = rhoij[i];
            //double si = szij [i];
            const Gauss::PairDerivs& A = pairs[i];
            double ci = cijs[i];
            double qi = ci*A.C;
            //for(int j=0; j<=i; j++){
            int j=1; {
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
                fPs[i].add(fp);    fPs[j].sub(fp);
                fSs[i] -= fs*A.s;  fSs[j] -= fs*B.s;
                fQs[i] += e*qj;    fQs[j] += e*qi; // ToDo : need to be made more stable ... different (qi,qj)
                E      += e*qij;
            }
        }
        
        // DEBUG
        if(DEBUG_iter==DEBUG_log_iter){
            for(int i=0; i<nqOrb; i++){
                printf( "[%i] e %g fx %g | x %g q %g \n", i,  fQs[i], fPs[i].x,     pairs[i].p.x, pairs[i].C*cijs[i] );
            };
        }
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
                cijs[ij] = cij;
                Gauss::PairDerivs& pair = pairs[ij];
                if(DEBUG_iter==DEBUG_log_iter) printf( "Exchange:fromRho[%i,%i] x(%g,%g) ", i, j, pi.x, pj.x );
                pair.backForce( Rij, cij, fQs[ij], fPs[ij], fSs[ij],  efpos[i], efpos[j], efsize[i], efsize[j] );
                ij++;
            }
        }
        //printf( " Exchange:FromRho DONE \n" );
        if(DEBUG_iter==DEBUG_log_iter)printf( "ExchangeOrbPair() E=%g \n", E  );
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

    void evalExchangeCorrelation(){
        // not sure how exactly to do that now - no using real space grid since it is slow
    }

    double eval(){
        double E=0;
        cleanForces();
        if( bEvalCoulomb  ) clearAuxDens();
        //projectOrbs( true );
        double Ek = projectOrbs( false );
        if( bEvalKinetic  ) E+=Ek;
        //reportCharges();
        if( bEvalPauli    ) E += evalPauli();
        if( bEvalExchange ) E += evalExchange();
        if( bEvalCoulomb  ) E += evalElectrostatICoulomb(); // repulsion between aux-density clouds => should not distinguish density terms here
        //printf( "eval() E=%g \n", E  );
        return E;
    }

    // ========== On Grid

    /// ToDo : Coulomb integral C{i,j}  on grid can be calculated like this:
    //  for i in MOs :
    //      wf_j = orb2grid( i, grid )  // potencial or density
    //      for j in MOs:
    //          Iij = 0
    //          for kbas in basis_decomposition[j]:
    //              for pos in grid:
    //                  Iij += kbas(pos)^2 * wf_j[pos]^2

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

    void orb2xsf( const GridShape& grid, int iorb, const char* fname ){
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

};

/// @}

#endif
