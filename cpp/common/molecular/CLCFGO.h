
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


    Vec3d DEBUG_dQdp;

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
    double* rhoEQ =0; /// coulomb energy

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
            _realloc( rhoEQ, nQtot );

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
            rhoEQ[i]=0;
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
        for(int i=0; i<nQtot; i++){ rhoEQ[i] = 0; };
    }

    inline void cleanForces(){
        for(int i=0; i<natom; i++){
            aforce[i] = Vec3dZero;
        }
        for(int i=0; i<nBas; i++){
            efpos [i] = Vec3dZero;
            efsize[i] = 0;
            efcoef[i] = 0;
        }

        //for(int i=0; i<nQtot; i++){
        //    rhofP[i] = Vec3dZero;
        //    rhofQ[i] = 0;
        //    rhofS[i] = 0;
        //}
        clearAuxDens();
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

    void evalPauli(){ // evaluate Energy components given by direct wave-function overlap ( below cutoff Rcut )
        double E = 0;
        for(int io=0; io<nOrb; io++){
            for(int jo=0; jo<nOrb; jo++){
                if( !checkOrbitalCutoff(io,jo) ) continue; // optimization, not to calculate needless interactions
                //interpolate( epos[j] - epos[i], Is[ityp][jtyp], f ); // ToDo: different basis function types
                E += evalPauli( io, jo );
            }
        }
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
        double DT=0; // kinetic energy change by orthognalization
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
            qcog.add_mul( pi, qii );
            /// ToDo: MUST USE PRODUCT OF GAUSSIANS !!!!   gaussProduct3D( double wi, const Vec3d& pi, double wj, const Vec3d& pj,  double& wij, Vec3d& pij ){
            Q      += qii;

            // DEBUG : testing Gaussian Transform
            Qs[ii]  = qii;
            Ps[ii]  = pi;
            Ss[ii]  = si*M_SQRT1_2;
            //Qs[ii]  = 0; // DEBUG

            //double fEki;
            //Ek += qii*addKineticGauss( si*M_SQRT2, fEki );
            //Ek += qii*Gauss::kinetic( si );
            double fr,fsi,fsj;
            Ek += qii*Gauss:: kinetic_s(  0.0, si, si,   fr, fsi, fsj );
            efsize[i]+= 2*fsi*qii;

            printf( "projectOrb[%i] [%i] q[%i] %g \n", io, i, ii, qii );

            //if(qii>1e-8)printf( "Kinetic[%i] %g | %g %g %g \n", i, fsi*ci*ci, fsi, ci, qii );

            //printf( "orb[%i|%i   ] s(%g):%g qii %g \n", io, i,  si, Ss[ii],  qii );
            ii++;
            //printf( "orb[%i|%i] s %g qii %g \n", io, i,  si*2,  qii );
            for(int j=i0; j<i; j++){
                Vec3d pj  = epos[j];
                Vec3d Rij = pj-pi;
                double r2 = Rij.norm2();
                if(r2>Rcut2) continue;

                double cj  = ecoef[j];
                double sj  = esize[j];
                // --- Evaluate Normalization, Kinetic & Pauli Energy
                double dSr, dSsi, dSsj;
                double dTr, dTsi, dTsj;
                const double resz = M_SQRT2; // TODO : PROBLEM !!!!!!!   getOverlapSGauss and getDeltaTGauss are made for density gaussians not wave-function gaussians => we need to rescale "sigma" (size)
                //double Sij  = getOverlapSGauss( r2, si*resz, sj*resz, dSr, dSsi, dSsj );
                double DTij = getDeltaTGauss  ( r2, si*resz, sj*resz, dTr, dTsi, dTsj ); // This is not normal kinetic energy this is change of kinetic energy due to orthogonalization

                //double Ekij = Gauss::kinetic(  r2, si, sj ) * 2; // TODO : <i|Lapalace|j> between the two gaussians
                double Ekij = Gauss:: kinetic_s(  r2, si, sj,   fr, fsi, fsj )*2; fr*=2; fsi*=2, fsj*=2;

                ///ToDo :   <i|Laplace|j> = Integral{ w1*(x^2 + y^2)*exp(w1*(x^2+y^2)) *exp(w2*((x+x0)^2+y^2)) }
                /// ToDo :  Need derivatives of Kinetic Overlap !!!!!
                // --- Project on auxuliary density functions
                Vec3d  pij;
                double sij;
                //double Cij = Gauss::product3D_s( si, pi, sj, pj, sij, pij );
                double Sij = Gauss::product3D_s_new( si, pi, sj, pj, sij, pij );
                double cij = ci *cj;
                //double qij = Sij*cij*2; // factor 2  because  Integral{(ci*fi + cj*fj)^2} = (ci^2)*<fi|fi> + (cj^2)*<fj|fj> + 2*ci*cj*<fi|fj>
                double qij = Sij*cij*2;


                //qij = 0; // DEBUG remove off-diagonal charge

                //printf( "DEBUG projectOrb[%i|%i,%i] r2 %g Tij %g \n", io, i,j, r2, DTij );
                //printf( "DEBUG projectOrb[%i|%i,%i] sij %g \n", io, i,j, sij );
                //printf( "orb[%i|%i,%i] r %g s(%g,%g):%g qS(%g,%g|%g):%g C %g \n", io, i,j, sqrt(r2),  si,sj,sij,  ci,cj,Sij,qij, Cij );
                qcog.add_mul( pij, qij );
                Q  +=   qij;
                DT += DTij*cij;
                Ek += Ekij*cij;

                /*
                // !!!!!!!!!!!!!!! DEBUG !!!!!!!!!!!!!!!!!!
                // !!!  Commented kinetc energy forces !!!!
                // !!!!!!!!!!!!!!! DEBUG !!!!!!!!!!!!!!!!!!
                // --- Derivatives ( i.e. Forces )
                Vec3d fij = Rij*(fr*cij);
                efpos [i].add( fij ); efpos[j].sub( fij );
                efsize[i]+= fsi*cij ; efsize[j]+= fsj*cij;
                efcoef[i]+= Ekij*cj ; efcoef[j]+= Ekij*ci;
                */

                //printf(  "Kinetic [%i,%i]  fsi,j %g %g  \n ", i,j, fsi*cij, fsj*cij  );

                // ToDo: MUST USE PRODUCT OF GAUSSIANS !!!!   gaussProduct3D( double wi, const Vec3d& pi, double wj, const Vec3d& pj,  double& wij, Vec3d& pij ){

                printf( "projectOrb[%i] [%i,%i] r %g q[%i] %g \n", io, i,j, r2, ii, qij );

                Qs[ii] = qij;
                Ps[ii] = pij;   // center of axuliary overlap density function in the middle between the two wavefunctions
                Ss[ii] = sij;
                ii++;
                // ToDo : Store qij to list of axuliary functions
            }
        }
        onq[io] = ii;
        oQs[io] = Q;
        if(bNormalize){
            double renorm  = sqrt(1./Q);
            double renorm2 = renorm*renorm;
            printf( "project orb[$i]: Q %g renorm %g renorm2 %g \n", io, Q, renorm, renorm2 );
            for(int i=i0; i<i0+perOrb; i++){ ecoef[i] *=renorm;  };
            for(int i= 0; i<ii       ; i++){ Qs   [i] *=renorm2; };
        }
        // ToDo: Renormalize also  rhos?
        return Ek;
    }

    void projectOrbs(bool bNormalize){   // project density of all orbitals onto axuliary charge representation ( charges, dipoles and axuliary functions )
        int nqOrb = perOrb*(perOrb+1)/2;
        int i0 =0;
        int ii0=0;
        for(int io=0; io<nOrb; io++){
            //int i0  = getOrbOffset(jo);
            //oQs[io] =
            projectOrb( io, odip[io], bNormalize );
            //projectOrb(  io, rhoP+ii0, rhoQ+ii0, rhoS+ii0, odip[io], true );
            //projectOrb(io, ecoefs+i0, erho+irho0, erhoP+irho0, odip[io] );
            ii0+=nqOrb;
            i0 +=perOrb;
        }
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
    void fromRhoDiag( int i, int ij, double& aij ){
        /// This function function cares about diagonal density terms rho_ii = wi * wi
        //Vec3d  pi  = epos [i];
        //double ci  = ecoef[i];
        //double si  = esize[i];
        Vec3d  Fpi = rhofP[ij];
        double Fqi = rhofQ[ij];
        double Fsi = rhofS[ij];
        //double Eqi = rhoEQ[ij];
        printf( "fromRho[%i]ii[%i] Fpi(%g,%g,%g) Fqi %g | Qi %g \n", i, ij, Fpi.x,Fpi.y,Fpi.z, rhoEQ[ij],    rhofQ[ij] );
        efpos [i].add( Fpi );
        efsize[i] += Fsi*aij*2;
    }


    //Vec3d fromRho( int i, int j, int ij ){
    void fromRho( int i, int j, int ij, double& Cij ){
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

        Vec3d  p;
        double s;
        double dSsi,dSsj;
        Vec3d  dXsi,dXsj;
        double dXxi,dXxj;
        double dCr;
        double dCsi,dCsj;
        Vec3d  dCdp;

        // NOTE: we must compute it twice
        // 1) in toRho() to obtain charges and positions
        // 2) here to obtain derivatives dXxi,dXxj,dCr
        // We cannot connect it because CoublombElement needs do (1)before and (2)after itself
        Cij = Gauss::product3D_s_deriv(
            si,   pi,
            sj,   pj,
            s ,   p ,
            dSsi, dSsj,
            dXsi, dXsj,
            dXxi, dXxj,
            dCsi, dCsj, dCr
        );

        Vec3d  Fpi = rhofP[ij];  //   fij = Rij*(-fr*qij);   ... from CoulombElement
        double Fqi = rhofQ[ij];
        double Fsi = rhofS[ij];
        double Eqi = rhoEQ[ij];

        double fsj = Fsi*dSsj + Fpi.dot( dXsj );
        double fsi = Fsi*dSsi + Fpi.dot( dXsi );
        Vec3d  fxi = Fpi*dXxi;  // ToDo : dSij/dxi == 0
        Vec3d  fxj = Fpi*dXxj;  //            dXxi == 0.5

        //fsi = ci*cj * ( fsi*aij + E*dCr );
        //fxi = ((Vec3d){1,1,1}) * dXxi ;

        dCdp = Rij*(-2*dCr*ci*cj);
        DEBUG_dQdp = dCdp;

        Vec3d Fq  = dCdp*Eqi;
        Vec3d Fxi = fxi + Fq;
        Vec3d Fxj = fxj - Fq;

        efpos [i].add( Fxi ); // TODO : Why 0.25 factor ? There is no reason for this !!!!!
        efpos [j].add( Fxj );
        //efpos [i].add( dCdp*Eqi ); // TODO : Why 0.25 factor ? There is no reason for this !!!!!
        //efpos [j].add( dCdp*Eqi );
        efsize[i] += fsi*Cij;
        efsize[j] += fsj*Cij;


        // Perhaps found a problem !!!!!
        //  rhoQ[ij]  = qi*qj*Sab*Scd is probably wrong !!!!!
        //  We are missing the j-orbital part  cj*Sj = cj*Scd
        // TODO : Due to summation we have to multiply it by Qj but not by Qi


        // --- Derivatives ( i.e. Forces )
        //printf( "fromRho r %g s %g E %g Fx %g fx %g  \n", sqrt(r2), s, Eqi,     );
        //printf( "fromRho r %g s %g | E %g e %g qij %g(%g) | F %g fx %g dQij %g \n", sqrt(r2), s, Eqi*rhoQ[ij],Eqi,rhoQ[ij],  Fxi.x,fxi.x,dCdp.x );
        //printf( "fromRho r %g s %g | E %g e %g qij %g(%g) | F %g fx %g dQij %g \n", sqrt(r2), s, Eqi*rhoQ[ij],Eqi,rhoQ[ij],Cij,  Fxi.x,fxi.x,dCdp.x );
        //printf( "fromRho r %g s %g | E %g e %g qij %g(%g) | fxi %g Fxi %g Fpi %g dQij %g \n", sqrt(r2), s, Eqi*rhoQ[ij],Eqi,rhoQ[ij],Cij,  fxi.x,Fxi.x,Fpi.x,dCdp.x );

        printf( "fromRho r %g  Eqi %g Cij %g | Fpi %g dXxi %g fxi %g Fxi %g \n", sqrt(r2), Eqi,Cij,  Fpi.x, dXxi, fxi.x, Fxi.x );

        //printf( "[%i,%i,%i] fxi %g Fpi %g dXxi %g \n",   i,j,ij,   fxi.x, Fpi, dXxi );
        //printf( "fsi, fsj, aij %g %g %g \n", fsi, fsj, aij );
        //printf( "fromRho[%i,%i][%i] Fpi(%g,%g,%g) Eqi %g Fqi \n", i, j, ij, Fpi.x,Fpi.y,Fpi.z, Eqi, Fqi);
        //printf( "fromRho[%i,%i][%i]  Q %g    fxi(%g,%g,%g) Eqi %g dCdp(%g,%g,%g) \n", i, j, ij,  1./rhoQ[ij],  fxi.x,fxi.y,fxi.z, Eqi, dCdp.x,dCdp.y,dCdp.z );
        //printf( "fromRho[%i,%i][%i] Eqi %g dCdp(%g,%g,%g) \n", i, j, ij, Eqi, dCdp.x,dCdp.y,dCdp.z );

        //dCsi*=-0.42;
        //dCsj*=-0.42;

        //return dCsi;
        //return Rij*(-2*dCr*ci*cj);
        //return dCr;
    }


    void assembleOrbForces_fromRho(int io ){
        // NOTE : This is less efficient but more mocular version of assembleOrbForces which use fromRho function per each element
        //        It is used for debuging since it is more modular, but after everything works original assembleOrbForces should be prefered
        //   line_Fana->ys[i]  = 0.5*solver.efpos[0].x + E_*line_dQi_ana->ys[i];
        //   Fana              = efpos                + EK*dQdp; ..... TODO
        int i0    = getOrbOffset(io);
        int irho0 = getRhoOffset(io);
        int ii    = irho0;
        double aij;
        double dCsi;
        double dCsj;
        Vec3d  dQdp;
        for(int i=i0; i<i0+perOrb; i++){
            // TODO : What to do with diagonal elements ? That is diagonal electron could rho_ii ( It should not change upon moving ? Or should only due to renormalization ? )
            //        Is the problem caused by renormalization ????????
            //fromRho( i, i, ii, aij, dCsi, dCsj, dQdp );
            fromRhoDiag( i, ii, aij );  // Why is this zero ????
            ii++;
            for(int j=i0; j<i; j++){
                //printf( "assembleOrbForces_fromRho[%i,%i][%i] \n", i, j, ii  );
                //fromRho( i, j, ii, aij, dCsi, dCsj, dQdp );
                fromRho( i, j, ii, aij );
                // ToDo: ad dQdp    //   Fana              = efpos                + EK*dQdp; ..... TODO
                // We need to copy EK from somewhere
                //dQdp
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

        fs *= qij*4;
        Vec3d fij = Rij*(-fr*qij);
        rhofP[i].add(fij);   rhofP[j].sub(fij);
        rhofS[i] -= fs*si;   rhofS[j] -= fs*sj; // Q: ??? Should not this be switched (i<->j)  rhofS[i] -= fs*sj instead of rhofS[i] -= fs*si ???
        //rhofP[i].add(Rij*(fr*qj));   rhofP[j].sub(Rij*(fr*qi));
        //rhofS[i] -= fs*si*qj;        rhofS[j] -= fs*sj*qi; // Q: ??? Should not this be switched (i<->j)  rhofS[i] -= fs*sj instead of rhofS[i] -= fs*si ???
        rhofQ[i] += E*qj;    rhofQ[j] += E*qi;  // ToDo : need to be made more stable ... different (qi,qj)
        rhoEQ[i] += E*qj;    rhoEQ[j] += E*qi;  // Coulombic energy per given density could (due to other density clouds)

        //printf( "CoublombElement r %g s %g E %g fr %g qij %g frq %g fij %g \n", r, s, E, fr, qij, frq, fij.x );
        //printf( "CoublombElement[%i,%i] E %g rhoEQij %g %g \n", i, j, E, rhoEQ[i], rhoEQ[j] );

        return E;
    }

    double CoulombOrbPair( int io, int jo ){
        /// Calculate Elecrostatic Energy&Force for projected density of orbitas @io and @jo
        int i0 = getRhoOffset(io);
        int j0 = getRhoOffset(jo);
        int nio = onq[io];
        int njo   = onq[jo];
        double Ecoul=0;
        //printf( "CoulombOrbPair[%i,%i] nV %i \n", io, jo, njo );
        for(int i=i0; i<i0+nio; i++){
            Vec3d  pi = rhoP[i];
            double qi = rhoQ[i];
            double si = rhoS[i];
            for(int j=j0; j<j0+njo; j++){
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
                //double Eqq  = CoulombGauss( r, s*2, fr, fs, qij );
                //fs*=4;

                double E  = Gauss::Coulomb( r, s*2, fr, fs );

                //printf( "CoulombOrbPair[%i,%i][%i,%i] qij %g(%g,%g) r %g E %g \n", io,jo, i,j, qij,qi,qj,  r, E );

                //printf(  " [%i,%i] q %g r %g E %g \n", i, j, qij, r, Eqq );

                fr *= qij;
                fs *= qij*4;

                // see    InteractionsGauss.h :: addCoulombGauss( const Vec3d& dR, double si, double sj, Vec3d& f, double& fsi, double& fsj, double qq ){

                // --- Derivatives (Forces)
                Vec3d fij = Rij*(-fr);
                rhofP[i].add(fij);   rhofP[j].sub(fij);
                rhofS[i] -= fs*si;   rhofS[j] -= fs*sj;
                rhofQ[i] += E*qj;    rhofQ[j] += E*qi; // ToDo : need to be made more stable ... different (qi,qj)
                rhoEQ[i] += E;       rhoEQ[j] += E;// Coulombic energy per given density could (due to other density clouds)

                //printf( "CoulombOrbPair[%i,%i] E %g rhoEQij %g %g \n", i, j, E, rhoEQ[i], rhoEQ[j] );

                //printf( "[%i,%i] E %g r %g sij %g(%g,%g) q %g(%g,%g) \n", i,j,  Eqq , r, s,si,sj, qi*qj,qi,qj );
                Ecoul      += E*qi*qj;
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
            for(int jo=0; jo<io; jo++){
                Vec3d dop = opos[jo] - opi;
                double r2 = dop.norm2();
                printf( "evalElectrostatICoulomb[%i,%i]  %g <? %g \n", io, jo, r2,RcutOrb2 );

                E += CoulombOrbPair( io, jo );

                /*
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
        return E;
    }

    void ExchangeOrbPair(int io, int jo){
        // It may be actually easier to evaluate Hartree-Fock like exchange, because it is linear
        // see:   https://en.wikipedia.org/wiki/Exchange_interaction#Exchange_of_spatial_coordinates
        // https://chemistry.stackexchange.com/questions/61176/what-is-the-exchange-interaction
        // => Exchange integral should be rather easy to compute, just by calculating overlap between the two functions   K = <rho_ij(r1)| 1/|r1-r2| |rho_ij(r2)>
        // scheme:
        //  1) for pair of orbital, construct auxiliary rho_ij (by expansion on auxuliary functions)
        //  2) callculate self-convolution of this auxuliary function

        const double R2Safe = sq(0.1);
        // we first project on axuliary density
        Vec3d  posij[nqOrb];
        double rhoij[nqOrb];
        double szij [nqOrb];
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
                if(r2>Rcut2) continue;
                double cj = ecoef[j];
                double sj = esize[j];

                double dSr, dSsi, dSsj;
                const double resz = M_SQRT2; // TODO : PROBLEM !!!!!!!   getOverlapSGauss and getDeltaTGauss are made for density gaussians not wave-function gaussians => we need to rescale "sigma" (size)
                double Sij = getOverlapSGauss( r2, si*resz, sj*resz, dSr, dSsi, dSsj );
                Vec3d  pij;
                double sij;
                double Cij = Gauss::product3D_s( si, pi, sj, pj, sij, pij );
                double cij = ci*cj;
                //double qij = Sij*cij*2; // factor 2  because  Integral{(ci*fi + cj*fj)^2} = (ci^2)*<fi|fi> + (cj^2)*<fj|fj> + 2*ci*cj*<fi|fj>
                double qij = Sij*cij*2;
                rhoij[ij]  = qij;           // ToDo: we may perhas evaluate it on even more extended axuliary basis
                posij[ij]  = pij;
                szij [ij]  = sij;
                ij++;
            }
        }
        // --- self-coulomb  of   rho_ij(r)
        double E = 0;
        for(int i=0; i<nqOrb; i++){
            Vec3d  pi = posij[i];
            double qi = rhoij[i];
            double si = szij [i];
            for(int j=0; j<i; j++){
                Vec3d Rij = posij[j] - pi;
                double r2 = Rij.norm2();
                double qj = rhoij[j];
                double sj = szij [j];

                double r  = sqrt(r2);
                double s2 = si*si + sj*sj;
                double s  = sqrt(s2);
                double fr,fs;
                double Cij = CoulombGauss( r, s*2, fr, fs, qi*qj );
                E += Cij; // ToDo : we should perhaps use some cutoff for large distances ?
                //E += qi*qj/sqrt(r2+R2Safe);  // ToDo  : we should perhaps calculate this more properly considering true shape of axuiliary density function
                //
            }
        }

    }

    void evalExchangeCorrelation(){
        // not sure how exactly to do that now - no using real space grid since it is slow
    }

    double eval(){
        double E=0;
        clearAuxDens();
        //projectOrbs( true );
        projectOrbs( false );
        reportCharges();
        //E += evalPauli();
        E += evalElectrostatICoulomb(); // repulsion between aux-density clouds => should not distinguish density terms here
        //for(int i=0; i<nOrb; i++) assembleOrbForces(i);
        for(int i=0; i<nOrb; i++) assembleOrbForces_fromRho(i); // DEBUG
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
