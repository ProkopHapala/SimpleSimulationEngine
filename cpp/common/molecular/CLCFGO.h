
#ifndef CLCFGO_h
#define CLCFGO_h

/*

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
        rho(r)=phi(r) = sum_i{ c_i chi_i(r)} sum_j{ c_j chi_j(r)}
 - The main problem is evaluation of exchange-correlation. Which is non-linear and cannot be easily experesed as combination of orbitals.
    - This can be done on a real-space grid, but that is very costly perhaps
    - We can use some approximation from Fireball ... e.g. McWeda approximation of Pavel Jelinek
    - We can oalso sort of fit some empirical approximation using machine learning etc.


    IDEA 1) - multiple function per one center. Each can have different cutoff.
       - Variation of coeficients is much faster than variation of position ( re-evaluation of integrals is not required )


 ## Make alternative with Gaussian Orbitals (instead of numerica basis)
    - It makes it possible to express auxuliary density functions exactly
        rho_ij(r)
           = Gauss(r-Ri) * Gauss(r-Rj)
           = exp( -wi*( (x-xi)^2 + (y-yi)^2 + (y-yi)^2 ) ) * exp( -wj*( (x-xj)^2 + (y-yj)^2 + (y-yi)^2 ) )
           = exp( -wi*( (x-xi)^2 + (y-yi)^2 + (y-yi)^2 ) -wj*( (x-xj)^2 + (y-yj)^2 + (y-yj)^2 )  )
           x : wi*(x-xi)^2 + wj*(x-xj)^2   =    (wi+wj)x^2   - 2*(wi*xi+wj*xj)x - (wi*xi^2+wj*xj^2)
        // http://www.tina-vision.net/docs/memos/2003-003.pdf
*/






#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"

#include "spline_hermite.h"
#include "GaussianBasis.h"
#include "Forces.h"

#include "Grid.h"
#include "AOIntegrals.h"


class CLCFGO{ public:

    double Rcut    =6.0;  // cutoff beyond which two basis functions chi has no overlap
    double RcutOrb =9.0;  // cutoff beoud which orbital (i.e. localized cluster of basis functions) has no overlap
    double Rcut2     =Rcut*Rcut;
    double RcutOrb2  =RcutOrb*RcutOrb;

    double   dsamp =0.1; // [A] sampling step for splines describing radial splines
    double  idsamp =1/dsamp;
    int nsamp     =0;   // number of samples in splines
    int nsampI    =0;
    int nsampMem  =0;
    int nsampIMem =0;

    int natypes  =0;

    int natom =0; // number of atoms (nuclei, ions)
    int perOrb=0; // number of spherical functions per orbital
    int nOrb  =0; // number of single-electron orbitals in system
    // this is evaluated automaticaly
    int nBas  =0; // number of basis functions
    int nqOrb =0; // number of charges (axuliary density elements) per orbital
    int nQtot =0; // total number of charge elements

    // atoms (ions)
    Vec3d*  apos   =0;  // positioon of atoms
    Vec3d*  aforce =0;  // positioon of atoms
    Vec3d*  aQs    =0;  // charge of atom
    int*    atype  =0;  // type of atom (in particular IKinetic pseudo-potential)

    // orbitals
    Vec3d*  opos =0;   // store positions for the whole orbital
    Vec3d*  odip =0;   // Axuliary array to store dipoles for the whole orbital
    double* oEs  =0;   // orbital energies
    int*    onq  =0;   // number of axuliary density functions per orbital

    // orbital basis functions
    Vec3d*  epos  =0; // position of spherical function for expansion of orbitals
    Vec3d*  epfs  =0; //   force acting on position of orbitals
    double* ecoefs=0; // expansion coefficient of expansion of given orbital
    double* ecfs  =0; //   force acting on combination coefficnet of orbitals
    // density approx
    //double  eQs  =0; // charge (occupation) of all orbitals is 1 (by definition)
    // denisty axuliary

    // denisty expansion on auxuliary basis
    double* erho  =0;  // temporary array to store density projection on pair overlap functions
    double* erhoE =0; // energy of density axuliary function
    Vec3d * erhoP =0; // position of density axuliary function
    Vec3d * erhoF =0; // force  on density axuliary function

    // ToDo: implement different basis function types
    //int* etypes = 0; // currently not used

    //PairInteractionTable** interactions;
    //PairInteractionTable** eInteractions;

    // Interaction tables   [itype][jtype][isample]
    //double*** IOverlap = 0;   // < chi_i(r) |    1    | chi_i(r) >  Overlap interactin between two  spherical bais functions
    //double*** IKs = 0;   // < chi_i(r) | Laplace | chi_i(r) >  Kinetic energy interactin between two spherical functions

    // ---- Radial Function Tables
    double* Wfs  = 0;  // wavefunction |Chi(r)>         radial wave function table
    double* Wf2s = 0;  // density      |Chi(r)>^2       denisty function from wavefunction
    double* Vfs  = 0;  // V(R) = (1/(r-R))*|Chi(r)>^2   potential generated by a wave function
    // ----- Integral Tables
    double* IOverlap = 0;  // Overlap  = < chi_i(r) |    1       | chi_j(r) >  Overlap interactin between two  spherical bais functions
    double* IKinetic = 0;  // Kinetic  = < chi_i(r) | Laplace    | chi_j(r) >  Kinetic energy interactin between two spherical functions
    double* ICoulomb = 0;  // Coulomb  = < chi_i(ri) | 1/|ri-rj| | chi_i(rj) >
    //double*** IPPWf2 = 0; // < rho(r) V(r) >  interactin between spherical_eletron_density_function and spherical electrostatic_potential_function
    double** PPs   = 0; // V(r) pseudopotential of atom
    double** IPPWf2 = 0; // < rho(r) V(r) >  interactin between spherical_eletron_density_function and spherical electrostatic_potential_function

    // ======= Functions

    void realloc( int natom_, int nOrb_, int perOrb_, int nsamp_, int natypes_ ){
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
            _realloc( onq , nOrb);

            // orbital basis functions
            _realloc( epos  , nBas );
            _realloc( ecoefs, nBas );
            //_realloc( eEs , nBas );
            _realloc( epfs,   nBas );
            _realloc( ecfs,   nBas );

            // --- density Axuliary
            _realloc( erho,  nQtot );
            _realloc( erhoP, nQtot );
            _realloc( erhoF, nQtot );
            _realloc( erhoE, nQtot );
        }
        if(nsamp != nsamp_){
            nsamp     = nsamp_;
            nsampI    = nsamp*2;
            nsampMem  = nsamp+3;
            nsampIMem = nsampI+3;
            _realloc( Wfs , nsampMem  );
            _realloc( Wf2s, nsampMem  );
            _realloc( Vfs , nsampMem  );
            _realloc( IOverlap , nsampIMem );
            _realloc( IKinetic , nsampIMem );
            _realloc( ICoulomb , nsampIMem );
            for(int i=0;i<natypes;i++){
                if(IPPWf2[i])delete [] IPPWf2[i];
                if(PPs  [i])delete [] PPs  [i];
            }
            if( natypes != natypes_ ){
                natypes = natypes_;
                _realloc( IPPWf2, natypes );
                _realloc( PPs  , natypes );
            }
            for(int i=0;i<natypes;i++){
                if(PPs   [i]) PPs   [i]=new double[nsampMem ];
                if(IPPWf2[i]) IPPWf2[i]=new double[nsampIMem];
            }
        }
    }

    inline int getOrbOffset(int iorb)const{ return iorb*nOrb;  }
    inline int getRhoOffset(int iorb)const{ return iorb*nqOrb; }

    inline void setRcut( double Rcut_ ){
        Rcut   = Rcut_;
        dsamp  = Rcut/nsamp;
        idsamp = 1/dsamp;
    }

    inline bool checkOrbitalCutoff(int i, int j)const{
        double r2 = (opos[i]-opos[j]).norm2();
        return r2<RcutOrb2;
    }

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

    //inline double evalShotRange( Vec3d Rij, double* Is, int i, int j ){
    inline double evalPauli( int io, int jo ){
        // Pauli Energy between two orbitals depend on IKinetic overlap  ~S^2
        //  See addPauliGauss() in common/molecular/InteractionsGauss.h
        int i0 = getOrbOffset(io);
        int j0 = getOrbOffset(jo);
         // ToDo : we may allow different number of orbitals per electron later (?)
         for(int i=i0; i<i0+perOrb; i++){
            Vec3d pi = epos[i];
            for(int j=j0; j<j0+perOrb; j++){
                // ToDo: Rcut may be read from
                Vec3d Rij = epos[j]-pi;
                double r2 = Rij.norm2();
                if(r2>Rcut2)continue; //{ fij = Vec3dZero; return 0; }

                double r = sqrt(r2);
                double S,dS;
                int    i   = (int)r;
                double x   =  r - i;
                Spline_Hermite::valdval( x, S, dS, IOverlap[i], IOverlap[i+1], IOverlap[i+2], IOverlap[i+3] ); // Overlap

                //Spline_Hermite::valdval( x, S, dS, IOverlap[i], IOverlap[i+1], IOverlap[i+2], IOverlap[i+3] ); // Overlap
                //Spline_Hermite::valdval( x, K, dK, IKs[i], IKs[i+1], IKs[i+2], IKs[i+3] ); // Kinetic energy
                //Pauli(); // Pauli is like S^2 ? ... that is not linear ... should we rather use density overlap ?


                double ci = ecoefs[i];
                double cj = ecoefs[j];
                double cij = ci*cj;

                double eP=0,fP=0,dcP=0;
                // Copied from addPauliGauss() in common/molecular/InteractionsGauss.h
                //double TfS = T*fS;
                //fsi +=         -(dTsi*eS + TfS*dSsi)*KRSrho.y;
                //fsj +=         -(dTsj*eS + TfS*dSsj)*KRSrho.y;
                //f.add_mul( dR, (dTr *eS + TfS*dSr )*KRSrho.x*KRSrho.x ); // second *KRSrho.x because dR is not multiplied

                Vec3d fij; fij.set_mul( Rij, cij*fP/r ); // ToDo:  evaluate fP properly

                // ToDo:  ecfs should be set as well
                ecfs[i] +=dcP;
                ecfs[j] +=dcP;
                epfs[i].add(fij);
                epfs[j].sub(fij);
                return eP*cij;
            }
        }
    }

    void evalPauli(){ // evaluate Energy components given by direct wave-function overlap ( below cutoff Rcut )
        for(int io=0; io<nOrb; io++){
            for(int jo=0; jo<nOrb; jo++){
                if( !checkOrbitalCutoff(io,jo) ) continue; // optimization, not to calculate needless interactions
                //interpolate( epos[j] - epos[i], Is[ityp][jtyp], f ); // ToDo: different basis function types
                evalPauli( io, jo );
            }
        }
    }

    double projectOrb(int io, double* coefs, double* rhos, Vec3d* erhoPs, Vec3d& dip){ // project orbital on axuliary density functions
        Spline_Hermite::Sampler<double> spline;
        int i0=getOrbOffset(io);
        double Q=0;
        double T=0; // kinetic energy
        dip         = Vec3dZero;
        Vec3d qcog  = Vec3dZero;
        Vec3d oqcog = opos[io];
        for(int i=i0; i<i0+perOrb; i++){
            Vec3d  pi  = epos[i];
            double ci  = coefs[i];
            double qii = ci*ci;
            qcog.add_mul(  pi, qii );
            *rhos  =qii; rhos++;
            *erhoPs=pi;  erhoPs++;
            Q += ci*ci;
            for(int j=i0; j<i; j++){
                Vec3d pj  = epos[j];
                Vec3d Rij = pj-pi;
                double r2 = Rij.norm2();
                if(r2>Rcut2) continue;
                double S,dS,  T,dT;
                //Spline_Hermite::valderiv( sqrt(r2), S, dS, IOverlap ); // Overlap
                double r   = sqrt(r2);
                int    i   = (int)r;
                double x   =  r - i;
                //Spline_Hermite::valdval( x, S, dS, IOverlap[i], IOverlap[i+1], IOverlap[i+2], IOverlap[i+3] ); // Overlap
                //Spline_Hermite::valdval( x, K, dK, IKinetic[i], IKinetic[i+1], IKinetic[i+2], IKinetic[i+3] ); // Kinetic energy
                spline.seek  ( sqrt(r2) );
                spline.preval( IOverlap ); S = spline.y(); dS = spline.dy();
                spline.preval( IKinetic ); T = spline.y(); dT = spline.dy();
                double cij = ci*coefs[j];
                double qij = S*cij;
                qcog.add_mul( pi+pj, 0.5*qij );
                Q +=   qij;
                T += T*cij;
                *rhos    += qij;
                *erhoPs   = (pi+pj)*0.5;   // center of axuliary overlap density function in the middle between the two wavefunctions
                rhos++; erhoPs++;
                // ToDo : Store qij to list of axuliary functions
            }
        }
        double renorm = sqrt(1./Q);
        for(int i=i0; i<i0+perOrb; i++){ coefs[i] *=renorm;  };
        // ToDo: Renormalize also  rhos?
        return Q;
    }

    double projectOrbs(){   // project density of all orbitals onto axuliary charge representation ( charges, dipoles and axuliary functions )
        int nqOrb = perOrb*(perOrb+1)/2;
        int i0=0;
        int irho0=0;
        for(int io=0; io<nOrb; io++){
            //int i0  = getOrbOffset(jo);
            //oQs[io] =
            projectOrb(io, ecoefs+i0, erho+irho0, erhoP+irho0, odip[io] );
            irho0+=nqOrb;
            i0+=perOrb;
        }
    }

    void CoulombOrbPair( int nrho, int nV, double* rhos, Vec3d* rhoPs, Vec3d* rhoFs,  double* Vs, Vec3d* VPs, Vec3d* VFs ){
        for(int i=0; i<nV; i++){
            Vec3d  pi  = rhoPs[i];
            double qi  = rhos [i];
            for(int j=0; j<nrho; j++){
                Vec3d  pj = rhoPs[j];
                Vec3d Rij = pj-pi;
                double r2 = Rij.norm2();
                if(r2>Rcut2) continue;
                double C,dC;
                double r   = sqrt(r2);
                int    i   = (int)r;
                double x   =  r - i;
                double qj  = rhos [j];
                Spline_Hermite::valdval( x, C, dC, ICoulomb[i], ICoulomb[i+1], ICoulomb[i+2], ICoulomb[i+3] );  // coulomb integrals
            }
        }
    }

    void evalElectrostatICoulomb( ){
        double r2safe = 0.1;
        double E = 0;
        for(int io=0; io<nOrb; io++){
            int i0 = getRhoOffset(io);
            Vec3d opi  = opos[io];
            Vec3d dipi = odip[io];
            int   nqi  = onq [io];
            for(int jo=0; jo<io; jo++){
                Vec3d dop = opos[jo] - opi;
                double r2 = dop.norm2();
                if( r2<RcutOrb2 ){
                //if( checkOrbitalCutoff( i, j ) ){
                    /// Calculate detailed short-range electrostatICoulomb
                    int j0 = getRhoOffset(jo);
                    // ToDo : nrho,nV may vary
                    CoulombOrbPair( nqi, onq[jo], erho+i0, erhoP+i0, erhoF+i0,  erho+j0, erhoP+j0, erhoF+j0 );
                }else{
                    double ir2 = 1/(r2+r2safe);
                    double ir  = sqrt(ir2);
                    // Calculate approximate long-range electrostatICoulomb (using charges and dipoles)
                    const Vec3d& dipj = odip[jo];
                    E += ir +                                   // charge-charge
                      + dop.dot(dipi) + dop.dot(dipj) *ir*ir2  // diple-charge
                      + dipi.dot(dipj);                        // dipole-dipole
                }
            }
        }
    }

    void ExchangeOrbPair(int io, int jo){
        const double R2Safe = sq(0.1);
        // we first project on axuliary density
        double rhoij[nqOrb];
        Vec3d  posij[nqOrb];
        int i0=io*perOrb;
        int j0=jo*perOrb;
        int ij=0;
        for(int i=i0; i<i0+perOrb; i++){
            Vec3d  pi  = epos  [i];
            double ci  = ecoefs[i];
            for(int j=j0; j<j0+perOrb; j++){
                Vec3d pj  = epos[j];
                Vec3d Rij = pj-pi;
                double r2 = Rij.norm2();
                if(r2>Rcut2) continue;
                double S,dS;
                //Spline_Hermite::valderiv( sqrt(r2), S, dS, IOverlap ); // Overlap
                double r   = sqrt(r2);
                int    i   = (int)r;
                double x   =  r - i;
                Spline_Hermite::valdval( x, S, dS, IOverlap[i], IOverlap[i+1], IOverlap[i+2], IOverlap[i+3] ); // Overlap
                double cij = ci*ecoefs[j];
                double qij = S*cij;
                rhoij[ij] = qij;           // ToDo: we may perhas evaluate it on even more extended axuliary basis
                posij[ij] = (pi+pj)*0.5;
                ij++;
            }
        }
        double E = 0;
        for(int i=0; i<nqOrb; i++){
            Vec3d  pi = posij[i];
            double qi = rhoij[i];
            for(int j=0; j<i; j++){
                Vec3d Rij = posij[j] - pi;
                double r2 = Rij.norm2();
                double qj = rhoij[j];
                E += qi*qj/sqrt(r2+R2Safe);  // ToDo  : we should perhaps calculate this more properly considering true shape of axuiliary density function
                                            // ToDo2 : what about coulomb constant ?????
            }
        }

    }

    void evalFockExchange(){
        // It may be actually easier to evaluate Hartree-Fock like exchange, because it is linear
        // see:   https://en.wikipedia.org/wiki/Exchange_interaction#Exchange_of_spatial_coordinates
        // https://chemistry.stackexchange.com/questions/61176/what-is-the-exchange-interaction
        // => Exchange integral should be rather easy to compute, just by calculating overlap between the two functions   K = <rho_ij(r1)| 1/|r1-r2| |rho_ij(r2)>
        // scheme:
        //  1) for pair of orbital, construct auxiliary rho_ij (by expansion on auxuliary functions)
        //  2) callculate self-convolution of this auxuliary function
    }


    void evalExchangeCorrelation(){
        // not sure how exactly to do that now - no using real space grid since it is slow
    }


    double eval(){
        projectOrbs();
        evalElectrostatICoulomb();
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
        Vec3d*  Ps = epos  +i0;
        double* Cs = ecoefs+i0;
        printf( "DEBUG orb2grid io %i i0 %i perOrb %i \n", io, i0, perOrb );
        for(int i=0; i<perOrb; i++){
            printf( "wf[%i] C(%e) P(%g,%g,%g)\n", i, Cs[i], Ps[i].x,Ps[i].y,Ps[i].z );
        }
        return evalOnGrid( gridShape, [&](int ig, const Vec3d& pos, double res){
            Spline_Hermite::Sampler<double> spline;
            double wfsum = 0.0;
            for(int i=0; i<perOrb; i++){
                Vec3d dR  = pos - Ps[i];
                double r2 = dR.norm2();
                //printf( "orb2grid[%i] pos(%g,%g,%g) epos(%g,%g,%g) \n", ig, pos.x, pos.y, pos.z,   epos[i].x, epos[i].y, epos[i].z );
                //printf( "orb2grid[%i] r2 %g  Rcut2 %g |  pos(%g,%g,%g) epos(%g,%g,%g) \n", ig, r2, Rcut2, pos.x, pos.y, pos.z,   epos[i].x, epos[i].y, epos[i].z );
                if(r2>=Rcut2){
                    buff[ig] = 0.0;
                    continue;
                }
                //double r   = sqrt(r2);
                //int    ix  = (int)r;
                //double x   =  r - ix;
                //double y   = Spline_Hermite::val( x, Wfs[ix], Wfs[ix+1], Wfs[ix+2], Wfs[ix+3] );
                spline.prepare( sqrt(r2)*idsamp, Wfs );
                wfsum += spline.y() * Cs[i];
                //wfsum += spline.y();
                //wfsum += exp( -r2 );
                //wfsum += exp( -r2 ) * Cs[i];
                //printf( "orb2grid[%i] r %g  y(r) %g  ci*y(r) %g \n", ig, sqrt(r2), spline.y(), ecoefs[i] );
            }
            buff[ig] = wfsum;
        });
    }

    // ========== EvalIntegrals

    void prepareWfs( bool bNormalize ){
        double sumQ = 0;
        double Qtot = 0;
        for(int i=0; i<nsampMem; i++){
            double r   = (i-1)*dsamp;
            double wf  = Wfs[i];
            double rho = wf*wf;
            if(i>0){
                //double rmid = r+dsamp*0.5;
                double rmid = r;
                //double rmid = r;
                sumQ       += rho*( (4*M_PI)*rmid*rmid )*dsamp;
                if(i==nsamp)Qtot=sumQ;
                //printf( "Wf[%i] %g \n", i, wf );
            }
            Wf2s[i]    = rho;       // density
            Vfs [i]    = sumQ/r; // potential of density basis function
        }
        Vfs[0] = Vfs[1]-(Vfs[2]-Vfs[1])*idsamp;
        if(bNormalize){
            double renorm =sqrt(1/Qtot);
            double renorm2=renorm*renorm;
            //printf( "Qtot %g Vspere %g \n", Qtot, (4/3.)*M_PI*(Rcut*Rcut*Rcut) );
            printf( "Qtot %g renorm %g \n", Qtot, renorm );
            for(int i=0; i<nsampMem; i++){
                Wfs [i] *= renorm;
                Wf2s[i] *= renorm2;
                Vfs [i] *= renorm2;
                printf( "Wf[%i] %g \n", i, Wfs[i] );
            }
        }
    }

    void prepareIntegralTables(){
        prepareWfs( true );
        integrateSK( nsampMem, 0, nsampI, dsamp, Rcut, dsamp, Wfs,     Wfs,  IOverlap+1, IKinetic+1 ); // Basis Overlap and Kinetic energy
        //integrateS ( nsampMem, 0, nsampI, dsamp, Rcut, dsamp, Wf2s,    Vfs,  ICoulomb+1  );          // Coulomb betwween density functions
        //integrateS ( nsampMem, 0, nsampI, dsamp, Rcut, dsamp, PPs[0],  Wf2s, IPPWf2[0]+1 );          // Coulomb with ion pseudopotential
    }

};

#endif
