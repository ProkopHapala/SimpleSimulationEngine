
#ifndef GaussianBasis_h
#define GaussianBasis_h

#include "fastmath.h"


namespace Gauss{

static int iDEBUG = 0;

// TODO : this should go elsewhere (physical constants or something)
const double const_hbar_SI      = 1.054571817e-34;    ///< [J.s]  #6.582119569e-16 # [eV/s]
const double const_Me_SI        = 9.10938356e-31;     ///< [kg]
const double const_e_SI         = 1.602176620898e-19; ///< [Coulomb]
const double const_eps0_SI      = 8.854187812813e-12; ///< [F.m = Coulomb/(Volt*m)]
const double const_eV_SI        = 1.602176620898e-19; ///< [J]
const double const_Angstroem_SI = 1.0e-10;

const double const_K_SI   =  const_hbar_SI*const_hbar_SI/const_Me_SI;
const double const_El_SI  =  const_e_SI*const_e_SI/(4.*M_PI*const_eps0_SI);
const double const_Ry_SI  = 0.5 * const_El_SI*const_El_SI/const_K_SI;

const double const_Ry_eV  = 13.6056925944;
const double const_El_eVA = const_El_SI/( const_e_SI*const_Angstroem_SI );
const double const_K_eVA  = (const_El_eVA*const_El_eVA)/(2*const_Ry_eV);
const double const_Ke_eVA = const_K_eVA*1.5;


#define _Gauss_Tr0(s)      const_K_eVA*-8.35249199525*s
#define _Gauss_dTr0_ds     const_K_eVA*-8.35249199525

#define _Gauss_sij_aux( si, sj ) \
    double si2  = si*si; \
    double sj2  = sj*sj; \
    double s2   = si2 + sj2; \
    double is2  = 1/s2; \
    double is4  = is2*is2; \

#define _Gauss_product(pi,pj,si,sj)    \
    double sqrt_is2 = sqrt(is2);       \
    double size_ij  = si*sj*sqrt_is2;            \
    Vec3d  pos_ij   = pj*(si2*is2) + pi*(sj2*is2); \

#define _Gauss_product_derivs(dSsi,dSsj,dXsi,dXsj,dXxi,dXxj){ \
    double inv3_2 = is2 * sqrt_is2; \
    dSsi   = sj*sj2*inv3_2; \
    dSsj   = si*si2*inv3_2; \
    dXsi   = dp*(-2*si*sj2*is4); \
    dXsj   = dp*( 2*sj*si2*is4); \
    dXxi   = sj2*is2; \
    dXxj   = si2*is2; \
} \

#define _Gauss_overlap( r2, si, sj ) \
    double sisj    = si*sj; \
    double inv_sisj= 1/sisj; \
    double g       = exp( -r2/(2*s2) );    \
    double S       =  (2*M_SQRT2) * g *  sisj*sisj*is2 * sqrt( inv_sisj*is2 ) ;  \
    double dS_dr   = -S * is2; \
    double inv_si  = sj*inv_sisj; \
    double inv_sj  = si*inv_sisj; \
    double S_s4    = S * is4; \
    double dS_dsi  = S_s4 * ( si2*r2 + 3*sj2*s2 ) * inv_si; \
    double dS_dsj  = S_s4 * ( sj2*r2 + 3*si2*s2 ) * inv_sj; \
    dS_dsi        -= 1.5*S * inv_si;    \
    dS_dsj        -= 1.5*S * inv_sj;    \

#define _Gauss_tau( r2, si, sj ) \
    double tau         = -(r2 - 3*s2)*is4; \
    double dTau_dr     = -2*is4; \
    double dTau_ds     =  2*( 2*r2 - 3*s2 )*is4*is2; \
    double dTau_dsi    =  si*dTau_ds; \
    double dTau_dsj    =  sj*dTau_ds; \

#define _Gauss_kinetic( r2, si, sj ) \
    double S_     = S  *const_K_eVA; \
    double tau_   = tau*const_K_eVA; \
    double T      = S_*tau; \
    double dT_dr  = S_*dTau_dr  + tau_*dS_dr;  \
    double dT_dsi = S_*dTau_dsi + tau_*dS_dsi; \
    double dT_dsj = S_*dTau_dsj + tau_*dS_dsj; \

/*
######  Product of Gaussians 1D

gij(x) =  exp( wi*(x-xi)^2 ) * exp( wj*(x-xj)^2 )

wi*( x^2 + 2*xi*x + xi^2 ) +  wj*( x^2 + 2*xj*x + xj^2 )
(wi+wj)*x^2  +   (wi*xi + wj*xj)*2*x   +   wi*x^2 + wj*x^2

wij*( x^2   +  2*x*xij +   xij^2 )   +    C

wij = (wi+wj)
wij*xij =   (wi*xi + wj*xj)

xij =    (wi*xi + wj*xj) / (wi + wj)

C =   wi*xi^2 + wj*xj^2     -    wij*xij^2
C =   wi*xi^2 + wj*xj^2     -   (wi*xi + wj*xj)^2 / (wi + wj)
C = wi*wj*(xi**2 - 2*xi*xj + xj**2)/(wi + wj)
C = wi*wj*(xi-xj)**2/(wi+wj)

Derived & Tested in Python Here:
/home/prokop/git/SimpleSimulationEngine/python/pyGaussAtom/GaussProduct.py
*/


inline double product1D_w( double wi,double xi,   double wj,double xj,   double& W, double& X  ){
    W           =  wi+wj;
    double wxi  =  wi*xi;
    double wxj  =  wj*xj;
    double wx   =  wxi + wxj;
    X           =  wx/W;
    double logC =  wxi*xi + wxj*xj - wx*X;
    //double C   = np.exp(-logC) * Ci * Cj
    return logC;
}

inline double product1D_s( double si,double xi,   double sj,double xj,   double& S, double& X  ){
    double si2   = si*si;
    double sj2   = sj*sj;
    double inv     = 1/(si2 + sj2);
    double sqrtinv = sqrt(inv);
    S          =  si*sj*sqrtinv;
    X          =  ( si2*xj + sj2*xi )*inv;
    //double logC =  wxi*xi + wxj*xj - wx*X;
    //double C   = np.exp(-logC) * Ci * Cj
    return 1;
}

inline double product1D_s_deriv( double si, double r2, double xi,   double sj,double xj,   double& S, double& X,    double& dSsi,double& dSsj,  double& dXsi,double& dXsj,  double& dXxi,double& dXxj ){
    double si2   = si*si;
    double sj2   = sj*sj;
    double s2   = si2 + sj2;
    double is2     = 1/(si2 + sj2);
    double sqrtis2 = sqrt(is2);
    S      =  si*sj*sqrtis2;
    X      =  ( si2*xj + sj2*xi )*is2;

    dSsi   = sj*sj2*sqrtis2*is2;
    dSsj   = si*si2*sqrtis2*is2;

    dXsi   = 4*si*xj;
    dXsj   = 4*sj*xi;

    dXxi   = 2*sj2;
    dXxj   = 2*si2;

   return 1;
}


inline double Coulomb( double r, double s, double& fr, double& fs ){
    // "s" is considered from rho ( density blobl; not wave function blob )
    // ToDo: maybe we can do without s=sqrt(s2) and r=sqrt(r2)
    //constexpr const double const_F2 = -2.*sqrt(2./M_PI);
    constexpr const double const_F2 = M_2_SQRTPI * M_SQRT2;
    //const double const_F2 = M_2_SQRTPI;
    // NOTE : there cannot be any non-constant prefactors because at large distance ( r -> +inf) it must converge to point-charge coulomb (1/r)*qi*qj 
    double ir   = 1./r; //(r+1.e-8);
    double is   = 1./s; //(s+1.e-8);
    double r_s  = r*is;
    double r_2s = M_SQRT1_2 * r_s; // This is for charge-density blobs (assuming si,sj comes from charge denisty)
    //double r_2s = r_s;
    //double r_2s = M_SQRT2   * r_s; // This is for wavefunction blobs (assuming si,sj comes from wavefunction)
    double e1   = ir * const_El_eVA;
    double e2   = erf(  r_2s      );            // ToDo : this should be possible to compute together !!!
    double g    = exp( -r_2s*r_2s ) * const_F2;
    double f1   = -e1*ir;
    double f2   = g*is*0.5;
    double e1f2 = e1*f2;
    fr          = (f1*e2 + e1f2)*ir      ;
    //printf( "r %g fr %g = (f1 %g * e2 %g )+(e1 %g *f2 %g) r_2s %g r %g s %g\n", r, fr, f1, e2, e1, f2, r_2s, r, s );
    fs          =          e1f2 *r_s * is;
    return e1 * e2;
}


inline double Coulomb( const Vec3d& Rij, double r2, double si, double sj, double qij, Vec3d& fp, double& fs ){
    constexpr const double R2SAFE = 1.0e-8;
    double r    = sqrt(r2 + R2SAFE);
    double s2   = si*si + sj*sj;
    double s    = sqrt(s2);
    double fr;
    // NOTE : there cannot be any non-constant prefactors because at large distance ( r -> +inf) it must converge to point-charge coulomb (1/r)*qi*qj 
    double e  = qij*Gauss::Coulomb( r, s, fr, fs ); // NOTE : remove s*2 ... hope it is fine ?
    //double e    = const_El_eVA*qij/sqrt( r*r - 0.01 );    // assymptotic limit for ( r -> +inf)
    //printf( "qij %g | s (%g,%g) r %g E %g  \n", qij, si, sj, r, e );
    //if(fabs(qij)>1e-16)printf( "Gauss::Coulomb() r %g s %g(%g,%g) q %g -> E %g fr %g fs %g \n", r, s,si,sj,  e, fr, fs  );
    // --- Derivatives (Forces)
    fp  = Rij*( fr * qij );   // use:   rhofP[i].add(fp);    rhofP[j].sub(fp);
    fs *=            qij;     // use:   rhofS[i] -= fs*si;   rhofS[j] -= fs*sj;
    return e;
}

inline double overlap( double r, double si, double sj,   double& fr, double& fsi, double& fsj ){
    // NOTE : this is for normalized gaussian
    // see pyMolecular/eFF_KineticAndOverlap.py
    const double C = 2*M_SQRT2; // 2^(3/2)

    _Gauss_sij_aux( si, sj )

    double r2    = r*r;

    double sij     = si*sj;
    double inv_sij = 1/sij;
    double inv_sqrt_sij = sqrt(inv_sij);
    double norm  = inv_sij*inv_sqrt_sij;     //  (si*sj)**(1.5)
    double g     = exp( -r2/(2*s2) );
    //double S   = const * norm * sij^3/( s2^1.5 ) * g
    double S     =  C * g *  sij*is2 * sqrt( inv_sij*is2 );
    fr           = -S / s2;
    fsi          =  S *  ( si2*r2 + 3*sj2*s2 )/(  si*s2*s2 );
    fsj          =  S *  ( sj2*r2 + 3*si2*s2 )/(  sj*s2*s2 );
    // corection by derivative of normalization
    fsi          = fsi - 1.5*S/si;
    fsj          = fsj - 1.5*S/sj;
    return S;
}

inline double tau_func( double r, double si, double sj, double& fr, double& fsi, double& fsj ){ 
    // Kinetic/Overlap Tij/Sij
    // // see pyMolecular/eFF_KineticAndOverlap.py
    double si2  = si*si;
    double sj2  = sj*sj;
    double s2   = si2 + sj2;
    double r2   = r*r;
    double is2  = 1/s2;
    double is4  = is2*is2; 

    double tau  = -(r2 - 3*s2)*is2*is2;
    fr          = -2*r*is4;
    double fs   =  2*si*( 2*r2 - 3*s2 )*is4*is2;
    fsi         =  fs*si;
    fsj         =  fs*sj;
    return tau;
}

inline double product3D_s_deriv(
    double si,    Vec3d   pi,
    double sj,    Vec3d   pj,
    double& sij,    Vec3d & pij,
    double& dSsi, double& dSsj,
    Vec3d & dXsi, Vec3d & dXsj,
    double& dXxi, double& dXxj,
    double& dCsi, double& dCsj, double& dCr
){
    double C;
    _Gauss_sij_aux( si, sj )
    Vec3d dp    = pi-pj;
    double r2   = dp.norm2();
    _Gauss_product(pi,pj,si,sj)
    _Gauss_product_derivs(dSsi,dSsj,dXsi,dXsj,dXxi,dXxj)
    _Gauss_overlap( r2, si, sj )
    C    = S; 
    dCsi = dS_dsi;
    dCsj = dS_dsj;
    dCr  = dS_dr;
    return C;
}

inline double overlap_s_deriv( // derived/simplified from product3D_s_deriv(
    double si,    Vec3d   pi,
    double sj,    Vec3d   pj,
    double& dCr, double& dCsi, double& dCsj
){
    double C;
    _Gauss_sij_aux( si, sj )
    Vec3d dp    = pi-pj;
    double r2   = dp.norm2();
    _Gauss_overlap( r2, si, sj )
    C    = S; 
    dCsi = dS_dsi;
    dCsj = dS_dsj;
    dCr  = dS_dr;
    return C;
}

/// ====== Product Blobs and its derivativs

struct PairInt{
    Vec3d  p;
    double si,sj;
    //double ci,cj;
    double S;

    //double fromDerivsT( double r2, double si, double sj ){
    //    _Gauss_overlap()
    //}
    //double fromDerivsS(){}
    inline void set( const Vec3d& p_, double si_, double sj_, double S_ ){ 
        p=p_; si=si_; sj=sj_; S=S_; 
        //printf( "PairInt::set() p(%g,%g,%g) s(%g,%g) S %g\n", p.x,p.y,p.z, si,sj,S  ); 
    }

    inline void applyForceScaled( double K, Vec3d& fpi, Vec3d& fpj, double& fsi, double& fsj )const{
        //Vec3d fp = K*p;
        fpi.add_mul( p,  K );
        fpj.add_mul( p, -K );
        fsi += si*K;
        fsj += sj*K;
    }
};

struct Blob{
    Vec3d  p;
    double s;
    double c;

    inline void setZero(){ p=Vec3dZero; s=0; c=0; };

    inline void set( const Vec3d& pos_, double size_, double charge_ ){ p=pos_; s=size_; c=charge_; }
    inline void add( const Vec3d& pos_, double size_, double charge_ ){ p.add(pos_); s+=size_; c+=charge_; }

    inline void applyForceScaled( double K, Vec3d& fp, double& fs, double& c_ )const{
        fp.add_mul( p,  K );
        fs +=       s * K  ;
        c_ +=       c * K  ;
    }
};

struct PairDeriv{
    double dSsi,dSsj;
    Vec3d  dXsi,dXsj;
    double dXxi,dXxj;
    double dCsi,dCsj,dCr;
};

inline double product3DDeriv( double si, Vec3d pi, double sj, Vec3d pj, Blob& out, PairDeriv& dOut ){
    double C = product3D_s_deriv( si, pi, sj, pj,
        out.s, out.p,
        dOut.dSsi, dOut.dSsj,
        dOut.dXsi, dOut.dXsj,
        dOut.dXxi, dOut.dXxj,
        dOut.dCsi, dOut.dCsj, dOut.dCr
    );
    out.c = C;
    return C;
}

inline void productBackForce(
    const Blob& fB, const PairDeriv& Ds,
    Vec3d Rij, double cij,
    Vec3d& Fxi, Vec3d& Fxj, double& fsi, double& fsj
){
    //if(DEBUG_iter==DEBUG_log_iter) printf( "cij %g dEdQ %g Fx %g \n", cij, dEdQ, Fp.x );
    fsi += ( fB.p.dot( Ds.dXsi ) + fB.s*Ds.dSsi + fB.c*Ds.dCsi*cij );
    fsj += ( fB.p.dot( Ds.dXsj ) + fB.s*Ds.dSsj + fB.c*Ds.dCsj*cij );
    Vec3d  Fq = Rij*(Ds.dCr*cij)*fB.c;
    Fxi.add( fB.p*Ds.dXxi + Fq );
    Fxj.add( fB.p*Ds.dXxj + Fq );
}




/// ====== this should be probably removed later

struct PairDerivs{
    double C;
    double s;
    Vec3d  p;
    double dSsi,dSsj;
    Vec3d  dXsi,dXsj;
    double dXxi,dXxj;
    double dCsi,dCsj,dCr;

    inline void get( double si, Vec3d pi,  double sj, Vec3d pj ){
        C = product3D_s_deriv( si, pi, sj, pj,
            s, p,
            dSsi, dSsj,
            dXsi, dXsj,
            dXxi, dXxj,
            dCsi, dCsj, dCr
        );
    }

    inline void backForce(
        Vec3d Rij, double cij, double dEdQ, Vec3d Fp, double Fs,
        Vec3d& Fxi, Vec3d& Fxj, double& fsi, double& fsj
    ){
        //if(DEBUG_iter==DEBUG_log_iter) printf( "cij %g dEdQ %g Fx %g \n", cij, dEdQ, Fp.x );
        fsi += ( Fp.dot( dXsi ) + Fs*dSsi + dEdQ*dCsi*cij );
        fsj += ( Fp.dot( dXsj ) + Fs*dSsj + dEdQ*dCsj*cij );
        Vec3d  Fq = Rij*(dCr*cij)*dEdQ;
        Fxi.add( Fp*dXxi + Fq );
        Fxj.add( Fp*dXxj + Fq );
    }

};






inline double product3D_s_new(
    double si,    Vec3d   pi,
    double sj,    Vec3d   pj,
    double& S,    Vec3d & p
){
    double si2   = si*si;
    double sj2   = sj*sj;
    double s2    = si2 + sj2;
    double is2   = 1/s2;
    double sqrtis2 = sqrt(is2);

    S      =  si*sj*sqrtis2;
    p      =  pj*(si2*is2) + pi*(sj2*is2);
    //X      =  ( si2*xj + sj2*xi )*inv;

    Vec3d dp = pi-pj;

    // --- constant

    double r2 = dp.norm2();

    double a2   = (si*sj)*is2;
    double a    = sqrt(a2);
    double e1   = a2*a*M_SQRT2;
    double e2   = exp( -r2*(is2*0.5) );

    // From  DensOverlapGauss_S() @  InteractionsGauss.h
    //double a    = 2.*(si*sj)*is2;
    //double a2   = a*a;
    //double Aa2  = amp*a2;
    //double e1   = Aa2*a;              // (2/(si/sj+si/sj))^3
    //double e2   = exp( -2*r2*is2 );   // exp( -2*r^2/(si^2+sj^2) )
    //printf( "r %g e1 %g e2 %g a2 %g is2 %g | si %g sj %g \n", sqrt(r2), e1, e2, a2, is2, si, sj );

    return e1 * e2;
}



inline double product1D_w_deriv( double wi,double xi,   double wj,double xj,   double& W, double& X,    double& dXdwi, double& dXdwj, double& dXdxi, double& dXdxj ){
    // derivatives in  /home/prokop/Dropbox/MyDevSW/Maxima/Gauss_Product.wxmx
    W          = wi+wj; // => dW/dwi = 1   ;    dW/dwj = 1   dW/dxi = 0   ;    dW/dxj = 0
    double wxi = wi*xi;
    double wxj = wj*xj;
    double wx  = wxi + wxj;
    X          = wx/W;
    dXdwi = xj/(wi*wi);
    dXdwj = xi/(wj*wj);
    dXdxi = 1/wj;
    dXdxj = 1/wi;
    //   X =  (wi*xi + wj*xj)/(wi*wj)
    double logC =  wxi*xi + wxj*xj - wx*X;
    //double C   = np.exp(-logC) * Ci * Cj
    return logC;
}



inline double product3D_w( double wi, const Vec3d& pi, double wj, const Vec3d& pj,  double& wij, Vec3d& pij ){
    double junk;
    double logC =
        + product1D_w( wi,pi.x,  wj,pj.x,   wij ,  pij.x );
        + product1D_w( wi,pi.y,  wj,pj.y,   junk,  pij.y );
        + product1D_w( wi,pi.z,  wj,pj.z,   junk,  pij.z );

    /*
    Vec3d wxi   =  pi*wi;
    Vec3d wxj   =  pj*wj;
    double W    =  wi+wj;
    double wx   =  wxi + wxj;
    double X    =  wx/W;
    double logC =  wxi*xi + wxj*xj - wx*X;
    */
    return logC;
}

// ToDo : needs derivatives of projection
inline double product3D_s( double si, const Vec3d& pi, double sj, const Vec3d& pj,  double& sij, Vec3d& pij ){
    double junk;
    double wi = 1/(2*si*si);
    double wj = 1/(2*sj*sj);
    double wij;
    double logC =
        + product1D_w( wi,pi.x,  wj,pj.x,   wij ,  pij.x );
        + product1D_w( wi,pi.y,  wj,pj.y,   junk,  pij.y );
        + product1D_w( wi,pi.z,  wj,pj.z,   junk,  pij.z );
    sij = 1/sqrt(2*wij);
    //printf( "DEBUG product3D_s sij %g(%g,%g) wij %g(%g,%g)\n", sij,si,sj,  wij,wi,wj );
    return logC;
}


// ToDo : needs derivatives of projection
inline double product3D_s_( double si, const Vec3d& pi, double sj, const Vec3d& pj,        double& sij, Vec3d& pij ){
    double junk;
    double wi = 1/(2*si*si);
    double wj = 1/(2*sj*sj);
    double wij;
    double logC =
        + product1D_s( wi,pi.x,  wj,pj.x,   wij ,  pij.x );
        + product1D_s( wi,pi.y,  wj,pj.y,   junk,  pij.y );
        + product1D_s( wi,pi.z,  wj,pj.z,   junk,  pij.z );
    sij = 1/sqrt(2*wij);
    return logC;
}


inline double norm3Dw( double w ){
    // https://math.stackexchange.com/questions/434629/3-d-generalization-of-the-gaussian-point-spread-function
    // n-dimansional Gaussian normalization C is :
    // G(x,y,z) = C*exp( -( x^2 + y^2 + z^2)/(2*sigma^2) ) =   C*exp( -w*( x^2 + y^2 + z^2) )
    // C  = 1/( sigma^n  (2*pi)^(n/2)  )
    // C3 = 1/( sigma^3 (2*pi)^(3/2) )
    // 1/sigma^2 = w
    // C3 = 1/( sqrt(2*sigma^2)^3/((2)^(3/2)) (2*pi)^(3/2) )
    // C3 = 1/( w^(3/2)/((2)^(3/2)) (2*pi)^(3/2) )
    // C3 = 1/( w*pi )^(3/2)
    double c2 = w/M_PI;
    double c  = sqrt( c2 );
    return c2*c;
}

inline double norm3Ds( double s ){
    const double sqrt_2pi = 2.50662827463;
    double c = 1/(sqrt_2pi*s);
    return c*c*c;
}

inline double sqnorm3Ds( double s ){
    const double sqrt_pi = 1.77245385091;
    // exp(-r2/(2*s^2))*exp(-r2/(2*s^2)) = exp(-r2/(s^2))
    // exp(-r2/(s^2)) = s*sqrt(pi)
    // ToDo: it may be better to store      1/sqrt(s) ... to avoid unnecessary sqrt() call
    double c = sqrt( 1/(sqrt_pi*s) ); // sqrt(sqrt(pi)*s)
    return c*c*c;
}

inline double sqnorm3Ds_deriv( double s, double& d ){
    const double C  =  0.42377720812; // pi^(-3/4)
    const double dC = -0.63566581218; // pi^(-3/4) * (3/2)
    double invs   = 1/s;
    double insqrt  = sqrt(invs);
    double N = invs*insqrt;   //   ()^(-3/2)
    d        = dC * N * invs; //   ()^(-5/2)
    return      C * N;
}

inline double sqnorm3Ds_sq( double s ){
    const double sqrt_pi = 1.77245385091;
    double c = 1/(sqrt_pi*s); // sqrt(sqrt(pi)*s)
    return c*c*c;
}

inline double  kinetic_r0   (double s){
    return 1.5*const_K_eVA/(s*s);
}

inline double  kinetic_r0_derivs(double s, double& dT_ds){
    double is = 1/s;
    double T = 1.5*const_K_eVA*is*is;
    dT_ds    = -T*is;
    return T;
}

//inline double uy3Ds( double r, double s ){ return norm3Ds(s)* exp( (r*r)/(-2*s*s) ); }
//inline double uy3Ds2( double r2, double s ){ return norm3Ds(s)* exp( r2/(-2*s*s) ); }
inline double bas3D_r2( double r2, double s ){ return sqnorm3Ds(s)* exp( r2/(-2*s*s) ); }
inline double rho3D_r2( double r2, double s ){ return   norm3Ds(s)* exp( r2/(-2*s*s) ); }

inline double LoG3D_r2( double r2, double s ){ double wr2=r2/(-2*s*s); return   norm3Ds(s)* ( (1+wr2) )*exp( wr2 ); }



inline double kinetic_w( double w ){
    // !!!!! NOT NORMALIZED GAUSSIAN
    //    -((3*pi^(3/2))/(2^(3/2) )   /sqrt(w))
    return -5.90610372965 * sqrt(w);
}


inline double kinetic_w(  double r2, double w1, double w2 ){
    // Derived in Maxima: see : /home/prokop/Dropbox/MyDevSW/Maxima/Gauss_Kinetic-Cleaned.wxmx
    //                                                              Gauss_Kinetic-Polar-.wxmx
    // Cartes : (2*%pi^(3/2)*w1*w2*(2*w1*w2*x1^2-3*w2-3*w1)*%e^(-(w1*w2*x1^2)/(w2+w1)))/(w2+w1)^(7/2)
    // Polar  : (2*%pi^(3/2)*w1*w2*(2*w1*w2*x1^2-3*w2-3*w1)*%e^(-(w1*w2*x1^2)/(w2+w1)))/(w2+w1)^(7/2)
    //   (2*pi^(3/2) *   w1*w2*    (2*w1*w2* x1^2 -3*(w2+w1) )*  exp( -(w1*w2*x1^2)/(w2+w1))   )   /(w2+w1)^(7/2)
    //                   wprod*    (2 *wprod* x*x - 3*wsum )  *  exp( -( wprod *x*x )/wsum  )
    double C         = 11.1366559937; //(pi^(3/2))*2
    double wsum      = w2+w1 ;
    double iwsum     = 1/wsum;
    double sqrtiwsum = sqrt(iwsum);
    double wprod     = w1*w2 ;
    //double x2        = x*x;
    double g         = exp( -( wprod * r2 )*iwsum );
    return C*wprod * ( 2*wprod*r2 - 3*wsum )*g*iwsum*iwsum*iwsum*sqrtiwsum;
}

inline double kinetic_s(  double r2, double si, double sj,   double& fr, double& fsi, double& fsj ){

 //  This is total kinetic energy ( not kinetic energy change upon renormalization, not normalized by overlap, not subtracted )
 // ToDo : Kinetic and Overlap share much of calculations => make sense to calculate them together in one function

 // Look also here  http://dx.doi.org/10.1016/j.mechmat.2015.02.008
 // (2^(3/2)*%pi^(3/2)*s1^3*s2^3*(  r2   -3*(s2^2+s1^2)*exp( -r2/(2*s2^2+2*s1^2) ) )/(s2^2+s1^2)^(7/2)
 //  2^(3/2)*pi^(3/2)    *si^3*sj^3*(  r2   -3*(sj^2+si^2)*exp(-r2/(2*sj^2+2*si^2)))     /     (s2^2+s1^2)^(7/2)
 //  2^(3/2)*pi^(3/2)    *sij^3*(  r2   -3* s2 *exp( -r2/(2*s2) ) ) /  s2^(7/2)
 const double C  = -15.7496099457 * const_K_eVA;  //    2^(3/2)*pi^(3/2)   // NOTE:  factor 2x was here for wrong reason
 double dCi,dCj;
 //C *= sqnorm3Ds(si) * sqnorm3Ds(sj);
 double Ci = sqnorm3Ds_deriv( si, dCi );
 double Cj = sqnorm3Ds_deriv( sj, dCj );
 double Cij = Ci*Cj;

 double sij   = si*sj;
 double si2   = si*si;
 double sj2   = sj*sj;
 double si4   = si2*si2;
 double sj4   = sj2*sj2;
 double si6   = si2*si4;
 double sj6   = sj2*sj4;
 double sij2  = sij*sij;
 double s2    = si2 + sj2;
 double invs2 = 1/s2;
 double invs4 = invs2*invs2;
 double invs  = sqrt( invs2 );
 double g     = exp( -r2*0.5*invs2 );
 double denom = invs4*invs2*invs;
 double comm  = C*sij2 * g * denom;
 //double poly     = sij*sij*sij  * ( r2  - 3*s2 );
 //double dpoly_si = 3*sij*sij*sj * ( r2  - 3*sj2 - 5*si2 );
 //double dpoly_sj = 3*sij*sij*sj * ( r2  - 3*si2 - 5*sj2 );

 //double T = const_K_eVA * ( 1.5*s2*isi2*isj2 -2.*( 3.*s2 - 2.*r2 )*is4 );
 //double T = const_K_eVA * ( -2.*( 3.*s2 - 2.*r2 )*is4 );

 double E = comm * sij * ( r2  -  3*s2 );
 comm*=invs2;
 fr       = comm * sij * ( r2  -  0.5*5*s2 );   // factor 0.5 here hellps but not fully
 comm*=invs2;
 fsi  = comm * sj  * (  si2*r2*r2  + ( 3*sj4  - 4*sij2 - 7*si4 )*r2 - 9*sj6 + (-12*sj2 + 3*si2)*sij2 + 6*si6 );
 fsj  = comm * si  * (  sj2*r2*r2  + ( 3*si4  - 4*sij2 - 7*sj4 )*r2 - 9*si6 + (-12*si2 + 3*sj2)*sij2 + 6*sj6 );

 //printf( "E %g fr %g fsi %g fsj %g \n", E, fr, fsi, fsj );

//double k=0.5;
//dCi*=k;
//dCj*=k;
 
 fsi = Cij*fsi + Cj*dCi*E;
 fsj = Cij*fsj + Ci*dCj*E;
 E *=Cij;
 fr*=Cij;

 //printf( "kinetic_s T %g r %g si %g sj %g ", E, sqrt(r2), si, sj );
 //printf( "Gauss::kinetic_s() E %g r %g s%g(%g,%g) \n", E, sqrt(r2), sij, si, sj );

 //double fsi  = C* (  *poly*g*denom +  ddenom
 //( s1 ^2*x1^4+3*s2^4*x1^2-4*s1^2*s2^2*x1^2-7*s1^4*x1^2-9*s2^6-12*s1^2*s2^4+3*s1^4*s2^2+6*s1^6)
 // si2*r2*r2 +    3*sj^4 * r2  -4*si^2*sj^2*r2  - 7*si^4*r2 -    9*sj^6 - 12*si^2*sj^4 + 3*si^4*sj^2 + 6*si^6  )
 // si2*r2*r2 +  ( 3*sj^4  -4*si^2*sj^2  - 7*si^4 )*r2  -    9*sj^6 - 12*si^2*sj^4 + 3*si^4*sj^2 + 6*si^6  )
 // si2*r2*r2 + ( 3*sj2*sj2  -4*si2*sj2  - 7*si2*si2 )*r2  -  9*sj2*sj2*sj2 - 12*si2*sj2*sj2 + 3*si2*si2*sj2 + si2*si2*si2  )
 // si2*r2*r2  + ( 3*sj4  -4*si4  - 7*si2*si2 )*r2 - 9*sj4*sj2 - 12*si2*sj4 + 3*si4*sj2 + si2*si4

 return   E;
}

/*
inline double kinetic_s(  double r2, double si, double& fsi ){
 //  from  Eq.   [Su, J. T. & Goddard, W. A. Excited Electron Dynamics Modeling of Warm Dense Matter. Phys. Rev. Lett. 99, 185003 (2007).]
 // (3/2)  * ( 1/si + 1/sj ) + 2*( 2*rij - 3*(si^2 + sj^2)  )/(si^2 + sj^2)^2
 // (3/2)  * ( 2/si ) - 2*(2*3*si^2  )/4*(si^2)^2
 // 3/si^2 - 3/(si^2)
 return   E;
}
*/


inline double kinetic( double s ){
    //     -(3*pi^(3/2)*s)/2
    //return -8.35249199525 * s * sqnorm3Ds(s) * sqnorm3Ds(s);
    return -8.35249199525 * s * sqnorm3Ds_sq(s);
}

inline double kinetic( double r2, double s1, double s2 ){
    return kinetic_w( r2, 1/(2*s1*s1), 1/(2*s2*s2) )  *  sqnorm3Ds(s1) * sqnorm3Ds(s2);
}


inline double kinetic( double r2, double s1, double s2, double& dr, double& ds1, double& ds2 ){
    /// ToDo :  Need derivatives of Kinetic Overlap !!!!!
    return kinetic_w( r2, 1/(2*s1*s1), 1/(2*s2*s2) )  *  sqnorm3Ds(s1) * sqnorm3Ds(s2);
}


/// Boys Function
//  https://chemistry.stackexchange.com/questions/41214/boys-function-for-gaussian-integrals-in-ab-initio-calculations
//  BoysF is used for Coulomb integrals like :
//  Vpq = Integral{   exp(-p*(r1)^2) * exp(-q*(r2)^2)  /|r1-r2| }
//  Vpq = 2*pi^(5/2) / ( p*q*sqrt(p+q) ) * BoysF( a*r^2 )
//  F0() is related to error-function erf()
//  F0(x) = sqrt(pi/(4x)) * erf( sqrt(x) )
//  Vpq = (pi/sqrt(p*q))^3 * erf( sqrt(a) * r ) / r
// => It makes sense to tabulate rather F0(x^2)

}

#endif



