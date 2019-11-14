
#ifndef EFF_h
#define EFF_h

#include "fastmath.h"
//#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Forces.h"



/*

eFF : Electron Force Field
---------------------------

[1] http://aip.scitation.org/doi/10.1063/1.3272671
The dynamics of highly excited electronic systems: Applications of the electron force field
Julius T. Su, William A. Goddard

[2] http://dx.doi.org/10.1016/j.mechmat.2015.02.008
Non-adiabatic dynamics modeling framework for materials in extreme conditions
Hai Xiao, Andrés Jaramillo-Botero, Patrick L. Theofanis, William A. Goddard,

NOTES:
-------

1) It seems that decrease of kinetic energy by sharing electron between two atoms is dominant contribution which makes formation of bonds fabourable
    * H2 Molecule perhaps cannot be stable without this contribution ( i.e. with fixed radius of electron blobs )

*/




/*
Erf approximation:
# Gaussian:    F = (x2-1)**2 / sqrtPi
# Erf          E = x*(1 + x2 * ( -0.66666666666 + 0.2*x2 ) ) * (2/(16.0/15.0))
*/

// ============= Constnts

//double wae = 1.0;
//double bAE = -5.0;
//double aAE = 20.0;

//#define QE -2.0
#define QE -1.0

const double const_hbar_SI      = 1.054571817e-34;    // [J.s]  #6.582119569e-16 # [eV/s]
const double const_Me_SI        = 9.10938356e-31;     // [kg]
const double const_e_SI         = 1.602176620898e-19; // [Coulomb]
const double const_eps0_SI      = 8.854187812813e-12; // [F.m = Coulomb/(Volt*m)]
const double const_eV_SI        = 1.602176620898e-19; // [J]
const double const_Angstroem_SI = 1.0e-10;

const double const_K_SI     =  const_hbar_SI*const_hbar_SI/const_Me_SI;
const double const_El_SI    =  const_e_SI*const_e_SI/(4.*M_PI*const_eps0_SI);
const double const_Ry_SI    = 0.5 * const_El_SI*const_El_SI/const_K_SI;

const double const_Ry_eV  = 13.6056925944;
const double const_El_eVA = const_El_SI/( const_e_SI*const_Angstroem_SI );
const double const_K_eVA  = (const_El_eVA*const_El_eVA)/(2*const_Ry_eV);
const double const_Ke_eVA = const_K_eVA*1.5;

// ToDo : Later properly
constexpr static const Vec3d default_eAbWs[] = {
{ 0.0,  0.0, 0.0},  // Q = 0 //
{ 0.0,  0.0, 0.25},  // Q = 1 // H
{10.0, -3.0, 0.25},  // Q = 2 // Be?
{10.0, -3.0, 0.25},  // Q = 3 // B
{10.0, -3.0, 0.25},  // Q = 4 // C
};

constexpr static const Vec3d default_aAbWs[] = {
{ 0.0,  0.0, 0.1},  // Q = 0 //
{ 0.0,  0.0, 0.25},  // Q = 1 // H
{20.0, -5.0, 0.25},  // Q = 2 // Be?
{20.0, -5.0, 0.25},  // Q = 3 // B
{20.0, -5.0, 0.25},  // Q = 4 // C
};


//double bAE = -6.0;
//double aAE = 100.0;




struct EFFAtomType{
    double Q; // nuncear charge
    double wa;
    double Aa; // prefactor
    double Ba; // exponential beta
    double we;
    double Ae; // prefactor
    double Be; // exponential beta
};

inline double frEFF_r2( double r2, double w1, double w2 ){
    return ( 1/(1 + w1*r2) - 1/(1 + w2*r2) )/( 1/w1 + 1/w2 );
}

inline void combineAbW( const Vec3d& abwi, const Vec3d& abwj, Vec3d& abw ){
    abw.x = abwi.x*abwj.x;
    abw.y = (abwi.y+abwj.y)*0.5;
    abw.z = abwi.z+abwj.z;
};

/*
inline double fr_ee(double r, double w){
    return 1/( r*r + w*w );
};

inline double fr_ae(double r, double w, double Q){
    return exp(5*r) - Q/( r*r + w*w );
};

inline double fr_aa(double r, double s){
    return 1/( r*r+ w*w );
};
*/

/*
inline double addPairEF_ee( const Vec3d& d, Vec3d& f, double w ){
    return 1/( r*r + w*w );
};

inline double addPairEF_ae( const Vec3d& d, Vec3d& f, double w, double q, double aPuli, double beta ){
    return exp(5*r) - Q/( r*r + w*w );
};

inline double addPairEF_aa( const Vec3d& d, Vec3d& f, double w, double qq, ){
    return qq/( r*r+ w*w );
};
*/

inline double addPairEF_expQ( const Vec3d& d, Vec3d& f, double w2, double qq, double bExp, double aExp ){
    double r2     = d.norm2();
    double invrw2 = 1./( r2 + w2 );
    double invrw  = sqrt(invrw2);
    double E      =  qq*invrw;
    double fr     = -qq*invrw2*invrw;
    if( bExp<0 ){
        double r      = sqrt( r2+R2SAFE );
        double Eexp  = aExp*exp( bExp*r );
        fr          += bExp*Eexp/r;
        E           += Eexp;
    }
    f.set_mul( d, fr );
    return E;
}

inline double addKineticGauss( double s, double& fs ){
    double is  = 1/s;
    double is2 = is*is*const_Ke_eVA;
    fs += -2.*is2*is;
    return is2;
}

inline double CoulombGauss( double r, double s, double& fr, double& fs, double qq ){
    // ToDo: maybe we can do without s=sqrt(s2) and r=sqrt(r2)
    constexpr const double const_F2 = -2.*sqrt(2./M_PI);

    double ir   = 1./r;
    double is   = 1./s;
    double r_s  = r*is;
    double r_2s = M_SQRT2 * r_s;
    double e1   = qq*ir * const_El_eVA;
    double e2   = erf( r_2s );
    double g    = exp( -r_2s*r_2s ) * const_F2;
    double f1   = -e1*ir;
    double f2   = g*is;
    double e1f2 = e1*f2;
    fr = (f1*e2 + e1f2)*ir;
    fs =          e1f2 *r_s * is;

    //printf( "CoulombGauss : e1 %g \n", e1 );
    //printf( "CoulombGauss : e2 %g \n", e2 );
    //printf( "CoulombGauss : f1 %g \n", f1 );
    //printf( "CoulombGauss : f2 %g \n", f2 );
    //printf( "CoulombGauss : fr %g \n", fr );
    //printf( "CoulombGauss : fs %g \n", fs );

    return e1 * e2;
}

inline double addCoulombGauss( const Vec3d& dR, Vec3d& f, double si, double& fsi, double qq ){

    double r    = dR.norm();

    double fr,fs;
    double E = CoulombGauss( r, si, fr, fs, qq );

    fsi += fs*si;
    f.add_mul( dR, fr );
    return E;
}

inline double addCoulombGauss( const Vec3d& dR, Vec3d& f, double si, double sj, double& fsi, double& fsj, double qq ){

    double s2   = si*si + sj*sj;
    double s    = sqrt(s2);
    double r    = dR.norm();

    double fs,fr;
    double E = CoulombGauss( r, s, fr, fs, qq );

    fsi += fs*si;
    fsj += fs*sj;
    f.add_mul( dR, fr );
    return E;

}

/*
inline double addCoulombGauss( const Vec3d& dR, Vec3d& f, double qq, double si, double sj, double& fsi, double& fsj ){

    // ToDo: maybe we can do without s=sqrt(s2) and r=sqrt(r2)

    constexpr const double const_F2 = -2.*sqrt(2./M_PI);

    double s2   = si*si + sj*sj;
    double s    = sqrt(s2);
    double r    = dR.norm();
    double ir   = 1./r;
    double is   = 1./s;
    double r_s  = r*is;
    double r_2s = M_SQRT2 * r_s;
    double e1   = qq*ir * const_El_eVA;
    double e2   = erf( r_2s );
    double f1   = -e1*ir;
    double g    = exp( -r_2s*r_2s );
    double f2   = const_F2*g*is;
    double fr   = f1*e2 + f2;
    f.add_mul( dR, fr*ir );

    double fs  = f2*is;
    fsi += fs*(si*is);
    fsj += fs*(sj*is);

    return const_El_eVA * e1 * e2;
}
*/

/*
//inline double getDeltaTGauss( const Vec3d& dR, Vec3d& f, double si, double sj, double& fsi, double& fsj ){
inline double getDeltaTGauss( double r2, double si, double sj,  double& fr, double& fsi, double& fsj ){
    double si2  = si*si;
    double sj2  = sj*sj;
    double isi2 = 1./si2;
    double isj2 = 1./sj2;
    double s2   = si2 + sj2;
    double is2  = 1./s2;
    double is4  = is2 * is2;
    //double r2   = dR.norm2();

    double B   =  4.*( 3.*s2 - 4.*r2 )*is4*is2;
    fsi = const_K_eVA*( -3.*isi2*isi2 + B*si )*si;
    fsj = const_K_eVA*( -3.*isj2*isj2 + B*sj )*sj;
    //f.add_mul( dR, const_K_eVA * (8.*is4) );
    fr  = const_K_eVA*(8.*is4);

    return const_K_eVA * ( 1.5*s2*isi2*isj2 - 2.*( 3.*s2 - 2.*r2 )*is4 );
}


//inline double getOverlapSGauss( const Vec3d& dR, Vec3d& f, double si, double sj, double& fsi, double& fsj ){
inline double getOverlapSGauss( double r2, double si, double sj, double& fr, double& fsi, double& fsj ){
    double si2  = si*si;
    double sj2  = sj*sj;
    double isi2 = 1./si2;
    double isj2 = 1./sj2;
    double s2   = si2 + sj2;
    double is2  = 1./s2;
    double is4  = is2 * is2;
    //double r2   = dR.norm2();

    double a2   = 2.*(si*sj)*is2;
    double a    = sqrt(a2);
    double e1   = a2*a;
    double e2   = exp( -r2*is2 );

    double f1   = 3.*a  * (si2-sj2)*is4;
    double f2   = 2.*e2 * r2*is4;

    fsi = e1*f2*si - e2*f1*sj;
    fsj = e1*f2*sj + e2*f1*si;
    //f.add_mul( dR, e1*e2*(-2.*is2) );
    fr  = e1*e2*(-2.*is2);
    return e1 * e2;
}
*/


//inline double getDeltaTGauss( const Vec3d& dR, Vec3d& f, double si, double sj, double& fsi, double& fsj ){
inline double getDeltaTGauss( double r2, double si, double sj,  double& dTr, double& dTsi, double& dTsj,
    double isi2, double isj2, double s2, double is2, double is4
  ){

    double B   =  4.*( 3.*s2 - 4.*r2 )*is4*is2;
    dTsi = const_K_eVA*( -3.*isi2*isi2 + B )*si;
    dTsj = const_K_eVA*( -3.*isj2*isj2 + B )*sj;
    //f.add_mul( dR, const_K_eVA * (8.*is4) );
    dTr  = const_K_eVA*(8.*is4);

    //printf( "getDeltaTGauss: e1 %g \n", 1.5*s2*isi2*isj2 );
    //printf( "getDeltaTGauss: e2 %g \n", -2.*( 3.*s2 - 2.*r2 )*is4 );
    //printf( "getDeltaTGauss: T  %g \n", const_K_eVA * ( 1.5*s2*isi2*isj2 -2.*( 3.*s2 - 2.*r2 )*is4 ) );
    //printf( "getDeltaTGauss: B  %g \n", B    );
    //printf( "getDeltaTGauss: dTsi %g \n", dTsi );
    //printf( "getDeltaTGauss: dTsj %g \n", dTsj );
    //printf( "getDeltaTGauss: dTr  %g \n", dTr  );
    //printf( "getDeltaTGauss: fr   %g \n", dTr*sqrt(r2) );

    return const_K_eVA * ( 1.5*s2*isi2*isj2 -2.*( 3.*s2 - 2.*r2 )*is4 );
}


//inline double getOverlapSGauss( const Vec3d& dR, Vec3d& f, double si, double sj, double& fsi, double& fsj ){
inline double getOverlapSGauss( double r2, double si, double sj, double& dSr, double& dSsi, double& dSsj,
    double si2, double sj2, double is2, double is4
){

    double a2   = 2.*(si*sj)*is2;
    double a    = sqrt(a2);
    double e1   = a2*a;
    double e2   = exp( -r2*is2 );

    double f1   = 3.*a  * (si2-sj2)*is4;
    double f2   = 2.*e2 * r2*is4;

    dSsi = e1*f2*si - e2*f1*sj;
    dSsj = e1*f2*sj + e2*f1*si;
    //f.add_mul( dR, e1*e2*(-2.*is2) );
    dSr  = e1*e2*(-2.*is2);

    //printf( "getDeltaSGauss :e1  %g \n", e1 );
    //printf( "getDeltaSGauss :e2 %g \n", e2 );
    //printf( "getDeltaSGauss :f1 %g \n", f1 );
    //printf( "getDeltaSGauss :f2 %g \n", f2 );
    //printf( "getDeltaSGauss :S    %g \n", e1*e2 );
    //printf( "getDeltaSGauss :dSsi %g \n", dSsi  );
    //printf( "getDeltaSGauss :dSsj %g \n", dSsj  );
    //printf( "getDeltaSGauss :dSr  %g \n", dSr   );
    //printf( "getDeltaSGauss :fr  %g \n",  dSr*sqrt(r2) );

    return e1 * e2;
}


inline double getDeltaTGauss( double r2, double si, double sj,  double& dTr, double& dTsi, double& dTsj ){
    double si2  = si*si;
    double sj2  = sj*sj;
    //double isi2 = 1./si2;
    //double isj2 = 1./sj2;
    double s2   = si2 + sj2;
    double is2  = 1./s2;
    double is4  = is2 * is2;
    //printf( "getDeltaTGauss: r2, si2, sj2, s2, s4 %g %g %g %g %g \n", r2, si2, sj2, s2, is4 );
    return getDeltaTGauss  ( r2, si, sj, dTr, dTsi, dTsj, 1./si2,1./sj2,s2,is2,is4 );
}

//inline double getOverlapSGauss( const Vec3d& dR, Vec3d& f, double si, double sj, double& fsi, double& fsj ){
inline double getOverlapSGauss( double r2, double si, double sj, double& dSr, double& dSsi, double& dSsj ){
    double si2  = si*si;
    double sj2  = sj*sj;
    //double isi2 = 1./si2;
    //double isj2 = 1./sj2;
    double s2   = si2 + sj2;
    double is2  = 1./s2;
    double is4  = is2 * is2;
    //printf( "getDeltaSGauss: r2, si2, sj2, s2, s4 %g %g %g %g %g \n", r2, si2, sj2, s2, is4 );
    return getOverlapSGauss( r2, si, sj, dSr, dSsi, dSsj, si2,sj2,is2,is4 );
}

inline double PauliSGauss_anti( double S, double& fS, double rho ){
    double S2    = S*S;
    double D     = 1./(1.+S2);
    double SDrho = rho*S*D;
    fS  = 2.*D*SDrho;
    return   S*SDrho;
}

inline double PauliSGauss_syn( double S, double& fS, double rho ){
    double S2   = S*S;
    double D    = 1./(1.+S2);
    double Dm   = 1./(1.-S2);
    //double S2 = S*S;
    double rhom = 1-rho;

    double SDm    = S*Dm;
    double SDrhom = S*D*rhom;

    fS   = 2.*( SDm*Dm  + SDrhom*D );
    return  S*( SDm     + SDrhom   );

}

inline double addPauliGauss( const Vec3d& dR, Vec3d& f, double si, double sj, double& fsi, double& fsj, bool anti, const Vec3d& KRSrho ){

    double r2 = dR.norm2();
    r2 *= KRSrho.x*KRSrho.x;
    si *= KRSrho.y;
    sj *= KRSrho.y;
    double r = sqrt(r2);

    double si2  = si*si;
    double sj2  = sj*sj;
    //double isi2 = 1./si2;
    //double isj2 = 1./sj2;
    double s2   = si2 + sj2;
    double is2  = 1./s2;
    double is4  = is2 * is2;

    double dTr,dTsi,dTsj;
    double dSr,dSsi,dSsj;
    double T = getDeltaTGauss  ( r2, si, sj, dTr, dTsi, dTsj, 1./si2, 1./sj2,s2,is2,is4 );
    double S = getOverlapSGauss( r2, si, sj, dSr, dSsi, dSsj,    si2,    sj2,is2,is4 );

    double eS,fS;
    if(anti){
        eS = PauliSGauss_anti( S, fS, KRSrho.z );
    }else{
        eS = PauliSGauss_anti( S, fS, KRSrho.z );
    }

    double TfS = T*fS;
    fsi +=         (dTsi*eS + TfS*dSsi)*KRSrho.y;
    fsj +=         (dTsj*eS + TfS*dSsj)*KRSrho.y;
    f.add_mul( dR, (dTr *eS + TfS*dSr )*KRSrho.x*KRSrho.x ); // second *KRSrho.x because dR is not multiplied

    //printf( "addPauliGauss r2, si, sj %g %g %g \n", r2, si, sj );
    //printf( "addPauliGauss T,dTr,dTsi,dTsj %g %g %g %g \n", T,dTr,dTsi,dTsj );
    //printf( "addPauliGauss S,dSr,dSsi,dSs  %g %g %g %g \n", S,dSr,dSsi,dSsj );
    //printf( "addPauliGauss S, eS, fS  %g %g %g \n", S, eS, fS  );
    //printf( "addPauliGauss fr1, fr2 fr %g %g %g %g \n", dTr*eS, TfS*dSr, (dTr *eS + TfS*dSr )*KRSrho.x, dR.z );
    //printf( "---------------------- \n" );

    return T * eS;

}


inline double interp_gx4(double r2, double y1, double y2 ){
    double c = (1-r2);
    c*=c;
    return c*y1 + (1-c)*y2;
}

class EFF{ public:
    //double dvmax = 0.1;
    //double dpmax = 0.1;
    /*
    double emass = 1.0;
    double see   = 0.5;
    double saa   = 0.1;
    double sa2  = saa*saa;
    double se2  = see*see;
    double sea2 = see*see;
    */

    constexpr static const Vec3d KRSrho = { 1.125, 0.9, 0.2 };
    //Vec3d KRSrho = { 1.125, 0.9, 0.2 };

    //double wee = 2.0;
    //double wae = 1.0;
    //double wee = 0.5;
    //double wee = 0.5;
    //double wae = 0.25;
    //double waa = 0.25;

    //double bEE     = -1.0;
    //double aEE     =  2.0;
    //double bEEpair = -1.0;
    //double aEEpair =  0.1;

    //double bAE = -6.0;
    //double aAE = 100.0;

    int ne=0,na=0,nDOFs=0;
    //int*   atype  =0;

    double * aQ     =0;
    Vec3d  * apos   =0;
    Vec3d  * aforce =0;

    Vec3d  * aAbWs  =0;
    Vec3d  * eAbWs  =0;

    //double * espin  =0;
    int8_t * espin  =0;
    Vec3d  * epos   =0;
    Vec3d  * eforce =0;
    double * esize  =0;
    double * fsize  =0;

    double* pDOFs =0;
    double* fDOFs =0;

    double Ek=0, Eee=0, EeePaul=0, Eae=0, Eaa=0;

void realloc(int na_, int ne_){
    na=na_; ne=ne_;
    nDOFs=na*3+ne*4;
    _realloc( pDOFs,   nDOFs);
    _realloc( fDOFs,   nDOFs);

    _realloc( aQ  ,   na);
    _realloc( aAbWs,  na);
    _realloc( eAbWs,  na);
    //_realloc( apos,   na);
    //_realloc( aforce, na);
    //_realloc( avel,   na);
    //_realloc( aorbs,  na);

    _realloc( espin, ne);
    //_realloc( epos   ,ne);
    //_realloc( eforce ,ne);
    //_realloc( evel   ,ne);

    apos   = (Vec3d*)pDOFs;
    aforce = (Vec3d*)fDOFs;

    epos   = (Vec3d*)(pDOFs + na*3);
    eforce = (Vec3d*)(fDOFs + na*3);
    esize  =          fDOFs + na*3 + ne*3;
    fsize  =          fDOFs + na*3 + ne*3;

}

void clearForce(){
    for(int i=0; i<nDOFs; i++){ fDOFs[i]=0; };
    //for(int i=0; i<nDOFs; i++){ fDOFs[i]=Vec3dZero; };
}

double evalKinetic(){
    Ek=0;
    for(int i=0; i<ne; i++){
        Ek += addKineticGauss( esize[i], esize[i] );
    }
    return Ek;
}

double evalEE(){
    Eee    =0;
    EeePaul=0;
    //double w2ee = wee*wee;
    const double qq = QE*QE;
    for(int i=0; i<ne; i++){
        Vec3d  pi    = epos[i];
        double si    = esize[i];
        int8_t spini = espin[i];
        double& fsi  = esize[i];
        Vec3d&  fpi  = eforce[i];
        for(int j=0; j<i; j++){
            Vec3d f=Vec3dZero;
            Vec3d   dR  = epos [j] - pi;
            double  sj  = esize[j];
            double& fsj = esize[i];
            Eee     += addCoulombGauss( dR, f, si, sj, fsi, fsj, qq );
            EeePaul += addPauliGauss  ( dR, f, si, sj, fsi, fsj, spini!=espin[j], KRSrho );
            //Eee += addPairEF_expQ( epos[j]-pi, f, w2ee, +1, 0, 0 );
            eforce[j].sub(f);
            fpi      .add(f);
            //glColor3f(1.0,0.0,0.0);
            //Draw3D::drawVecInPos( f*-1, epos[j] );
            //Draw3D::drawVecInPos( f   , pi      );
        }
    }
    return Eee+EeePaul;
}

double evalAE(){
    Eae=0;
    //double see2 = see*see;
    //double saa2 = saa*saa;
    //double invSae = 1/( see*see + saa*saa );
    //double w2ae = wae*wae;
    for(int i=0; i<na; i++){
        Vec3d  pi   = apos[i];
        double qqi  = aQ[i]*QE;
        Vec3d  abwi = eAbWs[i];
        for(int j=0; j<ne; j++){
            Vec3d f=Vec3dZero;
            //Eae += addPairEF_expQ( epos[j]-pi, f, abwi.z, qi*QE, abwi.y, abwi.x );
            //printf(  "a[%i]e[%i] r %g\n", i, j, (epos[j]-pi).norm() );
            Vec3d   dR  = epos [j] - pi;
            double fs_junk;
            Eee               += addCoulombGauss( dR, f, esize[j], fsize[j], qqi );
            if(qqi<-1.00001) EeePaul += addPauliGauss  ( dR, f, esize[j], abwi.z, fsize[j], fs_junk, false, KRSrho );

            eforce[j].sub(f);
            aforce[i].add(f);
            //glColor3f(1.0,0.0,1.0);
            //Draw3D::drawVecInPos( f*-1, epos[j] );
            //Draw3D::drawVecInPos( f   , pi      );
        }
    }
    return Eae;
}

double evalAA(){
    if( i_DEBUG>0 ) printf( "evalAA \n" );
    Eaa=0;
    //double invSaa = 1/(saa*saa);
    //double w2aa = waa*waa;
    for(int i=0; i<na; i++){
        Vec3d  pi   = apos[i];
        double qi   = aQ[i];
        Vec3d  abwi = aAbWs[i];
        for(int j=0; j<i; j++){
            Vec3d f;
            Vec3d  abw;
            combineAbW( abwi, aAbWs[j], abw );
            //if( (i_DEBUG>0) && (1==qi==aQ[j]) ){ printf( " abw(H-H): %i,%i A %g B %g w %g \n", i,j, abw.x, abw.y, abw.z ); }
            Eaa += addPairEF_expQ( apos[j]-pi, f, abw.z, qi*aQ[j], abw.y, abw.x );
            //   ToDo : Pauli Repulsion of core electrons ?????
            aforce[j].sub(f);
            aforce[i].add(f);
            //glColor3f(1.0,0.0,0.0);
            //Draw3D::drawVecInPos( f*-1, apos[j] );
            //Draw3D::drawVecInPos( f   , pi      );
        }
    }
    return Eaa;
}

double eval(){
    return
    evalKinetic()
    + evalEE()
    + evalAE()
    + evalAA()
    ;
}


void move_GD(double dt){
    for(int i=0;i<nDOFs;i++){
        //pDOFs[i].add_mul(fDOFs[i],dt);
        pDOFs[i] += fDOFs[i] * dt;
    }
}

void autoAbWs( const Vec3d * AAs, const Vec3d * AEs ){
    for(int i=0;i<na;i++){
        int ityp =(int)round(aQ[i]);
        aAbWs[i] = AAs[ityp];
        eAbWs[i] = AEs[ityp];
        printf( "atom[%i] typ %i aAbW(%g,%g,%g) eAbW(%g,%g,%g) \n", i, ityp, aAbWs[i].x, aAbWs[i].y, aAbWs[i].z,  eAbWs[i].x, eAbWs[i].y, eAbWs[i].z );
    }
}

bool loadFromFile_xyz( char const* filename ){
    printf(" filename: >>%s<< \n", filename );
    FILE * pFile;
    pFile = fopen (filename,"r");
    int ntot;
    fscanf (pFile, " %i \n", &ntot );
    fscanf (pFile, " %i %i\n", &na, &ne );
    printf("na %i ne %i \n", na, ne );
    realloc( na, ne );
    char buf[1024];
    //int ntot=na+ne;

    double Qasum = 0.0;

    int ia=0,ie=0;
    //int ia=0,ie=nDOFs-1;
    //int ia=0;
    //int ie=nDOFs-1;
    for (int i=0; i<ntot; i++){
        double x,y,z;
        int e;
        fgets( buf, 256, pFile); //printf( ">%s<\n", buf );
        int nw = sscanf (buf, " %i %lf %lf %lf", &e, &x, &y, &z );
        if( e<0){
            epos[ie]=(Vec3d){x,y,z};
            if     (e==-1){ espin[ie]= 1; }
            else if(e==-2){ espin[ie]=-1; }
            esize[ie] = 1.0;
            ie++;
            printf( " e[%i] ", ie );
        }else{
            apos[ia]=(Vec3d){x,y,z};
            aQ  [ia]=e;  // change this later
            //aAbws[ia] = default_aAbWs[e];
            //eAbws[ia] = default_eAbWs[e];
            Qasum += e;
            ia++;
            printf( " a[%i] ", ia );
        };
        printf( " %i %f %f %f  \n", e, x,y,z );
    }
    clearForce();
    printf( "Qtot = %g (%g - 2*%i) \n",  Qasum - ne, Qasum, ne );
    fclose (pFile);
    return 0;
}

};

#endif
