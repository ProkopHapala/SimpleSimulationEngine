
#ifndef EFF_old_h
#define EFF_old_h

/// @file
/// @ingroup eFF

#include "fastmath.h"
//#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Forces.h"


#include "InteractionsGauss_old.h"


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

Params
            a             b             c           d           e           s-core
Al     0.486000       1.049000      0.207000                              1.660000
Si     0.320852       2.283269      0.814857                              1.691398
C     22.721015       0.728733      1.103199     17.695345    6.693621    0.621427
O     25.080199       0.331574      1.276183     12.910142    3.189333    0.167813

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

// ToDo : Later properly
constexpr static const Vec3d default_eAbWs[] = {
{ 0.0,  0.0, 0.0},  // Q = 0 //
{ 0.0,  0.0, 0.01}, // Q = 1 // H
{ 2.0, -3.0, 0.1},  // Q = 2 // Be?
{ 2.0, -3.0, 0.1},  // Q = 3 // B
{ 2.0, -3.0, 0.1},  // Q = 4 // C
};

constexpr static const Vec3d default_aAbWs[] = {
{ 0.0,  0.0, 0.1},  // Q = 0 //
{ 0.0,  0.0, 0.01},  // Q = 1 // H
{ 1.0, -5.0, 0.1},  // Q = 2 // Be?
{ 1.0, -5.0, 0.1},  // Q = 3 // B
{ 1.0, -5.0, 0.1},  // Q = 4 // C
};

constexpr static const  double default_EPCs[] = {
// s-core       a                b               c             d             e
 0.621427,   22.721015,       0.728733,      1.103199,     17.695345,    6.693621,  // C
 0.167813,   25.080199,       0.331574,      1.276183,     12.910142,    3.189333,  // O
 1.660000,    0.486000,       1.049000,      0.207000,     -1       ,   -1       ,  // Al
 1.691398,    0.320852,       2.283269,      0.814857,     -1       ,  - 1       ,  // Si
};

// Reformulate    d=b  ,   e=c
/*
constexpr static const default_EPCs[] = {
// s-core        a               b              c             d               e
 0.621427,   22.721015,      17.695345,      6.693621,     ,    ,  // C
 0.167813,   25.080199,      12.910142,      3.189333,     ,    ,  // O
 1.660000,    0.486000,       1.049000,      0.207000,     -1       ,   -1       ,  // Al
 1.691398,    0.320852,       2.283269,      0.814857,     -1       ,  - 1       ,  // Si
};
*/

struct EFFAtomType{
    int ne;       // number of valence electrons
    double Rcore; // core electron radius
    double Acore; // core electron repulsion (depend on number of core electrons)
    // ToDo : maybe distinct core radius for pi-electrons (?)
    //double RcoreP; // core electron radius
    //double AcoreP; // core electron repulsion (depend on number of core electrons)
};


EFFAtomType EFFparams[11]{
    {0,0.,0.}, // void ... (or electron?)
    {1,0.,0.}, // H
    {2,0.,0.}, // He // only core ?

    {1,0.1,2.0}, // Li
    {2,0.1,2.0}, // Be
    {3,0.1,2.0}, // B
    {4,0.1,2.0}, // C
    {5,0.1,2.0}, // N
    {6,0.1,2.0}, // O
    {7,0.1,2.0}, // F
    {8,0.1,2.0}  // Ne // only core ?
};





inline double frEFF_r2( double r2, double w1, double w2 ){
    return ( 1/(1 + w1*r2) - 1/(1 + w2*r2) )/( 1/w1 + 1/w2 );
}

inline void combineAbW( const Vec3d& abwi, const Vec3d& abwj, Vec3d& abw ){
    abw.x = abwi.x*abwj.x;
    abw.y = (abwi.y+abwj.y)*0.5;
    abw.z = abwi.z+abwj.z;
};

inline double addPairEF_expQ( const Vec3d& d, Vec3d& f, double w2, double qq, double bExp, double aExp ){
    double r2     = d.norm2();
    double invrw2 = 1./( r2 + w2 );
    double invrw  = sqrt(invrw2);
    double E      =  const_El_eVA*qq*invrw;
    double fr     = -E*invrw2;
    if( bExp<0 ){
        double r      = sqrt( r2+R2SAFE );
        double Eexp  = aExp*exp( bExp*r );
        fr          += bExp*Eexp/r;
        E           += Eexp;
    }
    f.set_mul( d, fr );
    return E;
}

inline double interp_gx4(double r2, double y1, double y2 ){
    double c = (1-r2);
    c*=c;
    return c*y1 + (1-c)*y2;
}

/// EFF solver
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

    int iPauliModel=0;

    constexpr static const double default_esize = 0.5;
    constexpr static const Vec3d KRSrho = { 1.125, 0.9, 0.2 }; ///< eFF universal parameters
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

    int ne=0,na=0,nDOFs=0; ///< number of electrons, atoms, degrees of freedom
    //int*   atype  =0;

    double * aQ     =0; ///< atomic charges
    Vec3d  * apos   =0; ///< atomic positions
    Vec3d  * aforce =0; ///< atomic forces

    Vec3d  * aAbWs  =0; ///< atomic   parameters (amplitude, decay, width)
    Vec3d  * eAbWs  =0; ///< electron parameters (amplitude, decay, width)

    //double * espin  =0;
    int8_t * espin  =0; ///< electron spins
    Vec3d  * epos   =0; ///< electron positions
    Vec3d  * eforce =0; ///< electron forces
    double * esize  =0; ///< electron size
    double * fsize  =0; ///< electron force on size

    double* pDOFs =0;  ///< buffer of degrees of freedom
    double* fDOFs =0;  ///< buffer of forces on degrees of freedom

    double Ek=0, Eee=0,Eae=0,Eaa=0, EeePaul=0,EaePaul=0; ///< different kinds of energy

void realloc(int na_, int ne_){
    na=na_; ne=ne_;
    nDOFs=na*3+ne*3 + ne;
    _realloc( pDOFs, nDOFs);
    _realloc( fDOFs, nDOFs);

    _realloc( aQ  ,  na);
    _realloc( aAbWs, na);
    _realloc( eAbWs, na);
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
    esize  =          pDOFs + na*3 + ne*3;
    fsize  =          fDOFs + na*3 + ne*3;

    //_realloc(apos  ,na );
    //_realloc(aforce,na );
    //_realloc(epos  ,ne );
    //_realloc(eforce,ne );
    //_realloc(esize ,ne );
    //_realloc(fsize ,ne );

}

void dealloc(){
    delete [] pDOFs;
    delete [] fDOFs;
    delete [] aQ   ;
    delete [] aAbWs;
    delete [] eAbWs;
    delete [] espin;
}

/// evaluate kinetic energy of each electron
double evalKinetic(){
    Ek=0;
    for(int i=0; i<ne; i++){
        Ek += addKineticGauss( esize[i], fsize[i] );
        //if( i_DEBUG>0 ) printf( "evalKinetic[%i] s %g -> f %g Ek %g \n", i, esize[i], fsize[i], Ek );
    }
    return Ek;
}

/// evaluate Electron-Electron forces
double evalEE(){
    Eee    =0;
    EeePaul=0;
    //double w2ee = wee*wee;
    const double qq = QE*QE;
    for(int i=0; i<ne; i++){
        const Vec3d    pi  = epos[i];
        //Vec3d&   fi  = eforce[i];
        const int8_t spini = espin[i];
        const double   si  = esize[i];
        double&       fsi  = fsize[i];
        for(int j=0; j<i; j++){
            Vec3d  f  = Vec3dZero;
            const Vec3d  dR = epos [j] - pi;
            const double sj = esize[j];
            double& fsj = fsize[j];
            Eee     += addCoulombGauss( dR, si, sj, f, fsi, fsj, qq );
            //EeePaul += addPauliGauss  ( dR, si, sj, f, fsi, fsj, spini!=espin[j], KRSrho );
            if( spini==espin[j] ) EeePaul += addPauliGauss  ( dR, si, sj, f, fsi, fsj, false, KRSrho );
            //if( spini==espin[j] ) EeePaul += addDensOverlapGauss_S( dR,si,sj, 1, f, fsi, fsj );
            //Eee += addPairEF_expQ( epos[j]-pi, f, w2ee, +1, 0, 0 );
            //if( i_DEBUG>0 ) printf( "evalEE[%i,%i] dR(%g,%g,%g) s(%g,%g) q %g  ->   f(%g,%g,%g) fs(%g,%g) \n", i,j, dR.x,dR.y,dR.z, si,sj, qq,   f.x,f.y,f.z, fsi,fsj );
            eforce[j].sub(f);
            eforce[i].add(f);

            //DEBUG_fa_aa[j].sub(f);
            //DEBUG_fa_aa[i].add(f);
            //if( (i==DEBUG_i)&&(j==DEBUG_j) ){
            //    glColor3f(1.0,0.0,0.0);
            //    Draw3D::drawVecInPos( f*-1., epos[j] );
            //    Draw3D::drawVecInPos( f   , pi      );
            //}
        }
    }
    //if( i_DEBUG>0 )  for(int j=0; j<ne; j++){  printf( "evalEE: esize[%i] %g f %g \n", j, esize[j], fsize[j] ); }
    return Eee+EeePaul;
}

/// evaluate Atom-Electron forces
double evalAE(){
    Eae    =0;
    EaePaul=0;
    //double see2 = see*see;
    //double saa2 = saa*saa;
    //double invSae = 1/( see*see + saa*saa );
    //double w2ae = wae*wae;
    for(int i=0; i<na; i++){
        const Vec3d  pi   = apos[i];
        const double qqi  = aQ[i]*QE;
        const Vec3d  abwi = eAbWs[i];
        for(int j=0; j<ne; j++){
            Vec3d f=Vec3dZero;
            const Vec3d   dR  = epos [j] - pi;
            const double  sj  = esize[j];
            double& fsj = fsize[j];
            double  fs_junk;
            //Eae += addPairEF_expQ( epos[j]-pi, f, abwi.z, qi*QE, abwi.y, abwi.x );
            Eae                      += addCoulombGauss      ( dR,sj,                 f, fsj, qqi     );     // correct
            if(qqi<-1.00001) EaePaul += addDensOverlapGauss_S( dR,sj, abwi.z, abwi.a, f, fsj, fs_junk );     // correct
            //if(qqi<-1.00001) EaePaul += addPauliGauss  ( dR, sj, abwi.z, f, fsj, fs_junk, false, KRSrho );     // correct

            //if( i_DEBUG>0 ) printf( "evalAE[%i,%i] dR(%g,%g,%g) s %g q %g  ->   f(%g,%g,%g) fs %g \n", i,j, dR.x,dR.y,dR.z, sj, qqi,   f.x,f.y,f.z, fsj );
            eforce[j].sub(f);
            aforce[i].add(f);

            //DEBUG_fe_ae[j].sub(f);
            //DEBUG_fa_ae[i].add(f);
            //if( (i==DEBUG_i)&&(j==DEBUG_j) ){
            //    glColor3f(1.0,1.0,1.0); Draw3D::drawLine    ( pi, epos[j]    );
            //    glColor3f(0.0,1.0,0.0); Draw3D::drawVecInPos( f*-1., epos[j] );
            //    glColor3f(1.0,0.0,1.0); Draw3D::drawVecInPos( f    , pi      );
            //}
        }
    }
    //if( i_DEBUG>0 )  for(int j=0; j<ne; j++){  printf( "evalAE: esize[%i] %g f %g \n", j, esize[j], fsize[j] ); }
    return Eae+EaePaul;
}

/// evaluate Atom-Atom forces
double evalAA(){
    if( i_DEBUG>0 ) printf( "evalAA \n" );
    Eaa=0;
    //double invSaa = 1/(saa*saa);
    //double w2aa = waa*waa;
    for(int i=0; i<na; i++){
        const Vec3d  pi   = apos[i];
        const double qi   = aQ[i];
        const Vec3d  abwi = aAbWs[i];
        for(int j=0; j<i; j++){
            Vec3d f = Vec3dZero; // HERE WAS THE ERROR (missing initialization !!!! )
            Vec3d  abw;
            combineAbW( abwi, aAbWs[j], abw );
            const Vec3d dR  = apos[j] - pi;
            //if( (i_DEBUG>0) && (1==qi==aQ[j]) ){ printf( " abw(H-H): %i,%i A %g B %g w %g \n", i,j, abw.x, abw.y, abw.z ); }
            //Eaa += addPairEF_expQ( apos[j]-pi, f, abw.z, qi*aQ[j], abw.y, abw.x );
            //Eaa += addPairEF_expQ( apos[j]-pi, f, abw.z, qi*aQ[j], abw.y, abw.x );
            Eaa +=  addAtomicForceQ( dR, f, qi*aQ[j] );
            //   ToDo : Pauli Repulsion of core electrons ?????
            aforce[j].sub(f);
            aforce[i].add(f);

            //DEBUG_fa_aa[j].sub(f);
            //DEBUG_fa_aa[i].add(f);
            //if( (i==DEBUG_i)&&(j==DEBUG_j) ){
            //    glColor3f(1.0,0.0,0.0);
            //    //Draw3D::drawVecInPos( dR*-1., apos[j] );
            //    //Draw3D::drawVecInPos( dR    , pi      );
            //    Draw3D::drawVecInPos( f*-1., apos[j] );
            //   Draw3D::drawVecInPos( f    , pi      );
            //}
        }
    }
    return Eaa;
}

/// evaluate full Electron Forcefild
double eval(){
    return
    evalKinetic()
    + evalEE()
    + evalAE()
    + evalAA()
    ;
}

void clearForce(){
    for(int i=0; i<nDOFs; i++){ fDOFs[i]=0; }
}

void clearForce_noAlias(){
    for(int i=0;i<na;i++){
        aforce[i]=Vec3dZero;
    }
    for(int i=0;i<ne;i++){
        eforce[i]=Vec3dZero;
        fsize[i]=0;
    }
}

void move_GD(double dt){
    double sum = 0;
    for(int i=0;i<nDOFs;i++){
        //pDOFs[i].add_mul(fDOFs[i],dt);
        sum += fDOFs[i];
        pDOFs[i] += fDOFs[i] * dt;
    }
    //printf( "move_GD sum %g ", sum );
}

void move_GD_noAlias(double dt){
    //Vec3d fe=Vec3dZero,fa=Vec3dZero;
    //double fs=0;
    for(int i=0;i<na;i++){
        apos[i].add_mul( aforce[i], dt );
        //fa.add( aforce[i] );
    }
    for(int i=0;i<ne;i++){
        epos [i].add_mul( eforce[i], dt );
        //fe.add( aforce[i] );
        esize[i] += fsize[i] * dt;
        //fs += fsize[i];
    }
    //printf( "fs %g fe(%g,%g,%g) fa(%g,%g,%g)\n", fs, fe.x,fe.y,fe.z, fa.x,fa.y,fa.z );
}

void autoAbWs( const Vec3d * AAs, const Vec3d * AEs ){
    for(int i=0;i<na;i++){
        int ityp =(int)round(aQ[i]);
        aAbWs[i] = AAs[ityp];
        eAbWs[i] = AEs[ityp];
        printf( "atom[%i] typ %i aAbW(%g,%g,%g) eAbW(%g,%g,%g) \n", i, ityp, aAbWs[i].x, aAbWs[i].y, aAbWs[i].z,  eAbWs[i].x, eAbWs[i].y, eAbWs[i].z );
    }
}


double atomsPotAtPoint( const Vec3d& pos, double s, double Q )const{
    double Eae    =0;
    double EaePaul=0;
    /*
    Vec3d fp; double fs;
    for(int i=0; i<na; i++){
        const Vec3d  dR   = pos-apos[i];
        const Quat4d aPar = aPars[i]; // { x=Q,y=sQ,z=sP,w=cP }
        const double qq   = aPar.x*Q;
        Vec3d  f;
        double fsi,fsj;
        //addCoulombGauss      ( const Vec3d& dR, double si, double sj,             Vec3d& f, double& fsi, double& fsj, double qq ){
        //addDensOverlapGauss_S( const Vec3d& dR, double si, double sj, double amp, Vec3d& f, double& fsi, double& fsj            ){
        Eae  += addCoulombGauss              ( dR, s, aPar.y,         f, fsi, fsj, qq );
        if( aPar.w>1e-8 ){
            EaePaul+= addDensOverlapGauss_S( dR, s, aPar.z, aPar.w, f, fsi, fsj         );
        }
    }
    */
    return Eae + EaePaul;
}
double* atomsPotAtPoints( int n, Vec3d* ps, double* out=0, double s=0.0, double Q=1.0 )const{
    if(out==0){ out = new double[n]; };
    for(int i=0; i<n; i++){
        out[i] = atomsPotAtPoint( ps[i], s, Q );
    }
    return out;
}

void info(){
    for(int i=0; i<na; i++){
        printf( "a[%i] p(%g,%g,%g) q %g eAbW(%g,%g,%g) aAbW(%g,%g,%g) \n", i, apos[i].x, apos[i].y, apos[i].z, aQ[i], eAbWs[i].z,eAbWs[i].z,eAbWs[i].z, aAbWs[i].z,aAbWs[i].z,aAbWs[i].z );
    }
    for(int i=0; i<ne; i++){
        printf( "e[%i] p(%g,%g,%g) sz %g s %i \n", i, epos[i].x, epos[i].y, epos[i].z, esize[i], espin[i] );
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
        double x,y,z,s;
        int e;
        fgets( buf, 256, pFile); //printf( ">%s<\n", buf );
        int nw = sscanf (buf, " %i %lf %lf %lf", &e, &x, &y, &z, &s );
        if( e<0){
            epos[ie]=(Vec3d){x,y,z};
            if     (e==-1){ espin[ie]= 1; }
            else if(e==-2){ espin[ie]=-1; }
            if( nw>4 ){ esize[ie]=s; }else{ esize[ie]=default_esize; }
            //esize[ie] = 1.0;
            //esize[ie] = 0.5;
            //esize[ie] = 0.25;
            ie++;
            //printf( " e[%i] ", ie );
        }else{
            apos[ia]=(Vec3d){x,y,z};
            aQ  [ia]=e;  // change this later
            //aAbws[ia] = default_aAbWs[e];
            //eAbws[ia] = default_eAbWs[e];
            Qasum += e;
            ia++;
            printf( " a[%i] ", ia );
        };
        //printf( " %i %f %f %f  \n", e, x,y,z );
    }
    clearForce();
    printf( "Qtot = %g (%g - 2*%i) \n",  Qasum - ne, Qasum, ne );
    fclose (pFile);
    return 0;
}


bool loadFromFile_fgo( char const* filename ){
    //printf(" filename: >>%s<< \n", filename );
    FILE * pFile;
    pFile = fopen (filename,"r");
    if( pFile == NULL ){
        printf("ERROR in eFF::loadFromFile_fgo(%s) : No such file !!! \n", filename );
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
    if(perOrb_!=1){ printf("ERROR in eFF::loadFromFile_fgo(%s) : perOrb must be =1 ( found %i instead) !!! \n", filename, perOrb_ );};
    if(bClosedShell) nOrb_*=2;
    realloc( natom_, nOrb_ );
    double Qasum = 0.0;
    for(int i=0; i<na; i++){
        double x,y,z;
        double Q,sQ,sP,cP;
        fgets( buff, nbuff, pFile);
        int nw = sscanf (buff, "%lf %lf %lf %lf %lf %lf %lf", &x, &y, &z, &Q, &sQ, &sP, &cP );
        //printf( "atom[%i] p(%g,%g,%g) Q %g sQ %g sP %g cP %g \n", i, x, y, z,    Q, sQ, sP, cP );
        Q=-Q;
        apos  [i]=(Vec3d){x,y,z};
        //aPars[i].set(Q,sQ,sP,cP);
        aQ  [i]=Q;

        Qasum += Q;
    }
    int nBasRead = ne;
    if( bClosedShell ) nBasRead/=2;
    for(int i=0; i<nBasRead; i++){
        double x,y,z;
        double s,c;
        int spin;
        fgets( buff, nbuff, pFile); // printf( "fgets: >%s<\n", buf );
        int nw = sscanf (buff, "%lf %lf %lf %lf %lf %i", &x, &y, &z,  &s, &c, &spin );
        epos [i]=(Vec3d){x,y,z};
        esize[i]=s;
        //ecoef[i]=c;
        //int io=i/perOrb;
        if( !bClosedShell ){ if(nw>5)espin[i]=spin; }else{ espin[i]=1; };
        //printf( "ebasis[%i,%i|%i] p(%g,%g,%g) s %g c %g spin %i | nw %i io %i \n", i/perOrb, i%perOrb,i, x, y, z,  s, c, spin,  nw, io  );
    }
    if( bClosedShell ){
        for(int i=0; i<nBasRead; i++){
            int j = i+nBasRead;
            epos [j]=epos[i];
            esize[j]=esize[i];
            //ecoef[j]=ecoef[i];
            espin[j]=-1;
        }
    }
    printf( "na %i ne %i \n", na, ne );
    //printf( "Qtot = %g (%g - 2*%i) \n",  Qasum - nOrb, Qasum, nOrb );
    fclose (pFile);
    return 0;
}

};

#endif
