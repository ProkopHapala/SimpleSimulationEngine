
#ifndef EFF_h
#define EFF_h

#include "fastmath.h"
//#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
//#include "Forces.h"

/*
Erf approximation:
# Gaussian:    F = (x2-1)**2 / sqrtPi
# Erf          E = x*(1 + x2 * ( -0.66666666666 + 0.2*x2 ) ) * (2/(16.0/15.0))
*/

//double wae = 1.0;
//double bAE = -5.0;
//double aAE = 20.0;


// ToDo : Later properly
constexpr static const Vec3d default_AbWs[] = {
{ 0.0,  0.0, 0.0},  // Q = 0 //
{ 5.0, -8.0, 0.5},  // Q = 1 // H
{20.0, -5.0, 1.0},  // Q = 2 // Be?
{20.0, -5.0, 1.0},  // Q = 3 // B
{20.0, -5.0, 1.0},  // Q = 4 // C
{20.0, -5.0, 1.0},  // Q = 5 // N
{20.0, -5.0, 1.0},  // Q = 5 // O
{20.0, -5.0, 1.0},  // Q = 5 // F
};


struct EFFAtomType{
    double Q; // nuncear charge
    double A; // prefactor
    double B; // exponential beta
};

inline double frEFF_r2( double r2, double w1, double w2 ){
    return ( 1/(1 + w1*r2) - 1/(1 + w2*r2) )/( 1/w1 + 1/w2 );
}

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
    double r      = sqrt( r2+R2SAFE );
    double invrw2 = 1/( r2 + w2 );
    double E      =  qq*sqrt(invrw2);
    double fr     = -qq*invrw2;
    if( bExp<0 ){
        double Eexp  = aExp*exp( bExp*r );
        fr          += bExp*Eexp;
        E           += Eexp;
    }
    f.set_mul( d, fr/r );
    return E;
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

    double wee = 1.0;
    double wae = 1.0;
    double waa = 0.1;

    double bAE = -5.0;
    double aAE = 20.0;

    int ne=0,na=0,nDOFs=0;
    //int*   atype  =0;

    double * aQ     =0;
    Vec3d  * apos   =0;
    Vec3d  * aforce =0;

    Vec3d  * AbWs   =0;

    //double * espin  =0;
    int8_t * espin  =0;
    Vec3d  * epos   =0;
    Vec3d  * eforce =0;

    Vec3d * pDOFs =0;
    Vec3d * fDOFs =0;

    double Eee=0,Eae=0,Eaa=0;

void realloc(int na_, int ne_){
    na=na_; ne=ne_;
    nDOFs=na+ne;
    _realloc( pDOFs,   nDOFs);
    _realloc( fDOFs,   nDOFs);

    _realloc( aQ  ,   na);
    _realloc( AbWs,   na);
    //_realloc( apos,   na);
    //_realloc( aforce, na);
    //_realloc( avel,   na);
    //_realloc( aorbs,  na);

    _realloc( espin  ,ne);
    //_realloc( epos   ,ne);
    //_realloc( eforce ,ne);
    //_realloc( evel   ,ne);

    apos   = pDOFs;
    aforce = fDOFs;

    epos   = pDOFs + na;
    eforce = fDOFs + na;

}

void clearForce(){
    for(int i=0; i<nDOFs; i++){ fDOFs[i]=Vec3dZero; };
}

double evalEE(){
    Eee=0;
    double w2ee = wee*wee;
    for(int i=0; i<ne; i++){
        Vec3d pi = epos[i];
        for(int j=0; j<i; j++){
            Vec3d f;
            Eee += addPairEF_expQ( epos[j]-pi, f, w2ee, +1, 0, 0 );
            eforce[j].sub(f);
            eforce[i].add(f);
            glColor3f(1.0,0.0,0.0);
            Draw3D::drawVecInPos( f*-1, epos[j] );
            Draw3D::drawVecInPos( f   , pi      );
        }
    }
    return Eee;
}

double evalAE(){
    Eae=0;
    //double see2 = see*see;
    //double saa2 = saa*saa;
    //double invSae = 1/( see*see + saa*saa );
    double w2ae = wae*wae;
    for(int i=0; i<na; i++){
        Vec3d pi = apos[i];
        double qi = aQ[i];
        for(int j=0; j<ne; j++){
            Vec3d f;
            Eae += addPairEF_expQ( epos[j]-pi, f, w2ae, -qi, bAE, aAE );
            //printf(  "a[%i]e[%i] r %g\n", i, j, (epos[j]-pi).norm() );
            eforce[j].sub(f);
            aforce[i].add(f);
            glColor3f(1.0,0.0,1.0);
            Draw3D::drawVecInPos( f*-1, epos[j] );
            Draw3D::drawVecInPos( f   , pi      );
        }
    }
    return Eae;
}

double evalAA(){
    Eaa=0;
    //double invSaa = 1/(saa*saa);
    double w2aa = waa*waa;
    for(int i=0; i<na; i++){
        Vec3d pi = apos[i];
        double qi = aQ[i];
        for(int j=0; j<i; j++){
            Vec3d f;
            Eaa += addPairEF_expQ( apos[j]-pi, f, w2aa, qi*aQ[j], 0, 0 );
            aforce[j].sub(f);
            aforce[i].add(f);
            glColor3f(1.0,0.0,0.0);
            Draw3D::drawVecInPos( f*-1, apos[j] );
            Draw3D::drawVecInPos( f   , pi      );
        }
    }
    return Eaa;
}

void move_GD(double dt){
    for(int i=0;i<nDOFs;i++){
        pDOFs[i].add_mul(fDOFs[i],dt);
    }
}

void autoAbWs(const Vec3d * params ){
    for(int i=0;i<na;i++){
        int ityp =(int)round(aQ[i]);
        printf( "atom[%i] typ %i AbW(%g,%g,%g) \n", i, ityp, params[ityp].x, params[ityp].y, params[ityp].z );
        AbWs[i] = params[ityp];
    }
}



bool loadFromFile_bas( char const* filename ){
    printf(" filename: >>%s<< \n", filename );
    FILE * pFile;
    pFile = fopen (filename,"r");
    fscanf (pFile, " %i %i\n", &na, &ne );
    printf("na %i ne %i \n", na, ne );
    realloc( na, ne );
    char buf[1024];
    int ntot=na+ne;
    int ia=0,ie=0;
    double Qasum = 0.0;
    for (int i=0; i<ntot; i++){
        double x,y,z;
        int e;
        fgets( buf, 256, pFile); //printf( ">%s<\n", buf );
        int nw = sscanf (buf, " %i %lf %lf %lf", &e, &x, &y, &z );
        if( e<0){
            epos[ie]=(Vec3d){x,y,z};
            if     (e==-1){ espin[ie]= 1; }
            else if(e==-2){ espin[ie]=-1; }
            ie++;
            printf( " e[%i] ", ie );
        }else{
            apos[ia]=(Vec3d){x,y,z};
            aQ  [ia]=e;  // change this later
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
