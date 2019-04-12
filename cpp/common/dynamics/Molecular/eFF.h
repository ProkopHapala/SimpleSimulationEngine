
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


struct EFFAtomType{
    double Q; // nuncear charge
    double A; // prefactor
    double B; // exponential beta
};

inline double frEFF_r2( double r2, double s1, double s2 ){
    return (1/(1 + s1*r2) - 1/(1 + s2*r2))/( 1/s1 + 1/s2 );
}

inline double fr_ee(double r, double s){ 
    return 1/(r*r+s*s);
};

inline double fr_ae(double r, double s){ 
    return exp(5*r) - 1/(r*r+s*s);
};

inline double fr_aa(double r, double s){
    return 1/(r*r+s*s);
};

class EFF{ public:

    //double dvmax = 0.1;
    //double dpmax = 0.1;

    double emass = 1.0;
    double see = 1.0;
    double saa = 0.1;

    double sa2  = saa*saa;
    double se2  = see*see;
    double sea2 = see*see;

    int ne=0,na=0;
    //int*   atype  =0;

    double* aQ    =0;
    Vec3d* apos   =0;
    Vec3d* aforce =0;
    Vec3d* avel   =0;

    double* espin =0;
    Vec3d* epos   =0;
    Vec3d* eforce =0;
    Vec3d* evel   =0;

    void realloc(int na_, int ne_){
        na=na_; ne=ne_;
        _realloc( aQ  ,   na);
        _realloc( apos,   na);
        _realloc( aforce, na);
        _realloc( avel,   na);

        _realloc( espin  ,ne);
        _realloc( epos   ,ne);
        _realloc( eforce ,ne);
        _realloc( evel   ,ne);
    }

    void clearForce(){
        for(int i=0; i<ne; i++){ eforce[i].set(0.0); };
        for(int i=0; i<na; i++){ aforce[i].set(0.0); };
    }

    void clearVel(){
        for(int i=0; i<ne; i++){ evel[i].set(0.0); };
        for(int i=0; i<na; i++){ avel[i].set(0.0); };
    }

    double getFmax(){
        double fm=0;
        for(int i=0; i<ne; i++){ fm=fmax(fm,eforce[i].norm2()); };
        for(int i=0; i<na; i++){ fm=fmax(fm,aforce[i].norm2()); };
        return fm;
    }

    double getVmax(){
        double vm=0;
        for(int i=0; i<ne; i++){ vm=fmax(vm,evel[i].norm2()); };
        for(int i=0; i<na; i++){ vm=fmax(vm,avel[i].norm2()); };
        return vm;
    }

    void forceEE(){
        double invSee = 1/(see*see);
        for(int i=0; i<ne; i++){
            Vec3d pi = epos[i];
            for(int j=0; j<i; j++){
                Vec3d  d = epos[j] - pi;
                double r = d.norm();
                //fr       = Aee*exp(Bee*r);
                //double fr = frEFF_r2( r*r, 1.0, invSee );
                double fr = 4/(r*r + se2);
                //printf(" fee %i %i %g %g \n", i, j, r, fr);
                Vec3d f = d*(fr/r);
                double fmax=0;
                eforce[i].sub(f);
                eforce[j].add(f);
            }
        }
    }

    void forceAE(){
        //double see2 = see*see;
        //double saa2 = saa*saa;
        double invSae = 1/( see*see + saa*saa );
        for(int i=0; i<na; i++){
            Vec3d pi = apos[i];
            double qi = aQ[i];
            for(int j=0; j<ne; j++){
                Vec3d  d = epos[j] - pi;
                double r = d.norm();
                //double fr = qi* frEFF_r2( r*r, 1.0, invSae );
                //double fr = 2*qi*(  -2/(r*r + se2) + 1/(r*r + sa2) );
                //double fr = 2*qi*(  -2/(r*r + sea2) + 1/(r*r + sa2) );
                double fr = 2*qi*(  -1/(r*r + sea2) + 4*exp(-3*r) );
                //printf(" fae %i %i %g %g \n", i, j, r, fr);
                Vec3d f = d*(fr/r);
                aforce[i].sub(f);
                eforce[j].add(f);
            }
        }
    }

    void forceAA(){
        double invSaa = 1/(saa*saa);
        for(int i=0; i<na; i++){
            Vec3d pi = apos[i];
            double qi = aQ[i];
            for(int j=0; j<i; j++){
                Vec3d  d  = apos[j] - pi;
                double r  = d.norm();
                //double fr = aQ[j]*qi * frEFF_r2( r*r, 1.0, invSaa );
                double fr = aQ[j]*qi/(r*r + sa2);
                //printf(" faa %i %i %g %g \n", i, j, r, fr);
                Vec3d f   = d*(fr/r);
                aforce[i].sub(f);
                aforce[j].add(f);
            }
        }
    }

    void move(double dt, double damp){
        //for(int isub=0; i<; i++)   // sub electron step ?
        for(int i=0;i<ne;i++){
            Vec3d& v = evel[i];
            v.mul(damp);
            v.add_mul( eforce[i], dt );
            epos[i].add_mul( v,  dt );
        }
        for(int i=0;i<na;i++){
            Vec3d& v = avel[i];
            v.mul(damp);
            v.add_mul( aforce[i], dt );
            apos[i].add_mul( v,  dt );
        }
    }

    double getSafeTimeStep(double dt_, double dvmax, double dpmax){
        double dt = dt_;
        double fm = getFmax();
        double vm = getVmax();
        //if( fm*dt > dvmax ) dt_=dvmax/fm;
        //if( vm*dt > dpmax ) dt_=dpmax/vm;
        dt=fmin( dt, dvmax/fm );
        dt=fmin( dt, dpmax/vm );
        //printf( "fmax %g vmax %g dt %g \n", fm, vm, dt );
        return dt;
    }

    int run(int nStep, double dt, double damp){
        for(int i=0;i<nStep; i++){
            clearForce();
            forceEE();
            forceAE();
            forceAA();
            double dt_ = getSafeTimeStep(dt, 0.1, 0.1 );
            move(dt_, damp);
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
        if( e==0){
            epos[ie]=(Vec3d){x,y,z};
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
    clearVel();
    clearForce();
    printf( "Qtot = %g (%g - 2*%i) \n",  Qasum - 2*ne, Qasum, ne );
    fclose (pFile);
    return 0;
}

};

#endif
