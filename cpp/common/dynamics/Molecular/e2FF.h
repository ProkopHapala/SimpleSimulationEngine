
#ifndef EFF_h
#define EFF_h

#include "fastmath.h"
//#include "Vec2.h"
#include "Vec3.h"
#include "Mat4.h"
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

inline double interp_gx4(double r2, double y1, double y2 ){
    double c = (1-r2);
    c*=c;
    return c*y1 + (1-c)*y2;
}

class E2FF{ public:

    double lee = 0.5;

    //double dvmax = 0.1;
    //double dpmax = 0.1;

    double emass = 1.0;
    double see   = 0.5;
    double saa   = 0.1;

    double sa2  = saa*saa;
    double se2  = see*see;
    double sea2 = see*see;

    int ne=0,na=0;
    //int*   atype  =0;

    double* aQ    =0;
    Vec3d* apos   =0;
    Mat4d* aorbs  =0;
    Vec3d* aforce =0;
    Vec3d* avel   =0;

    //double* espin =0;
    Vec3d* epos   =0;
    Vec3d* edir   =0;
    Vec3d* eforce =0;
    Vec3d* evel   =0;

    void realloc(int na_, int ne_){
        na=na_; ne=ne_;
        _realloc( aQ  ,   na);
        _realloc( apos,   na);
        _realloc( aorbs,  na);
        _realloc( aforce, na);
        _realloc( avel,   na);

        //_realloc( espin  ,ne);
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

    void constrain_e_pairs(){
        printf("constrain_e_pairs \n");
        for(int i=0; i<ne; i+=2){
            Vec3d p1 = epos[i  ];
            Vec3d p2 = epos[i+1];
            Vec3d d  = p1-p2;
            double fac = 0.5*lee/( 0.0001 + d.norm() );
            printf("fac %g\n", fac );
            d.mul(fac );
            epos[i  ].sub(d);
            epos[i+1].add(d);
        }
    }

    void force_e_pairs(double K){
        //printf("constrain_e_pairs \n");
        for(int i=0; i<ne; i+=2){
            Vec3d p1 = epos[i  ];
            Vec3d p2 = epos[i+1];
            Vec3d d  = p1-p2;
            double l = d.norm();
            d.mul( K*(l-lee)/l );
            eforce[i  ].sub(d);
            eforce[i+1].add(d);
        }
    }

    void ortho_step(){
        for(int i=0; i<na; i++){
            Mat4d& m = aorbs[i];
            //double aa = m.a.dot(m.b);
            double ab = m.a.dot(m.b)*0.5;
            double ac = m.a.dot(m.b)*0.5;
            double ad = m.a.dot(m.b)*0.5;
            //double ba = m.b.dot(m.b);
            //double bb = m.b.dot(m.b);
            double bc = m.b.dot(m.b)*0.5;
            double bd = m.b.dot(m.b)*0.5;
            //double ca = m.c.dot(m.b);
            //double cb = m.c.dot(m.b);
            //double cc = m.c.dot(m.b);
            double cd = m.c.dot(m.b)*0.5;
            //double da = m.b.dot(m.b);
            //double db = m.b.dot(m.b);
            //double dc = m.b.dot(m.b);
            //double dd = m.c.dot(m.b);
            m.a = m.a - m.b*ab - m.c*ac - m.d*ad;
            m.b = m.b - m.a*ab - m.c*bc - m.d*bd;
            m.c = m.c - m.a*ac - m.b*bc - m.d*cd;
            m.d = m.d - m.a*ad - m.b*bd - m.c*cd;
            m.a.normalize_taylor3();
            m.b.normalize_taylor3();
            m.c.normalize_taylor3();
            m.d.normalize_taylor3();
            //outproject( const QUAT& q );
        }
    }


    void forceEE(){
        double invSee = 1/(see*see);
        for(int i=0; i<ne; i++){
            Vec3d pi = epos[i];
            for(int j=0; j<i; j++){
                Vec3d  d = epos[j] - pi;
                double r2 = d.norm2();
                double r  = sqrt(r2);
                //fr       = Aee*exp(Bee*r);
                //double fr = frEFF_r2( r*r, 1.0, invSee );
                double fr = 4/(r*r + se2);
                //double fr = 4*interp_gx4(r2, 4.0*r-2.0, 1.0/r2 );
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
                double r2 = d.norm2();
                double r  = sqrt(r2);
                //double fr = qi* frEFF_r2( r*r, 1.0, invSae );
                //double fr = 2*qi*(  -2/(r*r + se2) + 1/(r*r + sa2) );
                //double fr = 2*qi*(  -2/(r*r + sea2) + 1/(r*r + sa2) );
                double fr = 2*qi*( -1/(r*r + 0.15) + 10*exp(-8*r*r) );
                //printf(" fae %i %i %g %g \n", i, j, r, fr);
                //double fr = interp_gx4(r2*4.0, 14.0, -1.0/r2 );

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
                double r2 = d.norm2();
                double r  = sqrt(r2);
                //double fr = aQ[j]*qi * frEFF_r2( r*r, 1.0, invSaa );
                double fr = aQ[j]*qi/(r*r + sa2);
                //double fr = aQ[j]*qi * interp_gx4(r2, 0, 1/r2 );
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

    double getSafeTimeStep(double dt_, double dpmax, double dvmax ){
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
            force_e_pairs(10.0);
            double dt_ = getSafeTimeStep(dt, 0.1, 10.0 );
            move(dt_, damp);
            //constrain_e_pairs();
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
