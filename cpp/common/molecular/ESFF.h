
#ifndef ESFF_h
#define ESFF_h


/*
Electron Spin Foce Field
*/

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Forces.h"
//#include "GridFF.h"

#define N_ATOM_ORB 4

constexpr double cCoulomb =  14.399644;

struct ESFFelement{
    double aaa; // amplitude atom-atom interaction
    double baa; // decay     atom-atom interaction

    double aea; // amplitude electron-atom interaction
    double bea; // decay     electron-atom interaction

    double q;
};

class ESFF{ public:

    double aee = 1.0; // decay of electron-election interaction
    double bee = 1.0; // amplitude of electron-election interaction
    //double cee = 1.0; // coulomb electron-election interaction strenght
    double wee = 1.0;

    ESFFelement * elems = 0;

    int nelec=0,natom=0;

    int*    aZs  = 0;  // ion   chargs
    Vec3d*  aps  = 0;      // atom position
    Vec3d*  faps = 0;     // atom forces

    double* ess  = 0;      // electron spins
    Vec3d*  eps  = 0;      // eletron  positions

    double* fess = 0;     // electron spin forces
    Vec3d*  feps = 0;     // eletron positions forces

    // ion-ion interaction
    double aaEF(){
        double Eaa = 0.0;
        for(int i=0; i<natom; i++ ){
            int iZ = aZs[i];
            ESFFelement& elemi = elems[iZ];
            double aea = elemi.aea;
            double bea = elemi.bea;
            double qi  = elemi.q;
            Vec3d  pi  = aps[i];

            for(int j=i; j<nelec; j++ ){
                int jZ = aZs[j];
                ESFFelement& elemj = elems[jZ];
                Vec3d  d  = aps[j] - pi;
                double r  = d.norm();
                double ir   = 1/r;
                double b    = (bea+elemj.bea);
                double es = aea*elemj.aea*exp(b*r);  // Pauli repulsion of electrons ( for opposite sins is zero )
                double ec = qi*elemj.q*cCoulomb*ir;     // Coulombic attraction of electron clouds

                Eaa      += es + ec;
                Vec3d  f  = d*( es*ir*b + ec*ir*ir );

                faps[i].sub( f );
                feps[j].add( f );

            }
        }
    }

    // electorn-ion interaction
    double eaEF(){
        double Eae = 0.0;
        for(int ia=0; ia<natom; ia++ ){
            int iZ = aZs[ia];
            ESFFelement& elem = elems[iZ];
            double aea = elem.aea;
            double bea = elem.bea;
            double qi  = elem.bea;
            Vec3d  pi  = aps[ia];

            for(int je=0; je<nelec; je++ ){
                Vec3d  d  = eps[je] - pi;
                double r  = d.norm();
                double ir = 1/r;

                double es = aea*exp(bea*r);  // Pauli repulsion of electrons ( for opposite sins is zero )
                double ec = qi*cCoulomb*ir;     // Coulombic attraction of electron clouds

                Eae      += es + ec;
                Vec3d  f  = d*( es*ir*bea + ec*ir*ir );
                faps[ia].sub( f );
                feps[je].add( f );
            }
        }
    }

    // electron-electron interaction
    double eeEF(){
        double Eee = 0.0;
        for(int i=0; i<nelec; i++ ){
            Vec3d  pi = eps[i];
            double si = ess[i];
            for(int j=i; j<nelec; j++ ){
                Vec3d d   = eps[j] - pi;
                double sj = ess[j];
                double r  = d.norm();
                double ir = 1/r;
                double S  = aee*exp(bee*r);  // overlap between electron clouds
                double es = (1+si*sj)*S;       // Pauli repulsion of electrons ( for opposite sins is zero )
                double ec = cCoulomb*ir;       // Coulombic repulsion of electron clouds
                Eee      += es+ec;
                Vec3d f   = d*( es*ir*bee + ec*ir*ir );

                feps[i].sub( f );
                feps[j].add( f );
                fess[i] = -S*sj;
                fess[j] = -S*si;
            }
        }
    }


    double moveSpinGD (double dt){
        for(int i=0; i<nelec; i++){
            double s = ess[i] + fess[i]*dt;
            if(s<-1.0) s=-1.0;
            if(s>+1.0) s=+1.0;
            ess[i]=s;
        }
    }

    double moveElecGD (double dt){ for(int i=0; i<nelec; i++){ eps[i].add_mul(feps[i],dt); } }
    double moveAtomsGD(double dt){ for(int i=0; i<natom; i++){ aps[i].add_mul(faps[i],dt); } }


};



#endif
