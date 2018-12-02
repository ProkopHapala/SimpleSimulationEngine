
#ifndef RARFF_h
#define RARFF_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Forces.h"

/*
Rigid Atom Reactive Force-field
===============================

Problem:  binding of several atoms to one center
Solution: bonds repel for opposite side (back-side)

Electron Pairs
 - To prevent many atoms

Passivation atoms
 - except "scafold atoms" (C,O,N) there are termination atoms (-H, Cl) wich r



Optimization:
    - We can separate rotation and movement of atoms, rotation should be much faster
    - This makes only sense if we use forcefield with cheap evaluation of torque for fixed distance
        - we can perhaps save sqrt() calculation


*/

#define N_BOND_MAX 4
#define R1SAFE    1e-8

static const double sp3_hs[] = {
-0.57735026919, -0.57735026919, -0.57735026919,
+0.57735026919, +0.57735026919, -0.57735026919,
-0.57735026919, +0.57735026919, +0.57735026919,
+0.57735026919, -0.57735026919, +0.57735026919
};

static const double sp2_hs[] = {
+1.000000000000, -0.00000000000,  0.00000000000,
-0.500000000000, +0.86602540378,  0.00000000000,
-0.500000000000, -0.86602540378,  0.00000000000,
 0.00000000000,   0.00000000000,  1.00000000000     // electron - can be repulsive
};

static const double sp1_hs[] = {
+1.000000000000,  0.00000000000,  0.00000000000,
-1.000000000000,  0.00000000000,  0.00000000000,
 0.000000000000,  1.00000000000,  0.00000000000,    // electron - can be repulsive
 0.000000000000,  0.00000000000,  1.00000000000     // electron - can be repulsive
};


inline void overlapFE(double r, double amp, double beta, double& e, double& fr ){
    //https://www.wolframalpha.com/input/?i=exp(b*x)*(1%2Bb*x%2B(b*x)%5E2%2F3)+derivative+by+x
    double x     = r*beta;
    double expar = amp*exp(-x);
    e  =  expar*(1 +   x + 0.33333333*x*x );
    fr = (expar*(6 + 5*x +            x*x )*beta*0.33333333)/r;
}

template<typename T>
void rotateVectors(int n, const Quat4TYPE<T>& qrot, Vec3TYPE<T>* h0s, Vec3TYPE<T>* hs ){
    Mat3TYPE<T> mrot;
    qrot.toMatrix(mrot);
    for( int j=0; j<n; j++ ){
        Vec3TYPE<T> h;
        //mrot.dot_to_T    ( h0s[j], h );
        mrot.dot_to      ( h0s[j], h );
        hs[j] = h;
        //ps[j].set_add_mul( pos, p_, r0 );
    }
}

struct RigidAtomType{
    int    nbond = 4;  // number bonds
    double rbond0 = 0.5;
    double acore =  4.0;
    double bcore = -0.7;
    double abond = -2.0;
    double bbond = -1.1;
    double c6    = -100.0;
    double R2vdW =  8.0;
    Vec3d* bh0s = (Vec3d*)sp3_hs;

    inline void combine(const RigidAtomType& a, const RigidAtomType& b ){
        nbond  = a.nbond;
        acore  = a.acore  * b.acore ;  // TODO
        bcore  = a.bcore  + b.bcore ;
        abond  = -(a.abond  * b.abond);
        bbond  = a.bbond  + b.bbond ;
        rbond0 = a.rbond0 + b.rbond0;

        c6    = -(a.c6  * b.c6);
        R2vdW = a.R2vdW + a.R2vdW;
    }

    void print(){
        printf( "nbond %i rbond0 %g\n", nbond, rbond0 );
        printf( "abond %g bbond  %g\n", abond, bbond );
        printf( "acore %g bcore  %g\n", acore, bcore );
        printf( "c6 %g r2vdW  %g\n", c6, R2vdW );
        //exit(0);
    }
};

struct RigidAtom{
    RigidAtomType* type = 0;
    Vec3d  pos;
    Quat4d qrot;
    //Quat4d frot;
    Vec3d torq;
    Vec3d force;

    Vec3d omega;
    Vec3d vel;

    inline void cleanAux(){  vel=Vec3dZero; omega=Vec3dZero; torq=Vec3dZero; force=Vec3dZero; };

    inline void cleanForceTorq(){ torq=Vec3dZero; force=Vec3dZero; };

    inline void setPose(const Vec3d& pos_, const Quat4d& qrot_ ){ pos=pos_; qrot=qrot_; };

    inline void moveRotGD(double dt){
        qrot.dRot_exact( dt, torq );  qrot.normalize();          // costly => debug
        //dRot_taylor2( dt, torq );   qrot.normalize_taylor3();  // fast => op
    }

    inline void movePosGD(double dt){ pos.add_mul(force,dt); }

    inline void moveMDdamp( double dt, double invRotMass, double damp ){

        vel  .set_lincomb( damp, vel,   dt,             force );
        omega.set_lincomb( damp, omega, dt*invRotMass,  torq  );

        pos.add_mul    ( vel, dt    );
        qrot.dRot_exact( dt,  omega );
        qrot.normalize();
    }

/*
    inline moveMDcos( double dt, double invRotMass, damp ){
        //calculate damp for the whole thing
        //double cvel = vel  .dot(force)/sqrt( vel.norm2()*vel.norm2() );
        //double crot = omega.dot(torq )/sqrt( vel.norm2()*vel.norm2() );
        vel  .set_lincomb( damp, dt,          , vel, force );
        omega.set_lincomb( damp, dt*invRotMass, vel, torq  );
        pos.add_mul(force,dt);
        qrot.dRot_exact( dt, torq );
        qrot.normalize();
    }
*/

};


class RARFF2{ public:

    double invRotMass = 2.0;

    int natom            =0;
    RigidAtomType* types =0;
    RigidAtom*     atoms =0;
    Vec3d*         hbonds=0;
    Vec3d*         fbonds=0;

    void realloc(int natom_){
        natom=natom_;
        _realloc(atoms,natom   );
        _realloc(hbonds,natom*4);
        _realloc(fbonds,natom*4);
    }

    inline double pairEF( int ia, int ja, int nbi, int nbj, RigidAtomType& type){
    //inline double pairEF( const Vec3d& dij, const RigidAtomType& type, Vec3d* bhs, Vec3d& force, Vec3d& torq ){

        //printf( "==================== pairEF : \n" );
        Vec3d  dij = atoms[ja].pos - atoms[ia].pos;
        double r2  = dij.norm2();
        double rij = sqrt( r2 + R2SAFE);
        Vec3d hij  = dij*(1/rij);

        //type.bcore  = -1.8;
        //type.acore  =  1.0;
        //type.rbond0 =  1.2;

        double expar = exp(type.bcore*(rij-type.rbond0) );
        double E     =    type.acore*expar*expar;
        double Eb    = -2*type.acore*expar;

        //E += Eb;
        //E = 0;


        //E += Eb;
        //double E = sq( 1-exp( -1.8*(rij-1.4) ) );
        //double Eb = 0;
        //return E;
        double fr    =  2*type.bcore* E ;
        double frb   =    type.bcore* Eb;

        Vec3d force=hij*fr;

        Vec3d* his = hbonds+ia*N_BOND_MAX;
        Vec3d* hjs = hbonds+ja*N_BOND_MAX;
        Vec3d* fis = fbonds+ia*N_BOND_MAX;
        Vec3d* fjs = fbonds+ja*N_BOND_MAX;

        //torq =Vec3dZero;
        //printf( "eij %g \n", eij );
        //printf( "hij (%g,%g,%g) \n", dij.x, dij.y, dij.z, hij.x, hij.y, hij.z );


        //printf("nbi %i nbj %i \n", nbi, nbj);
        for(int ib=0; ib<nbi; ib++){
        //for(int ib=0; ib<1; ib++){
            const Vec3d& hi = his[ib];
            Vec3d& fi       = fis[ib];
            double ci       = hij.dot( hi );   // ci = <hi|hij>
            //printf( "ci  %g \n", ci );
            if(ci<0) continue;

            for(int jb=0; jb<nbj; jb++){
            //for(int jb=0; jb<1; jb++){
                const Vec3d& hj = hjs[jb];
                double cj       = hij.dot( hj );  // cj  = <hj|hij>
                double cij      = hi .dot( hj );  // cij = <hj|hi>

                //printf( "  ia(%i,%i) ib(%i,%i) nb(%i,%i)  c i,j,ij %g %g %g \n", ia, ja, ib, jb,  nbi,nbj,   ci, cj, cij );

                if( (cj>0)||(cij>0) ) continue;

                Vec3d& fj = fjs[jb];

                double cc  = ci*cj*cij;
                double cc2 = cc*cc;
                double e   = cc2*cc2;
                double de  = 4*cc2*cc;

                //printf( "e %g Eb %g \n", e, Eb );

                E += e * Eb;

                // derivative by dij
                Vec3d f; // d_E/d_hij
                force.x = ( cj*hi.x + ci*hj.x )*cij*de*Eb +    hij.x*frb*e;
                force.y = ( cj*hi.y + ci*hj.y )*cij*de*Eb +    hij.y*frb*e;
                force.z = ( cj*hi.z + ci*hj.z )*cij*de*Eb +    hij.z*frb*e;

                fi.x = ( cij*cj*hij.x + ci*cj*hj.x )*de*Eb;
                fi.x = ( cij*cj*hij.y + ci*cj*hj.y )*de*Eb;
                fi.x = ( cij*cj*hij.z + ci*cj*hj.z )*de*Eb;

                fj.x = ( cij*ci*hij.x + ci*cj*hi.x )*de*Eb;
                fj.x = ( cij*ci*hij.y + ci*cj*hi.y )*de*Eb;
                fj.x = ( cij*ci*hij.z + ci*cj*hi.z )*de*Eb;

            }
        }


        //printf( "E %g \n", E );
        //exit(0);
        return E;
    }

    double interEF(){
        //Vec3d bps[N_BOND_MAX];
        Vec3d bhs[N_BOND_MAX];
        double E = 0;

        for(int i=0; i<natom; i++){
            RigidAtom&     atomi = atoms[i];
            RigidAtomType& typei = *atomi.type;
            //int            nbi   = typei.nbond;
            Vec3d           pi   = atomi.pos;
            for(int j=0; j<natom; j++){
                RigidAtom&    atomj = atoms[j];
                RigidAtomType pairType;
                pairType.combine( typei, *atomj.type );
                E += pairEF( i, j, atomi.type->nbond, atomj.type->nbond, pairType );
            }
        }
    }

    void cleanAtomForce(){ for(int i=0; i<natom; i++){ atoms[i].cleanForceTorq(); } }

    double projectBonds(){
        for(int i=0; i<natom; i++){
            rotateVectors( N_BOND_MAX, atoms[i].qrot, atoms[i].type->bh0s, hbonds + i*N_BOND_MAX );
        }
    }

    void getCos( double& cosdRot, double& cosdPos ){
        double crot=0, cpos=0, r2f=0, r2tq=0, r2vel=0, r2omg=0;
        for(int i=0; i<natom; i++){
            RigidAtom& atomi = atoms[i];
            crot  += atomi.torq .dot(atomi.omega);
            r2tq  += atomi.torq .norm2();
            r2omg += atomi.omega.norm2();

            cpos  += atomi.force.dot(atomi.vel);
            r2f   += atomi.force.norm2();
            r2vel += atomi.vel  .norm2();
        }
        cosdRot = crot/sqrt(r2tq*r2omg);
        cosdPos = cpos/sqrt(r2f *r2vel);
        //return  + ;
    }


    void move(double dt){
        for(int i=0; i<natom; i++){
            atoms[i].moveRotGD(dt*invRotMass);
            atoms[i].movePosGD(dt);
        }
    }

    void moveMDdamp(double dt, double damp){
        for(int i=0; i<natom; i++){
           atoms[i].moveMDdamp( dt, invRotMass, damp);
        }
    }

};



#endif
