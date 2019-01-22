
#ifndef RARFFarr_h
#define RARFFarr_h

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


Simple reactive force-field for mix of sp2 and sp3 hybridized atoms. Energy is based on Morse potential where the attractive term is multiplied by product of cosines between oriented dangling-bonds  sticking out  of atoms (the white sticks). 

More specifically for two atoms with positions pi,pj (dij=pi-pj, rij=|dij|, hij = dij/rij )

energy is computed as 
E = exp( -2a(rij-r0)) - 2 K exp( -a(rij-r0) )

where rij is distance  K = (ci.cj.cj) where ci=dot(hi,hij), cj=dot(hj,hij), cij=dot(hi,hj) and where hi, hj are normalized vector in direction of bonds

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

    double rbond0 =  0.5;
    double aMorse =  4.0;
    double bMorse = -0.7;

    double c6    = -15.0;
    double R2vdW =  8.0;
    Vec3d* bh0s = (Vec3d*)sp3_hs;

    inline void combine(const RigidAtomType& a, const RigidAtomType& b ){
        nbond   = a.nbond;
        aMorse  = a.aMorse * b.aMorse;  // TODO
        bMorse  = a.bMorse + b.bMorse;
        rbond0  = a.rbond0 + b.rbond0;

        c6    = -(a.c6  * b.c6);
        R2vdW = a.R2vdW + a.R2vdW;
    }

    void print(){
        printf( "nbond  %i rbond0 %g\n", nbond, rbond0 );
        printf( "aMorse %g bMorse %g\n", aMorse, bMorse );
        printf( "c6     %g r2vdW  %g\n", c6, R2vdW );
        //exit(0);
    }
};

/*
struct RigidAtom{
    RigidAtomType* type = 0;
    Vec3d  pos;
    Quat4d qrot;

    Vec3d force;
    Vec3d torq;

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

};
*/


class RARFF2arr{ public:

    double invRotMass = 2.0;

    int natom            =0;
    //RigidAtomType* types =0;
    //RigidAtom*     atoms =0;

    // atom properties
    RigidAtomType** types = 0;
    Vec3d*  poss   = 0;
    Quat4d* qrots  = 0;
    Vec3d*  forces = 0;
    Vec3d*  torqs  = 0;
    Vec3d*  omegas = 0;
    Vec3d*  vels   = 0;

    // aux
    double*        ebonds=0;
    Vec3d*         hbonds=0;
    Vec3d*         fbonds=0;

    //double F2pos=0;
    //double F2rot=0;

    void realloc(int natom_){
        natom=natom_;
        //_realloc(atoms,natom   );
        _realloc(types  ,natom);
        _realloc(poss   ,natom);
        _realloc(qrots  ,natom);
        _realloc(forces ,natom);
        _realloc(torqs  ,natom);
        _realloc(omegas ,natom);
        _realloc(vels   ,natom);
        _realloc(ebonds ,natom*N_BOND_MAX);
        _realloc(hbonds ,natom*N_BOND_MAX);
        _realloc(fbonds ,natom*N_BOND_MAX);
    }

    inline double pairEF( int ia, int ja, int nbi, int nbj, RigidAtomType& type){

        Vec3d  dij = poss[ja] - poss[ia];
        double r2  = dij.norm2() + R2SAFE;
        double rij = sqrt( r2 );
        Vec3d hij  = dij*(1/rij);

        double ir3 = 1/(rij*r2);
        double xx  = dij.x*dij.x;
        double yy  = dij.y*dij.y;
        double zz  = dij.z*dij.z;

        double dxy=-dij.x*dij.y*ir3;
        double dxz=-dij.x*dij.z*ir3;
        double dyz=-dij.y*dij.z*ir3;
        double dxx=(yy+zz)*ir3;
        double dyy=(xx+zz)*ir3;
        double dzz=(xx+yy)*ir3;

        double ir2vdW = 1/(r2 + type.R2vdW);
        double evdW   =  type.c6*ir2vdW*ir2vdW*ir2vdW;

        double expar = exp( type.bMorse*(rij-type.rbond0) );
        double E     =    type.aMorse*expar*expar;
        double Eb    = -2*type.aMorse*expar;

        E +=evdW; 
        double fr    =  2*type.bMorse* E +   6*evdW*ir2vdW;
        double frb   =    type.bMorse* Eb;

        Vec3d force=hij*fr;

        Vec3d* his = hbonds+ia*N_BOND_MAX;
        Vec3d* hjs = hbonds+ja*N_BOND_MAX;
        Vec3d* fis = fbonds+ia*N_BOND_MAX;
        Vec3d* fjs = fbonds+ja*N_BOND_MAX;

        double* eis = ebonds+ia*N_BOND_MAX;
        double* ejs = ebonds+ja*N_BOND_MAX;

        for(int ib=0; ib<nbi; ib++){
            const Vec3d& hi = his[ib];
            Vec3d& fi       = fis[ib];
            double ci       = hij.dot( hi );   // ci = <hi|hij>
            if(ci<0) continue;

            for(int jb=0; jb<nbj; jb++){
                const Vec3d& hj = hjs[jb];
                double cj       = hij.dot( hj );  // cj  = <hj|hij>
                double cij      = hi .dot( hj );  // cij = <hj|hi>

                if( (cj>0)||(cij>0) ) continue;

                Vec3d& fj = fjs[jb];

                double cc  = ci*cj*cij;
                double cc2 = cc*cc;
                double e   = cc2*cc2;
                double de  = 4*cc2*cc;

                double eEb = e * Eb;
                eis[ib]+=eEb*0.5;
                ejs[jb]+=eEb*0.5;
                E += eEb;

                force.x += Eb*de*cij*( cj*( hi.x*dxx + hi.y*dxy + hi.z*dxz )     +    ci*( hj.x*dxx + hj.y*dxy + hj.z*dxz )   )  + hij.x*frb*e;
                force.y += Eb*de*cij*( cj*( hi.x*dxy + hi.y*dyy + hi.z*dyz )     +    ci*( hj.x*dxy + hj.y*dyy + hj.z*dyz )   )  + hij.y*frb*e;
                force.z += Eb*de*cij*( cj*( hi.x*dxz + hi.y*dyz + hi.z*dzz )     +    ci*( hj.x*dxz + hj.y*dyz + hj.z*dzz )   )  + hij.z*frb*e;

                fi.x -= ( cij*cj*hij.x + ci*cj*hj.x )*de*Eb;
                fi.y -= ( cij*cj*hij.y + ci*cj*hj.y )*de*Eb;
                fi.z -= ( cij*cj*hij.z + ci*cj*hj.z )*de*Eb;

                fj.x -= ( cij*ci*hij.x + ci*cj*hi.x )*de*Eb;
                fj.y -= ( cij*ci*hij.y + ci*cj*hi.y )*de*Eb;
                fj.z -= ( cij*ci*hij.z + ci*cj*hi.z )*de*Eb;

            }
        }
        forces[ja].sub(force);
        forces[ia].add(force);
        return E;
    }

    double interEF(){
        Vec3d bhs[N_BOND_MAX];
        double E = 0;
        for(int i=0; i<natom; i++){
            RigidAtomType& typei = *types[i];
            Vec3d           pi   = poss[i];
            for(int j=0; j<natom; j++){
                RigidAtomType pairType;
                pairType.combine( typei, *types[j] );
                E += pairEF( i, j, typei.nbond, types[j]->nbond, pairType );
            }
        }
        return E;
    }

    void cleanAux(){
        for(int i=0; i<natom; i++){  
            vels  [i]=Vec3dZero; 
            omegas[i]=Vec3dZero; 
            torqs [i]=Vec3dZero; 
            forces[i]=Vec3dZero;
        }
    }

    void cleanAtomForce(){ 
        for(int i=0; i<natom; i++){  
            //atoms[i].cleanForceTorq(); 
            torqs [i]=Vec3dZero; 
            forces[i]=Vec3dZero;
        }
        int nb   = natom*N_BOND_MAX;    for(int i=0; i<nb;   i++){ ebonds[i]=0; }
        int nval = natom*N_BOND_MAX*3;  for(int i=0; i<nval; i++){ ((double*)fbonds)[i]=0;}
    }

    double evalF2rot(){ double F2=0; for(int i=0; i<natom; i++){ F2+=torqs [i].norm2(); }; return F2; }
    double evalF2pos(){ double F2=0; for(int i=0; i<natom; i++){ F2+=forces[i].norm2(); }; return F2; }

    double projectBonds(){
        for(int i=0; i<natom; i++){
            rotateVectors( N_BOND_MAX, qrots[i], types[i]->bh0s, hbonds + i*N_BOND_MAX );
        }
    }

    void evalTorques(){
        for(int ia=0; ia<natom; ia++){
            //RigidAtom& atomi = atoms[ia];
            //Vec3d torq = Vec3dZero;
            int nbi =  types[ia]->nbond;
            for(int ib=0; ib<nbi; ib++){
                int io = 4*ia+ib;
                fbonds[io].makeOrthoU(hbonds[io]);
                //printf( "ia %i ib %i f(%g,%g,%g)\n", ia, ib,  fbonds[io].x,fbonds[io].y,fbonds[io].z );
                torqs[ia].add_cross( hbonds[io], fbonds[io] );
            }
        }
    }

    void move(double dt){
        for(int i=0; i<natom; i++){
            //atoms[i].moveRotGD(dt*invRotMass);
            //atoms[i].movePosGD(dt);
            qrots[i].dRot_exact( dt, torqs[i] );  qrots[i].normalize();          // costly => debug
            //dRot_taylor2( dt, torqs[i] );   qrots[i].normalize_taylor3();  // fast => op
            poss[i].add_mul(forces[i],dt);
        }
    }

    void moveMDdamp(double dt, double damp){
        for(int i=0; i<natom; i++){
            //atoms[i].moveMDdamp( dt, invRotMass, damp);
            vels  [i].set_lincomb( damp, vels  [i], dt,            forces[i] );
            omegas[i].set_lincomb( damp, omegas[i], dt*invRotMass, torqs [i] );
            poss  [i].add_mul    ( vels[i], dt    );
            qrots [i].dRot_exact ( dt,  omegas[i] );
            qrots [i].normalize();
        }
    }

};

#endif
