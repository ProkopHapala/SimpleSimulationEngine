
#ifndef FTRFF_h
#define FTRFF_h

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

*/

#define N_BOND_MAX 4

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

template<typename T>
void rotateVectors(int n, const Quat4TYPE<T>& qrot, Vec3TYPE<T>* h0s, Vec3TYPE<T>* hs ){
    Mat3TYPE<T> mrot; 
    qrot.toMatrix(mrot);
    for( int j=0; j<n; j++ ){
        Vec3TYPE<T> h;
        mrot.dot_to_T    ( h0s[j], h );
        hs[j] = h;
        //ps[j].set_add_mul( pos, p_, r0 );
    }
}

struct RigidAtomType{
    int    nbond = 0;  // number bonds
    Vec3d* bh0s = (Vec3d*)sp3_hs; 
    double rbond0;
    double acore;
    double bcore;
    double abond; 
    double bbond;
};

struct RigidAtom{
    RigidAtomType* type = 0;
    Vec3d  pos;
    Quat4d qrot;
    //Quat4d frot;
    Vec3d torq;
    Vec3d force;

    inline void moveRotGD(double dt){   
        qrot.dRot_exact( dt, torq );  qrot.normalize();          // costly => debug
        //dRot_taylor2( dt, torq );   qrot.normalize_taylor3();  // fast => op
    }

    inline void movePosGD(double dt){ pos.add_mul(force,dt); }

};

class RARFF{ public:

    int natom    =0;
    RigidAtomType* types=0;
    RigidAtom*     atoms=0;

    void interEF(){
        //Vec3d bps[N_BOND_MAX];
        Vec3d bhs[N_BOND_MAX];

        for(int i=0; i<natom; i++){
            RigidAtom&     atomi = atoms[i];
            RigidAtomType& typei = *atomi.type;
            int            nbi   = typei.nbond;
            Vec3d           pi   = atomi.pos;

            rotateVectors<double>(N_BOND_MAX, atomi.qrot, typei.bh0s, bhs );

            for(int j=0; j<natom; j++){
                RigidAtom&     atomj = atoms[j];
                RigidAtomType& typej = *atomi.type;

                double acore = typei.acore * typej.acore ;  // TODO
                double bcore = typei.bcore + typej.bcore ;
                double abond = typei.abond * typej.abond ; 
                double bbond = typei.bbond + typej.bbond ;

                // core-core interactions
                Vec3d  dij = atomj.pos - pi;
                double rij = dij.norm();

                double eij = acore*exp(bcore*rij);
                Vec3d  fij = dij*( eij*bcore );

                atomi.force.sub(fij);
                atomj.force.add(fij);

                double rbond0 = typei.rbond0 + typej.rbond0;

                // bonds interactions
                for(int ib=0; ib<nbi; ib++){

                        Vec3d  db = bhs[ib]*rbond0;

                        Vec3d  d  = dij + db;
                        double r  = d.norm();

                        double e = abond*exp(bbond*r);
                        Vec3d f  = d*( e*bbond );

                        //atomi.fbond[ib].add_mul(f,r); // TODO: think how this is transfered to the atom ?
                        atomi.torq .add_cross(f,db);
                        atomi.force.add(f);
                        atomj.force.sub(f);

                }
            }

        }

    }

};



#endif
