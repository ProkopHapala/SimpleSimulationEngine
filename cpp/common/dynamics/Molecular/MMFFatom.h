
#ifndef MMFFatom_h
#define MMFFatom_h

#include "fastmath.h"
#include "Vec3.h"

#define max_bonds 4

// forcefield parameters

Vec2i  * bond2atom;
double * bond_0;  // [A]
double * bond_k;  // [eV/A] ?

Vec2i  * ang2bond;
Vec2i  * ang2atom;
Vec3d  * ang_0; // [1]
double * ang_k; // [eV/A^2]

Vec3i  * tors2bond;
Quat4i * tors2atom;
double * tors_0; // [1]
double * tors_k; // [eV/A^2]

#define NBONDMAX 4;

class AtomTopo{
    int nbond;
    int ibond[NBONDMAX];
}


// molecular gemeotry
Vec3d  * apos;   // atomic position
double * lbond;  // bond lengths
Vec3d  * hbond;   // normalized bond unitary vectors
Vec3d  * aforce; 


void eval_bonds(){
    for(int ib=0; ib<nbonds; ib++){
        Vec2i iat = bond2atom[i];
        //Vec3d p1 = apos[iat.x];
        //Vec3d p2 = apos[iat.y];
        Vec3d dp; dp.set_sub( apos[iat.y], apos[iat.x] );
        double l = dp.norm();
        dp.mul(1.0d/l);
        rbond [i] = l;
        hbond [i] = dp;
        dp.mul( (l-bond_0[i])*bond_k[0] );
        aforce[i].add( dp );  
        aforce[i].add( dp );  
    }
}



#endif 
