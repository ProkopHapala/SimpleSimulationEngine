
#ifndef MMFF_h
#define MMFF_h

#include "fastmath.h"
#include "Vec3.h"

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

void eval_acos(){
    // simpler and faster angular forces cos(theta)
    // v      = -k*cos(a234) = k*dot(h23,h34)
    // dv/dr2 = ( h34 - dot(h23,h34)*h23 )/r23 
    // dv/dr4 = ( dot(h23,h34)*h34 - h23 )/r34
    for(int it=0; it<ntors; it++){
        Vec2i ib = angle2bond[ig];
        Vec3d h1 = hbond[ib.x];
        Vec3d h2 = hbond[ib.y];
        double c = h1.dot(h2);
        
        Vec3d hf1,hf2; // unitary vectors of force
        hf1 = h1 - h2*c;
        hf2 = h2 - h1*c;

        double fang = ang_k[ig] * c;
        
        Vec3i ia = ang2atom[ig]; 
        
        hf1.mul( fang/lbond[ib.x] );
        hf2.mul( fand/lbond[ib.y] );
        
        aforce[ia.y].add( hf1     );  
        aforce[ia.z].add( hf2     );
        aforce[ia.x].sub( hf1+hf2 );
        //aforce[ia.x] ????
        // TODO : zero moment condition 
    }
}


void eval_angles(){
    for(int it=0; it<ntors; it++){
        Vec2i ib = angle2bond[ig];
        Vec3d h1 = hbond[ib.x];
        Vec3d h2 = hbond[ib.y];
        Vec2d cs;
        cs.x = h1.dot(h2);
        cs.y = sqrt(1-cs.x*cs.x);
        
        Vec3d hf1,hf2; // unitary vectors of force
        hf1 = (h1 - h2*cs.x)/cs.y;
        hf2 = (h2 - h1*cs.x)/cs.y;
        
        cs.mul_cmplx( ang_0[ig] );
        //E = 0.5*ang_k[ig]*cs.x*cs.x;
        double fang = ang_k[ig] * cs.x * cs.y;
        
        Vec3i ia = ang2atom[ig]; 
        
        hf1.mul( fang/lbond[ib.x] );
        hf2.mul( fand/lbond[ib.y] );
        
        aforce[ia.y].add( hf1     );  
        aforce[ia.z].add( hf2     );
        aforce[ia.x].sub( hf1+hf2 );
        // TODO : check zero moment condition ... probably OK since pivot atom is in r=0 (COG, no torque by definition)
    }
}

// torsions should be special - for whole group - more atoms
void eval_torsion(){
    for(int it=0; it<ntors; it++){
        Vec3i ib = angle2bond[ig];
        Vec3d h1 = hbond[ib.x];
        Vec3d h2 = hbond[ib.y];
        Vec3d ax = hbond[ib.z];
        
        // --- outproject axis
        Vec3d h1.add_mul( ax, ax.dot(h1) );
        Vec3d h2.add_mul( ax, ax.dot(h2) );
        
        Vec2d cs;
        cs.x = h1.dot(h2);
        cs.y = sqrt(1-cs.x*cs.x);
        
        Vec3d hf1,hf2; // unitary vectors of force
        hf1 = (h1 - h2*cs.x)/cs.y;
        hf2 = (h2 - h1*cs.x)/cs.y;
        
        cs.mul_cmplx( ang_0[ig] );
        
        // TODO : torsion should be rather some fourier series or spline ? not just one rotation angle (?)
        double fang = ang_k[ig] * cs.x * cs.y;
        
        Vec3i ia = ang2atom[ig]; 
        
        hf1.mul( fang/lbond[ib.x] );
        hf2.mul( fand/lbond[ib.y] );
        
        aforce[ia.y].add( hf1 );  
        aforce[ia.z].add( hf2 );
        aforce[ia.x].sub( hf1 ); // is this correct ?
        aforce[ia.w].sub( hf2 ); // is this correct ?
        //aforce[ia.x] ????
        // TODO : zero moment condition 
    }
}

#endif 
