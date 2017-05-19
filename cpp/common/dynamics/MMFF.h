
#ifndef MMFF_h
#define MMFF_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"

// forcefield parameters

class MMFF{
    public:
    int  natoms=0, nbonds=0, nang=0, ntors=0;

    Vec2i  * bond2atom = NULL;
    double * bond_0    = NULL;  // [A]
    double * bond_k    = NULL;  // [eV/A] ?

    Vec2i  * ang2bond  = NULL;
    Vec3i  * ang2atom  = NULL;
    Vec2d  * ang_0     = NULL; // [1]
    double * ang_k     = NULL; // [eV/A^2]

    Vec3i  * tors2bond = NULL;
    Quat4i * tors2atom = NULL;
    Vec2d  * tors_0    = NULL; // [1]
    double * tors_k    = NULL; // [eV/A^2]

    // molecular gemeotry
    Vec3d  * apos   = NULL;   // atomic position
    double * lbond  = NULL;  // bond lengths
    Vec3d  * hbond  = NULL;   // normalized bond unitary vectors
    Vec3d  * aforce = NULL;


void allocate( int natoms_, int nbonds_, int nang_, int ntors_ ){
    natoms=natoms_; nbonds=nbonds_; nang=nang_; ntors=ntors_;
    if(apos     ==NULL) apos      = new Vec3d [natoms];
    if(aforce   ==NULL) aforce    = new Vec3d [natoms];

    if(lbond    ==NULL) lbond     = new double[nbonds];
    if(hbond    ==NULL) hbond     = new Vec3d [nbonds];

    if(bond2atom==NULL) bond2atom = new Vec2i [nbonds];
    if(bond_0   ==NULL) bond_0    = new double[nbonds];
    if(bond_k   ==NULL) bond_k    = new double[nbonds];

    if(ang2bond ==NULL) ang2bond  = new Vec2i [nang];
    if(ang2atom ==NULL) ang2atom  = new Vec3i [nang];
    if(ang_0    ==NULL) ang_0     = new Vec2d [nang];
    if(ang_k    ==NULL) ang_k     = new double[nang];

    /*
    tors2bond = new Vec3i [ntors];
    tors2atom = new Quat4i[ntors];
    tors_0    = new double[ntors];
    tors_k    = new double[ntors];
    */
}


void ang_b2a(){
    for(int i=0; i<nang; i++){
        Vec2i ib = ang2bond[i];
        int a00 = bond2atom[ib.x].x;
        int a01 = bond2atom[ib.x].y;
        int a10 = bond2atom[ib.y].x;
        int a11 = bond2atom[ib.y].y;
        /*
        if     ( a00 == a10 ){ ang2atom [i] = (Vec3i){ a01, a11, a00 }; }
        else if( a00 == a11 ){ ang2atom [i] = (Vec3i){ a01, a10, a00 }; }
        else if( a01 == a10 ){ ang2atom [i] = (Vec3i){ a00, a11, a01 }; }
        else if( a01 == a11 ){ ang2atom [i] = (Vec3i){ a00, a01, a01 }; }
        */
        if     ( a00 == a10 ){ ang2atom [i] = (Vec3i){ a01, a11, a00 }; }
        else if( a00 == a11 ){ ang2atom [i] = (Vec3i){ a01, a10, a00 }; ang2bond[i].x *=-1; }
        else if( a01 == a10 ){ ang2atom [i] = (Vec3i){ a00, a11, a01 }; ang2bond[i].y *=-1; }
        else if( a01 == a11 ){ ang2atom [i] = (Vec3i){ a00, a01, a01 }; ang2bond[i].x*=-1; ang2bond[i].y*=-1; }
    }
}


void eval_bonds(){
    for(int ib=0; ib<nbonds; ib++){
        Vec2i iat = bond2atom[ib];
        Vec3d dp; dp.set_sub( apos[iat.y], apos[iat.x] );
        //printf( "%i %i (%g,%g,%g) (%g,%g,%g) \n", iat.x, iat.y, apos[iat.x].x, apos[iat.x].y, apos[iat.x].z, apos[iat.y].x, apos[iat.y].y, apos[iat.y].z );
        double l = dp.norm();
        dp.mul(1.0d/l);
        //printf( " %i (%i,%i) (%g,%g,%g) %g \n", ib, iat.x, iat.y, dp.x, dp.y, dp.z, l );

        lbond [ib] = l;
        hbond [ib] = dp;

        dp.mul( (l-bond_0[ib])*bond_k[ib] );

        aforce[iat.x].add( dp );
        aforce[iat.y].sub( dp );
    }
    //exit(0);
}

void eval_angcos(){
    // simpler and faster angular forces cos(theta)
    // v      = -k*cos(a234) = k*dot(h23,h34)
    // dv/dr2 = ( h34 - dot(h23,h34)*h23 )/r23
    // dv/dr4 = ( dot(h23,h34)*h34 - h23 )/r34
    for(int ig=0; ig<nang; ig++){
        Vec2i ib = ang2bond[ig];
        //Vec3d h1 = hbond[ib.x];
        //Vec3d h2 = hbond[ib.y];
        Vec3d h1,h2;
        if(ib.x<0){ ib.x=-ib.x; h1 = hbond[ib.x]; h1.mul(-1.0d); }else{ h1 = hbond[ib.x]; };
        if(ib.y<0){ ib.y=-ib.y; h2 = hbond[ib.y]; h2.mul(-1.0d); }else{ h2 = hbond[ib.y]; };



        double c = h1.dot(h2);

        //printf( " %i (%i,%i) (%g,%g,%g) (%g,%g,%g) \n", ib, ib.x,ib.y, h1.x,h1.y,h1.z, h2.x,h2.y,h2.z );

        Vec3d hf1,hf2; // unitary vectors of force
        hf1 = h2 - h1*c;
        hf2 = h1 - h2*c;

        Vec3i ia = ang2atom[ig];

        double fang = ang_k[ig] * c;

        hf1.mul( fang/lbond[ib.x] );
        hf2.mul( fang/lbond[ib.y] );

        glColor3f(0.0f,1.0f,0.0f);
        Draw3D::drawVecInPos( h1, apos[ia.z] );
        Draw3D::drawVecInPos( h2, apos[ia.z] );

        glColor3f(1.0f,0.0f,0.0f);
        Draw3D::drawVecInPos( hf1, apos[ia.x] );
        Draw3D::drawVecInPos( hf2, apos[ia.y] );

        aforce[ia.x].add( hf1     );
        aforce[ia.y].add( hf2     );
        aforce[ia.z].sub( hf1+hf2 );
        //aforce[ia.x] ????
        // TODO : zero moment condition
    }
}

void eval_angles(){
    for(int ig=0; ig<nang; ig++){
        Vec2i ib = ang2bond[ig];
        Vec3d h1 = hbond[ib.x];
        Vec3d h2 = hbond[ib.y];
        Vec2d cs;
        cs.x = h1.dot(h2);
        cs.y = sqrt(1-cs.x*cs.x);

        Vec3d hf1,hf2; // unitary vectors of force
        double inv_sa = 1.0d/cs.y;
        hf1 = (h2 - h1*cs.x); // *inv_sa -> moved to fang for faster evaluation
        hf2 = (h1 - h2*cs.x);

        cs.mul_cmplx( ang_0[ig] );
        //E = 0.5*ang_k[ig]*cs.x*cs.x;
        double fang = (ang_k[ig] * cs.x * cs.y)*inv_sa;

        Vec3i ia = ang2atom[ig];

        hf1.mul( fang/lbond[ib.x] );
        hf2.mul( fang/lbond[ib.y] );

        aforce[ia.y].add( hf1     );
        aforce[ia.z].add( hf2     );
        aforce[ia.x].sub( hf1+hf2 );
        // TODO : check zero moment condition ... probably OK since pivot atom is in r=0 (COG, no torque by definition)
    }
}

// torsions should be special - for whole group - more atoms
void eval_torsion(){
    for(int it=0; it<ntors; it++){
        Vec3i ib = tors2bond[it];
        Vec3d h1 = hbond[ib.x];
        Vec3d h2 = hbond[ib.y];
        Vec3d ax = hbond[ib.z];

        // --- outproject axis
        h1.add_mul( ax, ax.dot(h1) );
        h2.add_mul( ax, ax.dot(h2) );

        Vec2d cs;
        cs.x = h1.dot(h2);
        cs.y = sqrt(1-cs.x*cs.x);

        double inv_sa = 1.0d/cs.y;
        Vec3d hf1,hf2; // unitary vectors of force
        hf1 = (h1 - h2*cs.x)*inv_sa;
        hf2 = (h2 - h1*cs.x)*inv_sa;

        cs.mul_cmplx( tors_0[it] );

        // TODO : torsion should be rather some fourier series or spline ? not just one rotation angle (?)
        double fang = ang_k[it] * cs.x * cs.y;

        Quat4i ia = tors2atom[it];

        hf1.mul( fang/lbond[ib.x] );
        hf2.mul( fang/lbond[ib.y] );

        aforce[ia.y].add( hf1 );
        aforce[ia.z].add( hf2 );
        aforce[ia.x].sub( hf1 ); // is this correct ?
        aforce[ia.w].sub( hf2 ); // is this correct ?
        //aforce[ia.x] ????
        // TODO : zero moment condition
    }
}

};

#endif
