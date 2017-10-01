
#ifndef MMFF_h
#define MMFF_h

#include <vector>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Forces.h"
#include "GridFF.h"

#include "integerOps.h"


#define SIGN_MASK 2147483648

//inline int invIndex( int i ){ return i^SIGN_MASK; }

// compied from RigidMolecule.h      ============>   TODO : make common lib of inter-molecular forcefields formulas

int pickParticle( int n, Vec3d * ps, Vec3d& ray0, Vec3d& hRay, double R ){
    double tmin =  1e+300;
    int imin    = -1;
    for(int i=0; i<n; i++){
        double ti = raySphere( ray0, hRay, R, ps[i] );
        if(ti<tmin){ imin=i; tmin=ti; }
    }
    return imin;
}

class MolType{ public:
    int nAtoms;
    Vec3d * apos   = NULL;
    int   * atypes = NULL;
    Vec3d * REQs   = NULL;
};

// ======================
// ====   MMFF
// ======================

class MMFF{ public:
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
    int    * atypes = NULL;
    Vec3d  * apos   = NULL;   // atomic position
    //Vec3d  * aLJq   = NULL;
    Vec3d  * aREQ   = NULL;
    Vec3d  * aPLQ   = NULL;   // this is used for grid-accelerated factorized potential
    double * lbond  = NULL;   // bond lengths
    Vec3d  * hbond  = NULL;   // normalized bond unitary vectors
    Vec3d  * aforce = NULL;

    std::vector<MolType> molTypes;

    int nFrag;  // rigid fragments
    //int *    imolTypes = NULL;
    Vec3d  ** fapos0s = NULL; // position of atoms in rigid fragment
    int    *  frag2a  = NULL; // start of the fragment in forcefield
    int    *  fragNa  = NULL; // lengh of the fragment
    double *  poses   = NULL; // rigd body pose of molecule (pos,qRot);
    double *  poseFs  = NULL; // rigd body pose of molecule (pos,qRot);

    //RigidSubstrate substrate;
    GridFF gridFF;


void allocate( int natoms_, int nbonds_, int nang_, int ntors_ ){
    natoms=natoms_; nbonds=nbonds_; nang=nang_; ntors=ntors_;
    printf( "MMFF::allocate natoms: %i  nbonds: %i  nang: %i ntors: %i \n", natoms, nbonds, nang, ntors );
    if(atypes   ==NULL) atypes    = new int   [natoms];
    if(apos     ==NULL) apos      = new Vec3d [natoms];
    if(aforce   ==NULL) aforce    = new Vec3d [natoms];
    if(aREQ     ==NULL) aREQ      = new Vec3d [natoms];

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

void allocFragment( int nFrag_ ){
    nFrag = nFrag_;
    printf( "MMFF::allocFragment nFrags: %i  nPosses: %i \n", nFrag, nFrag*8 );
    //imolTypes = new int[nFrag];
    frag2a    = new int   [nFrag];    // start of the fragment in forcefield
    fragNa    = new int   [nFrag];    // lengh of the fragment
    poses     = new double[nFrag*8];  // rigd body pose of molecule (pos,qRot);
    poseFs    = new double[nFrag*8];  // rigd body pose of molecule (pos,qRot);
    fapos0s   = new Vec3d*[nFrag];
}

void deallocate(){
    delete [] apos;     delete [] aforce;   delete [] aREQ;
    delete [] lbond;    delete [] hbond;    delete [] bond2atom; delete [] bond_0; delete [] bond_k;
    delete [] ang2bond; delete [] ang2atom; delete [] ang_0;     delete [] ang_k;
    if(aPLQ) delete [] aPLQ;
}

int pickBond( Vec3d& ray0, Vec3d& hRay, double R ){
    double dist_min =  R;
    int    imin     = -1;
    for(int ib=0; ib<nbonds; ib++){
        Vec2i iat = bond2atom[ib];
        double t1,t2;
        double dist = rayLine( ray0, hRay, apos[iat.x], hbond[ib], t1, t2 );
        //glColor3f(0.0f,0.0f,0.0f); Draw3D::drawLine( ray0+(hRay*t1), apos[iat.x]+(hbond[ib]*t2) );
        //printf( "%i %g %g %g \n", ib, t1, t2, dist );
        if( (dist<dist_min) && (t2>0) && (t2<lbond[ib]) ){
            imin=ib; dist_min=dist;
        }
    }
    return imin;
}

void genPLQ(){
    if(aPLQ==NULL) aPLQ = new Vec3d[natoms];
    for(int i=0; i<natoms; i++){
        aPLQ[i] = REQ2PLQ( aREQ[i], gridFF.alpha );
        printf( "genPLQ %i (%g,%g,%g)->(%g,%g,%g) \n", i, aREQ[i].x, aREQ[i].y, aREQ[i].z,   aPLQ[i].x, aPLQ[i].y, aPLQ[i].z );
    }
    //exit(0);
}

void translate( Vec3d dpos){ for(int i=0; i<natoms; i++) apos[i].add(dpos); }

void frags2atoms(){
    // : fragment pose -> atomic position
    //printf( ">>> MMFF::frags2atoms %i \n", nFrag );
    for(int ifrag=0; ifrag<nFrag; ifrag++){
        int im8 = ifrag<<3;
        Vec3d  pos = *((Vec3d* )(poses+im8  ));
        Quat4d rot = *((Quat4d*)(poses+im8+4));
        Mat3d T; rot.toMatrix(T);
        int ia = frag2a[ifrag];
        int na = fragNa[ifrag];
        //printf( "ia na %i %i \n", ia, na );
        //MolType& mtyp = mTypes[imolTypes[i]];
        //Vec3d * m_apos = mTypes[imolTypes[ifrag]].apos;
        Vec3d * m_apos = fapos0s[ifrag];
        for( int j=0; j<na; j++ ){
            Vec3d Tp;
            T.dot_to_T( m_apos[j], Tp );
            apos[ia].set_add( pos, Tp );
            //printf( "%i %i  (%g,%g,%g) (%g,%g,%g) \n", ifrag, j,  m_apos[j].x, m_apos[j].y, m_apos[j].z,   Tp.x, Tp.y, Tp.z  );
            ia++;
        }
    }
    //exit(0);
}

void cleanPoseTemps(){
    for( int ifrag=0; ifrag<nFrag; ifrag++ ){
        int im8 = ifrag<<3;
        double * poseF_i = poseFs+im8;
        ((Quat4d*)(poseF_i  ))->set(0.0);
        ((Quat4d*)(poseF_i+4))->set(0.0);
    }
}

void checkPoseNormal(){
    for( int ifrag=0; ifrag<nFrag; ifrag++ ){
        int im8 = ifrag<<3;
        double * pose_i = poseFs+im8;
        ((Quat4d*)(pose_i+4))->checkNormalized(1e-4);
    }
}

void aforce2frags(){
    // : atomic force -> force on fragment pose
    for(int ifrag=0; ifrag<nFrag; ifrag++){
        int ia = frag2a[ifrag];
        int na = fragNa[ifrag];
        Vec3d * m_apos   = fapos0s[ifrag];
        int im8 = ifrag<<3;
        double * pose_i  = poses +im8;
        double * poseF_i = poseFs+im8;
        for( int j=0; j<na; j++ ){
            ((Quat4d*)(pose_i+4))->addForceFromPoint( m_apos[j], aforce[ia], *((Quat4d*)(poseF_i+4)) );
            ((Vec3d *)(poseF_i)) ->add( aforce[ia] );
            ia++;
        }
    }
}

void RBodyForce(){
    // : Force between two rigid body molecules
    for(int ifrag=0; ifrag<nFrag; ifrag++){
        int ia = frag2a[ifrag];
        int ni = fragNa[ifrag];
        Vec3d * m_apos = fapos0s[ifrag];
        for(int iatom=0; iatom<ni; iatom++){
            int ig = ifrag+iatom;
            Vec3d REQi = aREQ[ig];
            Vec3d pi   = apos[ig];
            Vec3d f; f.set(0.0f);
            for(int jfrag=0; jfrag<nFrag; jfrag++){
                int ja = frag2a[ifrag];
                int nj = fragNa[ifrag];
                Vec3d * m_apos = fapos0s[ifrag];
                if(jfrag==ifrag){ continue; }
                for(int jatom=0; jatom<nj; jatom++){
                    int jg = jfrag+jatom;
                    Vec3d REQj = aREQ[jg];
                    Vec3d pj   = apos[jg];
                    addAtomicForceLJQ( pj-pi, f, REQi.x+REQj.x, REQi.y*REQj.y, REQi.z*REQj.z );
                }
            }
            //aforce[ig] = f; // Not necessary
            int im8 = ifrag<<3;
            double * pose_i  = poses +im8;
            double * poseF_i = poseFs+im8;
            Vec3d  * m_apos  = fapos0s[ifrag];
            ((Quat4d*)(pose_i+4))->addForceFromPoint( m_apos[iatom], f, *((Quat4d*)(poseF_i+4)) );
            ((Vec3d *)(poseF_i)) ->add( f );
        }
    }
    //exit(0);
}



void ang_b2a(){
    for(int i=0; i<nang; i++){
        Vec2i ib = ang2bond[i];
        Vec2i b1,b2;
        b1 = bond2atom[ib.x];
        b2 = bond2atom[ib.y];
        if     ( b1.x == b2.x ){ ang2atom[i] = (Vec3i){ b1.y, b2.y, b1.x }; }
        else if( b1.x == b2.y ){ ang2atom[i] = (Vec3i){ b1.y, b2.x, b1.x }; ang2bond[i].y|=SIGN_MASK; }
        else if( b1.y == b2.x ){ ang2atom[i] = (Vec3i){ b1.x, b2.y, b1.y }; ang2bond[i].x|=SIGN_MASK; }
        else if( b1.y == b2.y ){ ang2atom[i] = (Vec3i){ b1.x, b2.x, b1.y }; ang2bond[i].x|=SIGN_MASK; ang2bond[i].y|=SIGN_MASK; }
        //printf( "%i (%i,%i)(%i,%i) (%i,%i,%i) (%i,%i) (%i,%i) \n", i, b1.x+1,b1.y+1, b2.x+1,b2.y+1, ang2atom[i].x+1,ang2atom[i].y+1,ang2atom[i].z+1,  ang2bond[i].x&0xFFFF, ang2bond[i].y&0xFFFF, ang2bond[i].x, ang2bond[i].y );
    }
}

void eval_bonds( const bool substract_LJq ){
    for(int ib=0; ib<nbonds; ib++){
        Vec2i iat = bond2atom[ib];
        Vec3d dp; dp.set_sub( apos[iat.y], apos[iat.x] );
        Vec3d f;
        //printf( "%i %i (%g,%g,%g) (%g,%g,%g) \n", iat.x, iat.y, apos[iat.x].x, apos[iat.x].y, apos[iat.x].z, apos[iat.y].x, apos[iat.y].y, apos[iat.y].z );
        double l = dp.norm();
        f.set_mul(dp,1.0d/l);
        //printf( " %i (%i,%i) (%g,%g,%g) %g \n", ib, iat.x, iat.y, dp.x, dp.y, dp.z, l );

        lbond [ib] = l;
        hbond [ib] = f;

        f.mul( (l-bond_0[ib])*bond_k[ib] );

        if( substract_LJq ){
            //addAtomicForceLJQ( dp, f, aREQ[iat.x].x+aREQ[iat.y].x, -aREQ[iat.x].y*aREQ[iat.y].y, aREQ[iat.x].z*aREQ[iat.y].z );
            addAtomicForceMorseQ( dp, f, aREQ[iat.x].x+aREQ[iat.y].x, -aREQ[iat.x].y*aREQ[iat.y].y, aREQ[iat.x].z*aREQ[iat.y].z, gridFF.alpha );
        }

        aforce[iat.x].add( f );
        aforce[iat.y].sub( f );
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
        Vec3i ia = ang2atom[ig];
        //Vec3d h1 = hbond[ib.x];
        //Vec3d h2 = hbond[ib.y];
        Vec3d h1,h2;

        if(ib.x&SIGN_MASK){ ib.x&=0xFFFF; h1 = hbond[ib.x]; h1.mul(-1.0d); }else{ h1 = hbond[ib.x]; };
        if(ib.y&SIGN_MASK){ ib.y&=0xFFFF; h2 = hbond[ib.y]; h2.mul(-1.0d); }else{ h2 = hbond[ib.y]; };

        /*
        Vec3d h1_,h2_;
        h1_ = apos[ia.x] - apos[ia.z]; h1_.normalize();
        h2_ = apos[ia.y] - apos[ia.z]; h2_.normalize();
        printf( " ig %i c1 %g c2 %g   (%i,%i) \n", ig, h1_.dot(h1), h2_.dot(h2),   ang2bond[ig].x,  ang2bond[ig].y );
        */

        double c = h1.dot(h2);
        //double s = sqrt(1-c*c);

        //printf( " %i (%i,%i) (%g,%g,%g) (%g,%g,%g) \n", ib, ib.x,ib.y, h1.x,h1.y,h1.z, h2.x,h2.y,h2.z );

        Vec3d hf1,hf2; // unitary vectors of force
        hf1 = h2 - h1*c;
        hf2 = h1 - h2*c;
        //hf1 = h1*c - h2;
        //hf2 = h2*c - h1;

        //printf("ig %i ib (%i,%i) \n", ig, ib.x, ib.y );

        double fang = -ang_k[ig]/(1.02-c);
        hf1.mul( fang/lbond[ib.x] );
        hf2.mul( fang/lbond[ib.y] );

        //printf("ia (%i,%i,%i)\n", ia.x, ia.y, ia.z );

        /*
        glColor3f(0.0f,1.0f,0.0f);
        Draw3D::drawVecInPos( h1*0.25, apos[ia.z] );
        Draw3D::drawVecInPos( h2*0.25, apos[ia.z] );

        glColor3f(1.0f,0.0f,0.0f);
        Draw3D::drawVecInPos( hf1, apos[ia.x] );
        Draw3D::drawVecInPos( hf2, apos[ia.y] );
        */

        aforce[ia.x].add( hf1     );
        aforce[ia.y].add( hf2     );
        aforce[ia.z].sub( hf1+hf2 );

        //printf("ig %i\n", ig );
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
        double fang = (ang_k[ig]*cs.x*cs.y)*inv_sa;

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

void eval_LJq_On2(){
    for(int i=0; i<natoms; i++){
        Vec3d ljq_i = aREQ[i];
        Vec3d pi    = apos[i];
        Vec3d f; f.set(0.0);
        for(int j=0; j<natoms; j++){
            if(i!=j){
                Vec3d& ljq_j = aREQ[j];
                double rij = ljq_i.x+ljq_j.x;
                double eij = ljq_i.y*ljq_j.y;
                double qq  = ljq_i.z*ljq_j.z;
                addAtomicForceLJQ( pi-apos[j], f, rij, -eij, qq );
            }
        }
        aforce[i].add(f);
    }
}

void eval_MorseQ_On2(){
    for(int i=0; i<natoms; i++){
        Vec3d REQi = aREQ[i];
        Vec3d pi   = apos[i];
        Vec3d f; f.set(0.0);
        for(int j=0; j<natoms; j++){
            if(i!=j){
                Vec3d& REQj = aREQ[j];
                //printf("(%i,%i)", i, j );
                addAtomicForceMorseQ( pi-apos[j], f, REQi.x+REQj.x, -REQi.y*REQj.y, REQi.z*REQj.z, gridFF.alpha );
            }
        }
        aforce[i].add(f);
    }
    //exit(0);
}

void eval_MorseQ_Frags(){
    // we can use in principle eval_MorseQ_On2, but it will make it less numerically precise due to adding repulasive force between non-bonded atoms
    // also this will allow for more optimizations in future (e.g. per-fragment bounding boxes)
    for(int ifrag=0; ifrag<nFrag; ifrag++){
        int ia = frag2a[ifrag];
        int na = ia+fragNa[ifrag];
        for( int i=ia; i<na; i++ ){
            Vec3d REQi = aREQ[i];
            Vec3d pi   = apos[i];
            Vec3d f; f.set(0.0);
            for(int jfrag=0; jfrag<nFrag; jfrag++){
                if( jfrag == ifrag ) continue;
                int ja = frag2a[jfrag];
                int ma = ja+fragNa[jfrag];
                for( int j=ja; j<ma; j++ ){
                    Vec3d& REQj = aREQ[j];
                    addAtomicForceMorseQ( pi-apos[j], f, REQi.x+REQj.x, -REQi.y*REQj.y, REQi.z*REQj.z, gridFF.alpha );
                }
            }
            aforce[i].add(f);
        }
    }
}

void eval_FFgrid(){
    for(int i=0; i<natoms; i++){
        gridFF.addForce( apos[i], aPLQ[i], aforce[i] );
    }
}

void printBondParams(){
    for( int i=0; i<nbonds; i++ ){
        printf( "%i (%i,%i) %g %g \n", i, bond2atom[i].x+1, bond2atom[i].y+1, bond_0[i], bond_k[i] );
    }
}

}; // MMFF

#endif
