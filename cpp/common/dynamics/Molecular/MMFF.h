
#ifndef MMFF_h
#define MMFF_h

#include <vector>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Forces.h"
#include "GridFF.h"

#include "geom3D.h"

#include "integerOps.h"


#define SIGN_MASK 2147483648



template<typename T>
void frag2atoms(const Vec3TYPE<T>& pos, const Quat4TYPE<T>& qrot, int n, Vec3TYPE<T>* apos0, Vec3TYPE<T>* apos ){
    Mat3TYPE<T> mrot; qrot.toMatrix(mrot);
    for( int j=0; j<n; j++ ){
        Vec3TYPE<T> Mp;
        mrot.dot_to_T( apos0[j], Mp );
        apos[j].set_add( pos, Mp );
        printf( "frag2atoms[%i]  (%g,%g,%g) (%g,%g,%g) \n", j,  apos0[j].x, apos0[j].y, apos0[j].z,   apos[j].x, apos[j].y, apos[j].z  );
        //printf( "%i %i  (%g,%g,%g) (%g,%g,%g) \n", ifrag, j,  m_apos[j].x, m_apos[j].y, m_apos[j].z,   Tp.x, Tp.y, Tp.z  );
    }
}






void transformAtomRange( int n, Vec3d* inpos, Vec3d* outpos, const Vec3d& p0, const Vec3d& p1, const Mat3d& rot ){
    Vec3d p;
    //printf( "p0: (%f,%f,%f) \n", p0.x, p0.y, p0.z );
    for( int i=0; i<n; i++ ){
        rot.dot_to( inpos[i]-p0, p );
        outpos[i] = p + p1;
    }
}

class AtomicManipulator{ public:
    int     npick   =0;
    int   * picked  =0;
    int   * counts  =0;
    Vec3d * goalPos =0;

    int    nMaxCount  = 50;
    double Rconv      = 0.5;
    double Fmax       = 1.0;
    bool   newPerAtom = false;

    bool  bRelative = false;
    Vec3d goalSpan  = (Vec3d){1.0,1.0,1.0};
    Vec3d anizo     = (Vec3d){1.0,1.0,1.0};

    int   nenabled=0;
    int *  enabled=0;


    int natoms;
    Vec3d* aforce;
    Vec3d* apos;

    void bindAtoms( int natoms_, Vec3d* apos_, Vec3d* aforce_ ){ natoms=natoms_; apos=apos_; aforce=aforce_; }

    void realloc( int n){
        npick = n;
        _realloc( picked,  npick );
        _realloc( goalPos, npick );
        _realloc( counts,  npick );
    };

    void genGoal(int i){
        // should we ensure that goals does not repeat ?
        int iatom;
        if( enabled ){ iatom = enabled[rand()%nenabled]; }else{ iatom  = rand()%natoms; };
        picked [i] = iatom;
        Vec3d p = (Vec3d){ randf(-goalSpan.x,goalSpan.x), randf(-goalSpan.y,goalSpan.y), randf(-goalSpan.z,goalSpan.z) };
        if(bRelative) p.add( goalPos[i] );
        goalPos[i] = p;
        counts [i] = 0;
    }

    void genGoals(){ for(int i=0; i<npick; i++){ genGoal(i); }; }

    inline bool forceToTarget( int i, Vec3d pos, double Rconv, double Fmax ){
        Vec3d d = pos - apos[i];
        d.mul(anizo);
        double r2 = d.norm2();
        if( r2 < Rconv*Rconv ){
            double k=Fmax/Rconv;
            aforce[i].add_mul( d, k );
            return true;
        }else{
            aforce[i].add_mul( d, Fmax/sqrt(r2) );
            return false;
        }
    }

    bool forceToGoals(){
        bool ret=true;
        for(int i=0; i<npick; i++ ){
            bool b = forceToTarget( picked[i], goalPos[i], Rconv, Fmax );
            ret &= b;
            counts[i]++;
            if( (b&&newPerAtom) || (counts[i]>nMaxCount) ){
                printf( "new manipulation goal %i \n", i );
                genGoal(i);
            }
        }
        return ret;
    }

};



/*
Chimera generuje topologie ... zdrojak C/python
*/

//inline int invIndex( int i ){ return i^SIGN_MASK; }

// compied from RigidMolecule.h      ============>   TODO : make common lib of inter-molecular forcefields formulas

int pickParticle( int n, Vec3d * ps, const Vec3d& ray0, const Vec3d& hRay, double R ){
    double tmin =  1e+300;
    int imin    = -1;
    for(int i=0; i<n; i++){
        double ti = raySphere( ray0, hRay, R, ps[i] );
        //printf( "atom %i t %g %g hRay(%g,%g,%g) ps(%g,%g,%g) \n", i, ti, tmin, hRay.x,hRay.y,hRay.z,  ps[i].x, ps[i].y, ps[i].z );
        //printf( "atom %i t %g %g %g ray0(%g,%g,%g) ps(%g,%g,%g) \n", i, ti, R, tmin, ray0.x,ray0.y,ray0.z,  ps[i].x, ps[i].y, ps[i].z );
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
    int    * atom2frag = NULL; // [natom], inverted frag2a
    Vec3d  * apos   = NULL;   // atomic position
    //Vec3d  * aLJq   = NULL;
    Vec3d  * aREQ   = NULL;
    Vec3d  * aPLQ   = NULL;   // this is used for grid-accelerated factorized potential
    double * lbond  = NULL;   // bond lengths
    Vec3d  * hbond  = NULL;   // normalized bond unitary vectors
    Vec3d  * aforce = NULL;

    std::vector<MolType> molTypes;

    int nFrag=0;  // rigid fragments
    //int *    imolTypes = NULL;
    Vec3d  ** fapos0s = NULL; // position of atoms in rigid fragment
    int    *  frag2a  = NULL; // start of the fragment in forcefield
    int    *  fragNa  = NULL; // lengh of the fragment
    double *  poses   = NULL; // rigd body pose of molecule (pos,qRot);
    double *  poseFs  = NULL; //
    //double *  poseVs  = NULL; //

    // WHAT IS THIS?  - perhaps dynamic moleculers for ConfSearch ?
    int nDyn=0;
    double * dynInvMass = NULL;
    double * dynPos     = NULL;
    double * dynVel     = NULL;
    double * dynForce   = NULL;

    //RigidSubstrate substrate;
    GridFF gridFF;

    Box    Collision_box   = (Box){(Vec3d){0.0,0.0,0.0},(Vec3d){10.0,10.0,10.0}};
    double Collision_Rsc   = 0.5;
    double Collision_F2max = 1.0;

void setCollisionRF( double Rsc ){
    Collision_Rsc   = Rsc;
    Collision_F2max = exp( 2*2*(3.5*(1-Rsc))*gridFF.alpha );
}


void allocate( int natoms_, int nbonds_, int nang_, int ntors_ ){
    natoms=natoms_; nbonds=nbonds_; nang=nang_; ntors=ntors_;
    printf( "MMFF::allocate natoms: %i  nbonds: %i  nang: %i ntors: %i \n", natoms, nbonds, nang, ntors );
    if(atypes   ==NULL) atypes    = new int   [natoms];
    if(apos     ==NULL) apos      = new Vec3d [natoms];
    if(aforce   ==NULL) aforce    = new Vec3d [natoms];
    if(aREQ     ==NULL) aREQ      = new Vec3d [natoms];
    if(atom2frag==NULL) atom2frag = new int   [natoms];

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
    //poseVs    = new double[nFrag*8];
    fapos0s   = new Vec3d*[nFrag];
}

void allocateDyn(){
    nDyn = getNDym();
    if(dynInvMass)delete [] dynInvMass; dynInvMass = new double[nDyn];
    if(dynPos)    delete [] dynPos;     dynPos     = new double[nDyn];
    if(dynVel)    delete [] dynVel;     dynVel     = new double[nDyn];
    if(dynForce)  delete [] dynForce;   dynForce   = new double[nDyn];
}

void deallocate(){
    delete [] apos;     delete [] aforce;   delete [] aREQ;
    delete [] lbond;    delete [] hbond;    delete [] bond2atom; delete [] bond_0; delete [] bond_k;
    delete [] ang2bond; delete [] ang2atom; delete [] ang_0;     delete [] ang_k;
    if(aPLQ) delete [] aPLQ;
}

int pickBond( const Vec3d& ray0, const Vec3d& hRay, double R ){
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

int getNDym(){
    int nfree=0;
    for(int i=0; i<natoms; i++){  if( 0>atom2frag[i] ){ nfree++;  } }
    return nFrag*8 + nfree*3;
}

// TODO :
// this can be probably optimized if use one array of this structure:
//       | frag | atomic_free | atomic_frag |
// then  | frag | atomic_free |               is used for dynamical variables
// and          | atomic_free | atomic_frag | is used for atom-wise forces
//  this however assumes that atomic and free forces are re-ordered

void initDyn(){
    // copy to dynamical variables
    int off = 8*nFrag;
    memcpy( dynPos,   poses, sizeof(double)*off );
    memcpy( dynForce, poseFs, sizeof(double)*off );
    for( int i=0; i<off; i++ ){ dynVel[i]=0; };
    Vec3d *pF = (Vec3d*)( dynForce + off );
    Vec3d *pP = (Vec3d*)( dynPos   + off );
    Vec3d *pV = (Vec3d*)( dynVel   + off );
    for(int ia=0; ia<natoms; ia++){
        if( 0>atom2frag[ia] ){
            (*pP) = apos[ia];
            (*pF) = aforce[ia];
            (*pV).set(0.0);
            pP++; pF++; pV++;
        }
    };
};

void toDym( const bool bPos ){
    // copy to dynamical variables
    int off = 8*nFrag;
    if(bPos) memcpy( dynPos, poses, sizeof(double)*off ); // needed if positions modified outside optimizer
    memcpy( dynForce, poseFs, sizeof(double)*off );
    Vec3d *pF = (Vec3d*)( dynForce + off );
    Vec3d *pP = (Vec3d*)( dynPos + off   );
    for(int ia=0; ia<natoms; ia++){
        if( 0>atom2frag[ia] ){
            *pF = aforce[ia]; pF++;
            if(bPos){ *pP=apos[ia]; pP++; };              // needed if positions modified outside optimizer
        }
    };
}

void fromDym(){
    // copy from dynamical variables
    int off = 8*nFrag;
    memcpy( poses, dynPos, sizeof(double)*off );
    Vec3d *pP = (Vec3d*)( dynPos + off );
    for(int ia=0; ia<natoms; ia++){
        if( 0>atom2frag[ia] ){
            apos[ia] = *pP;
            pP++;
        }
    };
}

void frags2atoms(){
    // : fragment pose -> atomic position
    //printf( ">>> MMFF::frags2atoms %i \n", nFrag );
    for(int ifrag=0; ifrag<nFrag; ifrag++){
        //frag2atoms( ifrag, apos+frag2a[ifrag] );
        int im8 = ifrag<<3;
        Vec3d   pos = *((Vec3d* )(poses+im8  ));
        Quat4d  rot = *((Quat4d*)(poses+im8+4));
        frag2atoms( pos, rot, fragNa[ifrag], fapos0s[ifrag], apos+frag2a[ifrag] );
    }
    //exit(0);
}

void cleanAtomForce(){
    for(int i=0; i<natoms; i++){ aforce[i].set(0.0); }
}

void cleanPoseTemps(){
    for( int ifrag=0; ifrag<nFrag; ifrag++ ){
        int im8 = ifrag<<3;
        double * poseF_i = poseFs+im8;
        ((Quat4d*)(poseF_i  ))->set(0.0);
        ((Quat4d*)(poseF_i+4))->set(0.0);
    }
}

void checkPoseUnitary(){
    for( int ifrag=0; ifrag<nFrag; ifrag++ ){
        int ioff = (ifrag<<3) + 4;
        //Quat4d& rot  = *(Quat4d*)(poses+ioff);
        //Quat4d& frot = *(Quat4d*)(poseFs+ioff);
        //Quat4d& vrot = *(Quat4d*)(poseVs+ioff);
        Quat4d& rot  = *(Quat4d*)(dynPos  +ioff);
        Quat4d& frot = *(Quat4d*)(dynForce+ioff);
        Quat4d& vrot = *(Quat4d*)(dynVel  +ioff);
        rot.checkNormalized(1e-4);
        double cdot;
        cdot = rot.dot( frot );  frot.add_mul( rot, -cdot );
        cdot = rot.dot( vrot );  vrot.add_mul( rot, -cdot );
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

bool getCollosion( int i0, int n, double scR, const bool bFColPauli, Vec3d* poss, Vec3d* forces, double* R2mins ){
    int i1 = i0+n;
    bool ret = false;
    const bool justAny = ( (forces==0)&&(R2mins==0) );
    if(forces) for(int i=0; i<n; i++){ forces[i].set(0.0); };
    //if(energys) for(int i=0; i<n; i++){ dforces[i].set(0); };
    if(R2mins)   for(int i=0; i<n; i++){ R2mins[i]=1e+300; };
    for( int i=0; i<natoms; i++ ){
        Vec3d  pi = apos[i];
        //double ri  = aREQ[i].x;
        Vec3d REQi =  aREQ[i];
        if( (i>=i0)&&(i<i1) ) continue;
        for(int j=0; j<n; j++){
            Vec3d  d  = poss[j] - pi;
            double r2 = d.norm2();
            Vec3d&  REQj = aREQ[i0+j];
            double R  = ( REQj.x + REQi.x )*scR;
            double R2 = R*R;
            if( r2 < R2 ){
                if(justAny){
                    return true;
                }else{
                    ret=true;
                    if(R2mins){ double r2m=R2mins[j]; if(r2<r2m)R2mins[j]=r2; };
                    if(forces){
                        if (bFColPauli){ addAtomicForceExp( d, forces[j], REQi.x+REQj.x, REQi.y*REQj.y, gridFF.alpha*2 ); }
                        else           { forces[j].add_mul( d, 1/r2 - 1/R2 );                                             }
                    };
                    //if(dforces){ Energy [j].add_mul( d, 1/r2 ) };
                }
            }
        }
    }
    return ret;
}

double getCollisionGrid( int i0, int n, Vec3d* poss, Vec3d* forces ){
    double F2max=0;
    Vec3d f,gpos;
    for(int i=0; i<n; i++){
        //Vec3d , gridFF.addForce( apos[j], aPLQ[j], aforce[j] );
        gridFF.grid.cartesian2grid(poss[i],gpos);
        f = interpolate3DvecWrap( gridFF.FFPauli, gridFF.grid.n, gpos );
        f.mul( aPLQ[i0+i].x );
        double fr2 = f.norm2();
        if( fr2>F2max ) F2max=fr2;
        if(forces){ forces[i].add( f ); }
    }
    return F2max;
}

//void translate( int i0, int n, Vec3d dpos           ){ ( int i=i0, i<(i0+n); i++ ) apos[i].add(dpos); }

bool tryPose( int i0, int n, Vec3d p0, Vec3d p1, Mat3d rot ){

    Vec3d poss[n]; // Stack Allocation - should we change it ?
    transformAtomRange( n, apos+i0, poss, p0, p1, rot );
    double f2 = getCollisionGrid( i0, n, poss, 0 );
    bool bOK = ( f2<Collision_F2max );
    if( bOK ) bOK = !getCollosion( i0, n, Collision_Rsc, false, poss, 0, 0 );
    //bOK=true;
    if( bOK ){ for(int i=0; i<n; i++) apos[i0+i] = poss[i]; }
}

bool tryFragPose( int ifrag, bool bRel, Vec3d pos, Quat4d qrot ){

    int i0 = frag2a[ifrag];
    int n  = fragNa[ifrag];
    Vec3d poss[n]; // Stack Allocation - should we change it ?

    int im8    = ifrag<<3;
    Vec3d*  fragPos = ((Vec3d* )(poses+im8  ));
    Quat4d* fragRot = ((Quat4d*)(poses+im8+4));
    //printf( "fragPos (%f,%f,%f) \n", fragPos->x, fragPos->y, fragPos->z );
    if( bRel ){
        pos.add  ( *fragPos );
        qrot.qmul( *fragRot ); // Not sure about direction
    }
    Mat3d T; qrot.toMatrix(T);

    Vec3d* m_apos = fapos0s[ifrag];
    for( int j=0; j<n; j++ ){
        Vec3d Tp;
        T.dot_to_T( m_apos[j], Tp );
        poss[j].set_add( pos, Tp );
    }

    //double f2 = getCollisionGrid( i0, n, poss, 0 );
    //bool bOK  = ( f2<Try_F2max );

    bool bOK      = Collision_box.pointIn( pos );
    if( bOK ) bOK = Collision_F2max>getCollisionGrid( i0, n, poss, 0 );
    if( bOK ) bOK = !getCollosion( i0, n, Collision_Rsc, false, poss, 0, 0 );
    //bOK=true;
    //if( bOK ){ for(int i=0; i<n; i++) apos[i0+i] = poss[i]; }
    if( bOK ){ *fragPos=pos;  *fragRot=qrot; }
    //printf( "fragPos (%f,%f,%f) \n", fragPos->x, fragPos->y, fragPos->z );
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

void eval_MorseQ_On2_fragAware(){
    for(int i=0; i<natoms; i++){
        Vec3d REQi = aREQ[i];
        Vec3d pi   = apos[i];
        Vec3d f; f.set(0.0);
        int ifrag = atom2frag[i];
        for(int j=0; j<natoms; j++){
            if( (ifrag>=0) && (ifrag==atom2frag[j]) ){
                //printf( " %i %i | %i %i \n", i, j, ifrag, atom2frag[j] );
                continue;
            }
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
            //printf( " (%i,%i) (%g,%g,%g)  \n",  ifrag, i, REQi.x, REQi.y, REQi.z );
            for(int jfrag=0; jfrag<nFrag; jfrag++){
                if( jfrag == ifrag ) continue;
                int ja = frag2a[jfrag];
                int ma = ja+fragNa[jfrag];
                for( int j=ja; j<ma; j++ ){
                    Vec3d& REQj = aREQ[j];
                    //printf( " (%i,%i) (%i,%i) ()  \n",  ifrag jfrag i j );
                    //addAtomicForceMorseQ( pi-apos[j], f, REQi.x+REQj.x, -REQi.y*REQj.y, 0, gridFF.alpha );
                    addAtomicForceMorseQ( pi-apos[j], f, REQi.x+REQj.x, -REQi.y*REQj.y, REQi.z*REQj.z, gridFF.alpha );
                }
            }
            aforce[i].add(f);
        }
    }
    //exit(0);
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

void printAtomInfo(){
    printf("MMFF::printAtomInfo : \n" );
    for(int i=0; i<natoms; i++){
        printf( "%i %i %i %f %f %f \n", i, atom2frag[i], atypes[i], aREQ[i].x, aREQ[i].y, aREQ[i].z );
    }
}

}; // MMFF

#endif
