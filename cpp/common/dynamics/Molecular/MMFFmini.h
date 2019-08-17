
#ifndef MMFFmini_h
#define MMFFmini_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"

//#include "integerOps.h"

#define SIGN_MASK 2147483648

// ======================
// ====   MMFF
// ======================

/*
inline void pairs2triple( const  Vec2i b1, const Vec2i b2, Vec3i& tri, bool& flip1, bool& flip1 ){
    if     ( b1.x == b2.x ){ tri.set( b1.y, b2.y, b1.x ); flip1=false; flip1=false; }
    else if( b1.x == b2.y ){ tri.set( b1.y, b2.x, b1.x ); flip1=false; flip1=true;  }
    else if( b1.y == b2.x ){ tri.set( b1.x, b2.y, b1.y ); flip1=true;  flip1=false; }
    else if( b1.y == b2.y ){ tri.set( b1.x, b2.x, b1.y ); flip1=true;  flip1=true;  }
}
*/

inline void pairs2triple( const  Vec2i b1, const Vec2i b2, Vec3i& tri, bool& flip1, bool& flip2 ){
    if     ( b1.x == b2.x ){ tri.set( b1.y, b1.x, b2.y ); flip1=false; flip1=false; }
    else if( b1.x == b2.y ){ tri.set( b1.y, b1.x, b2.x ); flip1=false; flip1=true;  }
    else if( b1.y == b2.x ){ tri.set( b1.x, b1.y, b2.y ); flip1=true;  flip1=false; }
    else if( b1.y == b2.y ){ tri.set( b1.x, b1.y, b2.x ); flip1=true;  flip1=true;  }
}

class MMFFmini{ public:
    int  natoms=0, nbonds=0, nang=0, ntors=0;

    // --- Parameters

    Vec2i  * bond2atom = 0;
    double * bond_l0   = 0;  // [A]
    double * bond_k    = 0;  // [eV/A] ?

    Vec2i  * ang2bond  = 0;
    Vec3i  * ang2atom  = 0;
    //Vec2d  * ang_0   = 0; // alpha_0
    Vec2d  * ang_cs0   = 0; // cos(a),sin(a)
    double * ang_k     = 0; // [eV/A^2]

    // ToDo:
    //  We don't care about dihedrals yet
    //  Later we perhaps implement it using additional DOF (pi-orbitals,electron pairs)
    Vec3i  * tors2bond = 0;
    Quat4i * tors2atom = 0;
    //Vec2d  * tors_0    = 0; // [1]
    int    * tors_n    = 0;
    double * tors_k    = 0; // [eV/A^2]

    //Vec3d * hpi = 0;  // NOTE : nstead of pi define dummy atom !!!!


    // --- State variables

    Vec3d  * apos     = 0;   // atomic position
    Vec3d  * aforce   = 0;   // atomic position

    // --- Axuliary variables

    double * lbond  = 0;   // bond lengths
    Vec3d  * hbond  = 0;   // normalized bond unitary vectors





void realloc( int natoms_, int nbonds_, int nang_, int ntors_ ){
    natoms=natoms_; nbonds=nbonds_; nang=nang_; ntors=ntors_;
    //printf( "MMFF::allocate natoms: %i  nbonds: %i  nang: %i ntors: %i \n", natoms, nbonds, nang, ntors );
    _realloc( apos      , natoms );
    _realloc( aforce    , natoms );

    _realloc( lbond     , nbonds );
    _realloc( hbond     , nbonds );

    _realloc( bond2atom , nbonds );
    _realloc( bond_l0   , nbonds );
    _realloc( bond_k    , nbonds );

    _realloc( ang2bond  , nang   );
    _realloc( ang2atom  , nang   );
    _realloc( ang_cs0   , nang   );
    _realloc( ang_k     , nang   );

    _realloc( tors2bond , ntors  );
    _realloc( tors2atom , ntors  );
    _realloc( tors_n    , ntors  );
    _realloc( tors_k    , ntors  );

}

void dealloc(){
    _dealloc( apos      );
    _dealloc( aforce    );

    _dealloc( lbond     );
    _dealloc( hbond     );

    _dealloc( bond2atom );
    _dealloc( bond_l0   );
    _dealloc( bond_k    );

    _dealloc( ang2bond  );
    _dealloc( ang2atom  );
    _dealloc( ang_cs0   );
    _dealloc( ang_k     );

    _dealloc( tors2bond );
    _dealloc( tors2atom );
    _dealloc( tors_n    );
    _dealloc( tors_k    );
}

inline void readSignedBond(int& i, Vec3d& h){ if(i&SIGN_MASK){ i&=0xFFFF; h = hbond[i]; h.mul(-1.0d); }else{ h = hbond[i]; }; };

void cleanAtomForce(){ for(int i=0; i<natoms; i++){ aforce[i].set(0.0); } }


// ============== Evaluation


double eval_bond(int ib){
    //printf( "bond %i\n", ib );
    Vec2i iat = bond2atom[ib];
    Vec3d f; f.set_sub( apos[iat.y], apos[iat.x] );
    //printf( "%i %i (%g,%g,%g) (%g,%g,%g) \n", iat.x, iat.y, apos[iat.x].x, apos[iat.x].y, apos[iat.x].z, apos[iat.y].x, apos[iat.y].y, apos[iat.y].z );
    double l = f.normalize();
    //printf( " %i (%i,%i) (%g,%g,%g) %g \n", ib, iat.x, iat.y, dp.x, dp.y, dp.z, l );
    lbond [ib] = l;
    hbond [ib] = f;
    const double k = bond_k[ib];
    double dl = (l-bond_l0[ib]);
    f.mul( dl*k );
    aforce[iat.x].add( f );
    aforce[iat.y].sub( f );
    return k*dl*dl;
}

double eval_angle(int ig){
    //if(ig!=0) return 0;
    // efficient angular forcefield of cos(a/2) which prevent vanishing derivative of cos(a)
    // E = K*cos(a/2) = K*|ha+hb|/2
    Vec2i ib = ang2bond[ig];
    Vec3i ia = ang2atom[ig];

    // -- read bond direction ( notice orientation sign A->B vs. A<-B )
    Vec3d ha,hb;
    //if(ib.a&SIGN_MASK){ ib.a&=0xFFFF; ha = hbond[ib.a]; ha.mul(-1.0d); }else{ ha = hbond[ib.a]; };
    //if(ib.b&SIGN_MASK){ ib.b&=0xFFFF; hb = hbond[ib.b]; hb.mul(-1.0d); }else{ hb = hbond[ib.b]; };
    readSignedBond(ib.a, ha);
    readSignedBond(ib.b, hb);
    //printf( "ang[%i] ib(%i,%i) (%g,%g,%g) (%g,%g,%g)\n", ig, ib.a, ib.b, ha.x,ha.y,ha.z,   hb.x,hb.y,hb.z );

    // instead of
    //ha = hbond[ib.a];
    //hb = hbond[ib.b];
    //Vec3d ia; bool iflip,jflip;
    //pairs2triple( bond2atom[ib.i], bond2atom[ib.j], ia, irev, jrev );
    //if(iflip){ha.mul(-1);};
    //if(jflip){hb.mul(-1);};


    // DEBUG
    //ha.set_sub(apos[ia.a],apos[ia.c]);  ha.normalize();
    //hb.set_sub(apos[ia.b],apos[ia.c]);  hb.normalize();

    // angular force
    Vec3d h; h.set_add( ha, hb );
    double c2 = h.norm2()*0.25d;               // cos(a/2) = |ha+hb|
    double s2 = 1-c2;
    //printf( "ang[%i] (%g,%g,%g) (%g,%g,%g) (%g,%g,%g) c2 %g s2 %g \n", ig, ha.x,ha.y,ha.z,  hb.x,hb.y,hb.z,  h.x,h.y,h.z,   c2, s2 );
    double c  = sqrt(c2);
    double s  = sqrt(s2);
    const Vec2d& cs0 = ang_cs0[ig];
    const double k   = ang_k  [ig];
    double E         = k*( cs0.x*c -cs0.y*s - 1);  // just for debug ?
    double fr        = k*( cs0.x*s +cs0.y*c    );

    //printf( "ang[%i] cs0(%g,%g) cs(%g,%g) %g  \n", cs0.x,cs0.y, c,s, fr/k );

    // project to per leaver
    c2 *=-2;
    double lw     = 2*s*c;       //    |h - 2*c2*a| =  1/(2*s*c) = 1/sin(a)
    double fra    = fr/(lbond[ib.a]*lw);
    double frb    = fr/(lbond[ib.b]*lw);
    Vec3d fa,fb;
    fa.set_lincomb( fra,h,  fra*c2,ha );  //fa = (h - 2*c2*a)*fr / ( la* |h - 2*c2*a| );
    fb.set_lincomb( frb,h,  frb*c2,hb );  //fb = (h - 2*c2*b)*fr / ( lb* |h - 2*c2*b| );

    //printf( "ang[%i] %g (%g,%g,%g) (%g,%g,%g)\n", ig,  c, fa.x,fa.y,fa.z,   fb.x,fb.y,fb.z );

    //glColor3f(1.0,0.0,0.0); double fsc=100.0;
    //Draw3D::drawVecInPos( fa*fsc, apos[ia.x] );
    //Draw3D::drawVecInPos( fb*fsc, apos[ia.z] );
    //Draw3D::drawVecInPos( (fa+fb)*-fsc, apos[ia.y] );

    // to atoms
    aforce[ia.x].add(fa); aforce[ia.y].sub(fa);
    aforce[ia.z].add(fb); aforce[ia.y].sub(fb);

    return E;
}

double eval_torsion(int it){
    // efficient angular forcefield of cos(a/2) which prevent vanishing derivative of cos(a)
    // E = K*cos(a/2) = K*|ha+hb|/2
    Vec3i  ib = tors2bond[it];
    Quat4i ia = tors2atom[it];

    // -- read bond direction ( notice orientation sign A->B vs. A<-B )
    Vec3d ha,hb,hab;
    //if(ib.a&SIGN_MASK){ ib.a&=0xFFFF; ha = hbond[ib.a]; ha.mul(-1.0d); }else{ ha = hbond[ib.a]; };
    //if(ib.b&SIGN_MASK){ ib.b&=0xFFFF; hb = hbond[ib.b]; hb.mul(-1.0d); }else{ hb = hbond[ib.b]; };
    //if(ib.c&SIGN_MASK){ ib.c&=0xFFFF; hc = hbond[ib.b]; hc.mul(-1.0d); }else{ hb = hbond[ib.b]; };
    readSignedBond(ib.x, ha);
    readSignedBond(ib.y, hab);
    readSignedBond(ib.z, hb);

    Vec3d fa,fb;
    fa.set_cross( ha, hab );  //fa.mul( 1/fa.norm2() );
    fb.set_cross( hb, hab );  //fb.mul( 1/fb.norm2() );

    //glColor3f(1.0,0.0,0.0);
    //Draw3D::drawVecInPos(fa, apos[ia.x] );
    //Draw3D::drawVecInPos(fb, apos[ia.w] );

    double ira = 1/fa.norm();
    double irb = 1/fb.norm();

    fa.mul(ira);
    fb.mul(irb);

    Vec3d u; u.set_cross(fa,fb);
    //glColor3f(1.0,0.0,1.0); Draw3D::drawVecInPos(u, apos[ia.x]);
    Vec2d cs{
        fa.dot(fb ), // cos(a) = <a|b>
        u .dot(hab), // sin(a) = |a,b|
    };
    Vec2d csn = cs;

    const int n = tors_n[it];
    for(int i=0; i<n-1; i++){
        csn.mul_cmplx(cs);
    }

    const double k = tors_k[it];
    double E   = k*csn.x;
    double fr  = k*csn.y;

    //glColor3f(1.0,0.0,1.0); Draw3D::drawVecInPos(hab.normalized()*E, apos[ia.x]);
    //glColor3f(1.0,1.0,1.0); Draw3D::drawVecInPos(ha.normalized()*E, apos[ia.x]);

    fa.mul(fr*ira);
    fb.mul(fr*irb);

    //glColor3f(1.0,0.0,0.0);  Draw3D::drawVecInPos(fa, apos[ia.x]); Draw3D::drawVecInPos(fb, apos[ia.w]);

    //printf("tors[%i] (%i,%i,%i,%i)\n", it, ia.x, ia.y, ia.z, ia.w );

    // to atoms
    aforce[ia.x].add(fa);
    aforce[ia.y].sub(fa);
    aforce[ia.z].sub(fb);
    aforce[ia.w].add(fb);

    return E;
}

double eval_bonds(){
    double E=0;
    for(int ib=0; ib<nbonds; ib++){  E+=eval_bond(ib); }
    return E;
}

double eval_angles(){
    double E=0;
    for(int ig=0; ig<nang; ig++){ E+= eval_angle(ig); }
    return E;
}

double eval_torsions(){
    double E=0;
    for(int it=0; it<ntors; it++){ E+= eval_torsion(it); }
    return E;
}

void eval(){
    cleanAtomForce();
    eval_bonds();
    eval_angles();
    eval_torsions();
};

// ============== Preparation


void angles_bond2atom(){
    for(int i=0; i<nang; i++){
        Vec2i ib = ang2bond[i];
        Vec2i b1,b2;
        b1 = bond2atom[ib.i];
        b2 = bond2atom[ib.j];
        bool flip1,flip2;
        pairs2triple( b1, b2, ang2atom[i], flip1, flip2 );
        if(!flip1){ ang2bond[i].i|=SIGN_MASK; };
        if( flip2){ ang2bond[i].j|=SIGN_MASK; };
        printf( "ang[%i] ((%i,%i)(%i,%i))->(%i,%i,%i) (%i,%i)==(%i,%i) \n", i, b1.x,b1.y, b2.x,b2.y, ang2atom[i].x,ang2atom[i].y,ang2atom[i].z,ang2bond[i].x&0xFFFF, ang2bond[i].y&0xFFFF,ang2bond[i].x,        ang2bond[i].y         );
    }
}

void torsions_bond2atom(){
    for(int i=0; i<ntors; i++){
        Vec3i ib = tors2bond[i];
        Vec2i b1,b2,b12;
        b1  = bond2atom[ib.x];
        b12 = bond2atom[ib.y];
        b2  = bond2atom[ib.z];
        bool flip1,flip2,flip12;
        pairs2triple( b1, b12, *(Vec3i*)(        &tors2atom[i]    ), flip1,  flip12 );
        pairs2triple( b12, b2, *(Vec3i*)(((int*)(&tors2atom[i]))+1), flip12, flip2  );
        if( flip1 ){ tors2bond[i].i|=SIGN_MASK; };
        if( flip12){ tors2bond[i].j|=SIGN_MASK; };
        if( flip2 ){ tors2bond[i].j|=SIGN_MASK; };
        printf( "tors[%i] ((%i,%i)(%i,%i)(%i,%i))->(%i,%i,%i,%i) (%i,%i,%i)==(%i,%i,%i) \n", i, b1.x,b1.y, b12.x,b12.y, b2.x,b2.y,tors2atom[i].x,tors2atom[i].y,tors2atom[i].z,tors2atom[i].w,tors2bond[i].x&0xFFFF, tors2bond[i].y&0xFFFF, tors2bond[i].z&0xFFFF,tors2bond[i].x,        tors2bond[i].y,        tors2bond[i].z         );
    }
}

void printBondParams(){
    for( int i=0; i<nbonds; i++ ){
        printf( "bond[%i] (%i,%i) l0 %g k %g \n", i, bond2atom[i].x+1, bond2atom[i].y+1, bond_l0[i], bond_k[i] );
    }
}

void printAngleParams(){
    for( int i=0; i<nbonds; i++ ){
        printf( "angle[%i] (%i,%i|%i) cs0(%g,%g) k %g \n", i, ang2atom[i].a+1, ang2atom[i].b+1, ang2atom[i].c+1, ang_cs0[i].x, ang_cs0[i].y, ang_k[i] );
    }
}

}; // MMFF

#endif
