
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
    //Vec3i  * tors2bond = 0;
    //Quat4i * tors2atom = 0;
    //Vec2d  * tors_0    = 0; // [1]
    //double * tors_k    = 0; // [eV/A^2]

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
}

void cleanAtomForce(){
    for(int i=0; i<natoms; i++){ aforce[i].set(0.0); }
}

void angles_bond2atom(){
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

double eval_bonds(){
    double E=0;
    for(int ib=0; ib<nbonds; ib++){
        //printf( "bond %i\n", ib );
        Vec2i iat = bond2atom[ib];
        Vec3d f; f.set_sub( apos[iat.y], apos[iat.x] );
        //printf( "%i %i (%g,%g,%g) (%g,%g,%g) \n", iat.x, iat.y, apos[iat.x].x, apos[iat.x].y, apos[iat.x].z, apos[iat.y].x, apos[iat.y].y, apos[iat.y].z );
        double l = f.normalize();
        //printf( " %i (%i,%i) (%g,%g,%g) %g \n", ib, iat.x, iat.y, dp.x, dp.y, dp.z, l );
        lbond [ib] = l;
        hbond [ib] = f;
        f.mul( (l-bond_l0[ib])*bond_k[ib] );
        aforce[iat.x].add( f );
        aforce[iat.y].sub( f );
    }
    //exit(0);
    return E;
}


double eval_angles(){
    // efficient angular forcefield of cos(a/2) which prevent vanishing derivative of cos(a)
    // E = K*cos(a/2) = K*|ha+hb|/2
    double E=0;
    for(int ig=0; ig<nang; ig++){
        Vec2i ib = ang2bond[ig];
        Vec3i ia = ang2atom[ig];

        // -- read bond direction ( notice orientation sign A->B vs. A<-B )
        Vec3d ha,hb;
        //if(ib.a&SIGN_MASK){ ib.a&=0xFFFF; ha = hbond[ib.a]; ha.mul(-1.0d); }else{ ha = hbond[ib.a]; };
        //if(ib.b&SIGN_MASK){ ib.b&=0xFFFF; hb = hbond[ib.b]; ha.mul(-1.0d); }else{ hb = hbond[ib.b]; };

        // DEBUG
        ha.set_sub(apos[ia.a],apos[ia.c]);  ha.normalize();
        hb.set_sub(apos[ia.b],apos[ia.c]);  hb.normalize();

        // angular force
        Vec3d h; h.set_add( ha, hb );
        double c2 = h.norm2()*0.25d;               // cos(a/2) = |ha+hb|
        double s2 = 1-c2;
        //printf( "ang[%i] (%g,%g,%g) (%g,%g,%g) (%g,%g,%g) c2 %g s2 %g \n", ig, ha.x,ha.y,ha.z,  hb.x,hb.y,hb.z,  h.x,h.y,h.z,   c2, s2 );
        double c  = sqrt(c2);
        double s  = sqrt(s2);
        const Vec2d& cs0 = ang_cs0[ig];
        const double k   = ang_k  [ig];
        //E        +=  k*( cs0.x*c -cs0.y*s - 1);  // just for debug ?
        double fr = k*( cs0.x*s +cs0.y*c    );

        //printf( "ang[%i] cs0(%g,%g) cs(%g,%g) %g  \n", cs0.x,cs0.y, c,s, fr/k );

        // project to per leaver
        c2 *=-2;
        double lw     = 2*s*c;       //    |h - 2*c2*a| =  1/(2*s*c) = 1/sin(a)
        double fra    = fr/(lbond[ib.a]*lw);
        double frb    = fr/(lbond[ib.b]*lw);
        Vec3d fa,fb;
        fa.set_lincomb( fra,h,  fra*c2,ha );  //fa = (h - 2*c2*a)*fr / ( la* |h - 2*c2*a| );
        fb.set_lincomb( frb,h,  frb*c2,hb );  //fb = (h - 2*c2*b)*fr / ( lb* |h - 2*c2*b| );

        //glColor3f(1.0,0.0,0.0);
        //Draw3D::drawVecInPos( fa, apos[ia.a] );
        //Draw3D::drawVecInPos( fb, apos[ia.b] );

        // to atoms
        aforce[ia.a].add(fa); aforce[ia.c].sub(fa);
        aforce[ia.b].add(fb); aforce[ia.c].sub(fb);

    }
    return E;
}

void eval(){
    cleanAtomForce();
    eval_bonds();
    eval_angles();
};

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
