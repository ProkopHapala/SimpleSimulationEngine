
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
    if     ( b1.y == b2.x ){ tri.set( b1.x, b1.y, b2.y ); flip1=false; flip2=false;  }
    else if( b1.y == b2.y ){ tri.set( b1.x, b1.y, b2.x ); flip1=false; flip2=true;   }
    else if( b1.x == b2.x ){ tri.set( b1.y, b1.x, b2.y ); flip1=true;  flip2=false;  }
    else if( b1.x == b2.y ){ tri.set( b1.y, b1.x, b2.x ); flip1=true;  flip2=true;   }
}


/*
inline void pairs2triple( const  Vec2i b1, const Vec2i b2, Vec3i& tri, bool& flip1, bool& flip2 ){
    if     ( b1.x == b2.x ){ tri.set( b1.y, b1.x, b2.y ); flip1=false; flip1=false; }
    else if( b1.x == b2.y ){ tri.set( b1.y, b1.x, b2.x ); flip1=false; flip1=true;  }
    else if( b1.y == b2.x ){ tri.set( b1.x, b1.y, b2.y ); flip1=true;  flip1=false; }
    else if( b1.y == b2.y ){ tri.set( b1.x, b1.y, b2.x ); flip1=true;  flip1=true;  }
}
*/

void sum(int n, Vec3d* ps, Vec3d& psum){ for(int i=0;i<n;i++){ psum.add(ps[i]); } };

void sumTroq(int n, Vec3d* fs, Vec3d* ps, const Vec3d& cog, const Vec3d& fav, Vec3d& torq){
    for(int i=0;i<n;i++){  torq.add_cross(ps[i]-cog,fs[i]-fav);  }
};

void checkForceInvariatns( int n, Vec3d* fs, Vec3d* ps, Vec3d& cog, Vec3d& fsum, Vec3d& torq ){
    cog=Vec3dZero;
    fsum=Vec3dZero;
    torq=Vec3dZero;
    double dw = 1.d/n;
    sum(n, ps, cog ); cog.mul(dw);
    sum(n, fs, fsum); //cog.mul(dw);
    sumTroq(n, fs, ps, cog, fsum*dw, torq );
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

inline void setAngleParam(int i, double a0, double k){
    a0=a0*0.5; // we store half angle
    ang_cs0[i].fromAngle(a0);
    //ang_cs0[i] = (Vec2d){cos(a0),sin(a0)};
    ang_k  [i]=k;
};
inline void setBondParam(int i, double l0, double k){ bond_l0[i] = l0; bond_k[i]  = k; };
inline void setTorsParam(int i, int     n, double k){ tors_n [i] =  n; tors_k[i]  = k; };

inline void readSignedBond(int& i, Vec3d& h){ if(i&SIGN_MASK){ i&=0xFFFF; h = hbond[i]; h.mul(-1.0d); }else{ h = hbond[i]; }; };

void cleanAtomForce(){ for(int i=0; i<natoms; i++){ aforce[i].set(0.0); } }

// ============== Evaluation


int i_DEBUG = 0;

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
    f.mul( dl*k*2 );
    aforce[iat.x].add( f );
    aforce[iat.y].sub( f );
    return k*dl*dl;
}

double eval_angle(int ig){
    //if(ig!=i_DEBUG) return 0;
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
    //Vec3d h; h.set_sub( hb, ha );
    double c2 = h.norm2()*0.25d;               // cos(a/2) = |ha+hb|
    double s2 = 1-c2;
    //printf( "ang[%i] (%g,%g,%g) (%g,%g,%g) (%g,%g,%g) c2 %g s2 %g \n", ig, ha.x,ha.y,ha.z,  hb.x,hb.y,hb.z,  h.x,h.y,h.z,   c2, s2 );
    double c  = sqrt(c2);
    double s  = sqrt(s2);
    //const Vec2d cs0 = ang_cs0[ig];
    Vec2d cs = ang_cs0[ig];
    const double k   = ang_k  [ig];

    cs.udiv_cmplx({c,s});
    double E         =  k*( 1 - cs.x );  // just for debug ?
    double fr        = -k*(     cs.y );

    //printf(  "E %g c %g \n", E, cs.x  );

    //printf( "ang[%i] cs(%g,%g)%g  rcs(%g,%g)%g    %g %g %g \n",  c,s, atan2(s,c)*180/M_PI, cs.x,cs.y,atan2(cs.y,cs.x)*180/M_PI,  fr/k, E, k );

    //double E         = k*( cs0.x*c -cs0.y*s - 1);  // just for debug ?
    //double fr        = k*( cs0.x*s +cs0.y*c    );
    //printf( "ang[%i] cs0(%g,%g)%g cs(%g,%g)%g    %g %g %g \n", cs0.x,cs0.y,atan2(cs0.y,cs0.x)*180/M_PI, c,s, atan2(s,c)*180/M_PI,   fr/k, E, k );

    // project to per leaver
    c2 *=-2;
    double lw     = 2*s*c;       //    |h - 2*c2*a| =  1/(2*s*c) = 1/sin(a)
    double fra    = fr/(lbond[ib.a]*lw);
    double frb    = fr/(lbond[ib.b]*lw);
    Vec3d fa,fb;
    fa.set_lincomb( fra,h,  fra*c2,ha );  //fa = (h - 2*c2*a)*fr / ( la* |h - 2*c2*a| );
    fb.set_lincomb( frb,h,  frb*c2,hb );  //fb = (h - 2*c2*b)*fr / ( lb* |h - 2*c2*b| );

    //printf( "ang[%i] %g (%g,%g,%g) (%g,%g,%g)\n", ig,  c, fa.x,fa.y,fa.z,   fb.x,fb.y,fb.z );

    /*
    //i_DEBUG = 9;
    if(ig==i_DEBUG){


        Vec2i ib_ = ang2bond[ig];
        Vec2i b1,b2;
        b1 = bond2atom[ib.i];
        b2 = bond2atom[ib.j];
        bool flip1,flip2;
        pairs2triple( b1, b2, ang2atom[ig], flip1, flip2 );

        //printf( "ang[%i] (%i,%i)  (%i,%i)%i  (%i,%i)%i   \n", ig,   ib_.j, ib_.j,   b1.i,b1.j,flip1,    b2.i,b2.j,flip2  );
        //printf( "ang[%i] cs(%g,%g)%g  rcs(%g,%g)%g    %g %g %g \n",  c,s, atan2(s,c)*180/M_PI, cs.x,cs.y,atan2(cs.y,cs.x)*180/M_PI,  fr/k, E, k );
        glColor3f(1.0,0.0,0.0); double fsc=100.0;
        Draw3D::drawVecInPos( fa*fsc, apos[ia.x] );
        Draw3D::drawVecInPos( fb*fsc, apos[ia.z] );
        //Draw3D::drawVecInPos( (fa+fb)*-fsc, apos[ia.y] );
        Draw3D::drawVecInPos( fa*-fsc, apos[ia.y] );
        Draw3D::drawVecInPos( fb*-fsc, apos[ia.y] );
        glColor3f(0.0,1.0,0.0);
        Draw3D::drawArrow( apos[ia.x],apos[ia.y], 0.1 );
        Draw3D::drawArrow( apos[ia.y],apos[ia.z], 0.1 );

        glColor3f(0.0,1.0,0.0); Draw3D::drawArrow( apos[ia.b],apos[ia.b]+h, 0.1 );
        glColor3f(1.0,0.0,0.0); Draw3D::drawArrow( apos[ia.b],apos[ia.b]+ha, 0.1 );
        glColor3f(0.0,0.0,1.0); Draw3D::drawArrow( apos[ia.b],apos[ia.b]+hb, 0.1 );

        //Draw3D::drawArrow( apos[ia.x],apos[ia.y], 0.1 );
        //Draw3D::drawArrow( apos[ia.y],apos[ia.z], 0.1 );
        //Draw3D::d
    }
    */


    // to atoms
    aforce[ia.x].add(fa); aforce[ia.y].sub(fa);
    aforce[ia.z].add(fb); aforce[ia.y].sub(fb);

    return E;
}



double eval_torsion(int it){

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

    double ca   = hab.dot(ha);
    double cb   = hab.dot(hb);
    double cab  = ha .dot(hb);
    double sa2  = (1-ca*ca);
    double sb2  = (1-cb*cb);
    double invs = 1/sqrt( sa2*sb2 );
    //double c    = ;  //  c = <  ha - <ha|hab>hab   | hb - <hb|hab>hab    >

    Vec2d cs,csn;
    cs.x = ( cab - ca*cb )*invs;
    cs.y = sqrt(1-cs.x*cs.x); // can we avoid this sqrt ?

    const int n = tors_n[it];
    for(int i=0; i<n-1; i++){
        csn.mul_cmplx(cs);
    }


    // check here : https://www.wolframalpha.com/input/?i=(x+%2B+isqrt(1-x%5E2))%5En+derivative+by+x

    const double k = tors_k[it];
    double E   =  k  *(1-csn.x);
    double dcn =  k*n*   csn.x;
    //double fr  =  k*n*   csn.y;

    //double c   = cos_func(ca,cb,cab);

    //printf( "<fa|fb> %g cT %g cS %g \n", cs.x, cT, cS );

    // derivatives to get forces

    double invs2 = invs*invs;
    dcn *= invs;
    double dcab  = dcn;                           // dc/dcab = dc/d<ha|hb>
    double dca   = (1-cb*cb)*(ca*cab - cb)*dcn;  // dc/dca  = dc/d<ha|hab>
    double dcb   = (1-ca*ca)*(cb*cab - ca)*dcn;  // dc/dca  = dc/d<hb|hab>

    Vec3d fa,fb,fab;

    fa =Vec3dZero;
    fb =Vec3dZero;
    fab=Vec3dZero;

    Mat3Sd J;

    J.from_dhat(ha);    // -- by ha
    J.mad_ddot(hab,fa, dca ); // dca /dha = d<ha|hab>/dha
    J.mad_ddot(hb ,fa, dcab); // dcab/dha = d<ha|hb> /dha

    J.from_dhat(hb);    // -- by hb
    J.mad_ddot(hab,fb, dcb ); // dcb /dhb = d<hb|hab>/dha
    J.mad_ddot(ha ,fb, dcab); // dcab/dhb = d<hb|ha> /dha

    J.from_dhat(hab);         // -- by hab
    J.mad_ddot(ha,fab, dca);  // dca/dhab = d<ha|hab>/dhab
    J.mad_ddot(hb,fab, dcb);  // dcb/dhab = d<hb|hab>/dhab
    // derivative cab = <ha|hb>

    fa .mul( 1/lbond[ib.a] );
    fb .mul( 1/lbond[ib.c] );
    fab.mul( 1/lbond[ib.b] );

    aforce[ia.x].sub(fa);
    aforce[ia.y].add(fa -fab);
    aforce[ia.z].add(fab-fb);
    aforce[ia.w].add(fb);

    /*
    {
        double fsc = 100.0;
        glColor3f(1.0,0.0,1.0);  Draw3D::drawVecInPos(hab.normalized()*E, apos[ia.x]);
        glColor3f(1.0,1.0,1.0);  Draw3D::drawVecInPos(ha.normalized()*E, apos[ia.x]);
        glColor3f(1.0,0.0,0.0);  Draw3D::drawVecInPos(fa*fsc, apos[ia.x]);
        glColor3f(1.0,0.0,0.0);  Draw3D::drawVecInPos(fa*fsc*-1, apos[ia.y]);
        glColor3f(1.0,0.0,0.0);  Draw3D::drawVecInPos(fb*fsc*-1, apos[ia.z]);
        glColor3f(1.0,0.0,0.0);  Draw3D::drawVecInPos(fb*fsc, apos[ia.w]);

    }
    */

    return E;
}




double _bak_eval_torsion(int it){
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
    double E   =  k  *(1-csn.x);
    double fr  =  k*n*   csn.y;


    fa.mul(fr*ira);
    fb.mul(fr*irb);
    /*
    float fsc = 1000.0;
    glColor3f(1.0,0.0,1.0);  Draw3D::drawVecInPos(hab.normalized()*E, apos[ia.x]);
    glColor3f(1.0,1.0,1.0);  Draw3D::drawVecInPos(ha.normalized()*E, apos[ia.x]);
    glColor3f(1.0,0.0,0.0);  Draw3D::drawVecInPos(fa*fsc, apos[ia.x]);
    glColor3f(1.0,0.0,0.0);  Draw3D::drawVecInPos(fa*fsc*-1, apos[ia.y]);
    glColor3f(1.0,0.0,0.0);  Draw3D::drawVecInPos(fb*fsc*-1, apos[ia.z]);
    glColor3f(1.0,0.0,0.0);  Draw3D::drawVecInPos(fb*fsc, apos[ia.w]);
    //printf("tors[%i] (%i,%i,%i,%i)\n", it, ia.x, ia.y, ia.z, ia.w );
    */
    /*
    // to atoms
    aforce[ia.x].add(fa);
    aforce[ia.y].sub((fa+fb)*0.5);
    aforce[ia.z].sub((fb+fa)*0.5);
    aforce[ia.w].add(fb);
    */
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

double eval(){
    cleanAtomForce();
    double Eb = eval_bonds();
    double Ea = eval_angles();
    double Et = eval_torsions();
    printf( "Eb %g Ea %g Et %g\n", Eb, Ea, Et );
    return Eb+Ea+Et;
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
        //printf( "ang[%i] ((%i,%i)(%i,%i))->(%i,%i,%i) (%i,%i)==(%i,%i) \n", i, b1.x,b1.y, b2.x,b2.y, ang2atom[i].x,ang2atom[i].y,ang2atom[i].z,ang2bond[i].x&0xFFFF, ang2bond[i].y&0xFFFF,ang2bond[i].x,        ang2bond[i].y         );
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
        //printf( "tors[%i] ((%i,%i)(%i,%i)(%i,%i))->(%i,%i,%i,%i) (%i,%i,%i)==(%i,%i,%i) \n", i, b1.x,b1.y, b12.x,b12.y, b2.x,b2.y,tors2atom[i].x,tors2atom[i].y,tors2atom[i].z,tors2atom[i].w,tors2bond[i].x&0xFFFF, tors2bond[i].y&0xFFFF, tors2bond[i].z&0xFFFF,tors2bond[i].x,        tors2bond[i].y,        tors2bond[i].z         );
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
