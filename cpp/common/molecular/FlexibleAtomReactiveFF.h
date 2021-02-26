
/*

Flexible Atom sp-hybridization forcefield
 - each atom has 4 atomic orbitals of 3 types { sigma-bond, election-pair, pi-pond }
 - there are two kinds of interaction "onsite" and "bond"
 - "onsite" interaction
    - sigma and epair try to maximize angle between them ( leadign to tetrahedral configuration )
    - pi orbitals try to orhogonalize
 - "bond" interaction
    - spring constant for bonds
    - pi-pi alignment

*/

#ifndef FlexibleAtomReactiveFF_h
#define FlexibleAtomReactiveFF_h

//#include <vector>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
//#include "quaternion.h"
#include "Forces.h"

//include "geom3D.h"
//#include "integerOps.h"
//#define SIGN_MASK 2147483648

#define N_BOND_MAX 4

constexpr static const Vec2d sp_cs0s[] = {
{   0.0d,0.0d}, // 0, 0
{   1.0d,0.0d}, // 1, sp0
{  -1.0d,0.0d}, // 2, sp1
{  -0.5d, 0.86602540378d}, // 3, sp2
{-1/3.0d, 0.94280904158d}, // 4, sp3
};

 // |f-fE|: 27.9227   fE(8.31683,-3.56436,-4.75248)   f(-14.405,6.17289,8.2313)


struct FlexibleAtomType{

    double rbond0 =  0.8;  // Rbond
    double aMorse =  4.0;  // EBond
    double bMorse = -0.7;  // kBond

    double vdW_c6  = -15.0; // c6
    double vdW_R  =  8.0;  // RcutVdW

     //double fHlen =  1.0;
    //double fHlen =  1.6;  // C-H bond lenght = fHlen*l0ae    e.g.: 1.1=(1.6*0.7)
    //double l0ae  =  0.7;  // atom-electron equilibrium length
    //double Kae   =  20.0; // atom-electron strenght

    double Wee  =  1.0;   // electon-electron width
    double Kee  = 10.0;  // electon-electron strenght
    double Kpp  = 10.0;  // onsite  p-p orthogonalization strenght

    double Kpi  =  2.0;   // offsite p-p alignment         strenght
    //double Kpe  = -15.0;  // onsite  p-s orthogonalization strenght



    void print(){
        //printf( "nbond  %i rbond0 %g\n", nbond, rbond0 );
        printf( "aMorse %g bMorse %g\n", aMorse, bMorse );
        printf( "c6     %g r2vdW  %g\n", vdW_c6, vdW_R  );
        //exit(0);
    }

};

struct FlexiblePairType{
    //int    nbond = 4;  // number bonds

    double rbond0 =  1.6;  // Rbond
    double aMorse =  4.0;  // EBond
    double bMorse = -1.0*2;  // kBond

    double Kpi    =   0.5;  // kPz
    //double vdW_R2 =  8.0;  // RcutVdW
    double vdW_c6 = -15.0; // c6
    double vdW_w6 =  pow6( 2.0 );

    inline void combine(const FlexibleAtomType& a, const FlexibleAtomType& b ){
        aMorse  = a.aMorse * b.aMorse;  // TODO
        bMorse  = a.bMorse + b.bMorse;
        rbond0  = a.rbond0 + b.rbond0;

        Kpi    = a.Kpi * b.Kpi;
        vdW_c6 = -(a.vdW_c6  * b.vdW_c6);
        //vdW_w6 = pow6( a.vdw_R + b.vdw_R );
        vdW_w6 = pow3( a.vdW_R * b.vdW_R * 2 ); // replace aritmetic average by geometric ?
    }

    void print(){
        //printf( "nbond  %i rbond0 %g\n", nbond, rbond0 );
        printf( "aMorse %g bMorse %g\n", aMorse, bMorse );
        printf( "c6     %g r2vdW  %g\n", vdW_c6, vdW_w6 );
        //exit(0);
    }

};



class FARFF{ public:

    bool   substract_LJq = true;
    double Eatoms=0,Epairs=0;

    FlexibleAtomType atype0;
    FlexiblePairType ptype0;

    int natom   = 0;
    int nbond   = 0; // number of all bonding orbitals
    int norb    = 0;
    int nDOF    = 0; // natoms + neps + npis

    //int ncap    = 0; // number of capping
    //int nepair  = 0; // number of electron pairs
    //int nporb   = 0; // number of electron pi orbitals

    //Vec2i  * bondIs = 0;
    //Vec2d  * bLKs   = 0;  // lenght[A],K [eV/A]

    Vec3d  * dofs   = 0; // degrees of freedom
    Vec3d  * fdofs  = 0; // degrees of freedom
    //double * edofs  = 0;

    Vec3d  * apos    = 0;      // atomic position // ALIAS
    Vec3d  * aforce  = 0;      // atomic forces   // ALIAS
    double * aenergy = 0;

    //Vec3d  * capREQs = 0;
    //Vec3d  * aREQs   = 0;
    //Vec3i  * aconf   = 0;     // nH, nsigna,nt=(nsigna+ne)     npi = (4 - vt)=(4 - (nsigma+ne))
    Vec3ui8  * aconf   = 0;

    //Quat4i * atom2bond = 0; // atom to bond map

    //double * blen    = 0;   // bond lengths
    Vec3d  * opos     = 0;     // normalized bond unitary vectors  // ALIAS
    Vec3d  * oforce   = 0;
    double * oenergy  = 0;
    //double * invRs   = 0;
    //Vec3d  * bforce  = 0;     // forces acting on it              // ALIAS
    //double * bls     = 0;

    //Vec3d * capPos   = 0;
    //Vec3d * capForce = 0;

    bool * ignoreAtoms = 0;

    void realloc(int natom_){
        natom=natom_;
        norb =natom*N_BOND_MAX;
        nDOF =natom + norb;
        printf(  "FARFF realoc na %i \n", natom );

        _realloc(oenergy,norb );

        _realloc(aconf  ,natom);
        _realloc(aenergy,natom);
        _realloc(ignoreAtoms,natom);

        _realloc(dofs  ,nDOF);
        _realloc(fdofs ,nDOF);

        apos    = dofs;
        aforce  = fdofs;

        opos    =  dofs+natom;
        oforce  = fdofs+natom;

        for(int i=0; i<natom; i++){ ignoreAtoms[i]=false; }

    }

// ======== Force Evaluation

void cleanForce(){
    for(int i=0;i<nDOF; i++){ fdofs  [i].set(0.); }
    for(int i=0;i<natom;i++){ aenergy[i]=0;       }
}

double evalAtom(int ia){
    //printf( "atom[%i] \n", ia );

    const FlexibleAtomType& type = atype0;
    const Vec3ui8& conf = aconf[ia];
    const Vec3d&   pa   = apos [ia];
    const int ih        = ia*N_BOND_MAX;
    Vec3d*  hs  = opos  +ih;
    Vec3d*  fs  = oforce+ih;

    //Vec3d fa    = Vec3dZero;
    double E    = 0;
    const double w2ee   = sq(atype0.Wee);
    int    nsigma = conf.a;

    Vec2d cs0 = sp_cs0s[nsigma];
    //printf( "cs0.x %g \n", cs0.x );
    //cs0.x = 0.8; nsigma=3;

    // -- normalize all orbital directions
    //for(int i=0; i<N_BOND_MAX; i++){
    //    //hs[i].normalize();
    //    hs[i].normalize_taylor3();
    //}

    //int j_DEBUG=0;
    // -- repulsion between sigma bonds
    for(int i=1; i<nsigma; i++){
        const Vec3d& hi = hs[i];
        Vec3d&       fi = fs[i];
        for(int j=0; j<i; j++){ // electron-electron
            E += evalCos2  (hi,hs[j],fi,fs[j],type.Kee,cs0.x);
            //E += evalCos2_o(hi,hs[j],fi,fs[j],type.Kee,cs0.x);

            //Mat3Sd Ji,Jj;   // This gives exactly the same result for normalized vectors
            //Ji.from_dhat(hi   ); Ji.mad_ddot(hs[j],fs[i], 2*type.Kee*c );
            //Jj.from_dhat(hs[j]); Jj.mad_ddot(hi   ,fs[j], 2*type.Kee*c );
            //E+=type.Kee*c*c;

        }
    }

    // -- orthogonalization with p-orbitals
    for(int i=nsigma; i<N_BOND_MAX; i++){
        const Vec3d& hi = hs[i];
        Vec3d&       fi = fs[i];
        for(int j=0; j<i; j++){
            E += evalCos2   ( hi, hs[j], fi, fs[j], type.Kpp, 0);
        }
    }

    return E;
}




double evalPair( int ia, int ja, FlexiblePairType& type){
//double evalPair( int ia, int ja, int nbi, int nbj ){
    const Vec3ui8& confi = aconf[ia];
    const Vec3ui8& confj = aconf[ja];
    int nbi = confi.a;
    int nbj = confj.a;

    Vec3d  hij; hij.set_sub( apos[ja], apos[ia] );   // = apos[ja] - apos[ia];

    double r2   = hij.norm2() + R2SAFE;
    double r    = sqrt( r2 );
    double invr = 1/r;
    hij.mul(invr); // = dij*(1/rij);

    double r4     = r2*r2;
    double r6     = r4*r2;
    double invVdW = 1/( r6 + type.vdW_w6 );
    double EvdW   = type.vdW_c6*invVdW;
    double fvdW   = -6*EvdW*invVdW*r4*r;

    // todo: replace this by other short range force ... so we can cutoff bonds
    double expar =    exp( type.bMorse*( r - type.rbond0 ) );
    double Esr   =    type.aMorse*expar*expar;
    double Eb    = -2*type.aMorse*expar;
    double fsr   =  2*type.bMorse*Esr;
    double frb   =    type.bMorse*Eb;

    double E = Esr + EvdW;

    Vec3d force; force.set_mul( hij, fsr + fvdW );
    //printf( "r[%i,%i] %g   | Rb %g   EvdW %g \n", ia,ja, r,  type.rbond0, EvdW );

    //DEBUG
    //force.set(.0);
    //E=Eb+Esr;
    //E+=Eb;
    //force.add_mul(hij, frb );
    //printf( "pair[%i,%i] r %g \n", ia, ja,  r );
    //aforce[ia].add(force); aforce[ja].sub(force); return E;



    Mat3Sd Jbond;
    Jbond.from_dhat(hij);

    const int ioff= ia*N_BOND_MAX;
    const int joff= ja*N_BOND_MAX;
    Vec3d* his = opos  +ioff;
    Vec3d* hjs = opos  +joff;
    Vec3d* fis = oforce+ioff;
    Vec3d* fjs = oforce+joff;

    double* eis = oenergy+ioff;
    double* ejs = oenergy+joff;


    //double ccut    = 0.8;
    //double invcut  = 1-ccut;
    //double invcut2 = invcut*invcut;


    for(int ib=0; ib<nbi; ib++){

        const Vec3d& hi = his[ib];
              Vec3d& fi = fis[ib];

        double       ci = hij.dot( hi );   // ci = <hi|hij>

        if(ci<=0) continue;
        //if(ci<ccut) continue;

        //bool capi = ( capis[ib] >= 0 );

        for(int jb=0; jb<nbj; jb++){

            const Vec3d& hj = hjs[jb];
                  Vec3d& fj = fjs[jb];

            double cj       = hij.dot( hj );  // cj  = <hj|hij>

            // avoid symmetric image of bond
            //double cij      = hi .dot( hj );  // cij = <hj|hi>
            //if( (cj>0)||(cij>0) ) continue;
            //double cc  = ci*cj*cij;

            if(cj>=0) continue; // avoid symmetric image of bond
            double cc  = ci*cj;

            double cc2 = cc*cc;
            double e   = cc2*cc2;
            double de  = 4*cc2*cc;

            //if(cc<ccut)continue;
            //double ccm = cc-ccut;
            //double e  = invcut2*ccm*ccm;
            //double de = 2*invcut2*ccm;

            // # = w*cc/(1+cc+w) =   w*(1+cc+w-(1+w))/(1+cc+w) = w - w*(1+w)/(1-cc+w)
            // # = w*cc/(1-cc+w) =  -w*(1-cc+w-(1+w))/(1-cc+w) = w + w*(1+w)/(1-cc+w)
            //const double wBond  = 0.05;
            //const double wwBond = wBond*(1+wBond);
            //double invcc   = 1/(1+cc+wBond);
            //double invccww = wwBond*invcc;
            //double e       = (invccww + wBond);
            //double de      = -invccww*invcc;


            double eEb = e*Eb*0.5d;
            eis[ib]+=eEb;
            ejs[jb]+=eEb;
            E += eEb+eEb;


            double deEb     =    Eb*de;
            //double deEbcicj = -deEb*ci*cj;
            double deEbinvr =  deEb*invr;

            //Vec3d fi_=Vec3dZero,fj_=Vec3dZero,force_=Vec3dZero,force__=Vec3dZero; // DEBUG

            //force.add_mul( hi, -deEbinvr*cj );
            //force.add_mul( hj, -deEbinvr*ci );
            Jbond.mad_ddot( hi, force, deEbinvr*cj );
            Jbond.mad_ddot( hj, force, deEbinvr*ci );
            force.add_mul(hij, e*frb);

            fi.add_mul( hij, -cj*deEb );
            fj.add_mul( hij, -ci*deEb );

            //force.add(force_);
        }
    }


    // ---------- Pi-Pi interaction
    //printf( "%i,%i nbi %i nbj %i \n", ia, ja, nbi, nbj );
    if( (nbi==3)&&(nbj==3) ){ // align pz-pz in sp2-sp2 pair
        const Vec3d& hi = his[3];
        const Vec3d& hj = hjs[3];
        Vec3d&       fi = fis[3];
        Vec3d&       fj = fjs[3];
        double c = hi.dot(hj);
        //double E    = Epz * ( cdot*cdot );
        double de = -2*c*type.Kpi*Eb;
        //printf( "cdot %g Epz %g de %g \n", cdot, type.Epz, de );
        fi.add_mul(hj,de);
        fj.add_mul(hi,de);
        //force      // TODO: proper derivatives of energy
    }

    aforce[ia].add(force);
    aforce[ja].sub(force);
    return E;
}

double evalPairs(){
    Epairs = 0;
    for(int i=0; i<natom; i++){
        if(ignoreAtoms[i])continue;
        //FlexiblePairType& typei = *types[i];    // TODO add pairs later
        for(int j=i+1; j<natom; j++){
            if(ignoreAtoms[j])continue;
            //FlexiblePairType pairType;
            //pairType.combine( typei, *types[j] );
            //E += pairEF( i, j, typei.nbond, types[j]->nbond, pairType );
            Epairs += evalPair( i, j, ptype0 );
        }
    }
    return Epairs;
}

inline double evalAtoms(){
    Eatoms=0;
    for(int i=0; i<natom; i++){
        if(ignoreAtoms[i])continue;
        Eatoms+=evalAtom(i);
    }
    return Eatoms;
}


void normalizeOrbs(){
    // ToDo : is this numerically stable? if normal forces are very hi ?
    for(int i=0; i<norb; i++){
        //opos[i].normalize();
        if(ignoreAtoms[i>>2])continue;
        opos[i].normalize_taylor3();
    }
}

void transferOrbRecoil(){
    // ToDo : is this numerically stable? if normal forces are very hi ?
    for(int i=0; i<norb; i++){
        //aforce[i>>2].sub( oforce[i] );
        int ia=i>>2;
        if(ignoreAtoms[ia])continue;
        Vec3d f = oforce[i]*-1;
        aforce[ia].add( f );

        //glColor3f(1,1,0); Draw3D::drawVecInPos( oforce[i], apos[ia] );
        if((ia==0)&&(i-(ia<<2)==0) ){ glColor3f(1,1,0); Draw3D::drawVecInPos( f, apos[ia] ); }
    }
}

void removeNormalForces(){
    // ToDo : is this numerically stable? if normal forces are very hi ?
    for(int i=0; i<norb; i++){
        if(ignoreAtoms[i>>2])continue;
        oforce[i].makeOrthoU(opos[i]);
        //Vec3d fT = opos[i]*opos[i].dot(oforce[i]);
        //oforce[i   ].sub( fT );
        //aforce[i>>2].add( fT );
    }
}

double eval(){
    //cleanForce();
    normalizeOrbs();
    Eatoms=evalAtoms(); // do this first to orthonormalize ?
    //transferOrbRecoil();
    Epairs=evalPairs();
    //transferOrbRecoil();
    removeNormalForces();
    return Eatoms + Epairs;
}

void moveGD(double dt, bool bAtom, bool bOrbital ){
    //for(int i=0; i<nDOF; i++){ dofs[i].add_mul( fdofs[i], dt ); }
    if(bAtom   )for(int i=0; i<natom; i++){
            if(ignoreAtoms[i])continue;
            apos[i].add_mul( aforce[i], dt );
    }
    if(bOrbital)for(int i=0; i<norb;  i++){
        if(ignoreAtoms[i>>2])continue;
        opos[i].add_mul( oforce[i], dt );
    }
}


// ============== BACKUP

/*
void evalAtom(int ia){
    const Vec3ui8& conf = aconf[ia];
    const Vec3d&   pa   = apos [ia];
    int ih      = ia<<2;
    Vec3d*  hs  = opos  +ih;
    Vec3d*  fs  = oforce+ih;
    //double* irs = invRs +ih;

    Vec3d fa    = Vec3dZero;
    //double W2ee = Wee*Wee;
    double W2ee = sq(atype0.Wee);

    //aforce[ia].set(0.); // DEBUG
    //return;

    for(int i=conf.b; i<N_BOND_MAX; i++){ hs[i].normalize(); } // TODO: this can be made perhaps more efficiently

    for(int i=0; i<N_BOND_MAX; i++){
        Vec3d&       fi = fs[i];

        // non-pi
        if( i<conf.b ){  // non pi
            const Vec3d& pi = hs[i];
            Vec3d d; double l,fr;

            for(int j=i+1; j<conf.b; j++){ // electron-electron
                d.set_sub( hs[j], pi);
                double r2   = d.norm2();
                double ir2  = 1/( r2 + W2ee );
                double fr   = atype0.Kee*ir2*ir2;
                d.mul(fr);

                // Drift because bonds are not updated
                fi   .add(d);
                fs[j].sub(d);

                //glColor3f(1.0,0.0,0.0); Draw3D::drawVecInPos( fja*0.5,apos[ia]+hj*0.5 );

            }

            if(conf.b<N_BOND_MAX){
                Vec3d hi   = pi-pa;
                Mat3Sd D;
                D.from_dhat(hi); // used only if pi-bond used

                double ir2 = 1/(hi.norm2() + 1e-16);
                double ir  = sqrt(ir2);
                hi.mul(ir);

                for(int j=conf.b; j<N_BOND_MAX; j++){  // pi-electron i=sigma j=pi
                    //if(i!=2)continue; // DEBUG
                    Vec3d fia,fja;
                    // Note: moving this here save D-mat calculation
                    // ToDo : proper dhat
                    const Vec3d& hj = hs[j];
                    double c   = hi.dot(hj);
                    double dfc = atype0.Kpe*2*c;

                    //printf( " %i(%i,%i) c %g (%g,%g,%g) (%g,%g,%g) \n", ia, i,j, c, hi.x,hi.y,hi.z,  hj.x,hj.y,hj.z );
                    D.dhat_dot( hj, fia );
                    fia.mul(dfc*ir2*ir);

                    //fia.set_mul(hj,dfc);
                    fja.set_mul(hi,dfc);
                    fja.add_mul(hj,-c*dfc);    //printf( "<fj,fj> %g \n", hj.dot(fja) );
                    //fia.set(0.);
                    //fja.set(1.0,0.0,0.0);
                    //printf( "%i(%i,%i) (%g,%g,%g) (%g,%g,%g) \n", ia, i,j, fia.x,fia.y,fia.z, fja.x,fja.y,fja.z );

                    //glColor3f(1.0,0.0,0.0); Draw3D::drawVecInPos( fja*0.5,apos[ia]+hj*0.5 );

                    fi.add(fia); fs[j].add(fja);
                    fa.sub(fia); //fa   .sub(fja);
                }
            }

            // NOTE: moved here so we can out-project radial force component

            // atom-electorn
            d.set_sub(pi,pa);
            if(true){ // remove radial force component from "fi"
            //if(false){ // remove radial force component from "fi"
                l  = d.normalize();
                fr = atype0.Kae*(l-atype0.l0ae);
                // outproject radial component
                Vec3d frad;
                frad.set_mul(d,d.dot(fi));
                fa.add(frad);
                fi.sub(frad);
            }else{
                l  = d.norm();
                fr = atype0.Kae*(l-atype0.l0ae)/l;
            }
            d.mul(fr);
            fi.sub(d);
            fa.add(d);

        }else{           // pi

            double c,dfc;
            Vec3d fia,fja;
            const Vec3d& hi = hs[i];

            for(int j=i+1; j<N_BOND_MAX; j++){ // pi-pi
                const Vec3d& hj = hs[j];
                c   = hi.dot(hj);
                dfc = atype0.Kpp*2*c;
                fia.set_mul(hj,dfc);   fia.add_mul(hi,-c*dfc);
                fja.set_mul(hi,dfc);   fja.add_mul(hj,-c*dfc);

                fi.add(fia); fs[j].add(fja);
                //fa.sub(fia); fa   .sub(fja);
            }

        }

    }

    aforce[ia].add(fa);

    //glColor3f(1.0,0.0,0.0); Draw3D::drawVecInPos( aforce[ia], pa );

    //aforce[ia].set(0.); // DEBUG

}

*/


/*

double evalPair( int ia, int ja, FlexiblePairType& type){
//double evalPair( int ia, int ja, int nbi, int nbj ){
    const Vec3ui8& confi = aconf[ia];
    const Vec3ui8& confj = aconf[ja];
    int nbi = confi.a;
    int nbj = confj.a;

    Vec3d  hij; hij.set_sub( apos[ja], apos[ia] );   // = apos[ja] - apos[ia];
    double r2   = hij.norm2() + R2SAFE;
    double r    = sqrt( r2 );
    double invr = 1/r;
    hij.mul(invr); // = dij*(1/rij);


    double dxy=-hij.x*hij.y*invr;
    double dxz=-hij.x*hij.z*invr;
    double dyz=-hij.y*hij.z*invr;
    double xx  = hij.x*hij.x;
    double yy  = hij.y*hij.y;
    double zz  = hij.z*hij.z;
    double dxx=(yy+zz)*invr;
    double dyy=(xx+zz)*invr;
    double dzz=(xx+yy)*invr;


    //double ir3 = 1/(r*r2);
    //double xx  = dij.x*dij.x;
    //double yy  = dij.y*dij.y;
    //double zz  = dij.z*dij.z;

    //double dxy=-dij.x*dij.y*ir3;
    //double dxz=-dij.x*dij.z*ir3;
    //double dyz=-dij.y*dij.z*ir3;
    //double dxx=(yy+zz)*ir3;
    //double dyy=(xx+zz)*ir3;
    //double dzz=(xx+yy)*ir3;


    //double ir2vdW = 1/(r2 + type.vdW_w6 );


    double r6     = r2*r2*r2;
    double invVdW = 1/( r6 + type.vdW_w6 );
    double EvdW   = type.vdW_c6*invVdW;
    double fvdW   = -6*EvdW*invVdW*r6;

    // todo: replace this by other short range force ... so we can cutoff bonds
    double expar =  exp( type.bMorse*( r - type.rbond0 ) );
    double Esr   =    type.aMorse*expar*expar;
    double Eb    = -2*type.aMorse*expar;
    double fsr   =  2*type.bMorse*Esr;
    double frb   =    type.bMorse*Eb;

    double E = Esr + EvdW;

    Vec3d force; force.set_mul( hij, fsr + fvdW );
    //printf( "r[%i,%i] %g   | Rb %g   EvdW %g \n", ia,ja, r,  type.rbond0, EvdW );

    //DEBUG
    //force.add_mul(hij, frb );
    //nbi=nbj=1;


    Mat3Sd Jbond;
    Jbond.from_dhat(hij);

    const int ioff= ia*N_BOND_MAX;
    const int joff= ja*N_BOND_MAX;
    Vec3d* his = opos  +ioff;
    Vec3d* hjs = opos  +joff;
    Vec3d* fis = oforce+ioff;
    Vec3d* fjs = oforce+joff;

    double* eis = oenergy+ioff;
    double* ejs = oenergy+joff;


    double ccut    = 0.8;
    double invcut  = 1-ccut;
    double invcut2 = invcut*invcut;


    for(int ib=0; ib<nbi; ib++){

        const Vec3d& hi = his[ib];
              Vec3d& fi = fis[ib];

        double       ci = hij.dot( hi );   // ci = <hi|hij>

        if(ci<0) continue;
        //if(ci<ccut) continue;

        //bool capi = ( capis[ib] >= 0 );

        for(int jb=0; jb<nbj; jb++){

            const Vec3d& hj = hjs[jb];
                  Vec3d& fj = fjs[jb];

            double cj       = hij.dot( hj );  // cj  = <hj|hij>
            double cij      = hi .dot( hj );  // cij = <hj|hi>

            if( (cj>0)||(cij>0) ) continue; // avoid symmetric image of bond

            // --- more-concentrated ToDo - use lorenz ?
            double cc  = ci*cj*cij;
            //double cc  = -ci*cj;


            //double cc2 = cc*cc;
            //double e   = cc2*cc2;
            //double de  = 4*cc2*cc;



            //if(cc<ccut)continue;
            //double ccm = cc-ccut;
            //double e  = invcut2*ccm*ccm;
            //double de = 2*invcut2*ccm;


            // # = w*cc/(1+cc+w) =   w*(1+cc+w-(1+w))/(1+cc+w) = w - w*(1+w)/(1-cc+w)
            // # = w*cc/(1-cc+w) =  -w*(1-cc+w-(1+w))/(1-cc+w) = w + w*(1+w)/(1-cc+w)
            const double wBond  = 0.1;
            const double wwBond = wBond*(1+wBond);
            double invcc   = 1/(1-cc+wBond);
            double invccww = wwBond*invcc;
            double e       = invccww - wBond;
            double de      = invccww*invcc;


            double eEb = e*Eb*0.5;
            eis[ib]+=eEb;
            ejs[jb]+=eEb;
            E += eEb+eEb;


            double deEb     =    Eb*de;
            double deEbcij  =  deEb*cij;
            double deEbcicj = -deEb*ci*cj;
            double deEbcijinvr = deEbcij*invr;

            //Jbond.mad_ddot( hi, force, deEbcijinvr*cj );
            //Jbond.mad_ddot( hj, force, deEbcijinvr*ci );
            //force.add_mul ( hij, frb*e );
            //fi.add_lincomb( -cj*deEbcij, hij,    deEbcicj,  hj );
            //fj.add_lincomb( -ci*deEbcij, hij,    deEbcicj,  hi );



            Vec3d fi_,fj_,force_; // DEBUG
            Jbond.mad_ddot( hi, force_, deEbcijinvr*cj );
            Jbond.mad_ddot( hj, force_, deEbcijinvr*ci );
            force_.add_mul ( hij, frb*e );
            fi_.set_lincomb( -cj*deEbcij, hij,    deEbcicj,  hj );
            fj_.set_lincomb( -ci*deEbcij, hij,    deEbcicj,  hi );
            // without orthogonalized force




            //printf( "a[%i,%i]o[%i,%i] cc %g c(%g,%g,%g) e %g  f(%g,%g,%g) \n", ia,ja,ib,jb, cc,  ci,cj,cij,  e,  force.x,force.y,force.z );

            printf( "a[%i,%i]o[%i,%i] cc %g c(%g,%g,%g) e %g \n", ia,ja,ib,jb, cc,  ci,cj,cij,  e,  force.x,force.y,force.z );


            force_.x = Eb*de*cij*( cj*( hi.x*dxx + hi.y*dxy + hi.z*dxz )     +    ci*( hj.x*dxx + hj.y*dxy + hj.z*dxz )   )  + hij.x*frb*e;
            force_.y = Eb*de*cij*( cj*( hi.x*dxy + hi.y*dyy + hi.z*dyz )     +    ci*( hj.x*dxy + hj.y*dyy + hj.z*dyz )   )  + hij.y*frb*e;
            force_.z = Eb*de*cij*( cj*( hi.x*dxz + hi.y*dyz + hi.z*dzz )     +    ci*( hj.x*dxz + hj.y*dyz + hj.z*dzz )   )  + hij.z*frb*e;

            fi_.x = -( cij*cj*hij.x + ci*cj*hj.x )*de*Eb;
            fi_.y = -( cij*cj*hij.y + ci*cj*hj.y )*de*Eb;
            fi_.z = -( cij*cj*hij.z + ci*cj*hj.z )*de*Eb;

            fj_.x = -( cij*ci*hij.x + ci*cj*hi.x )*de*Eb;
            fj_.y = -( cij*ci*hij.y + ci*cj*hi.y )*de*Eb;
            fj_.z = -( cij*ci*hij.z + ci*cj*hi.z )*de*Eb;


            fi_.makeOrthoU(hi);
            fj_.makeOrthoU(hj);
            fi.add(fi_); fj.add(fj_);
            //aforce[ia].add(fi_); aforce[ja].add(fj_);
            force.add(force_);
            //if((ia==0)&&(ja==1)){
            if(ia==0){
                glColor3f(1,0,0);
                Draw3D::drawVecInPos( force_, apos[ia] );
                glColor3f(1,0,1);
                Draw3D::drawVecInPos( fi_, apos[ia]+hi );
                Draw3D::drawVecInPos( fj_, apos[ja]+hj );
                //Draw3D::drawVecInPos( hij*(-cj*deEbcij) + hj*(deEbcicj), apos[ia]+hi );
                //Draw3D::drawVecInPos( hij*(-ci*deEbcij) + hi*(deEbcicj), apos[ia]+hj );
            }

        }
    }

*/


}; // FFsp

#endif
