
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

#ifndef FFsp_h
#define FFsp_h

//#include <vector>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
//#include "quaternion.h"
//#include "Forces.h"

//include "geom3D.h"
//#include "integerOps.h"
//#define SIGN_MASK 2147483648

#define N_BOND_MAX 4

struct Mat3Sd{ // symmetric 3x3 matrix

    double xx,yy,zz,xy,xz,yz;

    inline void from_dhat(const Vec3d& h){
        //double ir  = irs[i];
        double hxx  = h.x*h.x;
        double hyy  = h.y*h.y;
        double hzz  = h.z*h.z;
        xy=-h.x*h.y;
        xz=-h.x*h.z;
        yz=-h.y*h.z;
        xx=(hyy+hzz);
        yy=(hxx+hzz);
        zz=(hxx+hyy);
    }

    inline void dhat_dot( const Vec3d& h, Vec3d& f ){
        f.x = h.x*xx + h.y*xy + h.z*xz;
        f.y = h.x*xy + h.y*yy + h.z*yz;
        f.z = h.x*xz + h.y*yz + h.z*zz;
    }

};


inline void addAtomicForceLJQ( const Vec3d& dp, Vec3d& f, double r0, double eps, double qq ){
    //Vec3f dp; dp.set_sub( p2, p1 );
    const double COULOMB_CONST = 14.3996448915d;  //  [V*A/e] = [ (eV/A) * A^2 /e^2]
    double ir2  = 1/( dp.norm2() + 1e+4 );
    double ir   = sqrt(ir2);
    double ir2_ = ir2*r0*r0;
    double ir6  = ir2_*ir2_*ir2_;
    double fr   = ( ( 1 - ir6 )*ir6*12*eps + ir*qq*-COULOMB_CONST )*ir2;
    f.add_mul( dp, fr );
}

class FspFF{ public:

    bool substract_LJq = true;

    double Wel  =  1.0;
    double le0  =  0.7;
    double Kea  =  1.0;

    double Kss  = -1.0;  // angular force between sigma-sigma
    double Kse  = -1.0;  // angular force between sigma-epair
    double Kee  = -1.0;  // angular force between epair-epair
    double Kp   = -1.0;  // angular force between pi-pi
    double Kpi  =  1.0;  // angular force pi bond

    int natom   = 0;
    int nbond   = 0; // number of all bonding orbitals
    int norb    = 0;
    int ncap    = 0;
    int nDOF    = 0; // natoms + neps + npis
    //int nsigmas = 0; // number of bonds between atoms
    //int npis    = 0; // number of electron pi orbitals
    //int neps    = 0; // number of electron pairs

    Vec2i  * bondIs = 0;
    Vec2d  * bLKs   = 0;  // lenght[A],K [eV/A]

    Vec3d  * dofs   = 0;  // degrees of freedom
    //Vec3d  * vdofs = 0;
    Vec3d  * fdofs  = 0; // degrees of freedom

    Vec3d  * apos   = 0;      // atomic position // ALIAS
    Vec3d  * aforce = 0;      // atomic forces   // ALIAS

    Vec3d  * capREQs = 0;
    Vec3d  * aREQs   = 0;
    Vec3i  * aconf   = 0;     // nH, nsigna,nt=(nsigna+ne)     npi = (4 - vt)=(4 - (nsigma+ne))
    //Quat4i * atom2bond = 0; // atom to bond map

    //double * blen    = 0;   // bond lengths
    Vec3d  * hdir    = 0;     // normalized bond unitary vectors  // ALIAS
    Vec3d  * hforce  = 0;
    double * invRs   = 0;
    //Vec3d  * bforce  = 0;     // forces acting on it              // ALIAS
    //double * bls     = 0;

    Vec3d * capPos   = 0;
    Vec3d * capForce = 0;

    void realloc(int natom_, int nbond_, int ncap_){
        natom=natom_;
        ncap =ncap_;
        nbond=nbond_;
        norb =natom*N_BOND_MAX;
        nDOF =natom + norb;
        //_realloc(atoms,natom   );
        _realloc(aconf  ,natom);
        _realloc(aREQs  ,natom);

        _realloc(bondIs ,nbond);
        _realloc(bLKs   ,nbond);

        _realloc(dofs  ,nDOF);
        _realloc(fdofs ,nDOF);

        _realloc(invRs ,norb );

        apos   = dofs;
        aforce = fdofs;

        hdir   =  dofs+natom;
        hforce = fdofs+natom;

        _realloc( capPos   ,ncap );
        _realloc( capForce ,ncap );
        _realloc( capREQs  ,ncap );
    }

void cleanForce(){
    for(int i=0; i<nDOF; i++){ fdofs[i].set(0.); }
}

void projectBondCenter(){
    for(int ib=0; ib<nbond; ib++){
        const Vec2i& b   = bondIs[ib];
        int ia = b.i>>2;
        int ja = b.j>>2;
        Vec3d ph =(apos[ja]+apos[ia])*0.5; // ToDo: optimize ?
        hdir[b.i]=ph;
        hdir[b.j]=ph;
    }
}

void evalBond(int ibond){
    const Vec2d& lk  = bLKs  [ibond];
    const Vec2i& b   = bondIs[ibond];
    int ia = b.i>>2;
    int ja = b.j>>2;

    //if(ibond>0) return;

    //printf( "bond %i atoms(%i,%i) hs(%i,%i) \n", ibond, ia, ja, b.i, b.j );
    Vec3d dp; dp.set_sub( apos[ja], apos[ia] );
    Vec3d f;
    double l     = dp.norm();
    double inv_l = 1.0d/l;
    f.set_mul(dp,inv_l);

    /*
    //lbond [ib] = l;
    hdir[b.i].set    (f   ); invRs[b.i] = inv_l;
    //hdir[b.j].set    (f   ); invRs[b.j] = inv_l;
    hdir[b.j].set_mul(f,-1); invRs[b.j] = inv_l;
    */

    /*
    Vec3d ph =(apos[ja]+apos[ia])*0.5; // ToDo: optimize ?
    hdir[b.i].set(ph);
    hdir[b.j].set(ph);
    */

    //aforce[ia].add(f);
    //aforce[ja].sub(f);
    //printf( "bond %i hs(%i,%i) dp(%g,%g,%g) %g \n", ibond, b.i, b.j, dp.x, dp.y, dp.z, l );
    //printf( "bond %i hs(%i,%i) h(%g,%g,%g) %g \n", ibond, b.i, b.j, f.x, f.y, f.z, l );

    f.mul( (l-lk.a)*lk.b );

    f.set(0.); // DEBUG : remove bond force


    aforce[ia].add(hforce[b.j]);
    aforce[ja].add(hforce[b.i]);
    //f.add(hforce[b.i]); f.add(hforce[b.j]);
    //f.add_mul(hforce[b.i],0.5); f.add_mul(hforce[b.j],0.5);
    hforce[b.i].set(0.); hforce[b.j].set(0.); // DEBUG : should not affect anything

    glColor3f(0.0,0.0,1.0); Draw3D::drawVecInPos( f, hdir[b.i] );

    /*
    if( substract_LJq ){
        addAtomicForceLJQ( dp, f, aREQs[ia].x+aREQs[ja].x, -aREQs[ia].y*aREQs[ja].y, aREQs[ia].z*aREQs[ja].z );
        //addAtomicForceMorseQ( dp, f, aREQ[iat.x].x+aREQ[iat.y].x, -aREQ[iat.x].y*aREQ[iat.y].y, aREQ[iat.x].z*aREQ[iat.y].z, gridFF.alpha );
    }
    */

    aforce[ia].add(f);
    aforce[ja].sub(f);
}

inline void evalBondPiForce(int ibond){
    // NOTE: this cannot be in bondForce, since some hdirs may not be updated yet
    const Vec2i& b = bondIs[ibond];
    int ia = b.i>>2;
    int ja = b.j>>2;
    if((aconf[ia].b==4)&&(aconf[ia].b==4)){
        int i = (ia<<2)+3;
        int j = (ja<<2)+3;
        Vec3d  hi  = hdir[i];
        Vec3d  hj  = hdir[j];
        double c   = hi.dot(hj);
        double dfc = Kpi*2*c;
        Vec3d  fia = hj*dfc;
        Vec3d  fja = hi*dfc;
        aforce[ia].add(fia); hforce[i].sub(fia);
        aforce[ja].add(fja); hforce[j].sub(fja);
    }
}

void evalAtom(int ia){

    const Vec3i& conf = aconf[ia];
    const Vec3d& pa   = apos [ia];
    int ih      = ia<<2;
    Vec3d*  hs  = hdir  +ih;
    Vec3d*  fs  = hforce+ih;
    double* irs = invRs +ih;

    Vec3d fa    = Vec3dZero;
    double W2el = Wel*Wel;

    //aforce[ia].set(0.); // DEBUG
    //return;

    for(int i=conf.b; i<N_BOND_MAX; i++){ hs[i].normalize(); } // TODO: this can be made perhaps more efficiently

    for(int i=0; i<N_BOND_MAX; i++){
        Vec3d&       fi = fs[i];

        // non-pi
        if( i<conf.b ){  // non pi
            const Vec3d& pi = hs[i];

            // atom-electorn
            Vec3d d;
            d.set_sub(pi,pa);
            double l  = d.norm();
            double fr = Kea*(l-le0)/l;
            d.mul(fr);
            fi.sub(d);
            fa.add(d);

            for(int j=i+1; j<conf.b; j++){ // electron-electron
                d.set_sub( hs[j], pi);
                double r2   = d.norm2();
                double ir2  = 1/( r2 + W2el );
                double fr   = Kee*ir2*ir2;
                d.mul(fr);

                // Drift because bonds are not updated
                fi   .add(d);
                fs[j].sub(d);

                //glColor3f(1.0,0.0,0.0); Draw3D::drawVecInPos( fja*0.5,apos[ia]+hj*0.5 );

            }
/*
            if(conf.b<N_BOND_MAX){
                Vec3d hi   = pi-pa;
                Mat3Sd D;
                D.from_dhat(hi); // used only if pi-bond used

                double ir2 = 1/(hi.norm2() + 1e-16);
                double ir  = sqrt(ir2);
                hi.mul(ir);

                for(int j=conf.b; j<N_BOND_MAX; j++){  // pi-electron
                    Vec3d fia,fja;
                    // Note: moving this here save D-mat calculation
                    // ToDo : proper dhat
                    const Vec3d& hj = hs[j];
                    double c   = hi.dot(hj);
                    double dfc = Kp*2*c;

                    //printf( " %i(%i,%i) c %g (%g,%g,%g) (%g,%g,%g) \n", ia, i,j, c, hi.x,hi.y,hi.z,  hj.x,hj.y,hj.z );
                    D.dhat_dot( hj, fia );
                    fia.mul(dfc*ir2*ir);

                    //fia.set_mul(hj,dfc);
                    //fja.set_mul(hi,dfc);

                    //glColor3f(1.0,0.0,0.0); Draw3D::drawVecInPos( fja*0.5,apos[ia]+hj*0.5 );

                    //fi.add(fia); fs[j].add(fja);
                    //fa.sub(fia); fa   .sub(fja);
                }
            }
*/
            if(i>=conf.c){ glColor3f(1.0,0.0,0.0); Draw3D::drawVecInPos( fi, pi ); }
        }else{           // pi
/*
            double c,dfc;
            Vec3d fia,fja;
            const Vec3d& hi = hs[i];

            for(int j=i+1; j<N_BOND_MAX; j++){ // pi-pi
                const Vec3d& hj = hs[j];
                c   = hi.dot(hj);
                dfc = Kp*2*c;
                fia.set_mul(hj,dfc);
                fja.set_mul(hi,dfc);

                //fi.add(fia); fs[j].add(fja);
                //fa.sub(fia); fa   .sub(fja);
            }
*/
        }

    }

    aforce[ia].add(fa);

    glColor3f(1.0,0.0,0.0); Draw3D::drawVecInPos( aforce[ia], pa );

    //aforce[ia].set(0.); // DEBUG

}


inline void evalBonds       (){ for(int i=0; i<nbond; i++){ evalBond(i);        } }
inline void evalBondPiForces(){ for(int i=0; i<nbond; i++){ evalBondPiForce(i); } }
inline void evalAtoms       (){ for(int i=0; i<natom; i++){ evalAtom(i);        } //exit(0);
}


void evalLJQs(){
    for(int i=0; i<natom; i++){
        const Vec3d& REQi = aREQs[i];
        const Vec3d& pi   = apos[i];
        Vec3d fi; fi.set(0.0);
        /*
        for(int j=0; j<natoms; j++){
            if(i!=j){
                Vec3d& ljq_j = aREQ[j];
                double rij = ljq_i.x+ljq_j.x;
                double eij = ljq_i.y*ljq_j.y;
                double qq  = ljq_i.z*ljq_j.z;
                addAtomicForceLJQ( pi-pos[j], f, rij, -eij, qq );
            }
        }
        force[i].add(f);
        */
        for(int j=0; j<i; j++){
            Vec3d  fij;fij.set(0.);
            Vec3d& REQj = aREQs[j];
            addAtomicForceLJQ( pi-apos[j], fij, REQi.x+REQj.x, -REQi.y*REQj.y, REQi.z*REQj.z );
            aforce[j].sub(fij);
            aforce[j].add(fij);
            //fi      .add(fij);
        }
        //force[i].add(fi);
    }
}

void evalLJQs(int n, const Vec3d* REQs, const Vec3d* ps, Vec3d* fs, const Vec3d& REQi, const Vec3d& pi, Vec3d& fi){
    for(int j=0; j<n; j++){    // atom-atom
        Vec3d  fij; fij.set(0.);
        const Vec3d& REQj = REQs[j];
        addAtomicForceLJQ( pi-ps[j], fij, REQi.x+REQj.x, -REQi.y*REQj.y, REQi.z*REQj.z );
        aforce[j].sub(fij);
        fi       .add(fij);
        //fi      .add(fij);
    }
}

void evalLJQs_H(){
    for(int i=0; i<natom; i++){ evalLJQs( natom,   aREQs,   apos,   aforce,  aREQs  [i],   apos[i],   aforce[i] ); }
    for(int i=0; i<natom; i++){ evalLJQs( ncap,  capREQs, capPos, capForce,  aREQs  [i], capPos[i], capForce[i] ); }
    for(int i=0; i<ncap;  i++){ evalLJQs( ncap,  capREQs, capPos, capForce,  capREQs[i], capPos[i], capForce[i] ); }
}

void evalForces(){
    cleanForce();
    projectBondCenter();
    //projectCaps();
    //evalLJQs();
    evalAtoms();
    evalBonds();
    //evalBondPiForces();
}

void moveGD(double dt){
    for(int i=0; i<nDOF; i++){ dofs[i].add_mul( fdofs[i], dt ); }
}

void setBondsAndHydrogens( Vec2i* bond2atom, Vec2i* Hps ){
    for(int ia=0; ia<natom; ia++){
        aconf[ia].a=0;
    }
    for(int ib=0; ib<nbond; ib++){
        const Vec2i& ba = bond2atom[ib];
        Vec3i& ci    = aconf[ba.i];
        Vec3i& cj    = aconf[ba.j];
        bondIs[ib].a = ba.i*N_BOND_MAX + ci.a; ci.a++;
        bondIs[ib].b = ba.j*N_BOND_MAX + cj.a; cj.a++;
        printf( "bond %i a2b(%i,%i) -> bIs(%i,%i) \n", ib, ba.i, ba.j, bondIs[ib].a, bondIs[ib].b );
    }
    for(int ia=0; ia<natom; ia++){
        const Vec2i& hp = Hps[ia];
        Vec3i& c = aconf[ia];
        // like this    [ sigmaBonds | hydrogens | epairs | pi ]
        c.c =c.a;      // end of bonds ( start of hydrogens )
        c.a+=hp.a; // end of sigma bonds (to atoms & hydrogens), not epairs
        c.b =4-hp.b;   // end of sigma&epairs is what lefts after assign pi bonds
    }
}



void guessOrbs(){
    // TODO: bonds should be before hydrogens, otherwise it is huge problem
    //cleanForce();
    //for(int i=0; i<nbond; i++){ evalBond(i);         }
    projectBondCenter();
    for(int ia=0; ia<natom; ia++){
        Vec3d  pa = apos[ia];
        Vec3i& c  = aconf[ia];
        Vec3d* hs = hdir + ia*N_BOND_MAX;
        int nb = c.c; // number of defined bonds
        Mat3d m;
        if      (nb==3){ // defined by 3 sigma bonds
            m.b.set_cross( hs[1]-hs[0], hs[2]-hs[0] );
            m.b.mul( -1/m.b.norm() );
            if(c.b==4){ // sp3 no-pi
                if( 0 < m.b.dot( hs[0]+hs[1]+hs[2]+pa*3 ) ){ m.b.mul(-1.); }
                hs[3]=pa+m.b;
            }else{
                hs[3]=m.b;
            }

        }else if(nb==2){ // defined by 2 sigma bonds
            m.fromCrossSafe( hs[0]-pa, hs[1]-pa );
            if      (c.b==4){ // -CH2- like sp3 no-pi
                const double cb = 0.81649658092; // sqrt(2/3)
                const double cc = 0.57735026919;  // sqrt(1/3)
                hs[nb  ] = pa+m.c*cc+m.b*cb;
                hs[nb+1] = pa+m.c*cc-m.b*cb;
            }else if(c.b==3){ // =CH- like  sp 1-pi
                hs[nb  ] = pa+m.c;
                hs[nb+1] = m.b;
                printf("like =CH- H(%g,%g,%g) pi(%g,%g,%g,) \n", hs[nb].x, hs[nb].y, hs[nb].z, hs[nb+1].x, hs[nb+1].y, hs[nb+1].z );
            }else{            // #C- like sp 2-pi
                hs[nb  ] = m.c;
                hs[nb+1] = m.b;
            }
        }else if(nb==1){
            m.c = hs[0]-pa; m.c.normalize();
            m.c.getSomeOrtho(m.b,m.a);
            if      (c.b==4){ // -CH3 like sp3 no-pi
                const double ca = 0.81649658092;  // sqrt(2/3)
                const double cb = 0.47140452079;  // sqrt(2/9)
                const double cc =-0.33333333333;  // 1/3
                hs[nb  ] = pa+m.c*cc + m.b*(cb*2) ;
                hs[nb+1] = pa+m.c*cc - m.b* cb    + m.a*ca;
                hs[nb+2] = pa+m.c*cc - m.b* cb    - m.a*ca;
            }else if(c.b==3){ // =CH2 like sp2 1-pi
                const double ca = 0.87758256189;  // 1/2
                const double cc =-0.5;            // sqrt(1/8)
                hs[nb  ] = pa+m.c*cc + m.a*ca;
                hs[nb+1] = pa+m.c*cc - m.a*ca;
                hs[nb+2] = m.b;
            }else{            // #CH sp  2-pi
                hs[nb  ] = pa+m.c*-1;
                hs[nb+1] = m.b;
                hs[nb+2] = m.a;
            }
        }else{
            printf( " WARRNING: atom %i not bonded to anything\n", ia );
        }
    }
}


///============= Backup
///============= Backup : Dirctions instead of positions
///============= Backup : Dirctions instead of positions
///============= Backup

/*
void bondRecoil(){
    for(int ib=0; ib<nbonds; ib++){
        const Vec2i& b   = bondIs[ibond];
        int ia = b.i>>2;
        int ja = b.j>>2;
        Vec3d f = hforce[b.i] + hforce[b.i];
        aforce[ia].add(f);
        aforce[ja].sub(f);
    }
}
*/

/*
void projectCaps(){
    int ic=0;
    for(int ia=0; ia<natom; ia++){
        const Vec3d& pi = apos[ia];
        //const Vec3d* hs = hdir[ia<<2];
        int ih = ia*N_BOND_MAX;
        for(int j=aconf[ia].c; j<aconf[ia].a; j++){ // all sigma bonds
            //capPos[ic] = pi + hs[j] * R;
            capPos[ic] = pi + hdir[ih+j];
            ic++;
        }
    }
};
*/

/*
void evalAtom(int ia){

// THIS DOES NOT WORK - cos = <hi|hj> does not work in 3D (e.g. does not form sp3 tetrahedron)

    const Vec3i& conf = aconf[ia];
    int ih      = ia<<2;
    Vec3d*  hs  = hdir  +ih;
    Vec3d*  fs  = hforce+ih;
    double* irs = invRs +ih;

    //double ivr3s[N_BOND_MAX];
    Mat3Sd Ds   [N_BOND_MAX];  // derivative determinants

    for(int i=0; i<conf.a; i++){ // all sigma bonds
        Ds[i].from_dhat(hs[i]);
    }

    // normalize  e-pairs and pi-bonds
    // ToDo : this can be done more efficiently later
    for(int i=conf.a; i<N_BOND_MAX; i++){
        hs[i].normalize();
    }

    aforce[ia].set(0.); // DEBUG
    //return;

    Vec3d fa = Vec3dZero;
    for(int i=0; i<N_BOND_MAX; i++){
        const Vec3d& hi = hs[i];
        Vec3d&       fi = fs[i];
        for(int j=i+1; j<N_BOND_MAX; j++){
            const Vec3d& hj = hs[j];
            Vec3d&       fj = fs[j];
            double c   = hi.dot(hj);
            Vec3d fia,fja;

            fia.set(0.); fja.set(0.); // DEBUG

            double dfc;
            if      (j<conf.a){ // bond-bond

                dfc = Kss/sqrt((1-c)*0.5);
                Ds[i].dhat_dot( hj, fia );  fia.mul(dfc*irs[i]);
                Ds[j].dhat_dot( hi, fja );  fja.mul(dfc*irs[j]);

            }else if(j<conf.b){ // epairs
                if(i<conf.a){   // bond-epair

                    //dfc = -Kse/sqrt((1-c)*0.5);
                    dfc = Kse;
                    //Ds[i].dhat_dot( hj, fia );  fia.mul(dfc*irs[i]);
                    fia.set_mul(hj,dfc);
                    fja.set_mul(hi,dfc);


                    //glColor3f(1.0,0.0,0.0); Draw3D::drawVecInPos( fia*0.5,apos[ia]+hi*0.5 );
                    glColor3f(1.0,0.0,0.0); Draw3D::drawVecInPos( fja*0.5,apos[ia]+hj*0.5 );

                }else{           // epair-epair
                    dfc = Kee/sqrt((1-c)*0.5);
                    //dfc = -Kee;
                    //printf( "%i (%i,%i) %g\n", ia, i,j, c );
                    fia.set_mul(hj,dfc);
                    fja.set_mul(hi,dfc);

                    glColor3f(0.0,0.0,1.0); Draw3D::drawVecInPos( fia*0.5,apos[ia]+hi*0.5 );
                    glColor3f(0.0,0.0,1.0); Draw3D::drawVecInPos( fja*0.5,apos[ia]+hj*0.5 );

                }
            }else{  // pi   - make orthogonal   =>   E = Kp * (<hi,hj>)^2

                if     (i<conf.a){   // pi-sigma
                    //dfc = Kp/sqrt((1+c)*0.5);
                    dfc = Kp*2*c;
                    Ds[i].dhat_dot( hj, fia );  fia.mul(dfc*irs[i]);
                    fja.set_mul(hi,dfc);
                //}else if(i<iconf.b){  // pi-epair
                //   Kpe
                }else{                // pi-pi   ( and e-pair)
                    dfc = Kp*2*c;
                    fia.set_mul(hj,dfc);
                    fja.set_mul(hi,dfc);
                }

            }

            //glColor3f(1.0,0.0,1.0); Draw3D::drawVecInPos( fj*0.5,apos[ia]+hj*0.5 );

            fi.add(fia); fj.add(fja);
            fa.sub(fia); fa.sub(fja);
        }

        glColor3f(1.0,0.0,1.0); Draw3D::drawVecInPos( fi*0.5,apos[ia]+hi*0.5 );
    }

    //aforce[ia].add(fa);

    aforce[ia].set(0.); // DEBUG

}
*/

/*
void guessOrbs(){
    // TODO: bonds should be before hydrogens, otherwise it is huge problem
    //cleanForce();
    for(int i=0; i<nbond; i++){ evalBond(i);         }
    for(int ia=0; ia<natom; ia++){
        Vec3i& c  = aconf[ia];
        Vec3d* hs = hdir + ia*N_BOND_MAX;
        int nb = c.c; // number of defined bonds
        Mat3d m;
        if      (nb==3){ // defined by 3 sigma bonds
            m.b.set_cross( hs[1]-hs[0], hs[2]-hs[0] );
            m.b.mul( -1/m.b.norm() );
            if(c.b==4){ // sp3 no-pi
                if( 0 < m.b.dot( hs[0]+hs[1]+hs[2]) ){ m.b.mul(-1.); }
            }
            hs[3]=m.b;
        }else if(nb==2){ // defined by 2 sigma bonds
            m.fromCrossSafe( hs[0], hs[1] );
            if      (c.b==4){ // -CH2- like sp3 no-pi
                const double cb = 0.81649658092; // sqrt(2/3)
                const double cc = 0.57735026919;  // sqrt(1/3)
                hs[nb  ] = m.c*cc+m.b*cb;
                hs[nb+1] = m.c*cc-m.b*cb;
            }else if(c.b==3){ // =CH- like  sp 1-pi
                hs[nb  ] = m.c;
                hs[nb+1] = m.b;
                printf("like =CH- H(%g,%g,%g) pi(%g,%g,%g,) \n", hs[nb].x, hs[nb].y, hs[nb].z, hs[nb+1].x, hs[nb+1].y, hs[nb+1].z );
            }else{            // #C- like sp 2-pi
                hs[nb  ] = m.c;
                hs[nb+1] = m.b;
            }
        }else if(nb==1){
            m.c = hs[0]; m.c.normalize();
            m.c.getSomeOrtho(m.b,m.a);
            if      (c.b==4){ // -CH3 like sp3 no-pi
                const double ca = 0.81649658092;  // sqrt(2/3)
                const double cb = 0.47140452079;  // sqrt(2/9)
                const double cc =-0.33333333333;  // 1/3
                hs[nb  ] = m.c*cc + m.b*(cb*2) ;
                hs[nb+1] = m.c*cc - m.b* cb    + m.a*ca;
                hs[nb+2] = m.c*cc - m.b* cb    - m.a*ca;
            }else if(c.b==3){ // =CH2 like sp2 1-pi
                const double ca = 0.87758256189;  // 1/2
                const double cc =-0.5;            // sqrt(1/8)
                hs[nb  ] = m.c*cc + m.a*ca;
                hs[nb+1] = m.c*cc - m.a*ca;
                hs[nb+2] = m.b;
            }else{            // #CH sp  2-pi
                hs[nb  ] = m.c*-1;
                hs[nb+1] = m.b;
                hs[nb+2] = m.a;
            }
        }else{
            printf( " WARRNING: atom %i not bonded to anything\n", ia );
        }
    }
}
*/

}; // FFsp

#endif
