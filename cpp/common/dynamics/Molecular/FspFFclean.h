
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

inline void combineREQ(const Vec3d& a, const Vec3d& b, Vec3d& out){
    out.a=a.a+b.a; // radius
    out.b=a.b*b.b; // epsilon
    out.c=a.c*b.c; // q*q
}

inline void addAtomicForceLJQ( const Vec3d& dp, Vec3d& f, const Vec3d& REQ ){
    //Vec3f dp; dp.set_sub( p2, p1 );
    const double COULOMB_CONST = 14.3996448915d;  //  [V*A/e] = [ (eV/A) * A^2 /e^2]
    double ir2  = 1/( dp.norm2() + 1e-4 );
    double ir   = sqrt(ir2);
    double ir2_ = ir2*REQ.a*REQ.a;
    double ir6  = ir2_*ir2_*ir2_;
    double fr   = ( ( 1 - ir6 )*ir6*12*REQ.b + ir*REQ.c*-COULOMB_CONST )*ir2;
    //printf( " (%g,%g,%g) r %g fr %g \n", dp.x,dp.y,dp.z, 1/ir, fr );
    f.add_mul( dp, fr );
}

class NBFF{ public:
// non bonded forcefield
    int n       = 0;
    int nmask   = 0;
    Vec3d* REQs = 0;
    Vec3d* ps   = 0;
    Vec3d* fs   = 0;
    Vec2i* pairMask = 0; // should be ordered ?

void realloc(int n_, int nmask_){
    n=n_;
    nmask=nmask_;
    _realloc(REQs,n);
    _realloc(ps  ,n);
    _realloc(fs  ,n);
    _realloc(pairMask ,nmask);
}

void setREQs(int i0,int i1, const Vec3d& REQ){
    for(int i=i0;i<i1;i++){ REQs[i]=REQ; }
}

void cleanForce(){
    for(int i=0; i<n; i++){ fs[i].set(0.); }
}

void evalLJQs(){
    const int N=n;
    for(int i=0; i<N; i++){
        Vec3d fi = Vec3dZero;
        Vec3d pi = ps[i];
        const Vec3d& REQi = REQs[i];
        for(int j=i+1; j<N; j++){    // atom-atom
            Vec3d fij = Vec3dZero;
            Vec3d REQij; combineREQ( REQs[j], REQi, REQij );
            addAtomicForceLJQ( ps[j]-pi, fij, REQij );
            fs[j].sub(fij);
            fi   .add(fij);
        }
        fs[i].add(fi);
    }
}

void evalLJQ_sortedMask(){
    //Vec2i* m=pairMask;
    int im=0;
    const int N=n;
    //printf( "N %i \n" );
    for(int i=0; i<N; i++){
        Vec3d fi = Vec3dZero;
        Vec3d pi = ps[i];
        const Vec3d& REQi = REQs[i];
        //printf( "-- ia %i REQi (%g,%g,%g) \n", i, REQi.x, REQi.y, REQi.z );
        for(int j=i+1; j<N; j++){    // atom-atom
            //printf( "(%i,%i) m[%i](%i,%i) ", i,j, im,pairMask[im].a, pairMask[im].b );
            if( (im<nmask)&&(i==pairMask[im].i)&&(j==pairMask[im].j) ){
                //printf( " masked \n" );
                im++; continue;
            }
            Vec3d fij = Vec3dZero;
            Vec3d REQij; combineREQ( REQs[j], REQi, REQij );
            //printf( "REQij (%g,%g,%g) \n", REQij.x, REQij.y, REQij.z );
            addAtomicForceLJQ( ps[j]-pi, fij, REQij );
            //printf( " fij(%g,%g,%g) \n", fij.x, fij.y, fij.z );
            //{   // DEBUG
            //    double c = -(ps[j]-pi).dot(fij);
            //    printf( " f[%i,%i] |f| %g |d| %g  \n", i,j, fij.norm()*signum(c), (ps[j]-pi).norm() );
            //    if(c>0)Draw3D::drawLine(pi,ps[j]);
            //}
            fs[j].sub(fij);
            fi   .add(fij);
        }
        fs[i].add(fi);
    }
}

};

class FspFF{ public:

    bool substract_LJq = true;

    //double fHlen =  1.0;
    double fHlen =  1.6;  // C-H bond lenght = fHlen*l0ae    e.g.: 1.1=(1.6*0.7)
    double l0ae  =  0.7;  // atom-electron equilibrium length
    double Kae   =  20.0;  // atom-electron strenght

    double Wee  =  1.0;  // electon-electron width
    double Kee  = -10.0;  // electon-electron strenght

    double Kpe  = -15.0;  // onsite  p-s orthogonalization strenght
    double Kpp  = -20.0;  // onsite  p-p orthogonalization strenght
    double Kpi  = -2.0;  // offsite p-p alignment         strenght

    int natom   = 0;
    int nbond   = 0; // number of all bonding orbitals
    int norb    = 0;
    int nDOF    = 0; // natoms + neps + npis

    int ncap   = 0; // number of capping
    int nepair = 0; // number of electron pairs
    int nporb  = 0; // number of electron pi orbitals

    Vec2i  * bondIs = 0;
    //Vec2d  * bLKs   = 0;  // lenght[A],K [eV/A]

    Vec3d  * dofs   = 0;  // degrees of freedom
    //Vec3d  * vdofs = 0;
    Vec3d  * fdofs  = 0; // degrees of freedom

    Vec3d  * apos   = 0;      // atomic position // ALIAS
    Vec3d  * aforce = 0;      // atomic forces   // ALIAS

    //Vec3d  * capREQs = 0;
    //Vec3d  * aREQs   = 0;
    //Vec3i  * aconf   = 0;     // nH, nsigna,nt=(nsigna+ne)     npi = (4 - vt)=(4 - (nsigma+ne))
    Vec3ui8  * aconf   = 0;

    //Quat4i * atom2bond = 0; // atom to bond map

    //double * blen    = 0;   // bond lengths
    Vec3d  * hdir    = 0;     // normalized bond unitary vectors  // ALIAS
    Vec3d  * hforce  = 0;
    //double * invRs   = 0;
    //Vec3d  * bforce  = 0;     // forces acting on it              // ALIAS
    //double * bls     = 0;

    //Vec3d * capPos   = 0;
    //Vec3d * capForce = 0;

    void realloc(int natom_, int nbond_){
        natom=natom_;
        nbond=nbond_;
        norb =natom*N_BOND_MAX;
        nDOF =natom + norb;

        _realloc(aconf  ,natom);

        _realloc(bondIs ,nbond);

        _realloc(dofs  ,nDOF);
        _realloc(fdofs ,nDOF);

        apos   = dofs;
        aforce = fdofs;

        hdir   =  dofs+natom;
        hforce = fdofs+natom;
    }

// ======== Force Evaluation

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

void transferBondForce(){
    for(int ib=0; ib<nbond; ib++){
        const Vec2i& b   = bondIs[ib];
        int ia = b.i>>2;
        int ja = b.j>>2;
        //aforce[ia].add(hforce[b.j]);
        //aforce[ja].add(hforce[b.i]);
        Vec3d f=(hforce[b.j]+hforce[b.i])*0.5;
        aforce[ia].add(f);
        aforce[ja].add(f);
        hforce[b.i].set(0.); hforce[b.j].set(0.); // DEBUG: does not realy matter
    }
}

inline void evalBondPiForce(int ibond){
    // NOTE: this cannot be in bondForce, since some hdirs may not be updated yet
    const Vec2i& b = bondIs[ibond];
    int ia = b.i>>2;
    int ja = b.j>>2;
    if((aconf[ia].b==3)&&(aconf[ja].b==3)){
        int i = (ia<<2)+3;
        int j = (ja<<2)+3;
        const Vec3d&  hi  = hdir[i];
        const Vec3d&  hj  = hdir[j];
        double c   = hi.dot(hj);
        double dfc = Kpi*-2*c;
        //printf( "pi-bond %i(%i,%i) %g \n", ibond, ia, ja,  c );
        Vec3d  fia; fia.set_mul(hj,dfc);  fia.add_mul(hi,-c*dfc);
        Vec3d  fja; fja.set_mul(hi,dfc);  fja.add_mul(hj,-c*dfc);
        hforce[i].add(fia); // aforce[ia].sub(fia);
        hforce[j].add(fja); // aforce[ja].sub(fja);
    }
}


void evalAtom(int ia){
    const Vec3ui8& conf = aconf[ia];
    const Vec3d& pa     = apos [ia];
    int ih      = ia<<2;
    Vec3d*  hs  = hdir  +ih;
    Vec3d*  fs  = hforce+ih;
    //double* irs = invRs +ih;

    Vec3d fa    = Vec3dZero;
    double W2ee = Wee*Wee;

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
                double fr   = Kee*ir2*ir2;
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
                    double dfc = Kpe*2*c;

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
                fr = Kae*(l-l0ae);
                // outproject radial component
                Vec3d frad;
                frad.set_mul(d,d.dot(fi));
                fa.add(frad);
                fi.sub(frad);
            }else{
                l  = d.norm();
                fr = Kae*(l-l0ae)/l;
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
                dfc = Kpp*2*c;
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

inline void evalBondPiForces(){ for(int i=0; i<nbond; i++){ evalBondPiForce(i); } }
inline void evalAtoms       (){
    //printf( ">>>BEGIN evalAtoms \n" );
    for(int i=0; i<natom; i++){ evalAtom(i);        }
    //printf( "<<<END   evalAtoms \n" );
}

void evalForces(){
    cleanForce();
    projectBondCenter();
    evalAtoms();
    transferBondForce();
    evalBondPiForces();
}

void moveGD(double dt){
    for(int i=0; i<nDOF; i++){ dofs[i].add_mul( fdofs[i], dt ); }
}


void checkForceTorque(Vec3d& cog,Vec3d& fsum,Vec3d& tqsum){
    cog  =Vec3dZero;
    fsum =Vec3dZero;
    tqsum=Vec3dZero;
    int n=0;
    for(int ia=0; ia<natom; ia++){
        cog .add(apos  [ia]);
        fsum.add(aforce[ia]);
        //printf( "fs[%i, ] %g   | (%g,%g,%g)  \n", ia,  aforce[ia].norm(), aforce[ia].x,aforce[ia].y,aforce[ia].z );
        n++;
        const Vec3ui8& c = aconf[ia];
        Vec3d* hs = hdir   + ia*N_BOND_MAX;
        Vec3d* fs = hforce + ia*N_BOND_MAX;
        for(int ih=c.c; ih<c.b; ih++){
            //printf( "fs[%i,%i] %g   | (%g,%g,%g)  \n", ia, ih, fs[ih].norm(), fs[ih].x,fs[ih].y,fs[ih].z );
            cog .add(hs[ih]);
            fsum.add(fs[ih]);
            n++;
        }
    }
    cog.mul(1./n);
    //Vec3d fav = fsum*(1./n);
    for(int ia=0; ia<natom; ia++){
        tqsum.add_cross(apos[ia]-cog,aforce[ia]);
        const Vec3ui8& c = aconf[ia];
        Vec3d* hs = hdir   + ia*N_BOND_MAX;
        Vec3d* fs = hforce + ia*N_BOND_MAX;
        for(int ih=c.c; ih<c.b; ih++){
            tqsum.add_cross(hs[ih]-cog,fs[ih]);
        }
    }
}

// ======== Paricle I/O

int outputParticlePositions( Vec3d* ps_ ){
    Vec3d* ps=ps_;
    for(int i=0; i<natom; i++){
        *ps=apos[i];
        ps++;
    }
    double mHlen = 1-fHlen;
    for(int ia=0; ia<natom; ia++){
        const Vec3ui8& c = aconf[ia];
        Vec3d* hs = hdir + ia*N_BOND_MAX;
        for(int ih=c.c; ih<c.a; ih++){
            //*ps=hs[ih];
            *ps=( hs[ih]*fHlen + apos[ia]*mHlen );
            ps++;
        }
    }
    return ps-ps_;
}

int inputParticleForces( Vec3d* fs_ ){
    Vec3d* fs=fs_;
    for(int i=0; i<natom; i++){
        aforce[i].add(*fs);
        fs++;
    }
    double mHlen = 1-fHlen;
    for(int ia=0; ia<natom; ia++){
        const Vec3ui8& c = aconf[ia];
        Vec3d* hfs = hforce + ia*N_BOND_MAX;
        for(int ih=c.c; ih<c.a; ih++){
            //hfs[ih].add(*fs);
            hfs   [ih].add_mul( *fs, fHlen );
            aforce[ia].add_mul( *fs, mHlen );
            fs++;
        }
    }
    return fs-fs_;
}

int outputParticleBonds( Vec2i* b2a_ ){
    Vec2i* b2a=b2a_;
    auto selectMin = [](int* arr, int n, int ilow){
        int amin=2147483647; // INT_MAX
        for(int i=0;i<n;i++){
            int ai = arr[i];
            if( (ai>ilow)&&(ai<amin)){amin=ai;}
        }
        return amin;
    };
    // bonds to neighbors
    int * nng = new int[natom];
    int * ngs = new int[natom*N_BOND_MAX];
    for(int ia=0; ia<natom; ia++){ nng[ia]=0; };
    for(int ib=0; ib<nbond; ib++){
        const Vec2i& b = bondIs[ib];
        int ia=b.i>>2;
        int ja=b.j>>2;
        //b.order();
        ngs[ia*N_BOND_MAX + nng[ia]]=ja; nng[ia]++;
        ngs[ja*N_BOND_MAX + nng[ja]]=ia; nng[ja]++;
        //printf( "ib %i (%i,%i) atoms(%i,%i) ng(%i,%i) \n", ib, b.i, b.j, ia,ja, nng[ia], nng[ja] );
    }
    //printf("----\n");
    // go over atom nieghbors
    int icap=0;
    for(int ia=0; ia<natom; ia++){
        // atom -> atom sorted
        //printf( "-- ia %i ng %i \n", ia, nng[ia] );
        int* ng = ngs + ia*N_BOND_MAX;
        int nh=nng[ia];
        int imin=-1;
        for(int ih=0;ih<nh;ih++){
            imin=selectMin(ng,nh,imin); // select sort - inefficient but the array is short
            //printf( "ng[%i] %i \n", ih, ng[ih] );
            //printf( "ng[%i] %i | imin %i | (%i,%i) \n", ih, ng[ih], imin, ia,imin );
            if(imin>ia){
                //printf("insert\n");
                *b2a={ia,imin};
                b2a++;
            }
        }

        // atom -> Hydrogen sorted
        //printf( "H-range (%i,%i) \n", aconf[ia].c, aconf[ia].a );
        const Vec3ui8& c = aconf[ia];
        for(int ih=c.c; ih<c.a; ih++){
            //*b2a={ia,natom+ia*N_BOND_MAX+ih};
            //printf( "ia %i ih %i icap %i (%i,%i) \n", ia, ih, icap,   ia,natom+icap );
            *b2a={ia,natom+icap};
            b2a++;
            icap++;
        }

    }
    delete [] nng;
    delete [] ngs;
    return b2a-b2a_;
}

// ======== Build Molecule

void setBondsAndHydrogens( Vec2i* bond2atom, Vec2i* Hps ){
    ncap  =0;
    nepair=0;
    nporb =0;
    for(int ia=0; ia<natom; ia++){
        aconf[ia].a=0;
    }
    for(int ib=0; ib<nbond; ib++){
        const Vec2i& ba = bond2atom[ib];
        Vec3ui8& ci    = aconf[ba.i];
        Vec3ui8& cj    = aconf[ba.j];
        bondIs[ib].a = ba.i*N_BOND_MAX + ci.a; ci.a++;
        bondIs[ib].b = ba.j*N_BOND_MAX + cj.a; cj.a++;
        //printf( "bond %i a2b(%i,%i) -> bIs(%i,%i) \n", ib, ba.i, ba.j, bondIs[ib].a, bondIs[ib].b );
    }
    for(int ia=0; ia<natom; ia++){
        const Vec2i& hp = Hps[ia];
        Vec3ui8& c = aconf[ia];
        // like this    [ sigmaBonds | hydrogens | epairs | pi ]
        c.c =c.a;      // end of bonds ( start of hydrogens )
        c.a+=hp.a; // end of sigma bonds (to atoms & hydrogens), not epairs
        c.b =4-hp.b;   // end of sigma&epairs is what lefts after assign pi bonds
        ncap  +=hp.a;
        nporb +=hp.b;
        nepair+=c.b-c.a;
    }
}

void guessOrbs(){
    // TODO: bonds should be before hydrogens, otherwise it is huge problem
    //cleanForce();
    //for(int i=0; i<nbond; i++){ evalBond(i);         }
    projectBondCenter();
    for(int ia=0; ia<natom; ia++){
        Vec3d  pa  = apos[ia];
        Vec3ui8& c = aconf[ia];
        Vec3d*  hs = hdir + ia*N_BOND_MAX;
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

}; // FFsp

#endif
