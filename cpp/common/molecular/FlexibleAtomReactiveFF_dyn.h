
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

#ifndef FlexibleAtomReactiveFF_dyn_h
#define FlexibleAtomReactiveFF_dyn_h

//#include <vector>

#include <vector>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
//#include "quaternion.h"
#include "Forces.h"

#include "FlexibleAtomReactiveFF.h"

struct FARFFAtom{
    FlexibleAtomType* type=0;
    Vec3d pos;
    Vec3d force;
    double energy;
    Vec3ui8 conf;
    Vec3d   opos   [N_BOND_MAX];     // normalized bond unitary vectors  // ALIAS
    Vec3d   oforce [N_BOND_MAX];
    double  oenergy[N_BOND_MAX];

    void cleanForce(){
        force=Vec3dZero;
        for(int i=0; i<N_BOND_MAX; i++){ oforce[i]=Vec3dZero; }
    };

    void normalizeOrbs(){
        for(int i=0; i<N_BOND_MAX; i++){ opos[i].normalize_taylor3(); }
    }

    void orbRecoil(){
        for(int i=0; i<N_BOND_MAX; i++){ force.sub( oforce[i] ); }
    }

    void removeNormalForces(){
        for(int i=0; i<N_BOND_MAX; i++){ oforce[i].makeOrthoU(opos[i]); }
    }

    void moveGD(double dt){
        pos.add_mul( force, dt );
        for(int i=0; i<N_BOND_MAX; i++){  opos[i].add_mul( oforce[i], dt ); }
    }

    double evalOnSite(){
        double E    = 0;
        const double w2ee   = sq(type->Wee);
        int    nsigma = conf.a;
        Vec2d cs0 = sp_cs0s[nsigma];
        //int j_DEBUG=0;
        // -- repulsion between sigma bonds
        for(int i=1; i<nsigma; i++){
            const Vec3d& hi = opos[i];
            Vec3d&       fi = oforce[i];
            for(int j=0; j<i; j++){ // electron-electron
                E += evalCos2  (hi,opos[j],fi,oforce[j],type->Kee,cs0.x);
            }
        }
        // -- orthogonalization with p-orbitals
        for(int i=nsigma; i<N_BOND_MAX; i++){
            const Vec3d& hi = opos  [i];
            Vec3d&       fi = oforce[i];
            for(int j=0; j<i; j++){
                E += evalCos2   ( hi, opos[j], fi, oforce[j], type->Kpp, 0);
            }
        }
        return E;
    }

    double evalInteraction( FARFFAtom& B, FlexiblePairType& type){
    //double evalPair( int ia, int ja, int nbi, int nbj ){
        //const Vec3ui8& confi = aconf[ia];
        //const Vec3ui8& confj = aconf[ja];
        int nbi =   conf.a;
        int nbj = B.conf.a;
        Vec3d  hij; hij.set_sub( B.pos, pos );   // = apos[ja] - apos[ia];
        double r2   = hij.norm2() + R2SAFE;
        double r    = sqrt( r2 );
        double invr = 1/r;
        hij.mul(invr); // = dij*(1/rij);
        double r4     = r2*r2;
        double r6     = r4*r2;
        double invVdW = 1/( r6 + type.vdW_w6 );
        double EvdW   = type.vdW_c6*invVdW;
        double fvdW   = -6*EvdW*invVdW*r4*r;
        double expar =    exp( type.bMorse*( r - type.rbond0 ) );
        double Esr   =    type.aMorse*expar*expar;
        double Eb    = -2*type.aMorse*expar;
        double fsr   =  2*type.bMorse*Esr;
        double frb   =    type.bMorse*Eb;
        double E = Esr + EvdW;
        Vec3d force_; force_.set_mul( hij, fsr + fvdW );
        Mat3Sd Jbond;
        Jbond.from_dhat(hij);
        for(int ib=0; ib<nbi; ib++){
            const Vec3d& hi = opos  [ib];
                  Vec3d& fi = oforce[ib];
            double       ci = hij.dot( hi );   // ci = <hi|hij>
            if(ci<=0) continue;
            for(int jb=0; jb<nbj; jb++){
                const Vec3d& hj = B.opos  [jb];
                      Vec3d& fj = B.oforce[jb];
                double cj       = hij.dot( hj );  // cj  = <hj|hij>
                if(cj>=0) continue; // avoid symmetric image of bond
                double cc  = ci*cj;
                double cc2 = cc*cc;
                double e   = cc2*cc2;
                double de  = 4*cc2*cc;
                double eEb = e*Eb*0.5d;
                  oenergy[ib]+=eEb;
                B.oenergy[jb]+=eEb;
                E += eEb+eEb;
                double deEb     =    Eb*de;
                double deEbinvr =  deEb*invr;
                Jbond.mad_ddot( hi, force_, deEbinvr*cj );
                Jbond.mad_ddot( hj, force_, deEbinvr*ci );
                force_.add_mul(hij, e*frb);
                fi.add_mul( hij, -cj*deEb );
                fj.add_mul( hij, -ci*deEb );
            }
        }
        // ---------- Pi-Pi interaction
        if( (nbi==3)&&(nbj==3) ){ // align pz-pz in sp2-sp2 pair
            const Vec3d& hi =   opos  [3];
            const Vec3d& hj = B.opos  [3];
            Vec3d&       fi =   oforce[3];
            Vec3d&       fj = B.oforce[3];
            double c = hi.dot(hj);
            double de = -2*c*type.Kpi*Eb;
            fi.add_mul(hj,de);
            fj.add_mul(hi,de);
            //force      // TODO: proper derivatives of energy
        }
        //aforce[ia].add(force);
        //aforce[ja].sub(force);
          force.add(force_);
        B.force.sub(force_);
        return E;
    }

};


class FARFF_dyn{ public:

    bool   substract_LJq = true;
    double Eatoms=0,Epairs=0;

    FlexibleAtomType atype0;
    FlexiblePairType ptype0;

    std::vector<FARFFAtom> atoms;

    //void realloc(int natom_){}

// ======== Force Evaluation

void cleanForce(){
    for(FARFFAtom& a:atoms){ a.cleanForce(); }
}

double evalPairs(){
    Epairs = 0;
    int natom = atoms.size();
    for(int i=0; i<natom; i++){
        FARFFAtom& A = atoms[i];
        for(int j=i+1; j<natom; j++){
            Epairs += A.evalInteraction( atoms[j], ptype0 );
        }
    }
    return Epairs;
}

inline double evalAtoms(){
    Eatoms=0;
    //for(int i=0; i<natom; i++){ Eatoms+=evalAtom(i);        }
    for(FARFFAtom& a:atoms){ a.evalOnSite(); }
    return Eatoms;
}

void normalizeOrbs(){ for(FARFFAtom& a:atoms){ a.normalizeOrbs(); } }
void transferOrbRecoil(){ for(FARFFAtom& a:atoms){ a.orbRecoil(); } }
void removeNormalForces(){ for(FARFFAtom& a:atoms){ a.removeNormalForces(); } }

inline int nDOFs(){ return atoms.size()*(1 + N_BOND_MAX ); }

inline void exportPoss( Vec3d* ps ){
    int j=0; int na=atoms.size();
    for(int i=0; i<na; i++){ ps[j]=atoms[i].pos; j++; };
    for(int i=0; i<na; i++){ for(int k=0; k<N_BOND_MAX; k++ ){ ps[j]=atoms[i].opos[k]; j++; } };
}

inline void importPoss( Vec3d* ps ){
    int j=0; int na=atoms.size();
    for(int i=0; i<na; i++){ atoms[i].pos=ps[j]; j++; };
    for(int i=0; i<na; i++){ for(int k=0; k<N_BOND_MAX; k++ ){ atoms[i].opos[k]=ps[j]; j++; } };
}

inline void exportForce( Vec3d* ps ){
    int j=0; int na=atoms.size();
    for(int i=0; i<na; i++){ ps[j]=atoms[i].force; j++; };
    for(int i=0; i<na; i++){ for(int k=0; k<N_BOND_MAX; k++ ){ ps[j]=atoms[i].oforce[k]; j++; } };
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

void moveGD(double dt ){
    for(FARFFAtom& a:atoms){ a.moveGD(dt); }
}


}; // FFsp

#endif
