
#ifndef FTRFF_h
#define FTRFF_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Forces.h"

/*

Fixed Type Reactive Force-field
===============================


force on orbitals
E = exp(b*rijk)

dijk = pi - pj - hjk*l0
lijk = |dijk|           = sqrt( (xi-xj-xjk*l0)^2 + (yi-yj-yjk*l0)^2 + (zi-zj-zjk*l0)^2 )

d_lij/d_xi  =   (xi-xj-xjk*l0)/lijk =    dijk.x/lijk
d_lij/d_xj  =   (xi-xj-xjk*l0)/lijk =    dijk.x/lijk
d_lij/d_xjk =   (xi-xj-xjk*l0)/lijk = l0*dijk.x/lijk



*/

#define N_BOND_MAX 4


inline double fBondRepel( double K, const Vec3d& ha, const Vec3d& hb, Vec3d& fa, Vec3d& fb ){
    Vec3d f = ha-hb;
    f.mul( K/f.norm2() );
    fa.add(f);
    fb.sub(f);
    return 0.0; // TODO: Energy
}

inline double fortho( double K, const Vec3d& ha, const Vec3d& hb, Vec3d& fa, Vec3d& fb ){
    double c = ha.dot(hb);
    fa.add_mul(ha,c);
    fb.add_mul(hb,c);
    return 0.0; // TODO: Energy
}



inline void expFE(double r, double beta, double amp, double& f, double& df ){ 
    double expar = amp*exp(r*beta);
    f  = expar; 
    df = expar*beta; 
}


inline void overlapFE(double r, double beta, double amp, double& f, double& df ){ 
    //https://www.wolframalpha.com/input/?i=exp(b*x)*(1%2Bb*x%2B(b*x)%5E2%2F3)+derivative+by+x
    double x = r*beta;
    double expar = amp*exp(x);
    f  = expar*(1 +   x + 0.33333333*x*x ); 
    df = expar*(6 + 5*x +            x*x )*beta*0.33333333; 
}



/*
inline fstright( double K, const Vec3d& ha, const Vec3d& hb, Vec3d& fa, Vec3d& fb ){
    double f = ha+hb;
    fa.add_mul(ha,c);
    fb.add_mul(hb,c);
}
*/

struct FTtype{
    enum kinds{ s=0, sp1=1, sp2=2, sp3=3, sp3O=4, sp2N=5 };
    int    nbond;  // number bonds
    double rbond0;
    double acore;
    double bcore;
    double abond; 
    double bbond;
};

struct FTAtom{
    //int ihyb         = 0; // 0:H,Cl   1=sp3   2=sp2  3=sp1  4=kink(-OH)   5=O   6=CN
    //FTtype* params = 0;
    FTtype* type = 0;
    Vec3d pos;
    Vec3d force;
    Vec3d hbond[N_BOND_MAX];
    Vec3d fbond[N_BOND_MAX];

    inline double forceToPlane(double K){
        /*
        Vec3d a = hbond[1]-hbond[0];
        Vec3d b = hbond[2]-hbond[0];
        Vec3d nr; nr.set_cross(a,b);
        */
        Vec3d nr; nr.add(hbond[0]); nr.add(hbond[1]); nr.add(hbond[2]);
        nr.normalize();
        hbond[3] = nr; // will be used later for dihedrals
        //K  = K/nr.norm();
        double f0 = nr.dot(hbond[0])*K; fbond[0].add_mul( nr, f0 );
        double f1 = nr.dot(hbond[1])*K; fbond[1].add_mul( nr, f1 );
        double f2 = nr.dot(hbond[2])*K; fbond[2].add_mul( nr, f2 );
        //force.add_mul( nr, -f2-f1-f2 );
        return 0.0; // TODO: Energy
    }

    inline double forceStright(double K){
        Vec3d nr; nr.add(hbond[0]); nr.add(hbond[1]);
        //nr.normalize();
        //hbond[3] = nr; // will be used later for dihedrals
        //K  = K/nr.norm();
        //double f0 = nr.dot(hbond[0])*K; fbond[0].add_mul( nr, f0 );
        //double f1 = nr.dot(hbond[1])*K; fbond[1].add_mul( nr, f1 );
        nr.mul(K);
        fbond[0].add( nr );
        fbond[1].add( nr );
        //force.add_mul( nr, -f2-f1 );
        return 0.0; // TODO: Energy
    }

    inline double fsp1( double Knr, double Kb ){
        fBondRepel( Kb, hbond[0],hbond[1], fbond[0],fbond[1] );
        fBondRepel( Kb, hbond[0],hbond[1], fbond[0],fbond[1] );
        fBondRepel( Kb, hbond[0],hbond[1], fbond[0],fbond[1] );
        forceToPlane( Knr);
        return 0.0; // TODO: Energy
    };

    inline double fsp2( double Knr, double Kb ){
        fBondRepel( Kb, hbond[0],hbond[1], fbond[0],fbond[1] );
        fBondRepel( Kb, hbond[0],hbond[1], fbond[0],fbond[1] );
        fBondRepel( Kb, hbond[0],hbond[1], fbond[0],fbond[1] );
        forceToPlane( Knr);
        return 0.0; // TODO: Energy
    };

    inline double fsp3( double Knr, double Kb ){
        fBondRepel( Kb, hbond[0],hbond[1], fbond[0],fbond[1] );
        fBondRepel( Kb, hbond[0],hbond[1], fbond[0],fbond[1] );
        fBondRepel( Kb, hbond[0],hbond[1], fbond[0],fbond[1] );
        fBondRepel( Kb, hbond[0],hbond[1], fbond[0],fbond[1] );
        return 0.0; // TODO: Energy
    };

};


class FTRFF{ public:

    int natom = 0;
    FTtype* types=0;
    FTAtom* atoms=0;

    void interEF(){
        
        for(int i=0; i<natom; i++){
            FTAtom& atomi = atoms[i];
            FTtype& typei = *atomi.type;
            int   nbi     = typei.nbond;
            Vec3d pi      = atomi.pos;

            for(int j=0; j<natom; j++){
                FTAtom& atomj = atoms[j];
                FTtype& typej = *atomj.type;
                int nbi       = typei.nbond;
                //double rbj = rbi; // TODO: correct this

                double acore = typei.acore * typej.acore ;  // TODO
                double bcore = typei.bcore + typej.bcore ;
                double abond = typei.abond * typej.abond ; 
                double bbond = typei.bbond + typej.bbond ;

                // core-core interactions
                Vec3d  dij = atomj.pos - pi;
                double rij = dij.norm();

                double eij = acore*exp(bcore*rij);
                Vec3d  fij = dij*( eij*bcore );

                atomi.force.sub(fij);
                atomj.force.add(fij);

                double rbond0 = typei.rbond0 + typej.rbond0;

                // bonds interactions
                for(int ib=0; ib<nbi; ib++){

                        Vec3d  d  = dij + atomi.hbond[ib]*rbond0;
                        double r  = d.norm();

                        double e = abond*exp(bbond*r);
                        Vec3d f  = d*( e*bbond );

                        atomi.fbond[ib].add_mul(f,r); // TODO: think how this is transfered to the atom ?
                        atomi.force.add(f);
                        atomj.force.sub(f);

                }
            }
        }

    }

};



#endif
