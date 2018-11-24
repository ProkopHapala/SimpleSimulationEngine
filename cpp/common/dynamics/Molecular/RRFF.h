
#ifndef RRFF_h
#define RRFF_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Forces.h"
//#include "GridFF.h"

struct RRFFtype{
    double R0;
    Vec2d  rot0;
    double dl0;
    double ef0;
};

struct RRFFelement{
    double Rmax = 3.0;
    double l0;
    double e0;
    int nkinds;
    RRFFtype *kinds[8];
};

class RRFFparams{ public:
    std::vector<RRFFelement> elems;
    std::vector<RRFFtype   > types;
};

struct NeighAtom{
    int iatom=1;
    double lbond;
    double expar; // exp(-a*r)
    Vec3d  hbond;
};

inline double getTripleEnergy( double lj, const Vec3d& hi, const Vec3d& hj, double l0, double aMorse, const Vec2d& rot0 ){
    //double ca = hat0.dot( hat );
    //double sa = sqrt(1-ca*ca);  // perhaps can be optimized ?
    Vec2d rot; rot.fromCos( hi.dot(hj) );
    rot.mul_cmplx (rot0);
    double Eang  = -rot.x + sq(rot.y); // Angular energy
    double expar = exp( aMorse*(lj-l0) );      //  TODO : this can be optimized by precalculation or "exapr"
    return ( expar*expar + Eang*2.0*expar );
}

class RRFF{ public:

    double Rcut   = 2.5;
    double aMorse = 1.8;


    RRFFparams * params = 0;

    const int nMaxKind = 4;
    int nNeighMax = 16;

    // atom positions
    int natoms=0, npairs=0;
    int*    aelems = 0;
    int*    atypes = 0;
    Vec3d*  apos   = 0;

    int*        nneighs = 0;
    NeighAtom * neighs  = 0;

    /*
    Vec2i*  pairs   = 0;
    Vec3d*  hpairs  = 0;
    double* lpairs  = 0;
    */

    void findPairs_On2( int nbmax ){
        //if(nbmax<0) nbmax=natoms**2;
        /*
        if( nbmax>0 ){
            _realloc(pairs,  nbmax);
            _realloc(hpairs, nbmax);
            _realloc(lpairs, nbmax);
        }
        */
        double R2cut = Rcut*Rcut;
        int ib = 0;
        for(int i=0; i<natoms; i++){
            Vec3d pi = apos[i];
            int nneigh = 0;
            int ing0 = i*nNeighMax;
            //for(int j=i; j<natoms; j++){ // TODO for simplicity we search all atoms not onlu J<i
            for(int j=0; j<natoms; j++){ // TODO for simplicity we search all atoms not onlu J<i
                if(j==i) continue;
                Vec3d  d = apos[j] - pi;
                double r2 = d.norm2();
                // TODO R2cut can be obtained from forcefield
                if(r2 > R2cut) continue;
                if( ib > nNeighMax ){
                    printf( "ERROR : atom %i has more than %i neighbors \n", i, nNeighMax );
                    exit(0);
                }
                NeighAtom& ng = neighs[ing0+ib];
                ng.iatom      = j;
                ng.lbond      = d.normalize();
                ng.hbond      = d;
                ib++;
                /*
                pairs [ib] = {i,j};
                lpairs[ib] = d.normalize();
                hpair [ib] = d;
                ib++;
                */
            }
            nneighs[i] = ib;
        }
    }

    /*
    void selectType( int i){
        nng = nneighs[i];
        NeighAtom* ngs = neighs + i*nNeighMax;

        RRFFelem& elem = params->elems[aelems[i]];
        int    Ns[nMaxKind];
        //double Es[nMaxKind];
        double Ss[nMaxKind]; // overlap or interaction strenght
        // initialize
        for(int k=0; k<elem.nkinds; k++){
            Ns[k] = 0.0;
            Ss[k] = 0.0;
        };

        // estimate energy for each type
        for(int ing=0; ing<nngs; ing++){ // loop over neighs
            double l  = ngs[ing].lbond;
            double Ss = exp(aMorse*(l-Ri-Rj));
            for(int k=0; k<elem.nkinds; k++){ // loop over types
                atype = ;
                // can calculate also angular forces ?
            }
        }

        // find type with minimum energy


        double rb_type = rb_type;
        double n = type.weight( ngs.lbond, Rtype );

    }
    */

    double minimizeAtomEnergy( int i ){
        int          nng  = nneighs[i];
        NeighAtom*   ngs  = neighs + i*nNeighMax;
        RRFFelement& elem = params->elems[aelems[i]];

        double ngs_l0[nNeighMax];
        double ngs_e0[nNeighMax];

        for(int j=0; j<nng; j++){
            RRFFelement& elemj = params->elems[aelems[j]];
            ngs_l0[i] = elemj.l0;
            ngs_e0[i] = elemj.e0;
        }

        // calculate energy for each kind, than take the lowest energy ( min() or softmax() )
        //   relative weight of force-fields is given by sum weighted by radial repulsive energy ( i.e. closest atoms has highest contribution )

        double Etot = 1e+300;

        for(int k=0; k<elem.nkinds; k++){
            double sumPi   = 0.0;
            double Ek      = 0.0;
            RRFFtype& kind = *elem.kinds[k];
            for(int i=0; i<nng; i++){
                Vec3d  hi = ngs[i].hbond;
                double li = elem.l0 + kind.dl0;
                double ei = elem.e0 * kind.ef0;
                double Ei = 0.0;
                for(int j=0; j<nng; j++){
                    if(j==i) continue;
                    double Eij = ei * ngs_l0[j];
                    double l0  = li + ngs_e0[j];
                    Ei += Eij * getTripleEnergy( ngs[j].lbond, hi, ngs[j].hbond, l0, aMorse, kind.rot0 );
                }
                double pi = ngs[i].expar; // contribution according to closenest
                Ek    += Ei * pi;
                sumPi += pi;
            }
            Ek /= sumPi; // normalize by weight
            if(Ek<Etot) Etot = Ek;
        }
        return Etot;
    }

};



#endif
