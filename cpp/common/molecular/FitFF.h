
#ifndef FitFF_h
#define FitFF_h

#include <vector>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
//#include "BasisFunc.h"

class AtomPair{ public:
    Vec2i  ij;
    double r;    // distance
    Vec3d  rhat; // normalized direction
    //RadialBasis* bas;
    int paramOffset;
};

void* relocate( int nbytes, char* arr ){
    void* arr_ = new char[nbytes];
    memcpy( arr_, arr, nbytes );
    delete [] arr;
    return arr_;
};


//template<typename FuncEF, typename FuncVar>
class FitStructure{ public:

    Mat3d  lvec;

    int natoms;
    int  * types = NULL;
    Vec3d* poss  = NULL;
    Vec3d* Fs    = NULL;
    Vec3d* Frefs = NULL;

    double E,Eref,Eerr;

    int npairs=0;
    AtomPair* pairs = NULL;
    //std::vector<AtomPair> pairs;

    double basisEF      ( double r, double* params, double& F ){ return 0.; };
    double basisVarDeriv( double r, double* params, double* dparams, double Ffactor, double Eerr ){ return 0.; };

    //int makePairList( double Rcut, Mat3d lvec, Vec3i npbc, int natoms, Vec3d* poss, int ntypes, int* types, int paramPerType, std::vector<AtomPair>& pairs ){
    int makePairs( double Rcut, Vec3i npbc, int ntypes, int paramPerType ){
        int n=0;
        double R2 = Rcut*Rcut;
        for(int i=0; i<natoms; i++){
            Vec3d pi  = poss[i];
            int itype = types[i];
            for(int ia=-npbc.a; ia<npbc.a; ia++){
                for(int ib=-npbc.b; ib<npbc.a; ib++){
                    for(int ic=-npbc.c; ic<npbc.a; ic++){
                        Vec3d pi_ = pi + lvec.a*ia + lvec.b*ib + lvec.c*ic;
                        for(int j=0; j<natoms; j++){
                            if( (i==j)&&(0==ia==ib==ic) ){
                                Vec3d d = poss[j] + pi_;
                                double r2 = d.norm2();
                                if( r2<R2 ){
                                    double r = sqrt(r2);
                                    d.mul(1/r);
                                    int jtype  = types[j];
                                    int ijtype = (jtype<itype)? itype*ntypes+jtype : jtype*ntypes+itype ;
                                    //pairs.push_back( {(Vec2i){i,j},r,d,ijtype*paramPerType} );
                                    pairs[n] = (AtomPair){(Vec2i){i,j},r,d,ijtype*paramPerType};
                                    n++;
                                }
                            }
                        }
                    }
                }
            }
        }
        npairs=n;
        pairs = (AtomPair*)relocate( npairs*sizeof(AtomPair), (char*)pairs );
        return n;
    }

    //template<typename Func>
    //double getPairwiseEF( int npairs, AtomPair* pairs, double * params, Vec3d* Fs, Func basisEF ){
    double getPairwiseEF( double * params, Vec3d* Fs ){
        int iold=0;
        double Fi = 0;
        double E  = 0;
        for( int ip=0; ip<npairs; ip++ ){
            AtomPair& pij = pairs[ip];
            int i  = pij.ij.a;
            if(i!=iold){ Fs[iold]=pij.rhat*Fi; iold=i; }
            E += basisEF( pij.r, params+pij.paramOffset, Fi );
        }
        return E;
    }

    //template<typename Func>
    //void getPairwiseVarDeriv( int npairs, AtomPair* pairs, double * params, double * dparams, double Eerr, Vec3d* Ferrs, Func basisVarDeriv ){
    void getPairwiseVarDeriv( double * params, double * dparams, double Eerr, Vec3d* Ferrs ){
        int iold=0;
        for( int ip=0; ip<npairs; ip++ ){
            AtomPair& pij = pairs[ip];
            int i  = pij.ij.a;
            double Ffactor = pij.rhat.dot( Ferrs[i] );
            basisVarDeriv( pij.r, params+pij.paramOffset, dparams+pij.paramOffset, Ffactor, Eerr );
        }
        //return E;
    }

    //template<typename FuncEF, typename FuncVar>
    // void getPairwiseVarDerivStep( int natoms, int npairs, AtomPair* pairs, double * params, double * dparams, double Eref, Vec3d* Frefs, Vec3d* Ferrs, FuncEF basisEF, FuncVar basisVarDeriv ){
    double getPairwiseVarDerivStep( double * params, double * dparams ){
        double E    = getPairwiseEF( params, Fs );
        Eerr = E - Eref;
        double F2sum=0.0;
        for(int i=0; i<natoms; i++){ Vec3d Fi=Fs[i]-Frefs[i]; F2sum+=Fi.norm2(); Fs[i]=Fi; }
        getPairwiseVarDeriv( params, dparams, Eerr, Fs );
        return F2sum;
    }

};


int replicatePBC( Mat3d lvec, Vec3i npbc, int natoms, int* types, Vec3d* poss, int*& types_, Vec3d*& poss_ ){
    int n=0;
    for(int i=0; i<natoms; i++){ types_[n]=types[i]; poss_[n]=poss[i]; n++; }
    for(int ia=-npbc.a; ia<npbc.a; ia++){
        for(int ib=-npbc.b; ib<npbc.a; ib++){
            for(int ic=-npbc.c; ic<npbc.a; ic++){
                if( ! (0==ia==ib==ic) ){
                   for(int i=0; i<natoms; i++){ types_[n]=types[i]; poss_[n]=poss[i]; n++; };
                }
            }
        }
    }
    return n;
}

/*
void GetVariational( int natoms, double Eerr, Vec3d * poss, Vec3d * Ferrs, int * atomTypes, RadialBasis* typeBasis, double** dparams ){
    for( int i=0; i<natoms; i++ ){
        int itype=atomTypes[i];
        // evaluate forces
        Vec3d pi = poss[i];
        //Vec3d Fi;


        // How to deal with energy ? ... take energy partitioning from previous step ?
        // evaluate fitting gradinet variational


        Vec3d Ferr_i = Ferrs[i];
        for( int j=0; j<natoms; j++ ){ // how to solve cell replication ?
            int jtype     = atomTypes[j];
            Vec3d dR      = poss[j] - pi;
            int ipairtype = itype*ntypes+jtype;
            double      *  dpairParams = dparams [ipairtype];
            PairBasis&  pairBasis = paramsOfTypes[ipairtype];
            //  dot( Ferr , dR ) ... somehow ... do chain derivs
            pairBasis.dParams( dR, dpairParams ); // somehow  do chain derivs
            //aforces[i]
        }
    }
}
*/

#endif
