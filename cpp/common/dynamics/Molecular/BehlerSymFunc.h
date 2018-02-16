#ifndef BehlerSymFunc_h
#define BehlerSymFunc_h

#include "fastmath.h"

#include "NN.h"

/*

Scheme:

* evalute sym function for 3 atoms:
  for all symFuncs of atom a:
    evaluate symFunction(a,b,c)
        store derivatives of Sf with respect to movement of neighbors
           dGdNeigh  // where dGdNeigh has size   nSf*nNeigh
    evaluate neural network
    evaluate atomic foces
        gather F =  dfdG[iSf] * dGdNeigh[iSf][ineigh]


*/

inline void dot3d( int n, Vec3d* d, double* s, Vec3d& out ){ for(int i=0; i<n; i++) out.add_mul(d[i],s[i]); }


inline void fcut( double x, double rScaling, double& y, double& dy ){
    x *= rScaling;
    y  = x*x*(3-x);
    dy = 6*x*(x-1)*rScaling;
}

class BehlerNNFF : public AnyShortRangeFF { public:

    NN nn;

    int nAtoms = 0;
    int nTypes=0;
    int nTypePairs=0;

    int mSfRad=0,nSfRad=0;
    int mSfAng=0,nSfAng=0;

    double Rcut     = 6.0;
    double rScaling = 1/Rcut;

    int nsf = 8;

    double*  SFR = NULL; //[nSfRadMax*nTypes];
    Vec3d*  dSFR = NULL; //[nSfRadMax*nTypes];

    double*  SFA = NULL; //[nSfAngMax*n];
    Vec3d*  dSFA = NULL; //[nSfAngMax*];

    double* atomEnerg = NULL;
    Vec3d*  atomForce = NULL;

    double * dfdG = NULL;

    int   *types = NULL; //   type of each atom

    //double* Es = NULL;
    //Vec3d * Fs = NULL;

    //double * nsfs = NULL;
    //double **params = NULL; // function params for each type combination

    //double * nsfsA   = NULL;
    //double **paramsA = NULL; // function params for each type combination

    //double *fvals  = NULL; // function values for each
    //double ***afvals  = NULL; // function values for each atom

    inline int type2i( int it0, int it1          ){ return nTypes*it0+it1; }
    inline int type3i( int it0, int it2, int it1 ){
        if(it2<it1){int itmp=it2; it2=it1; it1=itmp; }; // it1 > it2
        return it2 + ( it1*(it1-1) + nTypes*(nTypes+1)*it0 )/2;
    }

    inline int typePair( int it1, int it2 ){
        if(it2<it1){int itmp=it2; it2=it1; it1=itmp; }; // it1 > it2
        return it2 + it1*(it1-1);
    }

    //inline int iSfNeigh3( ){}
    inline int iSfNeigh3( int iTypePair, int  ){

    }

    int allocateSF( int nTypes_, int mSfRad_, int mSfAng_ ){
        mSfRad     = mSfRad_;
        mSfAng     = mSfAng_;
        nTypes     = nTypes_;
        nTypePairs = nTypes*(nTypes+1)/2;
        nSfRad = mSfRad*nTypes;
        nSfAng = mSfAng*nTypePairs;
        SFR    = new double[nSfAng];
        SFA    = new double[nSfAng];
        dSFR   = new Vec3d [nSfAng];
        dSFA   = new Vec3d [nSfAng*nAtoms];
    }

    /*
    double symfunc_R( int iatom, int ipar, double r1, const Vec3d& hi ){
        int nsf    = nsfs  [ipar];
        double* tp = params[ipar];
        for(int i=0; i<nsf; i++ ){
            double y,dy;
            symfunc_GaussCos( r1, Rmax, tp[0], tp[1], y, dy );
            tp += perSf;
        }
    };

    double symfunc_Ang( int iatom, int ipar, double ca, double r1, double r2, const Vec3d& h1, const Vec3d& h2 ){
        int nsf    = nsfsA  [ipar];
        double* tp = paramsA[ipar];
        for(int i=0; i<nsf; i++ ){
            double y,dy;
            //symfunc_GaussCosAng( ca,r1,r2, Rmax, tpi[0], tpi[i], y, dy_ca, dy_r1, dy_r2 );
            tp += perSf;
        }
    };
    */


    inline double symfunc_R( int iatom, int ipar, double r, const Vec3d& h ){
        //int npow   = nsfs  [ipar];
        //double* tp = params[ipar];

        double fc,dfc;
        fcut( r, rScaling, fc, dfc );
        double y=1,dy=1;
        for(int i=0; i<nSfAng; i++ ){
            dy *= y*dfc;
            y  *= y*fc;
        }
    };

    inline double symfunc_Ang( int iatom, int ipar, double ra, double rb, const Vec3d& ha, const Vec3d& hb, Vec3d * dsf ){
        //int nsf    = nsfsA  [ipar];
        //double* tp = paramsA[ipar];

        double fca, dfca; fcut( ra, rScaling, fca, dfca );
        double fcb, dfcb; fcut( rb, rScaling, fcb, dfcb );

        //  in FOrtran:
        //  sf = sf + fa*faN1*pfc
        //  pfc     = 0.5*pfc*zeta
        //  df_drij = faN1*( pfc*( rinvijik - costheta/r2ij ) + fa*pfcik*pfcjk*pdfcij/rij )  !  df/drij
        //  df_drik = faN1*( pfc*( rinvijik - costheta/r2ik ) + fa*pfcij*pfcjk*pdfcik/rik )  !  df/drik
        //  df_drjk = faN1*( pfc*( rinvijik                 ) - fa*pfcij*pfcik*pdfcjk/rjk )  !  df/drjk    ! is it correct? costheta independent of rjk ?

        //  without ab
        //df_drij = faN1*( pfc*( rinvijik - costheta/r2ij ) + fa*pdfcij* pfcik/rij )  !  df/drij
        //df_drik = faN1*( pfc*( rinvijik - costheta/r2ik ) + fa* pfcij*pdfcik/rik )  !  df/drik
        //df_drjk = faN1*( pfc*( rinvijik                 )                        )  !  why    rinvijik  ?


        double invRa   = 1/ra;
        double invRb   = 1/rb;
        double invRaRb = invRa*invRb;

        double ca     = ha.dot(hb);
        double fa     = 0.5+0.5*ca;
        double cfc    = 0.5*fca*fcb;
        double fa_da = cfc*( invRaRb - ca*ra*ra );
        double fa_db = cfc*( invRaRb - ca*rb*rb );
        double fa_dc = cfc*( invRaRb            ); // + fa*pfcij*pfcik*pdfcjk/rjk;
        double fc    = fca*fcb;
        double fc_da = fa*dfca* fcb*invRa;
        double fc_db = fa* fca*dfcb*invRb;

        Vec3d hc = ha - hb; // FIXE: this works if they are not normalized

        double fN=1;
        Vec3d * dsfa = dsf+mSfAng;
        Vec3d * dsfb = dsf+mSfAng*2;
        for(int i=0; i<nsf; i++ ){
            double f_da = fN*( fa_da*i + fc_da );
            double f_db = fN*( fa_db*i + fc_db );
            double f_dc = fN*( fa_db*i         );
            fN  *= fa;
            double f    = fN * fc;

            // --- dGdR
            Vec3d da = ha*f_da;
            Vec3d db = hb*f_db;
            Vec3d dc = hc*f_dc;
            dsf [i] = da+db;
            dsfa[i] = da+dc; // TODO check the signs
            dsfb[i] = db-dc;

            //
            // ?????  OVER WHAT THEY SUM ???????
            //
            //xptr[Nnum] += p1 * dxij + p2 * dxik; // summing? over what index ?
            //yptr[Nnum] += p1 * dyij + p2 * dyik;
            //zptr[Nnum] += p1 * dzij + p2 * dzik;
            //xptr[j] -= p1 * dxij + p3 * dxjk;
            //yptr[j] -= p1 * dyij + p3 * dyjk;
            //zptr[j] -= p1 * dzij + p3 * dzjk;
            //xptr[k] -= p2 * dxik - p3 * dxjk;
            //yptr[k] -= p2 * dyik - p3 * dyjk;
            //zptr[k] -= p2 * dzik - p3 * dzjk;

        }
    };

    void assembleAngForce( int nbonds, int kAtom, int* neigh_is){
        int type0 = types[ kAtom ];
        //int iAng3 = 0;
        Vec3d* dSFAi = dSFA;
        for(int ib=0; ib<nbonds; ib++){     // neigh i
            int  iatom = neigh_is[ib];
            int  type1 = types[ iatom ];
            for(int jb=0; jb<nbonds; jb++){ // neigh j
                if(ib==jb) continue;
                int jatom = neigh_is[ib];
                int type2 = types[ jatom];
                int iTypePair = typePair( type1, type2 );
                //int iAng3     = iAng*3;
                // Note:
                //  - dSFA  [ iAng, {k|i|j} , mSfAng ]
                //  - dfdGi [ iTypePair,      mSfAng ]
                double* dfdGi = dfdG+(iTypePair*mSfAng);
                dot3d( mSfAng, dSFAi, dfdGi, atomForce[kAtom] );  dSFAi += mSfAng;
                dot3d( mSfAng, dSFAi, dfdGi, atomForce[iatom] );  dSFAi += mSfAng;
                dot3d( mSfAng, dSFAi, dfdGi, atomForce[jatom] );  dSFAi += mSfAng;
                //iAng3+=3;
            }
        }
    }


    virtual void processAtomNeighbors( int nbonds, int iatom, Vec3d* hbonds, double* lbonds, int* neigh_is){
        int type0 = types[ iatom ];
        Vec3d* dSFAi = dSFA;
        for(int ib=0; ib<nbonds; ib++){
            Vec3d  hi    = hbonds[ib];
            double ri    = lbonds[ib];
            int    type1 = types[ neigh_is[ib] ];
            int ipar = type2i( type0, type1 );
            symfunc_R( iatom, ipar, ri, hi );
            for(int jb=0; jb<nbonds; jb++){
                if(ib==jb) continue;
                const Vec3d& hj = hbonds[jb];
                //double ca = hi.dot( hj );
                int    type2 = types[ neigh_is[ib] ];
                int iparA    = type3i( type0, type1, type2);

                symfunc_Ang( iatom, iparA, ri, lbonds[jb], hi, hj, dSFAi );
                dSFAi += 3*mSfAng;
            }
        }
    }

    double getEnergy(){

    }

    void getForces(){

    }

};

#endif
