#ifndef BehlerSymFunc_h
#define BehlerSymFunc_h

inline void symfunc_GaussCos( double x, double xcut, double eta, double xshift, double& y, double& dy ){
    double a   = M_PI*x/xcut;
    double ca  = 0.5+cos(a)*0.5;
    double dca = sin(a)*-0.5;
    double x_  = x-xshift;
    double  g  =  exp(-eta*x_*x_);
    double dg  = -eta*x_*g;
     y = ca*g;
    dy = dca*g + ca*dg;
}

inline void symfunc_GaussCosAng( double ca, double r1, double r2, double xcut, int npow, double eta, double& y, double& dy_ca, double& dy_r1, double& dy_r2 ){
/*
    double a   = M_PI*x/xcut;
    double ca  = 0.5+cos(a)*0.5;
    double dca = sin(a)*-0.5;
    double x_  = x-xshift;
    double  g  =  exp(-eta*x_*x_);
    double dg  = -eta*x_*g;
     y = ca*g;
    dy = dca*g + ca*dg;
*/
}

class BehlerNNFF : public AnyShortRangeFF { public:

    double Rmax = 6.0;

    int perSf  = 2;
    int natoms = 0;
    int ntypes = 0;

    int   *types = NULL;

    //double* Es = NULL;
    //Vec3d * Fs = NULL;

    double * nsfs = NULL;
    double **params = NULL; // function params for each type combination

    double * nsfsA   = NULL;
    double **paramsA = NULL; // function params for each type combination

    //double *fvals  = NULL; // function values for each
    //double ***afvals  = NULL; // function values for each atom

    inline int type2i( int it0, int it1          ){ return ntypes*it0+it1; }
    inline int type3i( int it0, int it2, int it1 ){
        if(it2<it1){int itmp=it2; it2=it1; it1=itmp; }; // it1 > it2
        return it2 + ( it1*(it1-1) + ntypes*(ntypes+1)*it0 )/2;
    }

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

    virtual void processAtomNeighbors( int nbonds, int iatom, Vec3d* hbonds, double* lbonds, int* neigh_is){
        int type0 = types[ iatom ];
        for(int ib=0; ib<nbonds; ib++){
            Vec3d  hi    = hbonds[ib];
            double ri    = lbonds[ib];
            int    type1 = types[ neigh_is[ib] ];
            int ipar = type2i( type0, type1 );
            symfunc_R( iatom, ipar, ri, hi );
            for(int jb=0; jb<nbonds; jb++){
                if(ib==jb) continue;
                const Vec3d& hj = hbonds[jb];
                double ca = hi.dot( hj );
                int    type2 = types[ neigh_is[ib] ];
                int iparA = type3i( type0, type1, type2);
                symfunc_Ang( iatom, iparA, ca, ri, lbonds[jb], hi, hj );
            }
        }
    }

};

#endif
