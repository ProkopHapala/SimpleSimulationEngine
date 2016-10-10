
#ifndef  SplineManager_h
#define  SplineManager_h

#include "spline_hermite.h"
#include "arrayAlgs.h"

#include <cstring>

class SplineManager{
    int n;
    int m;
    double *   ts;
    double **  CPs;
    double ** dCPs; // if explicit derivatives are avaible

    double getIt( int ip, double t,       int m_, const int * which, double * val, double * dval, double * ddval ){
        double t0  = ts[ip  ];
        double ut  = t-t0;
        double t1  = ts[ip+1];
        double dt  = t1-t0;
        double dtm = t0-ts[ip-1];
        double dtp = ts[ip+2]-t1;
        for( int im=0; im<m_; im++){
            double p0,p1,dp0,dp1;
            double * cpis = CPs[im];
            if(dCPs[im]){ // if explicit derivs avaible

            };
            // ToDo 1 - should be non uniform => rescale derivatives etc.
            // ToDo 2 - maybe using using hermite basis function will be faster
            if(  val) val[im] = Spline_Hermite::dval<double>( t,    cpis[ip-1], cpis[ip], cpis[ip+1], cpis[ip+2] );
            if( dval) val[im] = Spline_Hermite::dval<double>( t,    cpis[ip-1], cpis[ip], cpis[ip+1], cpis[ip+2] );
            if(ddval) val[im] = Spline_Hermite::dval<double>( t,    cpis[ip-1], cpis[ip], cpis[ip+1], cpis[ip+2] );
        };
        return ut;
    }

    int get( double t, int i0, int m_, const int * which, double * val, double * dval, double * ddval ){
        int ip = binSearch<double>( t, i0, ts );
        getIt( ip, t,  m_, which, val, dval, ddval );
        return ip;
    }

    int getUniform( double tstart, double tend, double dt, int m_, const int * which, double ** val, double ** dval, double ** ddval ){
        double T  = tend - tstart;
        int    nt = T/dt;
        int ip = binSearch<double>( tstart, 0, ts );
        double t = tstart;
        for(int i=0; i<nt; i++ ){
            getIt( i, t, m_, which, val[i], dval[i], ddval[i] ); // this would not work, arrays should be transposed
            t += dt;
            while( t>ts[ip] ) ip++;
        }
    }

    int removePoint( int i ){
        for( int im=0; im<m; im++){
            double * old = CPs[im];
            if(old){
                double *dest = new double[n-1];
                std::memcpy(old    ,dest  ,(i  )*sizeof(double) );
                std::memcpy(old+i+1,dest+i,(n-i)*sizeof(double));
                CPs[im] = dest;
                delete old;
            }
            old = dCPs[im];
            if(old){
                double *dest = new double[n-1];
                std::memcpy(old    ,dest  ,i  );
                std::memcpy(old+i+1,dest+i,n-i);
                dCPs[im] = dest;
                delete old;
            }
        }
        n--;
    }

    int insertPoint( double t, int m_, const int * which, double * val, double * dval ){
        // if val[which] or dval[which] not avaible it will be computed from interpolation

    }

};

#endif



