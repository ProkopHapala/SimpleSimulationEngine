
#ifndef  SplineManager_h
#define  SplineManager_h

#include "spline_hermite.h"
#include "arrayAlgs.h"

#include <cstring>

class SplineManager{
    public:
    int n;
    int m;
    double *   ts;
    double **  CPs;
    double ** dCPs; // if explicit derivatives are avaible

    void allocate( int m_, int n_, bool derivs = false){
        n=n_; m=m_;
        ts =new double[n];
        CPs =new double*[m];  for(int i=0;i<m;i++){  CPs[i]=new double[n];  }
        dCPs=new double*[n];  for(int i=0;i<m;i++){ if(derivs){ dCPs[i]=new double[n];}else{ dCPs[i]=NULL; } }
    }

    void deallocate(int n_, int m_){
        for(int i=0;i<m;i++){
            if( CPs[i]) delete  CPs[i];
            if(dCPs[i]) delete dCPs[i];
        }
        if(ts) delete ts;
        delete  CPs;
        delete dCPs;
    }

    inline double inferDeriv( int ip, int im ){
        double * dCPi = dCPs[im];
        double dt   =   ts[ip+1] -   ts[ip-1];
        double dval = dCPi[ip+1] - dCPi[ip-1];
        return dval/dt;
    }

    inline double getPointDeriv( int ip, int im ){
        double * dCPi = dCPs[im];
        if( dCPi ){ return dCPi[ip]; }else{ return inferDeriv(ip,im ); };
    }

    double evalIt( int ip, double t,       int m_, const int * which, double * val, double * dval, double * ddval ){
        double t0     = ts[ip  ];
        double t1     = ts[ip+1];
        double dt     = t1-t0;
        double invdt  = 1/dt;
        double u      = (t-t0)*invdt;
        double invdtm = 1/( t0-ts[ip-1] );
        double invdtp = 1/( ts[ip+2]-t1 );

        double r0,r1,rm,rp;
        double v0,v1,vm,vp;
        double a0,a1,am,ap;
        if( val  ){ Spline_Hermite::  basis<double>( u,  r0,r1,rm,rp ); rm*=dt;    rp*=dt;    }
        if( dval ){ Spline_Hermite:: dbasis<double>( u,  v0,v1,vm,vp ); v0*=invdt; v1*=invdt; }
        if( ddval){ Spline_Hermite::ddbasis<double>( u,  a0,a1,am,ap ); double invdt2=invdt*invdt; a0*=invdt2; a1*=invdt2; am*=invdt; ap*=invdt; }

        //printf("ip %i t1 %3.3f t2 %3.3f u %3.3f %3.3f %3.3f \n", ip, t0,t1,u,invdtm,invdtp  );
        //if( val  )  printf("r0,r1,rm,rp %3.3f %3.3f %3.3f %3.3f ", r0,r1,rm,rp );
        //if( dval )  printf("a0,a1,am,ap %3.3f %3.3f %3.3f %3.3f ", v0,v1,vm,vp );
        //if( dval )  printf("a0,a1,am,ap %3.3f %3.3f %3.3f %3.3f ", a0,a1,am,ap );
        //printf("\n");

        //printf("pointers %i %i %i \n", val, dval, ddval  );

        for( int im=0; im<m_; im++){
            int j = which[im];
            double p0,p1,d0,d1;
            double * CPi = CPs[j];
            p0 = CPi[ip  ];
            p1 = CPi[ip+1];

            double * dCPi = dCPs[j];
            if( dCPi ){ // if explicit derivs avaible
                d0 = dCPi[ip  ];
                d1 = dCPi[ip+1];
            }else{
                d0 = ( p1        - CPi[ip-1] )*invdtm;
                d1 = ( CPi[ip+2] - p0        )*invdtp;
            };

            if(    val )   val[im] =  p0*r0 + p1*r1 + d0*rm + d1*rp;
            if(   dval )  dval[im] =  p0*v0 + p1*v1 + d0*vm + d1*vp;
            if(  ddval ) ddval[im] =  p0*a0 + p1*a1 + d0*am + d1*ap;

            // ToDo 1 - should be non uniform => rescale derivatives etc.
            // ToDo 2 - maybe using using hermite basis function will be faster
            //if(  val) val[im] = Spline_Hermite::dval<double>( t,    cpis[ip-1], cpis[ip], cpis[ip+1], cpis[ip+2] );
            //if( dval) val[im] = Spline_Hermite::dval<double>( t,    cpis[ip-1], cpis[ip], cpis[ip+1], cpis[ip+2] );
            //if(ddval) val[im] = Spline_Hermite::dval<double>( t,    cpis[ip-1], cpis[ip], cpis[ip+1], cpis[ip+2] );

        };

        return t-t0;
    }

    int eval( double t, int i0, int m_, const int * which, double * val, double * dval, double * ddval ){
        int ip = binSearchFrom<double>( t, (n-i0), ts+i0 ) + i0;
        evalIt( ip, t,  m_, which, val, dval, ddval );
        return ip;
    }

    int evalUniform( double tstart, double tend, double dt, int m_, const int * which, double * val, double * dval, double * ddval ){
        double T  = tend - tstart;
        int    nt = T/dt;
        int ip = binSearchFrom<double>( tstart, n, ts );

        double t = tstart;
        for(int i=0; i<nt; i++ ){
            //printf(" %i %i %f \n", i, ip, t);
            evalIt( ip, t, m_, which, val, dval, ddval ); // this would not work, arrays should be transposed
            t += dt;
            if(  val)   val+=m_;
            if( dval)  dval+=m_;
            if(ddval) ddval+=m_;
            while( t>ts[ip+1] ) ip++;
        }
        return nt;
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
        return -1; // just place holder
    }

};

#endif



