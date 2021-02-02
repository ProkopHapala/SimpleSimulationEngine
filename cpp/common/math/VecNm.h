
#ifndef  VecNm_h
#define  VecNm_h

#include "VecN.h"

//#include <vector>
//#include <functional>
//#include <type_traits>
//#include <utility>
//#include "macroUtils.h"
//#include "fastmath.h"

template<typename T>
struct VecNm{
    int n,m;
    T* data;
    VecNm(int n_,int m_, T* data_):n(n_),m(m_),data(data_){};

    template<typename Func>
    void iter_N( Func func, T* out, int imin=0, int imax=-1){
        if(imin<0) imin=n+imin;
        if(imax<0) imax=n+imax;
        for(int i=imin; i<=imax; i++){
            //printf( "[%i] (%i,%i)\n", i, imin, imax );
            int j0=i*m;
            for(int j=0; j<m;j++){
                func( data[j0+j], out[j] );
                //out[j] = fmax( out[j], data[j0+j] );
            }
        }
    }

    template<typename Func>
    void iter_Nmat( Func func, T* M, int imin=0, int imax=-1 ){
        if(imin<0) imin=n+imin;
        if(imax<0) imax=n+imax;
        T tmpIn[m];
        for(int i=imin; i<=imax; i++){
            int j0=i*m;
            for(int j=0; j<m;j++){ tmpIn[j]=data[j0+j]; }
            func( m, tmpIn, data+j0, M );
        }
    }

    void apply_mat( T* M, int imin=0, int imax=-1){
        iter_Nmat( [](int m, T* vout, const T* vin, const T* M){
            for(int i=0; i<m; i++){ vout[i] = VecN::dot(vin,M[i*m]); }
        }, M, imin, imax );
    }

    void outer( T* M, int imin=0, int imax=-1){
        iter_Nmat(  [](int m, const T* a, const T* b, T* M){
            for(int i=0; i<m; i++){
                int j0=i*m; double ai=a[i];
                for(int j=0; j<m; j++){
                    M[j0+j] = ai*b[j];
                }
            }
        }, M, imin, imax );
    }

    void max( T* out, int imin=0, int imax=-1 ){ iter_N( [](T x,T& y){ y=fmax(y,x); }, out, imin, imax ); }
    void min( T* out, int imin=0, int imax=-1 ){ iter_N( [](T x,T& y){ y=fmin(y,x); }, out, imin, imax ); }
    void sum( T* out, int imin=0, int imax=-1 ){ iter_N( [](T x,T& y){ y+=x;        }, out, imin, imax ); }

    void shift( T* vec, int imin=0, int imax=-1 ){ iter_N( [](T& x,T y){ x+=y; }, vec, imin, imax ); }
    void scale( T* vec, int imin=0, int imax=-1 ){ iter_N( [](T& x,T y){ x*=y; }, vec, imin, imax ); }

    void average( T* out, int imin=0, int imax=-1 ){
        sum_N( out, imin, imax );
        double w=1/m;
        for(int j=0; j<m;j++){ out[j]*=w; };
    }

};

#endif





