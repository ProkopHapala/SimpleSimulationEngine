
#ifndef  NDarray_h
#define  NDarray_h

#include <vector>;


#include <functional>;


#include "fastmath.h"

namespace NDarray{

template<typename T> void o_max  (T ai,      T amax){ amax=((ai>amax)?ai:amax); };
template<typename T> void o_min  (T ai,      T amin){ amin=((ai<amin)?ai:amin); };
//template T _abs  (T ai,      T amax){ amax=ai>amax?ai:amax; };
//template T _normMax(T ai,      T amax){ sum =((ai>0)ai?ai:amax; };
template<typename T> void  o_norm1  (T ai,      T sum ){ sum+=((ai>=0)?ai:-ai); };
template<typename T> void  o_norm2  (T ai,      T sum ){ sum+=ai*ai; };
template<typename T> void  o_dot    (T ai,T bi, T sum ){ sum+=ai*bi; };

void print(double f){ printf("%g ",f); };

// =========================================
// =============== Tuple
// =========================================

template <typename T,int N>
struct Tuple{
    constexpr static const int nDim = N;
    T items[N];

    Tuple() = default;

    void print()const{
        printf("(");
        for(int i=0;i<N;i++){ items[i].print(); };
        printf(")");
    }
};

// =========================================
// =============== Range
// =========================================

struct Range{
    int from=0,to=0;
    Range() = default;
    Range(          int to_):from(0    ),to(to_){};
    Range(int from_,int to_):from(from_),to(to_){};
    void print()const{printf("[%i:%i]\n", from, to );}

    const Range& operator-(){ to=-to; return *this; }

};

// =========================================
// =============== Scan
// =========================================

struct Scan : Range {
    //Range range;
    int   by=1;
    Scan() = default;
    //operator Scan(int i) = delete;
    //Scan(){by=};
    Scan(const Range& r_        ):Range(r_){};
    Scan(const Range& r_,int by_):Range(r_),by(by_){};
    //Scan(          int to_        ):Range(to_      ),by(1){};
    Scan(int from_,int to_        ):Range(from_,to_),by(1){};
    Scan(int from_,int to_,int by_):Range(from_,to_),by(by_){};
    //Scan(int from_,int to_,int by_)::range(form_,to_){by=by_;};
    void print()const{printf("[%i:%i::%i]\n", from, to, by );}

    Scan operator-(){ to=-to; return *this; }
};

// =========================================
// =============== ScanN
// =========================================

template<int NDIM>
struct ScanN : public Tuple<Scan,NDIM> {
    //constexpr static const int nDim = NDIM;
    //static const int nDim = NDIM;
    ScanN()=default;
    //ScanN(const ScanN<NDIM-1>& o,Scan s){ };
    //ScanN<NDIM>& operator |(const Scan& s){ return ScanN<T,NDIM+1>{*this,s}; }

    void print()const{
        //printf("print nDim %i \n",nDim);
        //printf("print scans[0].to %i \n",scans[0].to);
        /*
        printf("[");
        for(int i=0;i<nDim;i++){
            if(i>0)printf("|");
            printf("%i:%i::%i", scans[i].from, scans[i].to, scans[i].by );
        }
        printf("]\n");
        */
        //printf("NDIM %i \n", NDIM);
        printf("[ ");
        for(int i=0;i<NDIM;i++){
            const Scan& s = (*this).items[i];
            //const Scan& s = ScanN<2>::items[i];
            //const Scan& s = Tuple<Scan,NDIM>::items[i];
            //const Scan& s = items[i];
            if(i>0)printf(" | ");
            printf("%i:%i::i", s.from, s.to, s.by );
        }
        printf(" ]\n");
    }
};

// https://stackoverflow.com/questions/23430296/how-do-i-create-a-user-defined-literal-for-a-signed-integer-type
Range operator"" _r (unsigned long long int to){ return Range((int)to); }
Scan  operator"" _s (unsigned long long int to){ return Scan (0,(int)to,1); }

inline Range operator+(const Range& r1,const Range& r2 ){ return Range(_min(r1.from,r2.from),_max(r1.to,r2.to)); }
inline Range operator-(const Range& r1,const Range& r2 ){ return Range(_max(r1.from,r2.to  ),_min(r1.from,r2.from)); }
inline Range operator*(const Range& r1,const Range& r2 ){ return Range(_max(r1.from,r2.from),_min(r1.from,r2.to));   }

inline Range operator<(int from, Range r){ r.from=from; return r; }
inline Scan  operator<(int from, Scan  s){ s.from=from; return s; }
inline Scan  operator/(const Range& r,int by){ return Scan{r,by}; }
inline Scan  operator/(int to,const Scan&  s){ return Scan{0,to,s.to}; }

inline ScanN<2> operator ,(const Scan& s0,const Scan& s1){
    //printf("  op Scan|Scan \n" );
    //return ScanN<2>{s1,s2};
    ScanN<2> s{};
    s.items[0]=s0;  //s0.print();
    s.items[1]=s1;  //s1.print();
    //s.print();
    return s;
    //return Tuple<Scan,2>{s1,s2};
}

template<int N>
inline ScanN<N+1> operator ,(ScanN<N> sn,const Scan& s){
    //printf("  ScanN<%i>|Scan \n", N );
    ScanN<N+1> sn_;
    for(int i=0;i<N;i++){ sn_.items[i]=sn.items[i]; };
    sn_.items[N]=s;
    //sn_.print();
    return sn_;
}

template<int N>
inline ScanN<N+1> operator ,(const Scan& s,ScanN<N> sn){
    //printf("  Scan|ScanN<%i>\n", N );
    ScanN<N+1> sn_;
    for(int i=0;i<N;i++){ sn_.items[i+1]=sn.items[i]; };
    sn_.items[0]=s;
    //sn_.print();
    return sn_;
}



// =========================================
// =============== VecN
// =========================================


template<typename T>
struct VecN{
    int n;
    T*  buff;

    void copy(int i0,int i1,T* from){ T* to=buff+i0; int n=i1-i0; for(int i=0;i<n;i++){ to[i]=from[i]; }; }

    VecN()=default;
    VecN(int n_){ n=n_; _realloc(buff,n); };
    VecN(std::initializer_list<T> ts):VecN(ts.size()){
        T* to=buff; for(const T& t: ts){ *to=t; to++; };
    }

    void print(){
        printf( "{ " );
        for(int i=0;i<n;i++){
            NDarray::print(buff[i]);
        };
        printf( " }" );
    }

    T apply_(int i0, int imax, int by, const std::function<void(T&,T&)>& func ){
        T val;
        for(int i=i0;i<imax;i+=by){ func(buff[i],val); };
        return val;
    }
    T apply_(const std::function<void(T&,T&)>& func ){ return apply_(0,n,1,func ); };

    template<typename Func>
    T apply(int i0, int imax, int by, Func func){
        T val;
        for(int i=i0;i<imax;i+=by){ val=func(buff[i],val); };
        return val;
    }
    template<typename Func> T apply(int i0,int imax,Func func){ return apply(i0,imax,1,func); };
    template<typename Func> T apply(Func func){ return apply(0,n,1,func); };
    template<typename Func> T apply(Scan  s,Func func){ return apply(s.from,s.to,s.by,func); };
    template<typename Func> T apply(Range r,Func func){ return apply(r.from,r.to,1   ,func); };

    template<typename Func>
    T apply(int i0, int imax, int by, Func func, VecN<T>& vec ){
        T val;
        for(int i=i0;i<imax;i+=by){ val=func(buff[i],vec.buff[i],val); };
        return val;
    }
    template<typename Func> T apply(int i0,int imax,Func func, VecN<T>& vec){ return apply(i0,imax,1,func,vec); };
    template<typename Func> T apply(Func func,                 VecN<T>& vec){ return apply(0,n,1,func,vec); };
    template<typename Func> T apply(Scan  s,Func func,         VecN<T>& vec){ return apply(s.from,s.to,s.by,func,vec); };
    template<typename Func> T apply(Range r,Func func,         VecN<T>& vec){ return apply(r.from,r.to,1   ,func,vec); };
    //void op(  ,[]){ };
    //apply;
};

template<typename T>
struct View: Scan{
    VecN<T>* vec;

    template<typename Func>
    T apply(Func func){
        T val;
        for(int i=from;i<to;i+=by){ val=func(vec.buff[i],val); }
        return val;
    }
};

template<typename T>
class View2{
    int n;
    View<T> view1;
    View<T> view2;

    template<typename Func>
    T apply(Func func){
        T val;
        T* buff1 = view1.vec.buff + view1.from;
        T* buff2 = view2.vec.buff + view2.from;
        for(int i=0;i<n;i++){
            val=func(buff1,buff2,val);
            buff1+=view1.by;
            buff1+=view2.by;
        }
        return val;
    }

};

template<typename T,typename Func>
T apply(VecN<T> v, Func func){
    T val;
    for(int i=0;i<v.n;i++){ val=func(v.buff[i],val); };
    return val;
}




template<typename T,int NDIM>
struct Tensor{
    constexpr static const int nDim = NDIM;
    int ns[];
    T* buff;

    bool isContinuous(){ return buff!=0; }



};




typedef VecN<double>  VecNd;






}; // namespace NDarray



#endif





