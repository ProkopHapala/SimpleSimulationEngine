
#ifndef  NumT_h
#define  NumT_h

#include <vector>

#include <functional>
//  --- template magig?
#include <type_traits>
#include <utility>

#include "macroUtils.h"
#include "fastmath.h"

namespace NumT{  // numerical Templates / Types

//_template_T

//template<typename T> void println(const T& t){t.print(); puts(""); };
_template_T void println(const T& t){t.print(); puts(""); };

//template<typename T> void  _dot(T ai,T bi,T sum){ sum+=ai*bi; };

_inline_T void _maximum(T ai,      T& amax){ amax=((ai>amax)?ai:amax); };
_inline_T void _minimum(T ai,      T& amin){ amin=((ai<amin)?ai:amin); };
_inline_T void _norm1  (T ai,      T& sum ){ sum+=((ai>=0)?ai:-ai); };
_inline_T void _norm2  (T ai,      T& sum ){ sum+=ai*ai; };
_inline_T void _dot    (T ai,T bi, T& sum ){ sum+=ai*bi; };

void print(double f){ printf("%g ",f); };

// =========================================
// =============== Tuple
// =========================================

//template <typename T,size_t N>
_template_TN
struct Tuple{
    constexpr static const int nDim = N;
    T items[N];

    //Tuple() = default;

    T operator[](size_t i){
        return items[i];
        //if(i<N){}else{};
    }

    void print()const{
        printf("(");
        _forN(i,N){ items[i].print(); };
        printf(")");
    }
};

/*
//https://stackoverflow.com/questions/57540299/initialize-array-inside-variadic-template-class/
//template<typename T, T...args>
//template<typename T, Args...args>
template<typename T, size_t N>
struct Tuple{
    //T data[sizeof...(args)];
    T data[N];
    //T data{args...};

    //template <typename ...Args>
    //Tuple_(T... args):data{args...}{};
    template<typename ...Args>
    Tuple(Args && ...args):data{std::forward<Args>(args)...}{};
    //Tuple_(const T& args...):data{args...}{};
    //Tuple_(T* ...args):data{args...}{};
};


template<typename ...Args>
Tuple(Args && ...args)
-> Tuple<std::common_type_t<std::remove_reference_t<Args>...>,sizeof...(args)>::Tuple(Args && ...args);
*/


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
    Scan(const Range& r_        ):Range(r_){};
    Scan(const Range& r_,int by_):Range(r_),by(by_){};
    Scan(int from_,int to_        ):Range(from_,to_),by(1){};
    Scan(int from_,int to_,int by_):Range(from_,to_),by(by_){};
    void print()const{printf("[%i:%i::%i]\n", from, to, by );}

    Scan operator-(){ to=-to; return *this; }
};

// =========================================
// =============== ScanN
// =========================================

//template<int NDIM>
_template_N
struct ScanN : public Tuple<Scan,N> {
    ScanN()=default;

    void print()const{
        printf("[ ");
        _forN(i,N){
            const Scan& s = (*this).items[i];
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

_template_N
inline ScanN<N+1> operator ,(ScanN<N> sn,const Scan& s){
    ScanN<N+1> sn_;
    _forN(i,N){ sn_.items[i]=sn.items[i]; };
    sn_.items[N]=s;
    return sn_;
}

_template_N
inline ScanN<N+1> operator ,(const Scan& s,ScanN<N> sn){
    ScanN<N+1> sn_;
    _forN(i,N){ sn_.items[i+1]=sn.items[i]; };
    sn_.items[0]=s;
    return sn_;
}

// =========================================
// =============== VecN
// =========================================

_template_T  using Func1 = std::function<T(const T&                  )>;
_template_T  using Func2 = std::function<T(const T&,const T&         )>;
_template_T  using Func3 = std::function<T(const T&,const T&,const T&)>;

_template_T
struct VecN{
    int n;
    T*  data;

    void copy(int i0,int i1,T*       from){ T* to=data+i0; int n=i1-i0;
        _forN(i,n){
            to[i]=from[i];
        };
    }
    void set (int i0,int i1,const T& from){ T* to=data+i0; int n=i1-i0; _forN(i,n){ to[i]=from; }; }

    VecN()=default;

    VecN(int n_){ n=n_; _realloc(data,n); };
    VecN(int n_,T t){ n=n_; _realloc(data,n); set(0,n,t); };
    VecN(int n_,T* ts){ n=n_; _realloc(data,n); copy(0,n,ts); };
    //VecN(const VecN& v){ n=v.n; _realloc(data,n); copy(0,n,v.data); };  // THIS IS IN DEFAULT CONSTRUCTOR


    inline T& operator[](size_t i){ return data[i]; };


    //https://tristanbrindle.com/posts/beware-copies-initializer-list

    /*
    template <typename T, typename... Args>
    std::vector<T> make_vector(Args&&... args){
        std::vector<T> vec;
        vec.reserve(sizeof...(Args));
        (vec.push_back(std::forward<Args>(args)), ...);
        return vec;
    }

    /*
    template <typename... Args>
    VecN(Args&&... args){
        n = (sizeof...(Args));
        //_realloc( data, n );
        //T* data_ = data;
        //( *(data++)=(std::forward<Args>(args)), ...);

        data = T[]{ ... };

        //using arr_t = int[];
        //(void) arr_t{0, (vec.push_back(std::forward<Args>(args)), 0)...};
        //(void) arr_t{0, ( *(data++)=(std::forward<Args>(args)),    , 0)...};
    }
    */

    VecN(const std::initializer_list<T>& ts):
    VecN(ts.size()){
        T* to=data; for(const T& t: ts){ *to=t; to++; };
    }

    void print()const{
        printf( "{ " );
        for(int i=0;i<n;i++){
            NumT::print(data[i]);
        };
        printf( " }" );
    }
    //void println()const{ print(); print("\n"); }

    using Func1T = Func1<T>;
    using Func2T = Func2<T>;
    using Func3T = Func3<T>;
    //using Func3T = Func3<T>;

    void scan(int i0, int imax, int by, const Func2T& func ){
        T val;
        for(int i=i0;i<imax;i+=by){ func(data[i],val); };
        return val;
    }
    void scan(const Func2T& func ){ return scan(0,n,1,func ); };

    void scan2(int i0, int imax, int by,const VecN<T>& v,const Func2T& func){
        for(int i=i0;i<imax;i+=by){ data[i]=func(data[i],v.data[i]); };
    }
    void scan2(const VecN<T>& v,const Func2T& func ){
        int n_ = _min(n,v.n);
        return scan2(0,n_,1,v,func);
    };



    /*
    T scan3(int i0, int imax, int by, VecN<T>& v1, VecN<T>& v2, Func2T& func){
        T val;
        for(int i=i0;i<imax;i+=by){ func(data[i],val); };
        return val;
    }
    T scan3(VecN<T>& v1, VecN<T>& v2,const Func3T& func ){
        int n_ = _min(n,_min(v1.n,v2.n));
        return scan3(0,n_,1,v1,v2,func);
    };
    */

    /*
    template<typename Func>
    T apply(int i0, int imax, int by, Func func){
        T val;
        for(int i=i0;i<imax;i+=by){ val=func(data[i],val); };
        return val;
    }
    template<typename Func> T apply(int i0,int imax,Func func){ return apply(i0,imax,1,func); };
    template<typename Func> T apply(Func func){ return apply(0,n,1,func); };
    template<typename Func> T apply(Scan  s,Func func){ return apply(s.from,s.to,s.by,func); };
    template<typename Func> T apply(Range r,Func func){ return apply(r.from,r.to,1   ,func); };
    */

    _template_Func
    T apply(int i0, int imax, int by, Func func, VecN<T>& vec ){
        T val;
        for(int i=i0;i<imax;i+=by){ val=func(data[i],vec.data[i],val); };
        return val;
    }
    _template_Func T apply(int i0,int imax,Func func, VecN<T>& vec){ return apply(i0,imax,1,func,vec); };
    _template_Func T apply(Func func,                 VecN<T>& vec){ return apply(0,n,1,func,vec); };
    _template_Func T apply(Scan  s,Func func,         VecN<T>& vec){ return apply(s.from,s.to,s.by,func,vec); };
    _template_Func T apply(Range r,Func func,         VecN<T>& vec){ return apply(r.from,r.to,1   ,func,vec); };
    //void op(  ,[]){ };
    //apply;
};

_template_T
struct View: Scan{
    VecN<T>* vec;

    template<typename Func>
    T apply(Func func){
        T val;
        for(int i=from;i<to;i+=by){ val=func(vec.data[i],val); }
        return val;
    }
};

_template_T
class View2{
    int n;
    View<T> view1;
    View<T> view2;

    //template<typename Func>
    //T apply(Func func){
    T apply(Func2<T> func){
        T val;
        T* data1 = view1.vec.data + view1.from;
        T* data2 = view2.vec.data + view2.from;
        for(int i=0;i<n;i++){
            val=func(data1,data2);
            data1+=view1.by;
            data1+=view2.by;
        }
        return val;
    }

};

//template<typename T,int NDIM>
_template_TN
struct Tensor{
    constexpr static const int nDim = N;
    int ns[];
    T* data;

    bool isContinuous(){ return data!=0; }



};

using VecNd = VecN<double>;
using VecNf = VecN<float>;
using VecNi = VecN<int>;

}; // namespace NDarray



#endif





