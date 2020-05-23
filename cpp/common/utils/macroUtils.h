
#ifndef  macroUtils_h
#define  macroUtils_h

#define SWAP( a, b, TYPE ) { TYPE t = a; a = b; b = t; }

//#define _max(a,b)      ((a>b)?a:b)
//#define _min(a,b)      ((a<b)?a:b)
//#define _abs(a)        ((a>0)?a:-a)
//#define _clamp(x,a,b)  max(a, min(b, x))
//#define _clamp(n, lower, upper) if (n < lower) n= lower; else if (n > upper) n= upper
//#define _clamp(a, lower, upper) ((a>lower)?((a<upper)?a:upper):lower)

#define _minit( i, x, imin, xmin )  if( x<xmin ){ xmin=x; imin=i; }
#define _maxit( i, x, imax, xmax )  if( x>xmax ){ xmax=x; imax=i; }

#define _setmin( xmin, x )  if( x<xmin ){ xmin=x; }
#define _setmax( xmax, x )  if( x>xmax ){ xmax=x; }

#define _circ_inc( i, n )   i++; if(i>=n) i=0;
#define _circ_dec( i, n )   i--; if(i< 0) i=n-1;

//#define _realloc(TYPE,arr,n){ if(var) delete [] arr; arr=new TYPE[n]; }

//#define BEGIN_WITH(x) { \
//    auto &_ = x;
//#define END_WITH() }

// ============= sorting

template<typename  T>
struct Buf{
    T*   data=0;
    bool own=false;

    void bind(T* data_){ data=data_; };
    void alloc(int n  ){ own=true; data=new T[n]; };
    //bool realloc(int n){ bool b=false; if(data){ delete [] data; b=true; }; data=new T[n]; return b; };
    //Arr(T* data_){own=true; };
    Buf() = default;
    Buf(int n   ):own{true} {data=new T[n]; };
    Buf(T* data_):own{false},data{data_}{};
    ~Buf(){ if(own&&data)delete[]data; }

    T& operator[](int i){
        //if (index >= size) { printf( "Array index [%i] out of bound [%i] \n", i, n ); exit(0); }
        return data[i];
    }
};

template<typename  T>
struct Arr{
    T*   data;
    int  n;
    bool own;

    //Arr(T* data_){own=true; };
    Arr() = default;
    Arr(int n_           ){ own=true;  n=n_; data=new T[n]; };
    Arr(int n_ , T* data_){ own=false; n=n_; data=data_;    };
    ~Arr(){ if(own&&data)delete[]data; }

    T& operator[](int i){
        //if (index >= size) { printf( "Array index [%i] out of bound [%i] \n", i, n ); exit(0); }
        return data[i];
    }
};


#define WITH(x) auto& _=x;

#define _forN(i,n)         for(int i=0 ;i<n;i++)
#define _for0N(i,i0,n)     for(int i=i0;i<n;i++)
#define _for0Nd(i,i0,n,d)  for(int i=i0;i<n;i+=d)

#define _template_Func   template<typename Func>
#define _template_T      template<typename T>
#define _template_N      template<size_t N>
#define _template_TN     template<typename T,size_t N>
#define _inline_T        template<typename T> inline


_inline_T void _swap (T& a, T& b) { T t=a; a=b; b=t; }

_inline_T const T& _min  (const T& a, const T& b){ return (a>b)?b:a; }
_inline_T const T& _max  (const T& a, const T& b){ return (a<b)?b:a; }
_inline_T const T& _clamp(const T& a, const T& amin, const T& amax){ return _max(amin,_min(amax,a)); }

_inline_T T        _abs  (const T& a ){ return !(a<0)?a:-a; }
_inline_T int      signum(T val)      { return (T(0) < val) - (val < T(0)); }

// ======= allocation

_inline_T bool _allocIfNull(T*& arr, int n){ if(arr==0){ arr=new T[n]; return true; } return false; }
_inline_T void _realloc(T*& arr, int n){ if(arr){ delete [] arr; } arr=new T[n]; }
_inline_T void _dealloc(T*& arr       ){ if(arr){ delete [] arr; } arr=0;        }
_inline_T bool _bindOrRealloc(int n, T* from, T*& arr ){ if(from){ arr=from; return false; }else{ _realloc(arr,n); return true; } }

_inline_T  bool _clone( int i0, int imax, T* from, T*& arr, int n){
    bool ret = _allocIfNull(arr,n);
    for(int i=i0; i<imax; i++){ arr[i]=from[i-i0]; } // use mem copy instead ?
    return ret;
}
_inline_T bool _set( int i0, int imax, const T& from, T*& arr, int n){
    bool ret = _allocIfNull(arr,n);
    for(int i=i0; i<imax; i++){ arr[i]=from; }
    return ret;
}


#endif
