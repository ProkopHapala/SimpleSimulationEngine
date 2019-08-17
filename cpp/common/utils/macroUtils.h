
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

template <typename T> inline const T& _min  (const T& a, const T& b) { return !(a>b)?a:b; }
template <typename T> inline const T& _max  (const T& a, const T& b) { return !(a<b)?a:b; }
template <typename T> inline const T& _clamp(const T& a, const T& amax, const T& amin){ return _max(amin,_min(amax,a)); }

template <typename T> inline const T& _abs  (const T& a ){ return !(a<0)?a:-a; }
template <typename T> inline int      signum(T val)      { return (T(0) < val) - (val < T(0)); }

// ======= allocation

template<typename T> inline void _allocIfNull(T*& arr, int n){ if(arr==0){ arr=new T[n];} }
template<typename T> inline void _realloc(T*& arr, int n){ if(arr){ delete [] arr;} arr=new T[n]; }
template<typename T> inline void _dealloc(T*& arr       ){ if(arr){ delete [] arr;} arr=0;        }
template<typename T> inline bool _bindOrRealloc(int n, T* from, T*& arr ){ if(from){arr=from;}else{_realloc(arr,n);} }

template<typename T> inline bool _clone( int i0, int imax, T* from, T*& arr, int n){
    _allocIfNull(arr,n);
    for(int i=i0; i<imax; i++){ arr[i]=from[i-i0]; } // use mem copy instead ?
}
template<typename T> inline bool _set( int i0, int imax, const T& from, T*& arr, int n){
    _allocIfNull(arr,n);
    for(int i=i0; i<imax; i++){ arr[i]=from; }
}


#endif
