
//#include "math.h"

#include "macros.h"

_inline TYPE METHOD(TPREF, sq( TYPE f ){ return f*f; }   )

_inline TYPE METHOD(TPREF, cos_xx_6 (TYPE xx){ return 1 +xx*(-0.3084251375340424 +xx*(0.01585247243656319 +xx*(-0.0003206691702035272 ))); } )
_inline TYPE METHOD(TPREF, cos_xx_8 (TYPE xx){ return 1 +xx*(-0.3084251375340424 +xx*(0.01585433742711407 +xx*(-0.0003259583464414167 +xx*(3.539800057252972e-06 )))); } )
_inline TYPE METHOD(TPREF, cos_xx_10(TYPE xx){ return 1 +xx*(-0.3084251375340424 +xx*(0.01585434422879243 +xx*(-0.0003259917751184507 +xx*(3.590575024064454e-06 +xx*(-2.430825694459951e-08 ))))); } )

_inline TYPE METHOD(TPREF, sin_xx_6 (TYPE xx){ return 0.7853981633974483 +xx*(-0.08074544390600076 +xx*(0.002490020287907521 +xx*(-3.596176288784333e-05 ))); } )
_inline TYPE METHOD(TPREF, sin_xx_8 (TYPE xx){ return 0.7853981633974483 +xx*(-0.08074551202387524 +xx*(0.002490393182719242 +xx*(-3.657233242102226e-05 +xx*(3.089664734016052e-07 )))); } )
_inline TYPE METHOD(TPREF, sin_xx_10(TYPE xx){ return 0.7853981633974483 +xx*(-0.08074551218802455 +xx*(0.002490394567124187 +xx*(-3.65761914008964e-05  +xx*(3.133376118159403e-07 +xx*(-1.736214417895343e-09 ))))); } )


_inline void  METHOD(TPREF, cos_sin( TYPE x_, TYPE* ca, TYPE* sa ){ )
    const TYPE inv2pi = 1/(3.14159265359/4);
    x_*=inv2pi;
    int ix = (int)x_;
    if(x_<0)ix--;
    TYPE x = x_ - ix;
    if(ix&1)x=x-1;
    TYPE xx = x*x;
    TYPE c,s;
    //c = cos_xx_6 (xx);
    c = METHOD(TPREF,cos_xx_8)(xx);
    //c = cos_xx_10(xx);
    s = x*METHOD(TPREF,sin_xx_6)(xx);
    //s = x*sin_xx_8 (xx);
    //s = x*sin_xx_10(xx);
    ix++;
    if(ix&2){ *ca=-s;*sa=c; }else{ *ca=c; *sa=s; };
    if(ix&4){ *ca=-*ca; *sa=-*sa; }
}


#undef  TYPE
#undef  TPREF