
#ifndef  gonioApprox_h
#define  gonioApprox_h

//   speed-accuracy sin/cos approx using Chebyshev
// https://stackoverflow.com/questions/345085/how-do-trigonometric-functions-work/345117#345117
// http://lolengine.net/blog/2011/12/21/better-function-approximations


// ========= Polar -> Cartesian ===========

template <class TYPE>
inline TYPE sin_taylor2( TYPE a ){
	constexpr TYPE c3 = 1.0d/6;
	constexpr TYPE c5 = 1.0d/120;
	TYPE a2 = a*a;
	return    a * ( 1 - a2*( c3 - c5*a2 ) );
}

template <class TYPE>
inline TYPE cos_taylor2( TYPE a ){
	constexpr TYPE c2 = 1.0d/2;
	constexpr TYPE c4 = 1.0d/24;
	TYPE a2 = a*a;
	return    1 - a2*( c2 - c4*a2 );
}

template <class TYPE>
inline void sincos_taylor2( TYPE a, TYPE& sa, TYPE& ca ){
	constexpr TYPE c2 = 1.0d/2;
	constexpr TYPE c3 = 1.0d/6;
	constexpr TYPE c4 = 1.0d/24;
	constexpr TYPE c5 = 1.0d/120;
	TYPE a2 = a*a;
	sa   = a * ( 1 - a2*( c3 - c5*a2 ) ) ;
	ca   =       1 - a2*( c2 - c4*a2 )   ;
}

template <class TYPE>
inline void sincosR2_taylor( TYPE r2, TYPE& sa, TYPE& ca ){
    constexpr TYPE c2 = -1.0d/2;
    constexpr TYPE c3 = -1.0d/6;
    constexpr TYPE c4 =  1.0d/24;
    constexpr TYPE c5 =  1.0d/120;
    constexpr TYPE c6 = -1.0d/720;
    //TYPE r2  = w.x*w.x + w.y*w.y + w.z*w.z;
    sa  =   1 + r2*( c3 + c5*r2 );
    ca  =  c2 + r2*( c4 + c6*r2 );
}

template <class TYPE>
inline void rot_csa( TYPE ca, TYPE sa, TYPE& ux, TYPE& uy ){
	double ux_;
	ux_ = ux*ca - uy*sa;
	uy  = ux*sa + uy*ca;
	ux  = ux_;
}

inline double cos_xx_6 (double xx){ return 1 +xx*(-0.3084251375340424 +xx*(0.01585247243656319 +xx*(-0.0003206691702035272 ))); }
inline double cos_xx_8 (double xx){ return 1 +xx*(-0.3084251375340424 +xx*(0.01585433742711407 +xx*(-0.0003259583464414167 +xx*(3.539800057252972e-06 )))); }
inline double cos_xx_10(double xx){ return 1 +xx*(-0.3084251375340424 +xx*(0.01585434422879243 +xx*(-0.0003259917751184507 +xx*(3.590575024064454e-06 +xx*(-2.430825694459951e-08 ))))); }

inline double sin_xx_6 (double xx){ return 0.7853981633974483 +xx*(-0.08074544390600076 +xx*(0.002490020287907521 +xx*(-3.596176288784333e-05 ))); }
inline double sin_xx_8 (double xx){ return 0.7853981633974483 +xx*(-0.08074551202387524 +xx*(0.002490393182719242 +xx*(-3.657233242102226e-05 +xx*(3.089664734016052e-07 )))); }
inline double sin_xx_10(double xx){ return 0.7853981633974483 +xx*(-0.08074551218802455 +xx*(0.002490394567124187 +xx*(-3.65761914008964e-05  +xx*(3.133376118159403e-07 +xx*(-1.736214417895343e-09 ))))); }

inline void cos_sin( double x_, double& ca, double& sa ){
    constexpr const double inv2pi = 1/(M_PI/4);
    x_*=inv2pi;
    int ix = (int)x_;
    if(x_<0)ix--;
    double x = x_ - ix;
    if(ix&1)x=x-1;
    double xx = x*x;
    double c,s;
    //c = cos_xx_6 (xx);
    c = cos_xx_8 (xx);
    //c = cos_xx_10(xx);
    s = x*sin_xx_6 (xx);
    //s = x*sin_xx_8 (xx);
    //s = x*sin_xx_10(xx);
    ix++;
    if(ix&2){ ca=-s;sa=c; }else{ ca=c; sa=s; };
    if(ix&4){ ca=-ca; sa=-sa; }
}

inline void cos_sin( double x_, double& ca, double& sa, const int order, const int nsplit ){
    constexpr const double inv2pi = 1/(M_PI/2);
    x_*=inv2pi;
    int ix = (int)x_;
    if(x_<0)ix--;
    double x = x_ - ix;
    if(ix&1)x=x-1;
    x/=(1<<nsplit);
    double xx = x*x;
    double CS[]{
        //-0.3084251375340424,     -0.08074551218828077,
         0.01585434424381549,     0.00249039457019272,
        -0.00032599188692739,    -3.657620418217724e-05,
         3.590860448591509e-06,   3.133616890378121e-07,
        -2.461136950494199e-08,  -1.7572476734434e-09,
         1.150115912797405e-10,   6.948453273886625e-12
    };
    double x2n=xx;
    double c=1                 -0.3084251375340424 *x2n;
    double s=0.7853981633974483-0.08074551218828077*x2n;
    for(int i=0; i<order; i+=2){
        x2n*=xx;
        c+=CS[i  ]*x2n;
        s+=CS[i+1]*x2n;
    }
    s*=x;
    rot_csa( c,s, c,s );
    for(int i=0;i<nsplit;i++){
        rot_csa( c,s, c,s );
    }
    if((ix+1)&2){c=-c;s=-s;};
    ca=c;
    sa=s;
}

inline double tan_xx_6 (double xx){ return 1 +xx*(0.3333333333333333 +xx*(0.1257172468409968 +xx*(0.08367259233981537 ))) ;};
inline double tan_xx_8 (double xx){ return 1 +xx*(0.3333333333333333 +xx*(0.1346369234100504 +xx*(0.04494732376211282 +xx*(0.04030214559905823 ))));};
inline double tan_xx_10(double xx){ return 1 +xx*(0.3333333333333333 +xx*(0.1331377015490709 +xx*(0.05602812424047391 +xx*(0.0144884610558432  +xx*(0.01919844564346994 )))));};
inline double tan_xx_12(double xx){ return 1 +xx*(0.3333333333333333 +xx*(0.1333599255952193 +xx*(0.05357641184260721 +xx*(0.02397340014368198 +xx*(0.003720346570791119 +xx*0.009082932867935635 )))));};

//inline double tan_xx_12(double xx){ return 1 +xx*(0.3333333333333333); };

inline double tan_xx_12_(double xx){
    double x4=xx*xx;
    return 1 + xx*( 0.3333333333333333 + x4*( 0.05357641184260721 +x4*(0.003720346570791119 )))
             + x4*( 0.1333599255952193 + x4*( 0.02397340014368198 +x4*(0.009082932867935635 )));
}

// ========= Cartesian -> Polar ===========

inline double atan2_xx_6 (double xx){ return -1 +xx*(0.3271414224903835 +xx*(-0.1584645896620068 +xx*(0.04626030273606059 ) ) ); };
//inline double atan2_xx_8 (double xx){ return -1 +xx*(0.331768825725     +xx*(-0.184940152398     +xx*(0.091121250024      +xx*(-0.0233480867489     ) ) ) ); };
inline double atan2_xx_8 (double xx){ return -1 +xx*(0.3319089228367821 +xx*(-0.1858750170176101 +xx*(0.09293719121053462 +xx*(-0.02441667691946986 ) ) ) ); };
inline double atan2_xx_10(double xx){ return -1 +xx*(0.3330176488472483 +xx*(-0.1956751542049746 +xx*(0.12132313880689250 +xx*(-0.05764138590306975 +xx*(0.01358448668618405 ) ) ) ) ); };
inline double atan2_xx_12(double xx){ return -1 +xx*(0.3332653548010540 +xx*(-0.1987725193225652 +xx*(0.13468905116342580 +xx*(-0.08358586640953759 +xx*(0.03683379230278266 +xx*(-0.007828996358246611 ) ) ) ) ) ); };
inline double atan2_xx_14(double xx){ return -1 +xx*(0.3333190269289065 +xx*(-0.1996711468724480 +xx*(0.14004322605623570 +xx*(-0.09876163152183830 +xx*(0.05901114389963034 +xx*(-0.02396583101041806 +xx*(0.004627201618174669 ) ) ) ) ) ) ); }; 

inline double atan_poly( double a ){
  return a * ( 1 + a * ( 0.00317436203193 + a * ( -0.362955638193 + a * ( 0.0920711110177 + a*( 0.105247322934 + a*-0.0521389943934 ) ) ) ) );
  /*
  final double c3 = -0.365791976855;
  final double c4 =  0.151190140253;
  return a * ( 1 + a * a * ( c3 + a * c4 ) );
  */
}

inline double atan2_a1( double y, double x ){
  //http://pubs.opengroup.org/onlinepubs/009695399/functions/atan2.html
  //Volkan SALMA
  double a, angle;
  double abs_y = fabs(y) + 1e-14;      // kludge to prevent 0/0 condition
  //double x_ = x + abs_y;
  //double y_ = x - abs_y;
  if ( x < 0 ){
    a     = ( x + abs_y ) / ( abs_y - x );
    //a     = -x_/y_;
    angle = 2.35619449019;
  }else{
    a     = ( x - abs_y ) / ( abs_y + x );
    //a     = y_/x_;
    angle = 0.78539816339;
  }
  double aa = a * a;
  //angle += a * ( -0.9817d + 0.1963d * aa );
  //angle +=  a * ( -1 + aa*( 0.326388646629 + aa*( -0.155559850719 + aa*0.0437730406925 ) ) );
  angle +=  a * ( -1 + aa*( 0.331768825725 + aa*( -0.184940152398 + aa*( 0.091121250024 -0.0233480867489*aa ) ) ) );
  return  ( y<0 )? -angle : angle ;
}

template<double (*xx_poly)(double)>
inline double atan2_t( double y, double x ){
  //http://pubs.opengroup.org/onlinepubs/009695399/functions/atan2.html
  //Volkan SALMA
  double a, angle;
  double abs_y = fabs(y) + 1e-14;      // kludge to prevent 0/0 condition
  if ( x < 0 ){
    a     = ( x + abs_y ) / ( abs_y - x );
    angle = 2.35619449019;
  }else{
    a     = ( x - abs_y ) / ( abs_y + x );
    angle = 0.78539816339;
  }
  angle += a * xx_poly(a*a);
  return  ( y<0 )? -angle : angle;
}

inline double atan2_a2( double y, double x ){
  double absx,absy;
  uint8_t kind=0;
  double a;
  if(    x > 0    ){ absx = x;   kind |= 2; }else{ absx=-x;  }
  if(    y > 0    ){ absy = y;   kind |= 4; }else{ absy=-y;  }
  if( absx > absy ){ a = absy / absx;       }else{ a = absx / absy; kind |= 1; }
  //double alfa = atan_poly( a );
  double alfa = a * ( 1 + a * ( 0.00317436203193 + a * ( -0.362955638193 + a * ( 0.0920711110177 + a*( 0.105247322934 + a*-0.0521389943934 ) ) ) ) );
  //double aa = a*a;
  //double alfa =  a * ( -1 + aa*( 0.331768825725 + aa*( -0.184940152398 + aa*( 0.091121250024 -0.0233480867489*aa ) ) ) );
  switch( kind ){
    case 0:  return -3.14159265359 + alfa;
    case 1:  return -1.57079632679 - alfa;
    case 2:  return                - alfa;
    case 3:  return -1.57079632679 + alfa;
    case 4:  return  3.14159265359 - alfa;
    case 5:  return  1.57079632679 + alfa;
    case 6:  return                  alfa;
    case 7:  return +1.57079632679 - alfa;
  }
  return 0;
}

inline double atan2_a3( double y, double x ){
  if( x > 0 ){
    if( y > 0 ){
      if( x > y ){
        double a = y / x;
        return atan_poly( a );
      }else{
        double a = x / y;
        return 1.57079632679 - atan_poly( a );
      }
    }else{
      if( x > -y ){
        double a = y / x;
        return -atan_poly( -a );
      }else{
        double a = x / y;
        return -1.57079632679 + atan_poly( -a );
      }
    }
  }else{
    if( y > 0 ){
      if( -x > y ){
        double a = y / x;
        return 3.14159265359 - atan_poly( -a );
      }else{
        double a = x / y;
        return 1.57079632679 + atan_poly( -a );
      }
    }else{
      if(  x < y ){
        double a = y / x;
        return -3.14159265359 + atan_poly( a );
      }else{
        double a = x / y;
        return -1.57079632679 - atan_poly( a );
      }
    }
  }
}


// from http://http.developer.nvidia.com/Cg/atan2.html

inline float atan2_nvidia( float y, float x ){
  float t0, t1, t2, t3, t4;

  t3 = _abs(x);
  t1 = _abs(y);
  t0 = _max(t3, t1);
  t1 = _min(t3, t1);
  t3 = 1 / t0;
  t3 = t1 * t3;

  t4 = t3 * t3;
  t0 =         - 0.013480470;
  t0 = t0 * t4 + 0.057477314;
  t0 = t0 * t4 - 0.121239071;
  t0 = t0 * t4 + 0.195635925;
  t0 = t0 * t4 - 0.332994597;
  t0 = t0 * t4 + 0.999995630;
  t3 = t0 * t3;

  t3 = ( _abs(y) > _abs(x) ) ? 1.570796327 - t3 : t3;
  t3 = (x < 0) ?  3.141592654 - t3 : t3;
  t3 = (y < 0) ? -t3 : t3;

  return t3;
}

inline float atan2_nvidia_fabs( float y, float x ){
  float t0, t1, t2, t3, t4;

  t3 = fabs(x);
  t1 = fabs(y);
  t0 = fmax(t3, t1);
  t1 = fmin(t3, t1);
  t3 = 1 / t0;
  t3 = t1 * t3;

  t4 = t3 * t3;
  t0 =         - 0.013480470;
  t0 = t0 * t4 + 0.057477314;
  t0 = t0 * t4 - 0.121239071;
  t0 = t0 * t4 + 0.195635925;
  t0 = t0 * t4 - 0.332994597;
  t0 = t0 * t4 + 0.999995630;
  t3 = t0 * t3;

  t3 = ( fabs(y) > fabs(x) ) ? 1.570796327 - t3 : t3;
  t3 = (x < 0) ?  3.141592654 - t3 : t3;
  t3 = (y < 0) ? -t3 : t3;

  return t3;
}


//////////////////////////////////////////
// ================= acos

    // acosf implementation:
    // https://github.com/bminor/glibc/blob/master/sysdeps/ieee754/flt-32/e_acosf.c

/*
inline float acos_logTable(float x){
    // 4194303   - mantisa mask
    int32_t i=*(int32_t*)&x;
    int32_t mantisa = i&4194303;
    int32_t exp     = i>>22;
    //table[exp];
}
*/

//dx = 0.059375000000
// generated by /home/prokop/Dropbox/MyDevSW/Python/_math/FunctionApproximation/generate_Tables.py
static const float TABLE_ACOS[66] = {2.824032224298, -0.190152182644, 2.669514068113, -0.130569545221, 2.552148856188, -0.106809015748, 2.452581381509, -0.093390057210, 2.363851039421, -0.084619296037, 2.282500230795, -0.078408703997, 2.206507612196, -0.073789885551, 2.134564472568, -0.070245696130, 2.065760358512, -0.067472680030, 1.999426185482, -0.065280533753, 1.935047823674, -0.063544084347, 1.872214770557, -0.062178221199, 1.810587726635, -0.061123905388, 1.749876975579, -0.060339964677, 1.689827206273, -0.059798120427, 1.630206268935, -0.059479937575, 1.570796326795, -0.059375000000, 1.511386384655, -0.059479937575, 1.451765447316, -0.059798120427, 1.391715678010, -0.060339964677, 1.331004926955, -0.061123905388, 1.269377883033, -0.062178221199, 1.206544829916, -0.063544084347, 1.142166468108, -0.065280533753, 1.075832295078, -0.067472680030, 1.007028181021, -0.070245696130, 0.935085041394, -0.073789885551, 0.859092422794, -0.078408703997, 0.777741614169, -0.084619296037, 0.689011272081, -0.093390057210, 0.589443797402, -0.106809015748, 0.472078585477, -0.130569545221, 0.317560429292, -0.190152182644};

inline float acos_table( float x ){
    constexpr float dx    = 0.059375f;
    constexpr float invdx = 1/dx;
    int sgn=1; float u = x*invdx;
    if(u<0){ sgn=-1; u=-u; }
    int   i  = (int)u;
    float du  = u - i;
    if(i<20){
        i<<1;
        float y0  = TABLE_ACOS[i  ];
        float y1  = TABLE_ACOS[i+1];
        float dy0 = TABLE_ACOS[i+2];
        float dy1 = TABLE_ACOS[i+3];
        float y01 = y0-y1;
        return  y0
          +du*(           dy0
          +du*( -3*y01 -2*dy0 - dy1
          +du*(  2*y01 +  dy0 + dy1 )));
    }else{
        return sqrt(1-du)*sgn;
    }
}


// acos() approximation by Sebastien Lagarde
// https://www.shadertoy.com/view/lsjXDc
//float sacos( float y ){
//    float x = abs( clamp(y,-1.0,1.0) );
//    float z = (-0.168577*x + 1.56723) * sqrt(1.0 - x);
//    return mix( 0.5*3.1415927, z, sign(y) );
//}

// http://http.developer.nvidia.com/Cg/acos.html
inline float acos_nVidia(float x) {
  float negate = float(x < 0);
  x = abs(x);
  float ret = -0.0187293;
  ret = ret * x;
  ret = ret + 0.0742610;
  ret = ret * x;
  ret = ret - 0.2121144;
  ret = ret * x;
  ret = ret + 1.5707288;
  ret = ret * sqrt(1.0-x);
  ret = ret - 2 * negate * ret;
  return negate * 3.14159265358979 + ret;
}

// sincos
// from here https://www.shadertoy.com/view/XsjXzt
/*
Vec2 sinCos_iq( in float x )
{
    vec2  f = abs(fract(x-vec2(0.25,0.0))-0.5);
    float h = abs(fract(x*4.0)-0.5);

    //     approx sin/cos             * approx renormalization (sqrt(2)Â·8/11)
    return (-1.0 + f*f*(24.0-32.0*f)) * (1.028519 - 0.0570379*h);
}
*/

#endif

