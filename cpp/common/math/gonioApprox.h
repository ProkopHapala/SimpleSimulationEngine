
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


inline sincos( double x, double& ca, double& sa, const int deg ){
    //Taylor Cos [ 1./479001600,  -1./3628800,   1./40320,  -1./720,  1./24,  -1./2, 1. ]
    //Taylor Sin [ 1./6227020800, -1./39916800,  1./362880, -1./5040, 1./120, -1./6, 1. ]
    constexpr const double inv2pi = 1/(2*M_PI);
    constexpr const double CS[]{ 1.0,1.0,  -1./2,-1./6,   1./24,1./120,  -1./720,-1./5040,  1./40320,1./362880,   -1./3628800,-1./39916800,   1./479001600,1./6227020800  };
    double s   = x_*inv2pi;
    double d   = s-(int)s;
    double d2  = d*d;
    double m   = 1-d;
    //double m2  = m*m;
    double d2n=d2;
    //double m2n=m2;
    double sa=0,ca=1;
    for(int i=0; i<deg; i+=2){
        ca += CS[i  ]*d2n;
        sa += CS[i+1]*d2n;
        d2n*=d2;
        //m2n*=m2;
    }
}

// ========= Cartesian -> Polar ===========

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
  double abs_y = fabs(y) + 1e-10f;      // kludge to prevent 0/0 condition
  if ( x < 0 ){
    a     = ( x + abs_y ) / ( abs_y - x );
    angle = 2.35619449019d;
  }else{
    a     = ( x - abs_y ) / ( x + abs_y );
    angle = 0.78539816339d;
  }
  double aa = a * a;
  //angle += a * ( -0.9817d + 0.1963d * aa );
  // angle +=  a * ( -1 + aa*( 0.326388646629 + aa*( -0.155559850719 + aa*0.0437730406925 ) ) );
  angle +=  a * ( -1 + aa*( 0.331768825725 + aa*( -0.184940152398 + aa*( 0.091121250024 -0.0233480867489*aa ) ) ) );
  return  ( y<0 )? -angle : angle ;
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


//////////////////////////////////////////
// ================= acos

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

