
#ifndef  PotentialFlow_h
#define  PotentialFlow_h


/*

References:
- short overview lif-indced-drag   http://aero.stanford.edu/reports/multop/multop.html

TODO:
 - Eliptical span-wise distribution of lift

 Integral( I(x)*(R x dl)/|R|^2 ) = x*sqrt(L^2-x^2)/(x^2 + a^2)
 = sqrt(L-x) - arcTgh(  )
 https://www.wolframalpha.com/input/?i=integral+sqrt(L%5E2-x%5E2)%2F(x%5E2%2Ba%5E2)+by+x

- solution of this eliptical integral is well described here ( 12.6
ELLIPTIC LIFT DISTRIBUTION page 18):
https://web.stanford.edu/~cantwell/AA200_Course_Material/AA200_Course_Notes/AA200_Ch_12_Wings_of_Finite_Span_Cantwell.pdf


Integration of linear lift distribution is quite simple
integral (x^2) / (x^2 + L^2) =  x - L/arctg(x/sqrt(L))


integral (x^3) / (x^2 + L^2)  = 0.5*(x^2 - L^2 log( L^2+x^2 ) )
https://www.wolframalpha.com/input/?i=integral+(x%5E3)+%2F+(x%5E2+%2B+L)

Any function can be inegrated when approximated by some polynominal expansion series (e.g. taylor expansion)

*/

inline Vec3d pointSource( Vec3d R ){
    double r2 = R.norm2();
    return R*(1/(r2*sqrt(r2)));
}

inline Vec3d sourceDipol( const Vec3d& R, const Quat4d& coefs ){
    // https://en.wikipedia.org/wiki/Dipole#Field_from_an_electric_dipole
    // F_dipole = e/(4pi*e0) ( 3<p|rhat>rhat-p )/|r|^3 - delta(r)
    // F_dipole = 3<p|r>r/|r|^5 - p/r^3
    double ir2 = 1/R.norm2();
    double ir3 = ir2*sqrt(ir2);
    //return R*(( ((Vec3d*)&coefs)->dot(R)*ir2  + coefs.w )*ir3) + p*ir3;
    return R*(( coefs.f.dot(R)*ir2  + coefs.e )*ir3) + coefs.f*ir3;
    //return R*(( coefs.f.dot(R)*ir2  + coefs.e )*ir3);
}

// http://s6.aeromech.usyd.edu.au/aerodynamics/index.php/sample-page/subsonic-aerofoil-and-wing-theory/3d-vortex-lattice-method/
// http://web.mit.edu/16.unified/www/SPRING/fluids/Spring2008/LectureNotes/f06.pdf

constexpr double VortexQuantum = 1.0/( 4*M_PI );

inline Vec3d dBiotSawart( Vec3d R, Vec3d dI ){
    // https://en.wikipedia.org/wiki/Biot%E2%80%93Savart_law
    // http://web.mit.edu/viz/EM/visualizations/coursenotes/modules/guide09.pdf
    //( Idl x r)/|r^3|  = ( dI x r )/|r^3|
    Vec3d dB;
    dB.set_cross( dI, R ); //printf( "(%f,%f,%f) (%f,%f,%f) (%f,%f,%f)\n", R.x,R.y,R.z, dI.x,dI.y,dI.z, dB.x,dB.y,dB.z );
    double r = R.norm();
    dB.mul( VortexQuantum/(r*r*r) );
    return dB;
}

inline Vec3d ILineFinite( Vec3d R, Vec3d hL, double l ){
    // https://en.wikipedia.org/wiki/Biot%E2%80%93Savart_law
    // http://web.mit.edu/viz/EM/visualizations/coursenotes/modules/guide09.pdf
    /*
    Vec3d dB;
    dB.set_cross( hL, R );    // actually sin(theta_1)*|R| ... can ve remove one sqrt() ?
    dB.normalize();
    double r1 = R.norm();
    double c1 = hL.dot( R )/r1;       // cos(Theta_1)
    double a  = r1*sqrt( 1-c1*c1 );   // minimum point-line distance
    R.add_mul( hL, l );
    double c2 = hL.dot( R )/R.norm(); // cos(Theta_2)
    dB.mul( (c2-c1) * VortexQuantum/a );
    */
    // Optimized:    since a = |hL x R|     B/a = (hL x R)/|hL x R|^2
    Vec3d B;
    B.set_cross( hL, R );
    double c1 = hL.dot( R )/R.norm();
    R.add_mul( hL, l );
    double c2 = hL.dot( R )/R.norm();
    //printf( "c1 %f c2 %f \n", c1, c2 );
    B.mul( (c2-c1)*VortexQuantum/B.norm2() ); // notice there is no "a"
    return B;
}

inline Vec3d ILineSemiInf( Vec3d R, Vec3d hL ){
    Vec3d B;
    B.set_cross( hL, R );
    double c1 = hL.dot( R )/R.norm();
    //printf( "(%f,%f,%f)  %f %f %f \n", B.x, B.y, B.z, R.norm(), c1, B.norm2() );
    //double out = (1-c1)/B.norm2();
    //printf( "%f %f   %f %f \n", R.x, R.y, c1, out );
    //B.mul( VortexQuantum*out );
    B.mul( VortexQuantum*(1-c1)/B.norm2() ); // notice there is no "a"
    return B;
}

inline Vec3d ILineSemiInfDecay( Vec3d R, Vec3d hL, double w2 ){
    Vec3d B;
    B.set_cross( hL, R );
    double cL = hL.dot( R );
    double c1 = cL/R.norm();
    B.mul( (1-c1)*VortexQuantum/( B.norm2() + cL*cL*w2 ) ); // Lorenzian Decay of vortex core
    return B;
}

inline void horseshoe( Vec3d& B, Vec3d R, Vec3d p0, Vec3d p1, Vec3d hDir, double strength ){
    B.add_mul( ILineSemiInf( R-p0, hDir ), -strength );
    B.add_mul( ILineSemiInf( R-p1, hDir ),  strength );
    hDir.set_sub(p0,p1);
    double l = hDir.normalize();
    B.add_mul( ILineFinite( R-p0, hDir, l ), strength );
    //printf( "R(%f,%f,%f) B(%f,%f,%f)\n", R.x, R.y, R.z,  B.x, B.y, B.z );
}

inline void horseshoeDecay( Vec3d& B, Vec3d R, Vec3d p0, Vec3d p1, Vec3d hDir, double strength, double w2 ){
    B.add_mul( ILineSemiInfDecay( R-p0, hDir, w2 ), -strength );
    B.add_mul( ILineSemiInfDecay( R-p1, hDir, w2 ),  strength );
    hDir.set_sub(p0,p1);
    double l = hDir.normalize();
    B.add_mul( ILineFinite( R-p0, hDir, l ), strength );
    //printf( "R(%f,%f,%f) B(%f,%f,%f)\n", R.x, R.y, R.z,  B.x, B.y, B.z );
}

inline Vec3d ISemiInfSheet( Vec3d R, Vec3d a, Vec3d b, double l ){
    // one line     ( 1-cos(theta) )/y
    // cos(theta) = x0/r = x0/sqrt(x0^2+y^2)
    // Integral_y  (1-x0*y/sqrt(x0^2 + y^2))/y  = log(sqrt(x0^+y^2)+x)
    Vec3d B; B.set_cross(a,b); B.normalize();
    double x  = a.dot(R); // along semiifinite line
    double y  = b.dot(R); //
    double x2 = x*x;
    double y_ = y + l;
    //printf( "x,y, %f %f %f \n", x,y,y_ );
    //B.mul( VortexQuantum* x*( log(sqrt(x2+y_*y_)+x) - log(sqrt(x2+y*y)+x) ) );
    B.mul( VortexQuantum* log( ( sqrt(x2+y_*y_)+x )/( sqrt(x2+y*y)+x ) ) );
    //double r1   = sqrt(x2+y*y  );
    //double r2   = sqrt(x2+y_*y_);
    //double out  = log( (r2+x)/(r1+x) );
    //printf( "%f,   %f,%f   %f,%f   %f \n", x,y,y_, r1,r2, out );
    //B.mul( VortexQuantum* out );
    return B;
};




// =================================================================
// ==================== Coulomb integrals ==========================
// =================================================================
//
// Calculated with sympy in :
// /home/prokop/Dropbox/MyDevSW/Python/_Physics/VortexLattice/basis_Coulomb_sympy.py
//
//#### rho:  1
//Fx  -1/sqrt(x**2 + y**2)
//Fy  x/(y*sqrt(x**2 + y**2))
//
//#### rho:  x
//Fx  -x/sqrt(x**2 + y**2) - log(y**2/x**2)/2 + log(1 + sqrt(x**2 + y**2)/x)
//Fy  -y/sqrt(x**2 + y**2)
//
//#### rho:  x**2
//Fx  (x**2 + 2*y**2)/sqrt(x**2 + y**2)
//Fy  -x*y/sqrt(x**2 + y**2) - y*log(y**2/x**2)/2 + y*log(1 + sqrt(x**2 + y**2)/x)
//
//#### rho:  x**3
//Fx  (x**3 + 3*x*y**2 - 3*y**2*sqrt(x**2 + y**2)*asinh(x/y))/(2*sqrt(x**2 + y**2))
//Fy  y*(x**2 + 2*y**2)/sqrt(x**2 + y**2)
// =============================================================
//#### rho:  1
//   r = sqrt(x**2 + y**2)
//Fx=  -1/r
//Fy=  (x/y)/r
//
//#### rho:  x
//Fx=  -x/r - log(  x/y * ( 1 + r/x ) )    =  -x/r   -
//Fy=  -y/r
//
//#### rho:  x**2
//Fx=  (x**2 + 2*y**2)/r
//Fy=  -x*y/r - y*log(y**2/x**2)/2 + y*log(1 + r/x)
//
//#### rho:  x**3
//Fx=  ( x*(x**3 + 3*y**2) - 3*y**2*r*log(  y/x  + r/y )    )/(2*r)
//Fy=    y*(x**2 + 2*y**2)/r

// use equality   log( y/x + r/x ) =====  asinh(y/x) ?



double sourceLineF_pow0( double x, double y, double& Fx, double& Fy ){
    double r  = sqrt(x*x+y*y);
    double ir = 1/r;
    Fx = -ir;
    Fy = (x/y)*ir;
}

double sourceLineF_pow1( double x, double y, double& Fx, double& Fy ){
    double r  = sqrt(x*x+y*y);
    double ir = 1/r;
    //Fx = -x*r - log(y*y/(x*x))/2 + log(1+r/x);
    //Fx = -x*r + log(x/y) + log(1+r/x);
    //Fx = -x*r + log( x/y * ( 1 + r/x ) );
    //Fx = -x*r + log( x/y + r/y );
    Fx = -x*r + asinh(x/y);
    Fy = -y*ir;
}

double sourceLineF_pow2( double x, double y, double& Fx, double& Fy ){
    double r  = sqrt(x*x+y*y);
    double ir = 1/r;
    Fx = x*x + 2*y*y)*ir
    Fy = y*( -x*r + asinh(x/y) );
}

double sourceLineF_pow3( double x, double y, double& Fx, double& Fy ){
    double r  = sqrt(x*x+y*y);
    double ir = 1/r;
    Fx = x*x + 2*y*y)*ir
    Fy = y*( -x*r + asinh(x/y) );
}

// ---------------- Vec3 functions

inline Vec3d sourceLine_const( Vec3d R, Vec3d hL, double L ){
    double Fx = (x*x + 2*y*y)*ir
    double Fy = y*( -x*r + asinh(x/y) );
}



#endif

