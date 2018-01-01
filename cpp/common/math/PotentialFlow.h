
#ifndef  PotentialFlow_h
#define  PotentialFlow_h

// http://s6.aeromech.usyd.edu.au/aerodynamics/index.php/sample-page/subsonic-aerofoil-and-wing-theory/3d-vortex-lattice-method/
// http://web.mit.edu/16.unified/www/SPRING/fluids/Spring2008/LectureNotes/f06.pdf

const double VortexQuantum = 1.0/( 4*M_PI );

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
    // Optimized:    since a = |hL x R|     dB/a = (hL x R)/|hL x R|^2
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
    B.mul( (1-c1)*VortexQuantum/B.norm2() ); // notice there is no "a"
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

inline void horseshoe( Vec3d& B, Vec3d R, Vec3d p0, Vec3d p1, Vec3d hDir, double strenght ){
    B.add_mul( ILineSemiInf( R-p0, hDir ), -strenght );
    B.add_mul( ILineSemiInf( R-p1, hDir ),  strenght );
    hDir.set_sub(p0,p1);
    double l = hDir.normalize();
    B.add_mul( ILineFinite( R-p0, hDir, l ), strenght );
    //printf( "R(%f,%f,%f) B(%f,%f,%f)\n", R.x, R.y, R.z,  B.x, B.y, B.z );
}

inline void horseshoeDecay( Vec3d& B, Vec3d R, Vec3d p0, Vec3d p1, Vec3d hDir, double strenght, double w2 ){
    B.add_mul( ILineSemiInfDecay( R-p0, hDir, w2 ), -strenght );
    B.add_mul( ILineSemiInfDecay( R-p1, hDir, w2 ),  strenght );
    hDir.set_sub(p0,p1);
    double l = hDir.normalize();
    B.add_mul( ILineFinite( R-p0, hDir, l ), strenght );
    //printf( "R(%f,%f,%f) B(%f,%f,%f)\n", R.x, R.y, R.z,  B.x, B.y, B.z );
}

#endif

