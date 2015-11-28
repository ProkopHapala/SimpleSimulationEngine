// https://en.wikipedia.org/wiki/Eccentric_anomaly
// https://en.wikipedia.org/wiki/True_anomaly
// https://en.wikipedia.org/wiki/Eccentric_anomaly
// http://www.castor2.ca/04_Propagation/04_True/index.html
// http://www.bogan.ca/orbits/kepler/e_anomly.html
// https://en.wikipedia.org/wiki/Kepler%27s_laws_of_planetary_motion


// find orbital elements from position
// http://space.stackexchange.com/questions/1904/how-to-programmatically-calculate-orbital-elements-using-position-velocity-vecto
/*
h=cross(r,v)
nhat=cross([0 0 1],h)

evec = ((mag(v)^2-mu/mag(r))*r-dot(r,v)*v)/mu
e = mag(evec)

energy = mag(v)^2/2-mu/mag(r)

if abs(e-1.0)>eps
   a = -mu/(2*energy)
   p = a*(1-e^2)
else
   p = mag(h)^2/mu
   a = inf

i = acos(h(3)/mag(h))

Omega = acos(n(1)/mag(n))

if n(2)<0
   Omega = 360-Omega

argp = acos(dot(n,evec)/(mag(n)*e))

if e(3)<0
   argp = 360-argp

nu = acos(dot(evec,r)/(e*mag(r))

if dot(r,v)<0
   nu = 360 - nu
*/

//https://en.wikipedia.org/wiki/Kepler%27s_laws_of_planetary_motion

inline double periodSqKepler ( double T_ref, double a_ref, double a ){	double aratio = a/a_ref; return aratio*aratio*aratio*T_ref*T_ref; }
inline double periodKepler   ( double T_ref, double a_ref, double a ){	return sqrt( periodSqKepler( T_ref, a_ref, a ) ); }

// https://en.wikipedia.org/wiki/Vis-viva_equation

inline double velocitySqVisVia( double GM, double a, double r ){ return GM*( 2/r - 1/a ); }

inline double velocityVisVia  ( double GM, double a, double r ){ return sqrt( velocitySqVisVia( GM, a, r ) ); }


inline double anomaly_MeanToEccentric( double M, double e, double maxErr ){
	// from http://alpheratz.net/dynamics/twobody/KeplerIterations_summary.pdf
	double E = M;
	for( int i=0; i<10; i++ ){
		double se = sin( E );
		double ce = cos( E );
		double A  = e*se;
		double B  = E - A - M;
		double C  = 1 - e*ce;
		double err  = E - B/( C - 0.5*A*(B/C) );
		if( fabs(err) > maxErr ){
			E = E - err; 
		}else{
			break;
		}
	}
	return E;
} 

/*
// probably not working, at leas for large angles
inline double anomaly_MeanToEccentric( double M, double e ){
	// from https://en.wikipedia.org/wiki/Orbital_mechanics
	double e2 = e*e;
	double e3 = e2*e;
	double e4 = e3*e;
	double a  = 1/(1-e);
	double a3 = a*a*a;
	double M2 = M*M;
	return a*M + M2*a3*6*(  -e + M2*a3*20*( ( 9*e2 + e ) + M2*a3*42*( -( 225*e3 + 54*e2 + e ) +  M2*a3*72*( 11025*e4 + 4131*e3 + 243*e2 + e ) ) ) );
} 
*/

inline double anomaly_EccentricToTrue( double E, double e ){
	double ce = cos( E );
	return  acos( ( ce - e ) / ( 1.0d - e*ce ) );
}

/*
inline void pos_on_orbit( double E, double e, const Vec3d& a, const Vec3d& b, Vec3d& p ){
	double ca = cos( E );
	double sa = sin( E );
	p.set_mul( a, ca - e );
	p.add_mul( b, sa     );
}
*/

inline void pos_on_orbit( double E, const Vec3d& a, const Vec3d& b, Vec3d& p ){
	p.set_mul( a, cos( E ) );
	p.add_mul( b, sin( E ) );
}

inline void getOrbitAxes( double incl, double Omega, double omega, Vec3d& ahat, Vec3d& bhat ){
	//	incl  ... inclination
	//  Omega ... Longitude of ascending node  
	//  omega ... Argument of periapsis ( Longitude of perihelion )
	double ci = cos( incl  );
	double si = sin( incl  );
	double co = cos( omega );
	double so = sin( omega );
	double cO = cos( Omega );
	double sO = sin( Omega );
	//Vec3d a,b;
	//a.set(  co, -so,     0 );
	//b.set(  so,  co,     0 );
	//a.set(  co, -so*ci, -so*si );
	//b.set(  so,  co*ci,  co*si );
	//a.set(  cO*co - sO*( -so*ci ) , cO*(-so*ci) + sO*co , -so*si );
	//b.set(  cO*so - sO*(  co*ci ) , cO*( co*ci) + sO*so ,  co*si );
	double soci = -so*ci;
	double coci =  co*ci;
	ahat.set(  cO*co - sO*soci , cO*soci + sO*co , -so*si );
	bhat.set(  cO*so - sO*coci , cO*coci + sO*so ,  co*si );
};


inline void getOrbitAxes( double incl, double Omega, double omega, double a, double e, Vec3d& avec,  Vec3d& bvec ){
	//	a     ... semimajor axis length
	//  e     ... eccentricity
	//	incl  ... inclination
	//  Omega ... Longitude of ascending node  
	//  omega ... Argument of periapsis ( Longitude of perihelion )
	getOrbitAxes( incl, Omega, omega, avec, bvec );
	avec.mul( a  );
	bvec.mul( a * sqrt( 1 - e*e )  );
};











