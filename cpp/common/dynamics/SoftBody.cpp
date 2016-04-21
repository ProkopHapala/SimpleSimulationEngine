
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include <SoftBody.h>

void SoftBody::bondFromType( int ib, int * bondTypes, const BondTypes& bondTypesBooks ){
	int itype = bondTypes[ib];
	int ib2 = ib<<1;
	int i   = bonds[ib2  ];
	int j   = bonds[ib2+1];
	Vec3d d; d.set_sub( points[i], points[j] );
	l0s[ib]			= d.norm();
	double area     = bondTypesBooks.area   [itype];
	double bondMass = bondTypesBooks.density[itype] * area * l0s[ib];
	kTens [ ib ]    = bondTypesBooks.kTens  [itype] * area;
	kPres [ ib ]    = bondTypesBooks.kPres  [itype] * area;
	//double tensStrength = bondTypes_tensStrength[i] * area;
	//double presStrength = bondTypes_presStrength[i] * area;
	mass[i] += bondMass*0.5;
	mass[j] += bondMass*0.5;
	//printf( " %i %i   %i   %i %i   %f %f %f %f  \n", ib, ib2, itype, i, j,   bondMass, l0s[ib], kTens[ib], kPres[ib]   );
}

void SoftBody::evalForces( const Vec3d& gravity, const Vec3d& airFlow ){
	//printf( "==============\n" );
	for( int i=0; i<npoints; i++ ){
		//forces[i].set(0.0);
		evalPointForce( i, gravity, airFlow );
		//printf( " %i  %f %f %f   %f %f %f \n", i, forces[i].x, forces[i].y, forces[i].z, mass[i], invMass[i], drag[i] );
	}
	for( int i=0; i<nbonds;  i++ ){	evalBondForce( i ); }
}

void SoftBody::applyConstrains(){
	for( int i=0; i<nfix;   i++ ){
		int ip = fix[i];
		forces    [ip].set(0.0);
		velocities[ip].set(0.0);
	}
}

void SoftBody::move( double dt, double damp ){
	for (int i=0; i<npoints; i++){
		velocities[i].mul( damp );
		velocities[i].add_mul( forces[i], invMass[i] * dt );
		points[i].  add_mul( velocities[i], dt );
	}
}

void SoftBody::draw( float forceScale ){
	glDisable (GL_LIGHTING);
	glBegin   (GL_LINES);
	for( int ib=0; ib<nbonds;  ib++ ){
		int ib2 = ib<<1;
		Vec3f pi,pj, d;
		convert( points[ bonds[ib2   ] ], pi );
		convert( points[ bonds[ib2+1 ] ], pj );
		d.set_sub( pi, pj );
		float l  = d.norm(); // this should be optimized
		float dl = ( l - l0s[ib] ) / l;
		if( dl > 0 ){
			float f = (float)(kTens[ib]*dl)*forceScale;
			glColor3f( f,0,0 );
		}else{
			float f = (float)(kPres[ib]*dl)*forceScale;
			glColor3f( 0,0,-f );
		}
		glVertex3d( pi.x, pi.y, pi.z );
		glVertex3d( pj.x, pj.y, pj.z );
		//printf( " %f %f %f   %f %f %f \n", pi.x, pi.y, pi.z,   pj.x, pj.y, pj.z  );
	}
	glEnd();
}

void SoftBody::init(
	int npoints_, int nbonds_, int nfix_,
	Vec3d  * points_, double * mass_, double * drag_,
	int    * bonds_, int * bondTypes, const BondTypes& bondTypeBooks,
	int    * fix_
){
	npoints = npoints_;  nbonds = nbonds_; nfix = nfix_;
	points  = points_; bonds   = bonds_; fix = fix_;
	if( mass_ != NULL ){ mass = mass_; }else{  mass=new double[npoints]; for(int i=0; i<npoints; i++ ){ mass[i]=0; } }
	if( drag_ != NULL ){ drag = drag_; }else{  drag=new double[npoints]; for(int i=0; i<npoints; i++ ){ drag[i]=0; } }
	velocities = new Vec3d[npoints];
	forces     = new Vec3d[npoints];
	invMass    = new double[ npoints ];

	l0s   = new double[nbonds];
	kTens = new double[nbonds];
	kPres = new double[nbonds];
	for (int i=0; i<nbonds; i++ ){ bondFromType( i, bondTypes, bondTypeBooks ); }
	for (int i=0; i<npoints; i++){
		invMass[i] = 1/mass[i]; velocities[i].set(0.0);
		// printf( " %f   %f %f %f \n", invMass[i], velocities[i].x, velocities[i].y, velocities[i].z  );
	}
}

void SoftBody::init(
	int npoints_, int nbonds_, int nfix_,
	Vec3d  * points_, double * mass_,    double * drag_,
	int    * bonds_,  double * kTens_,   double * kPres_,   double * l0s_,
	int    * fix_
){
	npoints = npoints_;  nbonds = nbonds_; nfix = nfix_;
	points = points_;    mass    = mass_;
	bonds   = bonds_;    kTens      = kTens_;      kPres = kPres_; l0s = l0s_;
	if( drag_ != NULL ){ drag = drag_; }else{  drag=new double[npoints]; for(int i=0; i<npoints; i++ ){ drag[i]=0; } }
	fix     = fix_;
	velocities = new Vec3d[npoints];
	forces     = new Vec3d[npoints];
	invMass    = new double[ npoints ];
	for (int i=0; i<npoints; i++){ invMass[i] = 1/mass[i]; velocities[i].set(0.0); }
}

