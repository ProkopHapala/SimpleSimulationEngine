
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include <SoftBody.h>

// ==== Dynamics

void SoftBody::evalForces( ){
	//printf( "==============\n" );
	for( int i=0; i<npoints; i++ ){
		//forces[i].set(0.0);
		evalPointForce( i, gravity, airFlow );
		//printf( " %i  %f %f %f   %f %f %f \n", i, forces[i].x, forces[i].y, forces[i].z, mass[i], invMass[i], drag[i] );
	}
	for( int i=0; i<nbonds;  i++ ){	addBondForce( bonds[i] ); }
}

void SoftBody::applyConstrains(){
	for( int i=0; i<nfix; i++ ){
		int ip = fix[i];
		if( (ip>0)&&(ip<npoints) ){
            forces    [ip].set(0.0);
            velocities[ip].set(0.0);
		}
	}
}

void SoftBody::move_LeapFrog( ){
	for (int i=0; i<npoints; i++){
		velocities[i].mul( damp );
		velocities[i].add_mul( forces[i], invMass[i] * dt );
		points[i].    add_mul( velocities[i], dt );
/*
		printf( " (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) %3.3f %3.3f %3.3f \n",
            points[i].x, points[i].y, points[i].z,
            velocities[i].x, velocities[i].y, velocities[i].z,
            forces[i].x, forces[i].y, forces[i].z,
            invMass[i] ,  dt,  damp
        );
*/
	}
}

void SoftBody::step( ){
    //for (int i=0; i<npoints; i++){ forces[i].set(0.0); }
    evalForces( );
    applyConstrains();
    move_LeapFrog( );
}

// ==== Setup

void SoftBody::allocate(	int npoints_, int nbonds_, int nfix_,  	Vec3d  * points_, double * mass_, double * drag_,  	int    * fix_ ){
	npoints = npoints_; nbonds = nbonds_; nfix = nfix_;

	if( points_ != NULL ){ points = points_; }else{  points=new Vec3d [npoints]; }
	if( mass_   != NULL ){ mass = mass_;     }else{  mass  =new double[npoints]; for(int i=0; i<npoints; i++ ){ mass[i]=0; } }
	if( drag_   != NULL ){ drag = drag_;     }else{  drag  =new double[npoints]; for(int i=0; i<npoints; i++ ){ drag[i]=0; } }
	velocities = new Vec3d [npoints ];
	forces     = new Vec3d [npoints ];
	invMass    = new double[npoints ];

    bonds  = new Bond[ nbonds ];

    if( fix_ != NULL ){ fix = fix_;  }else{ fix = new int[nfix]; for(int i=0; i<nfix; i++ ){ fix[i]=-1; } }

    gravity.set(0.0,-9.81,0.0);
    airFlow.set(0.0,0.0,0.0);
}

void SoftBody::deallocateAll( ){
    if( points     != NULL ) delete points;       //printf( "DEBUG delete points\n" );
    if( velocities != NULL ) delete velocities;   //printf( "DEBUG delete velocities\n" );
    if( forces     != NULL ) delete forces ;      //printf( "DEBUG delete forces\n" );
    if( mass       != NULL ) delete mass;         //printf( "DEBUG delete mass\n" );
    if( drag       != NULL ) delete drag;         //printf( "DEBUG delete drag\n" );
    if( invMass    != NULL ) delete invMass;      //printf( "DEBUG delete invMass\n" );
    if( bonds      != NULL ) delete bonds;        //printf( "DEBUG delete bonds\n" );
    if( fix        != NULL ) delete fix;          //printf( "DEBUG delete fix\n" );
}

void SoftBody::prepareBonds( bool l0_fromPos ){
    for( int i=0; i<nbonds; i++ ){
        Bond& bond = bonds[i];
        if( l0_fromPos ){  bond.l0 = getBondLength( bond.i, bond.j ); }
        double dmass = bond.getMass() * 0.5;
        double ddrag = bond.getDrag() * 0.5;
        mass[bond.i] += dmass;
        mass[bond.j] += dmass;
        drag[bond.i] += ddrag;
        drag[bond.j] += ddrag;
        //printf( " bond %i %i %e %e %e \n", i, bond.type.id, bond.type.linearDensity, bond.type.kPress, bond.type.kTens );
        //printf( " bond %i %f %f  %f \n", i, bond.l0, dmass, ddrag );
    }
}

void SoftBody::preparePoints(  bool clearVelocity, double constDrag, double constMass ){
    for (int i=0; i<npoints; i++){
        if( clearVelocity ) velocities[i].set(0.0);
        if( constDrag > 0 ) drag[i] = constDrag;
        if( constMass > 0 ) mass[i] = constMass;
		invMass[i] = 1/mass[i];
		//printf( " poins %i %e %e   (%e,%e,%e) \n", i, mass[i], invMass[i], velocities[i].x, velocities[i].y, velocities[i].z  );
	}
}



/*
void allocate_RVF          ( ){  }
void deallocate_RVF        ( ){ delete pos; delete vel; delete force; }
void allocate_pointParams  ( ){ mass = new double[npoints]; drag = new double[npoints]; invMass = new double[npoints]; }
void deallocate_pointParams( ){ delete mass; delete drag; delete invMass; }
void init_Axuliary         ( ){ for( int i=0; i<npoints; i++ ){ invMass = 1/mass; } }

void allocate_Bonds        ( ){ bonds = new Bonds[nbonds]; }
void deallocate_Bonds      ( ){ delete bonds; }

void init( int npoints_, int nbonds_, int nfix_, bool own_points_, bool own_mass_ ){
    own_points = own_points_;
    own_mass   = own_mass_;
    own_fix    = own_fix_;
    npoints = npoints_;
    nbonds  = nbonds_;
    nfix    = nfix_;
    if(own_points){
        pos = new Vec3d[npoints]; vel = new Vec3d[npoints]; force = new Vec3d[npoints];
    }


}
*/


/*
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

*/

/*
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
*/


/*

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

*/

