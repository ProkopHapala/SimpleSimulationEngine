

/*

kinetic viscosity of watter   1.002e-3 [m^2/s]   1.004e-6 [m^2/s]			0.1 * 1 / 1.004e-6 = 100000
kinetic viscosity of air      1.983e-5 [m^2/s]   1.568e-5 [m^2/s]           1.0 * 5 / 1.568e-5 = 
=> Re = 100,000 ... 300,000

CD limits ( Re = 10,000 )
0.47 sphere
1.28 flat plate perpendicular to flow (3D)


*/

class KinematicBody2D {
	public:
	double phi;
	Vec2d  rot;
	Vec2d  pos;
	inline void setAngle( double phi_ ){  phi = phi_; rot.fromAngle( phi_ ); }
};

const double SAFETY_v = 1e-6;

class AeroSurf2D : public KinematicBody2D {
	public:
	double area;  // [m^2]
	double CD0    = 0.02;  
	double dCD    = 0.9;  
	double dCDS   = 0.9;  
	double dCL    = 6.28;
	double dCLS   = 2.82743338823;
	double sStall = 0.16;
	double wStall = 0.08;
	
	int      nPolar;
	double * cDs;
	double * cLs;

	bool loadPolar( char const* filename ){
		printf(" loading polar from : %s \n", filename );
		FILE * pFile;
		pFile = fopen (filename,"r");
  		fscanf (pFile, " %i", &nPolar);
		//printf("natoms %i \n", natoms );
		cDs = new double[ nPolar ];
		cLs = new double[ nPolar ];
		for (int i=0; i<nPolar; i++){
  			int nw = fscanf (pFile, " %i %lf %lf %lf %lf", &cDs[i], &cLs[i] );
		}
  		fclose (pFile);
  		return 0;
	}

	void plot_polar( double x0, double y0, double fscale, double phi0 ){
		glColor3f( 0.8f, 0.2f, 0.2f );	drawPolarFunc( x0, y0, fscale, nPolar, phi0, cDs );
		glColor3f( 0.2f, 0.2f, 0.8f );  drawPolarFunc( x0, y0, fscale, nPolar, phi0, cLs );
	}


	inline void polarModel( double ca, double sa, double& CD, double& CL ){
		//double abs_sa = fabs( sa ); 
		double abs_sa = (sa>0)?sa:-sa;
		double wS     = trashold_cub( abs_sa, sStall, sStall+wStall );
		double mS     = 1 - wS; 
		CD     = CD0 + ( mS*dCD*abs_sa + wS*dCDS        ) * abs_sa;
		//CL     =       ( mS*dCL        + wS*dCLS*abs_ca ) * sa;
		if( ca <0 ){ ca=-ca; sa=-sa; };
		CL     =       ( mS*dCL        + wS*dCLS*ca ) * sa;
	}


	inline void assertAeroForce( RigidBody2D& platform, const Vec2d& vel, double density ){
		Vec2d gpos, grot, vhat, force;

		grot  .set_mul_cmplx( platform.rot, rot );
		gpos  .set_mul_cmplx( platform.rot, pos );

		vhat.x = vel.x + gpos.y * platform.omega;
		vhat.y = vel.y - gpos.x * platform.omega;

		double vr2 = vhat.norm2(); 
		double ivr = 1/sqrt( vr2 + SAFETY_v ); 
		vhat.mul( ivr ); 


		double sa = vhat.dot  ( grot );
		double ca = vhat.cross( grot );

		//double sa = vhat.cross  ( grot );
		//double ca = vhat.dot    ( grot );

		//double CD = CD0 + ddCD*sa*sa;
		//double CL = dCL*sa*ca*ca;

		double CD,CL;
		polarModel( ca, sa, CD, CL );

		double prefactor = vr2 * density * area;
		force.x = prefactor*( CD*vhat.x - CL*vhat.y );  
		force.y = prefactor*( CD*vhat.y + CL*vhat.x );

		//printf( " sa, ca, CD, CL %f %f %f %f vhat %f %f grot %f %f \n", sa, ca, CD, CL,  vhat.x, vhat.y, grot.x, grot.y );
		//printf( " sa ca  %f %f CD CL %f %f omega |v| |f| %f %f %f \n", sa, ca, CD, CL, platform.omega, sqrt(vr2), force.norm() );

		//printf( " force %f %f \n",  force.x, force.y );
		//printf( "   vhat %g %g   force %g %g   ca sa %g %g    CD CL %g %g \n", vhat.x, vhat.y,  force.x, force.y,   ca, sa,    CD, CL );
		printf( "   vhat %g %g   force %g %g   CD CL %g %g  prefactor %g %g %g %g \n", vhat.x, vhat.y,  force.x, force.y,   CD, CL,   prefactor,  vr2, density, area );

		platform.apply_force( force, gpos );

		gpos.add( platform.pos );
		glColor3f( 0.8f, 0.1f, 0.1f ); drawVecInPos( force*0.01, gpos );
		glColor3f( 0.1f, 0.1f, 0.8f ); drawVecInPos(  vhat, gpos );
		glColor3f( 0.1f, 0.1f, 0.1f ); drawVecInPos(  {CD,CL}, (Vec2d){0.0d,0.0d} );

		
	}

	virtual void draw( RigidBody2D& platform ){
		Vec2d gpos, grot;  		
		grot  .set_mul_cmplx( platform.rot, rot );
		gpos  .set_mul_cmplx( platform.rot, pos );
		gpos.add( platform.pos );
		//printf( " platform.rot %f %f grot %f %f  gpos %f %f \n",  platform.rot.x, platform.rot.y, grot.x, grot.y, gpos.x, gpos.y );
		//drawPointCross( gpos, 0.1 ); drawVecInPos( grot, gpos );

		float lperp = 0.1;
		float llong = 0.5;
		glBegin(GL_LINES);
			glVertex3f( (float)( gpos.x-grot.x*lperp), (float)(gpos.y-grot.y*lperp), 1 );   glVertex3f( (float)(gpos.x+grot.x*lperp), (gpos.y+grot.y*lperp), 1 );
			glVertex3f( (float)( gpos.x-grot.y*llong), (float)(gpos.y+grot.x*llong), 1 );   glVertex3f( (float)(gpos.x+grot.y*llong), (gpos.y-grot.x*llong), 1 );
		glEnd();
	}

};
