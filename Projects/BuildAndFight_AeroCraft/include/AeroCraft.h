
// ==== Classes Forward declaration

class AeroSurface; 
class AeroCraft;

// ==== Classes

class AeroSurface{
	public:
	AeroCraft *craft;
	
	Vec3d C;   // Aerodynamic coefficients in each direction
	Vec3d lpos;
	Mat3d lrot;

	// Methods
	inline void  applyForceSimple( const Vec3d& vair0 );
	void render( );
};

class AeroCraft : public RigidBody {
	public:
	AeroSurface wingLeft,wingRight, elevator, rudder;
	// Methods
	inline void applyAeroForces( const Vec3d& vwind );
	void render();
};

// ==== Method implementation


///////////////////////////////////////////
//             AeroSurface
///////////////////////////////////////////

inline void  AeroSurface::applyForceSimple( const Vec3d& vair0 ){
	Mat3d grot;  grot.set_mmul( lrot, craft->rotMat );
	Vec3d gdpos; craft->rotMat.dot_to( lpos, gdpos );

	Vec3d uair; 
	uair.set_cross( gdpos, craft->omega );
	//uair.set(0);
	uair += vair0;
	double vrair2  = uair.norm2();
	if( vrair2 >0 ){
		double vrair = sqrt(vrair2);
		uair *= (1/vrair);		

		// plane force		
		double ca = grot.a.dot( uair );
		double cb = grot.b.dot( uair );
		double cc = grot.c.dot( uair ); 
		
		Vec3d force;
		//force.set( uair*(C.y*vrair2) );
		//force.set( grot.b*(C.y*cb*vrair2) );
		force.set( grot.a*(C.x*ca*vrair2) + grot.b*(C.y*cb*vrair2) + grot.c*(C.z*cc*vrair2) );

		//printf( "vrair %f \n", vrair );
		//printVec( uair ); printf("uair\n");
		//printVec( uair ); printf("force\n");

		craft->applyForce( force, gdpos );

		//drawMatInPos( craft->rotMat, craft->pos );
		drawMatInPos( grot, craft->pos + gdpos );

		//glColor3f( 0.0f,0.0f,0.0f );
		//drawLine( craft->pos + gdpos,craft->pos + gdpos + (grot.b*5));
		glColor3f( 0.9f,0.0f,0.9f );
		drawLine( craft->pos + gdpos,craft->pos + gdpos + (force*0.1));
		//glColor3f( 0.0f,0.5f,0.0f );
		//drawLine( craft->pos + gdpos,craft->pos + gdpos + (uair*2));

	}
};

void AeroSurface::render( ){
	//drawLine( const Vec3d& p1, const Vec3d& p2 );
	double sz = 3*sqrt( C.z );
	glEnable (GL_LIGHTING);
	glColor3f( 0.5f,0.5f,0.5f );
	glBegin  (GL_QUADS);
		glNormal3d( lrot.b.x, lrot.b.y, lrot.b.z );	          	     
		glVertex3d( lpos.x-sz*lrot.a.x, lpos.y-sz*lrot.a.y, lpos.z-sz*lrot.a.z ); 
		glVertex3d( lpos.x-sz*lrot.c.x, lpos.y-sz*lrot.c.y, lpos.z-sz*lrot.c.z ); 
		glVertex3d( lpos.x+sz*lrot.a.x, lpos.y+sz*lrot.a.y, lpos.z+sz*lrot.a.z );
		glVertex3d( lpos.x+sz*lrot.c.x, lpos.y+sz*lrot.c.y, lpos.z+sz*lrot.c.z ); 
	glEnd();		
};


///////////////////////////////////////////
//              AeroCraft 
///////////////////////////////////////////

inline void AeroCraft::applyAeroForces( const Vec3d& vwind ){
	Vec3d vair = vwind - vel;
	wingLeft .applyForceSimple( vair );
	wingRight.applyForceSimple( vair );
	rudder   .applyForceSimple( vair );
	elevator .applyForceSimple( vair );
};

void AeroCraft::render(){
	glPushMatrix();
	float glmat[16];
	toGLMat( pos, rotMat, glmat);
	glMultMatrixf( glmat );
	glDisable (GL_LIGHTING);
	glColor3f( 0.3f,0.3f,0.3f );
	drawLine( {0,0,0},wingLeft.lpos);
	drawLine( {0,0,0},wingRight.lpos);
	drawLine( {0,0,0},elevator.lpos);
	wingLeft.render();
	wingRight.render();
	rudder.render();
	elevator.render();
	glPopMatrix();
};
