
// according to 
// An Introduction to Physically Based Modeling:
// Rigid Body Simulation I—Unconstrained Rigid
// Body Dynamics
// David Baraff
// https://www.cs.cmu.edu/~baraff/sigcourse/notesd1.pdf


const double GRAV_CONTS = 6.67384e-11;

class PointBody;

///////////////////////////
//   CLASS :   PointBody
///////////////////////////

class PointBody{
	public:
	// parameters
	double	mass;
	// auxiliary parameters
	double 	invMass; 
	// State variables 
	Vec3d pos;
	Vec3d vel;
	// auxiliary variables
	Vec3d accel;

	inline void move_PointBody( double dt ){
		vel.add_mul( accel, dt );
		pos.add_mul( vel,   dt );
		//printf( "dt: %f force: ",dt ); printVec( force ); printf( " vel: " ); printVec( vel ); printf( " pos: " ); printVec( pos ); printf( "\n" );
	};

	inline void clean_temp( ){  accel.set(0.0); }

	virtual void evalForce()   { accel.set( 0,-9.81f,0 ); }
	virtual void move(double dt){ move_PointBody(dt);      }

	virtual void render(){
		drawPoint( pos );
	}
};

///////////////////////////
//   CLASS :   SpaceBody
///////////////////////////

class SpaceBody : public PointBody {
	public:
	double albedo;
	double radius;
	double surface;
	
	Vec3d color;
	int LOD[3];
	int path_list;
	int pathAcc_list;

	static const int nMaxCenters = 4;
	int nCenters                 = 0;
	SpaceBody* centers[ nMaxCenters ];

	// trajectory
	int nPath;
	double * path_t;
	Vec3d  * path_pos;
	Vec3d  * path_vel;
	Vec3d  * path_acc;

	void init(){
		path_list    = 0;
		pathAcc_list = 0;
		nPath        = 0;
		path_t   = NULL;
		path_pos = NULL;
		path_vel = NULL;
		path_acc = NULL;
	}

	virtual void evalForce()   {
		for (int i=0; i<nCenters; i++) {
			//printf( ":" );
			Vec3d r12;
			r12.set( centers[i]->pos - pos );
			double ir2 = 1.0d/r12.norm2( );
			double fr =  ir2*sqrt(ir2)*GRAV_CONTS*centers[i]->mass;
			accel.add_mul( r12, fr ); 
		}
	}

	virtual void render( ){
		Vec3d pos_  = (pos - VIEW_CENTER)*(1/VIEW_ZOOM);
		float rad_  = radius / VIEW_ZOOM;
		float rpix  = rad_ * HEIGHT; 
		if ( rpix < 1 ){
	    	if ( rpix > 0.01 ){ glDisable(GL_LIGHTING); glColor( color*rpix ); glBegin(GL_POINTS);  glVertex( pos_ ); glEnd();   }
		} else {
			glEnable (GL_LIGHTING);
			glColor( color        );
			glPushMatrix (        );
			glTranslate( pos_ );
			glScale    ( rad_ );
			if       ( rpix < 20   ){ glCallList( LOD[0] ); } 
			else if  ( rpix < 100  ){ glCallList( LOD[1] ); } 
			else                    { glCallList( LOD[2] ); }
			glPopMatrix();
		}
	}


	int deallocatePath(){
		delete [] path_t;
		delete [] path_pos;
		delete [] path_vel;
		delete [] path_acc;
	}

	int allocatePath( int n ){
		if( path_t != NULL ) deallocatePath();
		nPath    = n;
		path_t   = new double[nPath];
		path_pos = new  Vec3d[nPath];
		path_vel = new  Vec3d[nPath];
		path_acc = new  Vec3d[nPath];
	}

	double makePath( double t0, double dt, int n, int nsub ){	
		allocatePath( n );
		Vec3d pos0; pos0.set(pos);
		Vec3d vel0; vel0.set(vel); 
		double t=t0;
		for ( int i=0; i<nPath; i++  ){
			path_t  [i] = t; 
			path_pos[i].set( pos );
			path_vel[i].set( vel );
			for ( int j=0; j<nsub; j++ ){
				clean_temp( );
				evalForce ( );
				if(j==0) path_acc[i].set( accel );
				move      ( dt );
			}
			t += dt;
		}
		pos.set(pos0);
		vel.set(vel0);
		return t;
	}


	int pathToList( ){	
		if( path_list > 0 ) {	glDeleteLists( path_list, 1 );	}		
		int ilist=glGenLists(1);
		glNewList( ilist, GL_COMPILE );
		glDisable ( GL_LIGHTING );
		glBegin(GL_LINE_STRIP);
		double izoom = 1.0/VIEW_ZOOM_DEFAULT;
		glVertex( pos*izoom );
		for ( int i=0; i<nPath; i++  ){		
			glVertex( path_pos[i] * izoom );
		}
		glEnd();
		glEndList();
		path_list = ilist;
		return ilist;
	}

	int pathToListAcc( double accScale ){
		if( pathAcc_list > 0 ) {	glDeleteLists( pathAcc_list, 1 );	}		
		int ilist=glGenLists(1);
		glNewList( ilist, GL_COMPILE );
		glDisable ( GL_LIGHTING );
		//glColor3f ( 0.6, 0.8, 0.9 );
		glBegin( GL_LINES );
		double izoom = 1.0/VIEW_ZOOM_DEFAULT;
		printf( " accScale izoom %g %g \n", accScale, izoom );
		for ( int i=0; i<nPath; i++  ){		
			Vec3d p1,p2;
			p1.set(  ( path_pos[i]                        ) * izoom );
			p2.set(  ( path_pos[i] + path_acc[i]*accScale ) * izoom );
			glVertex( p1 );
			glVertex( p2 );
			//printf( " %g %g %g   %g %g %g \n", path_acc[i].x,path_acc[i].y,path_acc[i].z,   path_acc[i].x*accScale*izoom, path_acc[i].y*accScale*izoom , path_acc[i].z*accScale*izoom  );
			printf( " %g %g %g   %g %g %g \n", p1.x, p1.y, p1.z,    p2.x, p2.y, p2.z );
		}
		glEnd();
		glEndList();
		pathAcc_list = ilist;
		return ilist;
	}



};



