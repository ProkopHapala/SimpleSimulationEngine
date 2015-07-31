
class PhysicalSystem{
	public:
	int nmols;
	MoleculeType ** molecules;
	Vec3d  *pos,*vpos,*fpos;
	Quat4d *rot,*vrot,*frot;
	int nptmp;
	Vec3d  *Tps_i, *Tps_j;
	Vec3d  *fs_i,  *fs_j;

	AtomTypes * atypeList;
	double    *C6s,*C12s;

	double surf_z0;
	double surf_zMin; 
	double surf_Emin;
	Vec3d  surf_hat;

	OptimizerDerivs * optimizer;

// ======== initialization

	void initParams( ){
		// optimization parameters
		int nparams = (3+4)*nmols; 
		optimizer   = new OptimizerFIRE( nparams, new double[nparams], new double[nparams], new double[nparams], NULL );
		//optimizer   = new OptimizerDerivs( nparams, new double[nparams], new double[nparams], new double[nparams], NULL );
		optimizer->dt = 0.05;
		optimizer->damping = 0.05;
		int irot = 3*nmols;
		pos   = (Vec3d* )&optimizer->pos[0];
		vpos  = (Vec3d* )&optimizer->vel[0];
		fpos  = (Vec3d* )&optimizer->force[0];
		rot   = (Quat4d*)&optimizer->pos[irot];
		vrot  = (Quat4d*)&optimizer->vel[irot]; 
		frot  = (Quat4d*)&optimizer->force[irot]; 
	}

	void initTPoints(){
		// point temporary
		nptmp = 0; 
		for (int i=0; i<nmols; i++){ if( molecules[i]->natoms > nptmp ) nptmp = molecules[i]->natoms; }
		Tps_i = new Vec3d[nptmp]; 
		Tps_j = new Vec3d[nptmp];
		fs_i  = new Vec3d[nptmp]; 
		fs_j  = new Vec3d[nptmp];
		atypeList = molecules[0]->typeList;
	}

	PhysicalSystem( char const* filename, MoleculeType * molTypeList ){
		printf(" filename: %s \n", filename );

		FILE * pFile;
		pFile = fopen (filename,"r");
  		fscanf (pFile, " %i", &nmols);
		printf("mols %i \n", nmols );
		printf( " PhysicalSystem 1 \n" );
		molecules = new MoleculeType*[nmols];
		initParams( );
		printf( " PhysicalSystem 2 \n" );

		for (int i=0; i<nmols; i++){
			int itype;
			Mat3d M;
  			fscanf (pFile, " %i %lf %lf %lf     %lf %lf %lf     %lf %lf %lf", &itype, &pos[i].x, &pos[i].y, &pos[i].z,    &M.ax, &M.ay, &M.az,    &M.bx, &M.by, &M.bz );
			//printf(  " %i %f %f %f     %f %f %f     %f %f %f \n", itype, pos[i].x, pos[i].y, pos[i].z,    M.ax, M.ay, M.az,    M.bx, M.by, M.bz );
			M.a.normalize();	
			M.b.add_mul( M.a, -M.a.dot( M.b ) );
			M.b.normalize();
			M.c.set_cross( M.a, M.b );
			rot[i].fromMatrix( M );
/*
			printf(" matrix:  \n");
			printVec(M.a);
			printVec(M.b);
			printVec(M.c);
			printf(" ortonormality:  \n");
			printf( " a.a b.b c.c %f %f %f \n", M.a.norm2(),  M.b.norm2(),  M.c.norm2() );
			printf( " a.b a.c b.c %f %f %f \n", M.a.dot(M.b), M.a.dot(M.c), M.b.dot(M.c) );
			printf( " det(M) %f \n", M.determinant() );
			printf(" matrix -> quat \n");
			rot[i].fromMatrix( M );
			printf(" q: %f %f %f %f  qnorm %f  \n", rot[i].x, rot[i].y, rot[i].z, rot[i].w,   rot[i].norm2() );
			printf(" quat -> matrix  \n");
			rot[i].toMatrix( M );
			printVec(M.a);
			printVec(M.b);
			printVec(M.c);
*/
			molecules[ i ] = &molTypeList[ itype-1 ]; 
			printf("=====\n");
		}
  		fclose (pFile);
		printf( "PhysicalSystem 3 \n" );
  		initTPoints();
		printf( "PhysicalSystem 4 \n" );
	}

	PhysicalSystem( int nmols_, MoleculeType ** molecules_ ){
		nmols = nmols_;
		molecules = molecules_;
		initParams ( );
		initTPoints( );
	}


	void makeFF( ){
		makeLJparams( atypeList->ntypes, atypeList->vdwRs, atypeList->vdwEs, C6s, C12s );
	}


// =========== Rotation optimization 

	void transformPoints( const Vec3d& pos, const Quat4d& rot, int npoints, Vec3d * points, Vec3d * Tpoints ){
		Mat3d T;
		//q.toMatrix_unitary2( T );
		rot.toMatrix( T);
		//printVec( T.a );
		//printVec( T.b );
		//printVec( T.c );
		//printf( " %f %f %f \n", T.a.dot(T.b), T.a.dot(T.c), T.b.dot(T.c) );
		//printf( " %f %f %f \n", T.a.dot(T.a), T.b.dot(T.b), T.c.dot(T.c) );
		for( int i=0; i<npoints; i++ ){ 
			Vec3d Tp;
			T.dot_to_T(   points[i],  Tp );  
			//T.dot_to(   points[i],  Tp ); 
			Tpoints[i].set_add( pos, Tp  ); 
		}
	}

	void forceFromPoints( int npoints, Vec3d * points, Vec3d * forces,  const Quat4d& q,  Vec3d& fp, Quat4d& fq ){
		//printf( "----------------" );
		for( int i=0; i<npoints; i++ ){ 
			q .addForceFromPoint( points[i], forces[i], fq ); 
			fp.add( forces[i] );
			//printf( " %i   %f %f %f   %f %f %f %f \n", i,   forces[i].x, forces[i].y, forces[i].z,    fp.x, fp.y, fp.z,  fq.x, fq.y, fq.z, fq.w  );
		}
	}

	void cleanPointForce( int npoints, Vec3d * forces ){	for( int i=0; i<npoints; i++ ){  forces[i].set(0.0);  }	}

	void assembleForces( ){
		for (int i=0; i<nmols; i++){
			MoleculeType * moli = molecules[i];
			int npi = moli->natoms;
			transformPoints( pos[i], rot[i], npi, moli->xyzs, Tps_i );
			cleanPointForce( moli->natoms, fs_i );
			for (int j=0; j<i; j++){
				MoleculeType * molj = molecules[j];
				int npj = molj->natoms;
				transformPoints( pos[j], rot[j], npj, molj->xyzs, Tps_j );
				//cleanPointForce( moli->natoms, fs_i );
				cleanPointForce( molj->natoms, fs_j );
				interMolForceLJE( 
					npi, moli->atypes, moli->Qs, Tps_i, fs_i, 
					npj, molj->atypes, molj->Qs, Tps_j, fs_j,
					atypeList->ntypes, C6s, C12s
				);
				forceFromPoints( npj, molj->xyzs, fs_j,  rot[j], fpos[j], frot[j] );
				//forceFromPoints( moli->natoms, moli->xyzs, fs_i,  rot[i],  fpos[i], frot[i] );
			}
			forceMolSurf( surf_z0, surf_zMin, surf_Emin, surf_hat, npi,  Tps_i, fs_i );
			forceFromPoints( npi, moli->xyzs, fs_i,  rot[i],  fpos[i], frot[i] );
		}
	}

	void rigidOptStep( ){
		for (int i=0; i<nmols; i++){ 
			fpos[i].set(0.0); frot[i].set(0.0,0.0,0.0,0.0);   // set all forces to zero
			rot[i].normalize();                               // keep quaternion normalized, otherwise unstable !!!
		}
		assembleForces( );
		for (int i=0; i<nmols; i++){ 
			double qfq = rot[i].dot( frot[i] );	
			frot[i].add_mul( rot[i], -qfq ); // out-project component which would harm unitarity of quaternion
		}  
		optimizer->move();
	}

// =========== view utils

	void draw(){
		for (int i=0; i<nmols; i++){
			if( molecules[i]->viewlist > 0 ){	
				Mat3d rotmat;
				//rot[i].toMatrix_unitary2( rotmat );
				//rot[i].toMatrix_unitary( rotmat );
				rot[i].toMatrix( rotmat );
				float glMat[4*4];
				toGLMat( pos[i], rotmat, glMat );
				glPushMatrix();
				glMultMatrixf( glMat );
				glCallList   ( molecules[i]->viewlist );
				glPopMatrix();
			}
		};
	}

};
