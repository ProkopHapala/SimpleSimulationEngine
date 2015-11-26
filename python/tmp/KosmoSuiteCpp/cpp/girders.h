
// girder composed of two elements
void makeGirder_Type1( 
	int n, double length, double width, const Vec3d& dir, const Vec3d& up, const Vec3d& side,
	int& npoints, int& nbonds, Vec3d*& points, int*& bonds, int*& bondTypes 
 ){

	double dl = length/(2*n); 
	int dnb = 2+4+4+4;
	int dnp = 4;
	nbonds  = n*dnb - (4+4);
	npoints = n*dnp;
	points     = new Vec3d[ dnp*n ];
	bonds      = new int  [ 2*nbonds ];
	bondTypes  = new int  [   nbonds ];

	int dnb2=2*dnb;
	int it  = 0;
	int i00 = 0;
	int b0  = 0;
	for (int i=0; i<n; i++){
		// points
		int i01=i00+1; int i10=i00+2; int i11=i00+3;
/*
		points[i00].set( -width, 0, dl*(2*i  ) );	
		points[i01].set( +width, 0, dl*(2*i  ) );	
		points[i10].set( 0, -width, dl*(2*i+1) );	
		points[i11].set( 0, +width, dl*(2*i+1) );
*/
		points[i00].set_lincomb( -width, 0, dl*(2*i  ), side, up, dir );	
		points[i01].set_lincomb( +width, 0, dl*(2*i  ), side, up, dir );	
		points[i10].set_lincomb( 0, -width, dl*(2*i+1), side, up, dir );	
		points[i11].set_lincomb( 0, +width, dl*(2*i+1), side, up, dir );
		// bonds
		// perperndicular
		bonds[ b0    ] = i00; bonds[ b0    +1 ] = i01;
		bonds[ b0+2  ] = i10; bonds[ b0+2  +1 ] = i11;
		for( int k=0; k<2; k++ ){ bondTypes[it]=1; it++; };
		// zigzag inside cell
		bonds[ b0+4  ] = i00; bonds[ b0+4  +1 ] = i10;
		bonds[ b0+6  ] = i00; bonds[ b0+6  +1 ] = i11;
		bonds[ b0+8  ] = i01; bonds[ b0+8  +1 ] = i10;
		bonds[ b0+10 ] = i01; bonds[ b0+10 +1 ] = i11;
		for( int k=0; k<4; k++ ){ bondTypes[it]=2; it++; };
		if( i<(n-1) ){
			// zigzag outside cell
			bonds[ b0+12 ] = i10; bonds[ b0+12 +1 ] = i00 + dnp;
			bonds[ b0+14 ] = i10; bonds[ b0+14 +1 ] = i01 + dnp;
			bonds[ b0+16 ] = i11; bonds[ b0+16 +1 ] = i00 + dnp;
			bonds[ b0+18 ] = i11; bonds[ b0+18 +1 ] = i01 + dnp;
			for( int k=0; k<4; k++ ){ bondTypes[it]=2; it++; };
			// longitudinatl
			bonds[ b0+20 ] = i00; bonds[ b0+20 +1 ] = i00 + dnp;
			bonds[ b0+22 ] = i01; bonds[ b0+22 +1 ] = i01 + dnp;
			bonds[ b0+24 ] = i10; bonds[ b0+24 +1 ] = i10 + dnp;
			bonds[ b0+26 ] = i11; bonds[ b0+26 +1 ] = i11 + dnp;
			for( int k=0; k<4; k++ ){ bondTypes[it]=0; it++; };
		}
		i00+=dnp;
		b0 +=dnb2; 
	}
}

/*

// more general version of the same ( for m=2 the same, for m=3 also interesting, for m>6 tube with thinner and thiner shell )
void makeGirder_TypeM( 
	int n, int m, double length, double r1, double r2, 
	const Vec3d& dir, const Vec3d& up, const Vec3d& side,
	int& npoints, int& nbonds, Vec3d*& points, int*& bonds, int*& bondTypes 
 ){

	double dl = length/(2*n); 
	int dnp = 2*m;
	int dnb = 2*( (m-1) +m+m+m );
	nbonds  = n*dnb - (m+m);
	npoints = n*dnp;
	points     = new Vec3d[ dnp*n ];
	bonds      = new int  [ 2*nbonds ];
	bondTypes  = new int  [   nbonds ];

	int dnb2=2*dnb;
	int it  = 0;
	int i00 = 0;
	int b0  = 0;

	double dphi = 2*M_PI/m;
	double sa = sin(dphi);
	double ca = cos(dphi);

	int b = 0;
	for (int i=1; i<n; i++){
		int pi0 = 2*m*i;
		for ( int j=0;i<m; j+=2 ){
			int p0 = pi0 + j*2;
			int p1 = p0  + 1;
			points[ p0   ].set_lincomb( -r1, 0, dl*(2*i  ), side, up, dir );	rot_csa<double>( ca, sa, ux, uy );
			points[ p0+1 ].set_lincomb( -r2, 0, dl*(2*i  ), side, up, dir );	rot_csa<double>( ca, sa, ux, uy );
			// equatorial
			int jw = 2*((j+1)%m);
			bonds [ b ] = p0; bonds[ b0+1 ] = pi0+jw;    b+=2;
			bonds [ b ] = p1; bonds[ b0+1 ] = pi0+jw+1;  b+=2;
			// zigzag in cell
			bonds [ b ] = p0; bonds[ b0+1 ] = p0-2*m;
			bonds [ b ] = p1; bonds[ b0+1 ] = p0-2*m;
			if( i<(n-1) ){
				// longitudinal
				bonds [ b ] = p0; bonds[ b0+1 ]     = p0+2*m;  b+=2;
				bonds [ b ] = p1; bonds[ b0+1 ]     = p1+2*m;  b+=2;
				// zigzag out cell
				bonds [ b  ] = i00; bonds[ b+1 ] = p0-2*m;
				bonds [ b  ] = i10; bonds[ b+1 ] = p0-2*m;
			}
		}
	}
}

*/


static Uint32 typeColors[]={ 0xFFFF0000, 0xFF008F00, 0xFF0000FF };

int drawGirder( int npoints, int nbonds, Vec3d* points, int* bonds, int* bondTypes ){
	int nvert=0;
	glDisable (GL_LIGHTING);
	glBegin   (GL_LINES);
	int ib=0;
	for( int i=0; i<nbonds; i++ ){
		setColorInt32( typeColors[bondTypes[i]] );
		Vec3f p1,p2;
		convert( points[bonds[ib  ]], p1 );	glVertex3f( p1.x, p1.y, p1.z );	nvert++;
		convert( points[bonds[ib+1]], p2 );	glVertex3f( p2.x, p2.y, p2.z ); nvert++;
		//printf( " %i    %i %i %i   %f %f %f    %f %f %f \n", i, bonds[ib], bonds[ib+1], bondTypes[i], p1.x, p1.y, p1.z,   p2.x,p2.y,p2.z  );
		ib+=2;
	}
	glEnd();
	return nvert;
}


































