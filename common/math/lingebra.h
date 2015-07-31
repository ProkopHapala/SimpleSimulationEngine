

// ===============================================
// ======                                  =======
// ======              vector              =======
// ======                                  =======
// ===============================================

inline double dot(	int n,  double* a, double* b ){
	double sum = 0;
	for (int i=0; i<n; i++ ){		sum+= a[i]*b[i];	} 
	return sum;
}

inline void set( int n, double  f,            double* out ){  	for (int i=0; i<n; i++ ){ out[i] = f;	      } }
inline void add( int n, double  f, double* b, double* out ){  	for (int i=0; i<n; i++ ){ out[i] = f+b[i];    } }
inline void mul( int n, double  f, double* b, double* out ){  	for (int i=0; i<n; i++ ){ out[i] = f*b[i];    } }

inline void set( int n, double* a,            double* out ){  	for (int i=0; i<n; i++ ){ out[i] = a[i];      } }
inline void add( int n, double* a, double* b, double* out ){  	for (int i=0; i<n; i++ ){ out[i] = a[i]+b[i]; } }
inline void sub( int n, double* a, double* b, double* out ){  	for (int i=0; i<n; i++ ){ out[i] = a[i]-b[i]; } }
inline void mul( int n, double* a, double* b, double* out ){  	for (int i=0; i<n; i++ ){ out[i] = a[i]*b[i]; } }
inline void div( int n, double* a, double* b, double* out ){  	for (int i=0; i<n; i++ ){ out[i] = a[i]/b[i]; } }
inline void fma( int n, double* a, double* b, double f, double* out ){ for(int i=0; i<n; i++) { out[i]=a[i]+f*b[i]; }  }



void random_vector ( int n, double xmin, double xmax, double * out ){
	double xrange = xmax - xmin;
	for (int i=0; i<n; i++ ){		out[i] = xmin + xrange*randf();	} 
}

void print_vector( int n, double * a ){
	for (int i=0; i<n; i++ ){	printf( "%f ", a[i] );	} 
	printf( "\n" );
}

// ===============================================
// ======                                  =======
// ======              Matrix              =======
// ======                                  =======
// ===============================================


// creates double** from any continuous memory block folowing *p
double ** from_continuous( int m, int n, double *p ){
	double ** A = new double*[m];
	for (int i=0; i<m; i++ ){ A[i] = &(p[i*n]); }
	return A;
}

// allocates new matrix 
double ** new_matrix( int m, int n ){
	double ** A = new double*[m];
	for (int i=0; i<m; i++ ){ A[i] = new double[n]; }
	return A;
}

// delete matrix
double ** delete_matrix( int m, double** A ){
	for (int i=0; i<m; i++ ){ delete A[i]; }
	delete [] A;
}


// transpose matrix
void transpose( int m, int n, double** A, double** TA ){
	for (int i=0; i<m; i++ ){
		for (int j=0; j<n; j++ ){
			TA[i][j] = A[j][i];
		} 
	} 
}

// dot product
void dot( int m, int n, double** A, double* x, double* out ){
	for (int i=0; i<m; i++ ){
		out[i] = dot( n, A[i], x );
	} 
}

void dotT( int m, int n, double** A, double* x, double* out ){
	for (int i=0; i<m; i++ ){
		double doti = 0;
		for (int j=0; j<n; j++ ){ doti += A[j][i] * x[j];	}
		out[i] = doti;
	} 
}

inline void set( int m, int n, double   f,             double** out ){  	for (int i=0; i<m; i++ ){ set( n, f,       out[i] );  } }
inline void add( int m, int n, double   f, double** B, double** out ){  	for (int i=0; i<m; i++ ){ add( n, f, B[i], out[i] );  } }
inline void mul( int m, int n, double   f, double** B, double** out ){  	for (int i=0; i<m; i++ ){ mul( n, f, B[i], out[i] );  } }

inline void set( int m, int n, double** A,             double** out ){  	for (int i=0; i<m; i++ ){ set( n, A[i],       out[i] );   } }
inline void add( int m, int n, double** A, double** B, double** out ){  	for (int i=0; i<m; i++ ){ add( n, A[i], B[i], out[i] );   } }
inline void sub( int m, int n, double** A, double** B, double** out ){  	for (int i=0; i<m; i++ ){ sub( n, A[i], B[i], out[i] );   } }
inline void mul( int m, int n, double** A, double** B, double** out ){  	for (int i=0; i<m; i++ ){ mul( n, A[i], B[i], out[i] );   } }
inline void div( int m, int n, double** A, double** B, double** out ){  	for (int i=0; i<m; i++ ){ div( n, A[i], B[i], out[i] );   } }



void mmul_ik_kj( int ni, int nj, int nk, double** A, double** B, double** out ){
	for(int i=0; i<ni; i++){
		for(int j=0; j<nj; j++){
			double dotij = 0;
			for(int k=0; k<nk; k++){	dotij += A[i][k] * B[k][j];	}
			out[i][j] = dotij;
		}
	}
}

void mmul_ik_jk( int ni, int nj, int nk, double** A, double** B, double** out ){
	for(int i=0; i<ni; i++){
		for(int j=0; j<nj; j++){
			double dotij = 0;
			for(int k=0; k<nk; k++){	dotij += A[i][k] * B[j][k];	}
			out[i][j] = dotij;
		}
	}
}

void mmul_ki_kj( int ni, int nj, int nk, double** A, double** B, double** out ){
	for(int i=0; i<ni; i++){
		for(int j=0; j<nj; j++){
			double dotij = 0;
			for(int k=0; k<nk; k++){	dotij += A[k][i] * B[k][j];	}
			out[i][j] = dotij;
		}
	}
}

void mmul_ki_jk( int ni, int nj, int nk, double** A, double** B, double** out ){
	for(int i=0; i<ni; i++){
		for(int j=0; j<nj; j++){
			double dotij = 0;
			for(int k=0; k<nk; k++){	dotij += A[k][i] * B[j][k];	}
			out[i][j] = dotij;
		}
	}
}



void random_matrix( int m, int n, double xmin, double xmax, double** out ){
	double xrange = xmax - xmin;
	for (int i=0; i<m; i++ ){  random_vector ( n, xmin, xmax, out[i] ); } 
}

void print_matrix( int m, int n, double ** A ){
	for (int i=0; i<m; i++ ){	print_vector( n, A[i] );	} 
}

// ===============================================
// ======                                  =======
// ======              Qadric              =======
// ======                                  =======
// ===============================================

//  Q =  Sum_i |a_i> k_i <a_i|
void makeQuadricFormMatrix( int m, int n, double * ks, double ** A, double ** Q ){
	for (int im=0; im<m; im++ ){
		double kim   = ks[im];
		double * Aim =  A[im]; 
		for (int i=0; i<n; i++ ){
			for (int j=0; j<i; j++ ){
				double Qij = kim*Aim[i]*Aim[j];
				Q[i][j] += Qij;
				Q[j][i] += Qij;
			}
			Q[i][i] += kim*Aim[i]*Aim[i];
		}
	}
}

//  y =  < x | Q.x >
double evalQudraticForm( int n, double* x, double** Q ){
	double y =  0;
	for (int i=0; i<n; i++ ){	y += x[i] * dot( n, Q[i], x );   } 
	return y;
}

//  y = Sum_i k_i < a_i | x >^2
double evalQudraticFormDirs( int m, int n, double* x, double* k, double** A ){
	double y =  0;
	for (int i=0; i<m; i++ ){
		double ri = dot( n, A[i], x );
		y += k[i]*ri*ri;
	} 
	return y;
}


// ===============================================
// ======                                  =======
// ======       Gauss Elimination          =======
// ======                                  =======
// ======     Linsolve, invertMatrix       =======
// ======                                  =======
// ===============================================


// GaussElimination
// Method to carry out the partial-pivoting Gaussian
// elimination.  Here index[] stores pivoting order.

void GaussElimination( int n, double ** A, double * c, int * index ) {

	// Initialize the index
	for (int i=0; i<n; ++i) index[i] = i;

	// Find the rescaling factors, one from each row
	for (int i=0; i<n; ++i) {
	  double c1 = 0;
	  for (int j=0; j<n; ++j) {
		double c0 = fabs( A[i][j] );
		if (c0 > c1) c1 = c0;
	  }
	  c[i] = c1;
	}

	// Search the pivoting element from each column
	int k = 0;
	for (int j=0; j<n-1; ++j) {
	  double pi1 = 0;
	  for (int i=j; i<n; ++i) {
		double pi0 = fabs( A[ index[i] ][j] );
		pi0 /= c[ index[i] ];
		if (pi0 > pi1) {
		  pi1 = pi0;
		  k = i;
		}
	  }

	  // Interchange rows according to the pivoting order
	  int itmp = index[j];
	  index[j] = index[k];
	  index[k] = itmp;
	  for (int i=j+1; i<n; ++i) {
		double pj = A[ index[i] ][ j ]/A[ index[j] ][j];

	   // Record pivoting ratios below the diagonal
		A[ index[i] ][j] = pj;

	   // Modify other elements accordingly
		for (int l=j+1; l<n; ++l)
		  A[ index[i] ][l] -= pj*A[ index[j] ][l];
	  }
	}
}

void linSolve_gauss( int n, double ** A, double * b, int * index, double * x ) {

	// Transform the matrix into an upper triangle
	GaussElimination( n, A, x, index);

	// Update the array b[i] with the ratios stored
	for(int i=0; i<n-1; ++i) {
		for(int j =i+1; j<n; ++j) {
			b[index[j]] -= A[index[j]][i]*b[index[i]];
		}
	}

	// Perform backward substitutions
	x[n-1] = b[index[n-1]]/A[index[n-1]][n-1];
	for (int i=n-2; i>=0; --i) {
		x[i] = b[index[i]];
		for (int j=i+1; j<n; ++j) {
			x[i] -= A[index[i]][j]*x[j];
		}
		x[i] /= A[index[i]][i];
	}
}

/*
public static double[][] InvertMatrix_Gauss( int n, double ** A, int * index, double ** invA ) {
	int n = A.length;
	double x[][] = new double[n][n];
	double b[][] = new double[n][n];
	int index[] = new int[n];
	for (int i=0; i<n; ++i) b[i][i] = 1;

	// Transform the matrix into an upper triangle
	gaussian(A, index);

	// Update the matrix b[i][j] with the ratios stored
	for (int i=0; i<n-1; ++i)
	  for (int j=i+1; j<n; ++j)
		for (int k=0; k<n; ++k)
		  b[index[j]][k] -= A[index[j]][i]*b[index[i]][k];

	// Perform backward substitutions
	for (int i=0; i<n; ++i) {
	  x[n-1][i] = b[index[n-1]][i]/A[index[n-1]][n-1];
	  for (int j=n-2; j>=0; --j) {
		x[j][i] = b[index[j]][i];
		for (int k=j+1; k<n; ++k) {
		  x[j][i] -= A[index[j]][k]*x[k][i];
		}
		x[j][i] /= A[index[j]][j];
	  }
	}
	return x;
}
*/


void linSolve_CG( int n, double ** A, double * b, double * x ){
	const int maxIters   = 10;
	const double maxErr2 = 1e-5;
	double *  r     = new double[n];
	double *  r2    = new double[n];
	double *  p     = new double[n];
	double *  Ap    = new double[n];
	dot( n, n, A, x, r );
	sub( n, b, r, r );
	set( n, r, p );
	double rho = dot(n, r,r);
	double alpha = 0;
	for ( int i =0; i<maxIters; i++) {
		dot( n, n, A, p, Ap);
		alpha = rho / dot(n, p, Ap);
		fma( n, x, p ,  alpha,   x );
		fma( n, r, Ap, -alpha,   r2 );
		double err2 = dot(n, r2,r2);
		//printf( " iter: %i  err2: %f |  alpha %f \n", i, err2,     alpha );
		printf( " iter: %i  err2: %f \n", i, err2 );
		if (err2 < maxErr2 ) break;
		double rho2 = dot(n, r2,r2);
		double beta = rho2 / rho;
		fma( n, r2, p, beta, p );
		rho = rho2;
		double * swap = r; r = r2; r2 = swap;
	}
	delete r;
	delete r2;
	delete p;
	delete Ap;
}

void linSolve_BCG( int n, int m, double ** A, double * b, double * x ){
	const int maxIters   = 10;
	const double maxErr2 = 1e-5;
	double *  r     = new double[n];
	double *  r2    = new double[n];
	double *  p     = new double[n];
	double *  Ap    = new double[n];
	dot ( n, m, A, x, r );
	dotT( m, n, A, x, r );
	sub( n, b, r, r );
	set( n, r, p );
	double rho = dot(n, r,r);
	double alpha = 0;
	for ( int i =0; i<maxIters; i++) {
		dot ( n, m, A, p, Ap);
		dotT( m, n, A, p, Ap);
		alpha = rho / dot(n, p, Ap);
		fma( n, x, p ,  alpha,   x );
		fma( n, r, Ap, -alpha,   r2 );
		double err2 = dot(n, r2,r2);
		printf( " iter: %i  err2: %f \n", i, err2 );
		if (err2 < maxErr2 ) break;
		double rho2 = dot(n, r2,r2);
		double beta = rho2 / rho;
		fma( n, r2, p, beta, p );
		rho = rho2;
		double * swap = r; r = r2; r2 = swap;
	}
	delete r;
	delete r2;
	delete p;
	delete Ap;
}





void leastSquareFit_Gauss( int n, int m, double ** A, double * b, double * x ){
	double ** AA    = new_matrix( n, n );
	double *  Ab    = new double[n];
	int    *  index = new int   [n];
	dotT( n, m, A, b, Ab );
	mmul_ki_kj( n, n, m, A, A, AA );
	//linSolve_gauss( n, AA, Ab, index, x );
	linSolve_CG( n, AA, Ab, x );
	print_vector( n, x );
	delete AA;
	delete Ab;
	delete index;
}




