
#ifndef SchroedingerLeapFrog_h
#define SchroedingerLeapFrog_h

//  Schroedinger eqution
//  dd Y    = (E-V) Y
//  d( dY ) = (E-V) Y

const double SCHR_CONST = 10;

inline double LJfunc( double ia ){
	double ia3 = ia*ia*ia;
	double ia6 = ia3*ia3;
	return ia6*ia6 - 2*ia6;
}

inline double potential( double x ){

	//double V = 3*( 2*LJfunc( 2/(x-4.0) ) + LJfunc( 2/(4.0+x) ) ) + 6;

	//double V = 5*( 2*LJfunc( 2/(x-4.0) ) ) + 10;
	//V = (x> -2.2)? V : 10.0 ;
	//V = ( V < 10.0)? V : 10.0;

	double V = -3*x;
	double x_ = x+4;
	//V += -20/(1+(x_*x_));
	V += -50.0*exp(-(1+x_*x_));

/*
	// quadric
	double V = x*x*x*x - 4*x*x + 4;
	V = ( V < 10.0)? V : 10.0;
*/
/*
	// box
	double V = -0.05;
	V = (x< 2.2)? V : 10.0 ;
	V = (x>-1.8)? V : 10.0 ;
*/
	return V;
}

inline void leapfrog_step( double dx, double V, double& y, double& dy  ){
	dy +=  y * SCHR_CONST * V * dx;
	y  += dy * dx;
}

void eval_potential( int n, double dx, double x0,
					 double * x_buff, double * V_buff
){
	double x  =  x0;
	for ( int i=0; i<n; i++ ){
		x_buff[i] = x;
		V_buff[i] = potential( x );
		//printf( " %i  :  %f %f \n", i, x, V_buff[i]  );
		x += dx;
	}
}

void integrate( int n, int i0, double dx, double y0, double dy0, double E, double* V_buff, double &yLeft, double &yRight ){
	// Right
	double y   =  y0;
	double dy  = dy0;
	for ( int i=i0+1; i<n; i++ ){
		double V = ( V_buff[i-1] - E );
		leapfrog_step( dx, V, y, dy  );
	}
	yRight = y;
	// Left
	dx = -dx;
	y   =  y0;
	dy  = dy0;
	for ( int i=i0-1; i>=0; i-- ){
		double V = ( V_buff[i+1] - E );
		leapfrog_step( dx, V, y, dy  );
	}
	yLeft = y;
}


void integrate_buff( int n, int i0, double dx, double y0, double dy0, double E,
				double* V_buff, double* y_buff, double* dy_buff
){
	// Right
	double y   =  y0;
	double dy  = dy0;
	 y_buff[i0] =  y;
	dy_buff[i0] = dy;
	for ( int i=i0+1; i<n; i++ ){
		double V = ( V_buff[i-1] - E );
		leapfrog_step( dx, V, y, dy  );
		 y_buff[i] =  y;
		dy_buff[i] = dy;
		//printf( " %i  : %f  %f %f \n", i, x, y, dy  );
	}
	// Left
	dx = -dx;
	y   =  y0;
	dy  = dy0;
	for ( int i=i0-1; i>=0; i-- ){
		double V = ( V_buff[i+1] - E );
		leapfrog_step( dx, V, y, dy  );
		 y_buff[i] =  y;
		dy_buff[i] = dy;
		//printf( " %i  : %f  %f %f \n", i, x, y, dy  );
	}
}

void evalResidua( int n, int i0, double dx, double E, double * V_buff, double& fb, double& ySqOut ){
	double yAL=0,yAR=0,yBL=0,yBR=0;
	integrate     ( n, i0, dx, 1.0, 0.0, E, V_buff, yAL, yAR );
	integrate     ( n, i0, dx, 0.0, 1.0, E, V_buff, yBL, yBR );
	double yABL    = yBL - yAL;
	double yABR    = yBR - yAR;
	fb  =  -( yABL*yAL + yABR*yAR ) / ( yABL*yABL + yABR*yABR );
	double fa  =  1.0-fb;
	double yL      =  fa*yAL + fb*yBL;
	double yR      =  fa*yAR + fb*yBR;
	ySqOut  =  yL*yL + yR*yR;
}


void scanResidua( int n, int i0, double dx,  int nE, double Emin, double Emax,
				  double * V_buff, double * E_buff, double * y2out_buff
){
	double dE = (Emax-Emin)/nE;
	for (int i=0; i<nE; i++){
		double E = i*dE;
		double y2out = 0, fb_best = 0;
		evalResidua( n, i0, dx, E, V_buff, fb_best, y2out );
		//y2out_buff[i] = y2out * 0.0000001;
		y2out_buff[i] = log( y2out )*0.05;
		E_buff[i] = E;
		printf( " E= %f y2out= %f  log(y2out) %f \n", E, y2out, y2out_buff[i]  );
	}
}

void lincomb( int i0, int i1, double fa, double fb, double * a, double * b, double * out ){
	for( int i=i0; i<i1; i++ ){ out[i] = fa*a[i] + fb*b[i]; }
}


double dot( int i0, int i1,	double * a, double * b ){
	double sum_ab = 0;
	for( int i=i0; i<i1; i++ ){ sum_ab += a[i] * b[i]; }
	return sum_ab;
}

#endif
