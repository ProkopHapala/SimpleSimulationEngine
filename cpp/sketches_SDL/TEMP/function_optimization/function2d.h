

typedef double (*Function1d)( double x);
typedef double (*DiffFunc1d)( double x, double& dx );


typedef double (*Function2d)( double x, double y);
typedef double (*DiffFunc2d)( double x, double y, double& dfdx, double& dfdy );


typedef double (*FunctionNd)( int n, double * xs );
typedef double (*DiffFuncNd)( int n, double * xs, double * dfs );



double sigmoideAbs( double x, double& dfdx ){
	double D  = 1/( 1 + fabs(x) );
	dfdx      = D*D; 
	return  x*D; 
}

double sigmoideSqrt( double x, double& dfdx ){
	double D2 = 1/(1 + x*x);
	double D  = sqrt( D2 );
	dfdx      = D*D2; 
	return  x*D; 
}

double lorenz( double x, double& dfdx ){
	double f = 1/(1 + x*x);
	dfdx     = 2*x*f*f; 
	return  f; 
}

double x2period( double x, double& dfdx ){
	int ix = (int)(x+1000.0)-1000;
	double dx = x-ix-0.5;
	if( ix&0x1 ){
		dfdx =   -2; 
		return 1-dx*dx;
	}else{
		dfdx =  2;
		return dx*dx-1;
	}
}



/*
int      badass_nkink;
int      badass_nfall;
double * badass_xkink;
double * badass_ykink;
double * badass_kkink;
double * badass_xfall;
double * badass_Efall;
double * badass_kfall;

double badassFunc( double x, double y ){

}
*/


double warp_x2period( double x, double y, double& dfdx, double& dfdy ){
	double cwarp = 0.5;
	for(int i=0; i<5; ){
		double dsx,dsy;
		double sx = x2period( x, dsx );
		double sy = x2period( y, dsy );
 		x = x + cwarp*sx;
		y = y + cwarp*sy;
	}
}
















double harmonic( double x, double y ){
	return x*x + y*y;
}

double rosenbrok( double x, double y ){
	double f = (x*x - y); 
	return f*f + x*x*0.1;
}

double sinValey( double x, double y ){
	double f = ( sin(x)  - y); 
	return f*f + x*x*0.1;
}

double cosValey( double x, double y ){
	double f = ( cos(4*x) - y); 
	return (0.2*f*f + x*x*0.05)*0.2;
}

double spiral( double x, double y ){
	double  phi = atan2 (y,x);
	double  r2  = x*x + y*y;
	double  r   = sqrt(r2);
	return  (1+sin( 2*M_PI*r + phi ))*0.1     + 0.05*r2;
}

/*
double mandelbort( double cX, double cY ){
	double x = 0, y = 0;
	int i=0;
	for (i=0; i<256; i++){
		double x_ =  x * x - y * y + cX;
        y         = 2.0 * x * y    + cY;
        x = x_;
		if( x*x + y*y > 4.0d ) break;
	}
	return i / 256.0;
}
*/

double mandelbort( double cX, double cY ){
	double x = 0, y = 0;
	int i=0;
	for (i=0; i<10; i++){
		double x_ =  x * x - y * y + cX;
        y         = 2.0 * x * y    + cY;
        x = x_;
	}
	double r = sqrt(x*x + y*y);
	//return 0.5*(sin(r)+1.0);
	return tanh(r);
}
