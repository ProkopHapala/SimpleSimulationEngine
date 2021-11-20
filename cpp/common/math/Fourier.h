
#ifndef  Fourier_h
#define  Fourier_h

#include <math.h>

//sincos_basis( ){}

// evaluate fourier coefficients of scalar function val=f(phi) given complex phase z=ca+i*sa=exp(i*phi)
inline void fourier_coef( double ca, double sa, double val, int n, double * coefs ){
	double cna=ca;
	double sna=sa;
	// should we use Vec2d ? perhaps not, to save unnecesary include
	coefs[0]+=cna*val;
	coefs[1]+=sna*val;
	for(int i=2; i<(n<<1); i+=2){
		double im=cna*ca-sna*sa;
		double re=cna*sa+sna*ca;
		coefs[i  ]+=re*val;
		coefs[i+1]+=im*val;
		cna=im; sna=re;
	}
}

// evaluate a scalar function val=f(phi) defined by fourier coefficients given complex phase z=ca+i*sa=exp(i*phi)
inline double fourier_eval( double ca, double sa, int n, double * coefs ){
	double cna=ca;
	double sna=sa;
	double result = cna*coefs[0] + sna*coefs[1];
	for(int i=2; i<(n<<1); i+=2){
		double im=cna*ca-sna*sa;
		double re=cna*sa+sna*ca;
		result += cna*coefs[i] + sna*coefs[i+1];
		cna=im; sna=re;
	}
	return result;
}

void genSinCosArray( int n, double * phis, double * ca, double * sa ){
	for( int i=0; i<n; i++ ){
		double phi = phis[i];
		ca[i]  = cos(phi);
		sa[i]  = sin(phi);
	}
}

void fourier_coef_array( int n, double * ca, double * sa, double * vals, int m, double * coefs ){
	for( int i=0; i<m; i++ ){ coefs[i]=0; }
	for( int i=0; i<n; i++ ){
		fourier_coef( ca[i], sa[i], vals[i], m, coefs );
	}
}

void fourier_eval_array( int n, double * ca, double * sa, double * vals, int m, double * coefs ){
	for( int i=0; i<n; i++ ){ vals[i] = fourier_eval( ca[i], sa[i], m, coefs ); }
}

// FFT/IFFT routine. (see pages 507-508 of Numerical Recipes in C)
// from https://www.google.cz/url?sa=t&rct=j&q=&esrc=s&source=web&cd=2&ved=0ahUKEwiymK34kIDSAhUrEJoKHbkIA10QFggjMAE&url=https%3A%2F%2Fsoftware.intel.com%2Fsites%2Fdefault%2Ffiles%2Fforum%2F392242%2Ffft.c&usg=AFQjCNHcT2Cx2XQYTuD_R2gUYFNme6vniQ&sig2=Tj6WyNS3SDbO5K4rSZoqTg
// data[2*i ]=Re data[2*i+1]=Im;   nn=2^pow;



// USE LIKE THIS:
//  Forward-FFT:   FFT( (double*)data, n,  1);    where n = 2^m
//  inverse-FFT:   FFT( (double*)data, n, -1);


#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
void FFT(double *data, unsigned long nn, int isign){
    unsigned long n,mmax,m,j,istep,i;
    double wtemp,wr,wpr,wpi,wi,theta;
    double tempr,tempi;

    n=nn << 1;
    j=1;
    for (i=1;i<n;i+=2) {
        if (j > i) {
            SWAP(data[j],data[i]);
            SWAP(data[j+1],data[i+1]);
        }
        m=n >> 1;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    mmax=2;
    while (n > mmax) {
        istep=mmax << 1;
        theta=isign*(2*M_PI/mmax);
        wtemp=sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi=sin(theta);
        wr=1.0;
        wi=0.0;
        for (m=1;m<mmax;m+=2) {
            for (i=m;i<=n;i+=istep) {
                j=i+mmax;
                tempr=wr*data[j]-wi*data[j+1];
                tempi=wr*data[j+1]+wi*data[j];
                data[j]=data[i]-tempr;
                data[j+1]=data[i+1]-tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr=(wtemp=wr)*wpr-wi*wpi+wr;
            wi=wi*wpr+wtemp*wpi+wi;
        }
        mmax=istep;
    }
    if( isign<0 ){
		double renorm=1./nn;
		for (int i=0; i<n; i+=2) {
			data[i  ]*=renorm;
			data[i+1]*=renorm;
			//ops++;
		}
	}
}
#undef SWAP



/*
int FFT(double * data, int nn, int isign){
    // see   https://stackoverflow.com/questions/2220879/c-numerical-recipies-fft
    //int n, mmax, m, j, istep, i;
    //double wtemp, wr, wpr, wpi, wi, theta;
    //double tempr, tempi;
    unsigned long n,mmax,m,j,istep,i;
    double wtemp,wr,wpr,wpi,wi,theta;
    double tempr,tempi;
    n=nn << 1;
    j=1;

	//int ops = 0;
    //int n = nn*2;
    int j = 1;
    //DEBUG
    for (int i = 1; i < n; i += 2) {
		if (j > i) {
			//double tempr;
			//tempr = data[j  ]; data[j  ] = data[i  ]; data[i  ] = tempr; // swap(i,j)
			//tempr = data[j+1]; data[j+1] = data[i+1]; data[i+1] = tempr;
			_swap(data[j  ],data[j  ]);
			_swap(data[j+1],data[j+1]);
			//ops++;
		}
		int m = n >> 1;
		while (m >= 2 && j > m) { j -= m; m >>= 1; }
		j += m;
    }
    int mmax = 2;
    //DEBUG
    while (n > mmax) {
		int istep = 2*mmax;
		double theta = (2*M_PI)/(isign*mmax);
		double wtemp = sin(0.5*theta);
		double wpr   = -2.0*wtemp*wtemp;
		double wpi   = sin(theta);
		double wr    = 1.0;
		double wi    = 0.0;
		//int iDEBUG = 0;
		for (int m = 1; m < mmax; m += 2) {
			for (int i = m; i <= n; i += istep) {
				j =i + mmax;
				double tempr = wr*data[j  ] - wi*data[j+1];
				double tempi = wr*data[j+1] + wi*data[j  ];
				data[j  ]    = data[i  ] - tempr;
				data[j+1]    = data[i+1] - tempi;
				data[i  ]   += tempr;
				data[i+1]   += tempi;
				//ops++;
				//iDEBUG++;
			}
			wr = (wtemp = wr)*wpr - wi*wpi + wr;
			wi = wi*wpr + wtemp*wpi + wi;
		}
		//printf( "iDEBUG %i | nn %i n %i mmax %i istep %i \n", iDEBUG, nn, n, mmax, istep );
		mmax = istep;
    }
    //DEBUG
	if( isign<0 ){
		double renorm=1./nn;
		for (int i=0; i<n; i+=2) {
			data[i  ]*=renorm;
			data[i+1]*=renorm;
			//ops++;
		}
	}
	return ops;
}
*/

void FFT_2D(int nnx, int nny, double * data, int isign, double* tmp = 0 ){
    bool bAllocTmp = (tmp == 0);
    int nx=(1<<nnx);
    int ny=(1<<nny);
    if(bAllocTmp){
        printf( "FFT_2D alloc tmp \n" );
        tmp = new double[nx*ny];
    }
    for(int iy=0; iy<ny; iy++)                          { FFT(data+iy*nx, nnx, isign);         }  // FFT along x-axis
    for(int ix=0; ix<nx; ix++)for(int iy=0; iy<ny; iy++){ tmp[ix*ny+iy] = data[iy*nx+ix];      }  // transpose data
    for(int ix=0; ix<ny; ix++)                          { FFT( tmp+ix*ny, nny, isign);         }  // FFT along y-axis
    if(bAllocTmp){
        for(int iy=0; iy<ny; iy++)for(int ix=0; ix<nx; ix++){ data[iy*nx+ix] = tmp[ix*ny+iy]; }  // transpose data back
        delete [] tmp;
    }
}

void FFT_3D(int nx, int ny, int nz, double * data, int isign){
    int nxy=nx*ny;
    int nxz=nx*nz;
    int nyz=ny*nz;
    double* tmp  = new double[nxy*nz];
    double* tmp2 = new double[nxy*nz];
    for(int iz=0; iz<nz; iz++)for(int iy=0; iy<ny; iy++)                          { FFT(data+ iz*nxy+iy*nx    , nx, isign);                 }  // FFT along x-axis     [x,y,z]
    for(int iz=0; iz<nz; iz++)for(int ix=0; ix<nx; ix++)for(int iy=0; iy<ny; iy++){     tmp  [iz*nxy+ix*ny+iy] = data[iz*nxy+iy*nx+ix]; }  // transpose data (x<->y)
    for(int iz=0; iz<nz; iz++)for(int ix=0; ix<ny; ix++)                          { FFT(tmp+  iz*nxy+ix*ny    , ny, isign);                 }  // FFT along y-axis     [y,x,z]
    for(int iy=0; iy<ny; iy++)for(int ix=0; ix<nx; ix++)for(int iz=0; iz<nz; iz++){     tmp2[ iy*nxz+ix*nz+iz] = tmp [iz*nxy+ix*ny+iy]; }  // transpose data (x<->y)
    for(int iy=0; iy<ny; iy++)for(int ix=0; ix<nx; ix++)                          { FFT(tmp2+ iy*nxz+ix*nz    , nx, isign);                 }  // FFT along y-axis     [z,x,y]
    for(int iz=0; iz<nz; iz++)for(int iy=0; iy<ny; iy++)for(int ix=0; ix<nx; ix++){ data[iz*nxy+iy*nx+ix]      = tmp2[iy*nxz+ix*nz+iz]; }  // transpose data back  [z,x,y] -> [x,y,z]
    delete [] tmp;
    delete [] tmp2;
}


#endif
