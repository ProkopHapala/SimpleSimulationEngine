
#ifndef  Lingebra_h
#define  Lingebra_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "fastmath.h"
#include "VecN.h"


void (dotFunc)(int nx, int ny, double * x, double * y );


//#define IJ(i,j) n*i+j

namespace Lingebra{

	// ==== function declarations

    //double ** from_continuous( int m, int n, double *p );
    double ** from_continuous( int m, int n, double *p, double ** A=0 );
	double ** new_matrix     ( int m, int n );
	void      delete_matrix  ( int m,        double** A );
	void transpose           ( int m, int n, double** A, double** TA );
	void dot                 ( int m, int n, double** A, double* x, double* out );
	void dotT                ( int m, int n, double** A, double* x, double* out );

    void symCopy             ( int m, double* A, bool toLower );


	void mmul_ik_kj( int ni, int nj, int nk, double** A, double** B, double** out );
	void mmul_ik_jk( int ni, int nj, int nk, double** A, double** B, double** out );
	void mmul_ki_kj( int ni, int nj, int nk, double** A, double** B, double** out );
	void mmul_ki_jk( int ni, int nj, int nk, double** A, double** B, double** out );

	void random_matrix( int m, int n, double xmin, double xmax, double** out );
	void print_matrix( int m, int n, double ** A );

	void   makeQuadricFormMatrix( int m, int n, double * ks, double ** A, double ** Q );
	double evalQudraticForm( int n, double* x, double** Q );
	double evalQudraticFormDirs( int m, int n, double* x, double* k, double** A );

	void GaussElimination    ( int n, double ** A, double * c, int * index );
	void linSolve_gauss      ( int n, double ** A, double * b, int * index, double * x );
	void linSolve_CG         ( int n, double ** A, double * b, double * x );
	void linSolve_BCG        ( int n, int m, double ** A, double * b, double * x );
	void leastSquareFit_Gauss( int n, int m, double ** A, double * b, double * x );

	int eig_Jacobi( int n, double* A, double* V, double* es, double tol, int nMaxIter  );
	double eig_Jacobi_init( int n, double* A, double* V, int* mjs, int& imax, int& jmax );
	void eig_Jacobi_step  ( int n, double* A, double* V, int* mjs, int& imax, int& jmax, double& vmax );

	// ==== inline function

	inline void set( int m, int n, double   f,             double** out ){  for (int i=0; i<m; i++ ){ VecN::set( n, f,       out[i] );  } }
	inline void add( int m, int n, double   f, double** B, double** out ){  for (int i=0; i<m; i++ ){ VecN::add( n, f, B[i], out[i] );  } }
	inline void mul( int m, int n, double   f, double** B, double** out ){  for (int i=0; i<m; i++ ){ VecN::mul( n, f, B[i], out[i] );  } }

	inline void set( int m, int n, double** A,             double** out ){  for (int i=0; i<m; i++ ){ VecN::set( n, A[i],       out[i] );   } }
	inline void add( int m, int n, double** A, double** B, double** out ){  for (int i=0; i<m; i++ ){ VecN::add( n, A[i], B[i], out[i] );   } }
	inline void sub( int m, int n, double** A, double** B, double** out ){  for (int i=0; i<m; i++ ){ VecN::sub( n, A[i], B[i], out[i] );   } }
	inline void mul( int m, int n, double** A, double** B, double** out ){  for (int i=0; i<m; i++ ){ VecN::mul( n, A[i], B[i], out[i] );   } }
	inline void div( int m, int n, double** A, double** B, double** out ){  for (int i=0; i<m; i++ ){ VecN::div( n, A[i], B[i], out[i] );   } }

	//template<void dotFunc(int nx, int ny, double* x, double* y)>
	template<typename F>
    void genLinSolve_CG(  int n, double * b, double * x , F dotFunc, int maxIters = 10 ){
        //const int    maxIters   = 10;
        const double    maxErr2 = 1e-5;
        double *  r     = new double[n];
        double *  r2    = new double[n];
        double *  p     = new double[n];
        double *  Ap    = new double[n];
        dotFunc( n, n, x, r );
        VecN::sub( n, b, r, r );
        VecN::set( n, r, p );
        double rho = VecN::dot(n, r,r);
        double alpha = 0;
        for ( int i =0; i<maxIters; i++) {
            dotFunc( n, n, p, Ap);

            //printf("p  "); VecN::print_vector(n, p);
            //printf("Ap "); VecN::print_vector(n, Ap);

            alpha = rho / VecN::dot(n, p, Ap);
            VecN::fma( n, x, p ,  alpha,   x );
            VecN::fma( n, r, Ap, -alpha,   r2 );
            double err2 = VecN::dot(n, r2,r2);
            //printf( " iter: %i  err2: %f |  alpha %f \n", i, err2,     alpha );
            //printf( " iter: %i  err2: %f \n", i, err2 );
            if (err2 < maxErr2 ) break;
            double rho2 = VecN::dot(n, r2,r2);
            double beta = rho2 / rho;
            VecN::fma( n, r2, p, beta, p );
            rho = rho2;
            double * swap = r; r = r2; r2 = swap;
        }
        delete [] r;
        delete [] r2;
        delete [] p;
        delete [] Ap;
    }

};


//template<void dotFunc(int nx, int ny, double* x, double* y)>
class LinSolver{ public:
    //static constexpr int       maxIters = 10;
    //static constexpr double    maxErr2  = 1e-5;
    int n;
    // temp
    int    istep=0;
    double *  r=0;
    double *  r2=0;
    double *  p=0;
    double *  Ap=0;
    double rho = 0;
    double alpha = 0;

    // to solve
    double* x = 0;
    double* b = 0;
    double* M = 0;

    virtual void dotFunc( int n, double * x, double * Ax )=0;

    void realloc(int n_, bool bx=false){
        if(n_!=n){
            n = n_;
            _realloc(r,n);
            _realloc(r2,n);
            _realloc(p,n);
            _realloc(Ap,n);
            if(bx){
                _realloc(x,n);
                _realloc(b,n);
            }
        }
    };

    void dealloc(){
        delete [] r;
        delete [] r2;
        delete [] p;
        delete [] Ap;
    };

    void setLinearProblem(int n_, double* x_, double* b_, double* M_ = 0 ){
        realloc(n_);
        x=x_; b=b_; M=M_;
    }

    double step_GD(double dt){
        dotFunc  ( n, x, r );
        VecN::sub( n, b, r, r );
        //VecN::add( n, b, r, r );
        VecN::fma( n, x, r, dt, x );
        return VecN::dot(n, r,r);
    }

    double step_CG(){
        // see https://en.wikipedia.org/wiki/Conjugate_gradient_method
        //printf( "LinSolver::step_CG %i \n", istep );
        if(istep==0){
            printf( "STEP1 START \n" );
            dotFunc  ( n, x, r );
            printf("r   "); VecN::print_vector(n, r);
            printf("b   "); VecN::print_vector(n, b);
            VecN::sub( n, b, r, r ); // r = b - A*x
            printf("r_  "); VecN::print_vector(n, r);
            VecN::set( n, r, p );    // p = r
            rho = VecN::dot(n, r,r);
            alpha = 0;
            printf( "STEP1 END rho %f alpha %f \n", rho, alpha );
        }else{
            double rho2 = VecN::dot(n, r2,r2);
            double beta = rho2 / rho;
            VecN::fma( n, r2, p, beta, p );
            rho = rho2;
            double * tmp = r; r = r2; r2 = tmp;
        }
        // NOTE : BCQ can be done if (A.T()*A) is applied instead of A in dotFunc
        printf("p  "); VecN::print_vector(n, p);
        dotFunc( n, p, Ap);
        printf("Ap "); VecN::print_vector(n, Ap);
        alpha = rho / VecN::dot(n, p, Ap);    // a  = <r|r>/<p|A|p>
        printf( "rho %f alpha %f \n", rho, alpha );
        VecN::fma( n, x, p ,  alpha,   x );   // x  = x - a*p
        VecN::fma( n, r, Ap, -alpha,   r2 );  // r2 = r - a*A|p>
        double err2 = VecN::dot(n, r2,r2);
        printf( "err2 %f \n", err2 );
        istep++;
        return err2;
        //printf( " iter: %i  err2: %f |  alpha %f \n", i, err2,     alpha );
        //printf( " iter: %i  err2: %f \n", i, err2 );
        //if (err2 < maxErr2 ) break;
    }

    // Kd  = np.dot(K,d)
    // dt  = np.dot(d,f) / np.dot(d,Kd)    # step length such that f is orthogonal to d 
    // x   = x + d *dt
    // f   = f - Kd*dt
    // f2  = np.dot(f,f)
    // d   = f + d*(f2/f2old)
    // return x, f, d, f2

    double step0_CG(){
        // f  = f0 - np.dot(K,x)
        // d  = f.copy()
        // f2 = np.dot(f,f)
        dotFunc  ( n, x, r );    // r = A*x
        VecN::sub( n, b, r, r ); // r = b - A*x
        VecN::set( n, r, p    ); // p = r
        double err2 = VecN::dot(n, r,r);
        rho = err2;
        //printf("[_]x :"); VecN::print_vector(n, x);
        //printf("[_]f0:"); VecN::print_vector(n, b);
        //printf("[_]f :"); VecN::print_vector(n, r);
        //printf("[_]d :"); VecN::print_vector(n, p);
        return err2;
    }

    double step_CG_simple(){
        // NOTE : BCQ can be done if (A.T()*A) is applied instead of A in dotFunc
        //printf("p  "); VecN::print_vector(n, p);
        dotFunc( n, p, Ap);   // Kd = K*d
        //printf("Ap "); VecN::print_vector(n, Ap);
        
        //printf("[%i]r :",istep-1); VecN::print_vector(n, r );
        //printf("[%i]p :",istep-1); VecN::print_vector(n, p );
        //printf("[%i]Ap:",istep-1); VecN::print_vector(n, Ap);

        //alpha = VecN::dot(n, r, p) / VecN::dot(n, p, Ap);
        alpha = rho / VecN::dot(n, p, Ap);   // dt  = dot(d,f) / dot(d,Kd)    # step length such that f is orthogonal to d
        //printf( "### CG_step %i dt=%g rho=%g \n", istep, alpha, rho );
        VecN::fma( n, x, p ,  alpha,   x );  // x = x + d  * dt;
        VecN::fma( n, r, Ap, -alpha,   r );  // f = f - Kd * dt;
        //printf("[%i]Kd:",istep-1); VecN::print_vector(n, Ap);
        //printf("[%i]x :",istep  ); VecN::print_vector(n, x);
        //printf("[%i]f :",istep  ); VecN::print_vector(n, r);
        double err2 = VecN::dot(n, r,r);
        //printf( "err2 %f \n", err2 );
        double beta = err2 / rho;
        VecN::fma( n, r, p, beta,   p );                       // d   = f + d*(f2/f2old)
        //printf("[%i]d :",istep); VecN::print_vector(n, p);
        istep++;
        rho = err2;
        return err2;
        //printf( " iter: %i  err2: %f |  alpha %f \n", i, err2,     alpha );
        //printf( " iter: %i  err2: %f \n", i, err2 );
        //if (err2 < maxErr2 ) break;
    }

    void solve_CG( int maxIters, double maxErr2 ){
        istep = 0;
        for ( int i =0; i<maxIters; i++) {
            if ( step_CG() < maxErr2 ) break;
        }
    }

    int solve_CG_( int maxIters, double maxErr2 ){
        istep = 0;
        step0_CG();
        int i=0;
        for ( i =0; i<maxIters; i++) {
            //printf( "[%i]========\n", istep );
            double err2 = step_CG_simple();
            //printf("CG[%i]x :",istep  ); VecN::print_vector(n, x);
            //printf( "[%i]err2 %g maxErr2 %g \n", istep, err2, maxErr2 );
            if ( err2< maxErr2 ){ 
                //printf("### CG converged at step %i \n", istep );
                break;
            }
        }
        return i;
    }

};


#endif

