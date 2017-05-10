#ifndef  TestFuncND_h
#define  TestFuncND_h

#include "Vec2.h"
#include "VecN.h"

class TestFuncND{
    public:
    int ndim  = 0;
    int nfreq = 0;
    double  freq0 = 1.0;
    double* stiffness = NULL;
    Vec2d*  coefs     = NULL;

    double evalX0( const Vec2d& cs, Vec2d * coefsi ){
        double x0 = 0.0d;
        Vec2d csn = cs;
        for(int k=0; k<nfreq; k++){
            //x0 += csn.x*coefsi[k].x + csn.y*coefsi[k].y;
            x0 += csn.dot( coefsi[k] );
            csn.mul_cmplx(cs);
        }
        return x0;
    }

    // TODO : xs can be translated and multiplied by unitary matrix to remove axial bias
    double eval( double * xs ){
        Vec2d cs; cs.fromAngle( xs[0]*freq0 );
        double E = 0.0d;
        for(int i=1; i<ndim; i++){
            double x0 = evalX0( cs, coefs + i*nfreq );
            double dx = xs[i]-x0;
            E += stiffness[i]*dx*dx;
        }
        return E;
    }

    void getNodeAt(double t, double * x0s){
        Vec2d cs; cs.fromAngle( t*freq0 );
        for(int i=1; i<ndim; i++){
            x0s[i] = evalX0( cs, coefs + i*nfreq );
        }
    }

    void setRandomStiffness( double kmin, double kmax ){ VecN::random_vector ( ndim, kmin, kmax, stiffness ); }
    void setRandomCoefs( double * scaling, double fsin, double fcos ){
        for( int i=0; i<ndim; i++ ){
            //VecN::random_vector ( nfreq*2, , , (double*)(coefs + i*nfreq) );
            Vec2d * coefsi = coefs + i*nfreq;
            for(int k=0; k<nfreq; k++){
                double sc = 2.0*scaling[k];
                coefsi[k].set( (randf()-0.5f)*sc*fcos, (randf()-0.5f)*sc*fsin );
            }
        }
    }

    void allocate(int ndim_, int nfreq_){
        ndim=ndim_; nfreq=nfreq_;
        stiffness = new double[ndim];
        coefs     = new Vec2d[ndim*nfreq];
    }

    void dealocate(){
        if(stiffness) delete [] stiffness;
        if(coefs)     delete [] coefs;
    }
};

#endif
