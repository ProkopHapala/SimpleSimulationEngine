
#ifndef  radial_splines_h
#define  radial_splines_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

//#include <Vec3.h>

class GeneralSpline{
    public:
    int n;
    double   m0,dm,invdm;
    double * coefs;

    inline double get_splineR2( double r2 ){
        double ir2 = 1/r2;
        double  r4 = r2*r2;
        //double  m  = ir2;
        //double  m_ = (m-m0)*invdm;
        //double  m_ = (m-m0);
        //int     i  = (n-1)-(int)(m_);
        //int     i4 = i<<2;
        //printf("%i %f %f\n", i, m, m_ );
        //printf("%i %f %f %f %f\n", i, coefs[i4], coefs[i4+1], coefs[i4+2], coefs[i4+3]  );
        //printf("%i %f %f %f %f\n", i, 1.0, r2, r4, ir2  );
        int     i4 = ((n-1)-int((ir2-m0)*invdm))<<2;
        return  coefs[i4] + coefs[i4+1]*r2 + coefs[i4+2]*r4 +  coefs[i4+3]*ir2;
    }

    inline int loadFromFile( const char * fname ){
        printf(" loading spline from: >>%s<<\n", fname );
        FILE * pFile;
        pFile = fopen(fname,"r");
        fscanf( pFile, " %i %lf %lf\n", &n, &dm, &m0 );
        printf(        " %i %lf %lf\n",  n,  dm,  m0 );
        invdm=1/dm;
        coefs = new double[4*n];
        for(int i=0; i<n; i++){
            int i4 = i<<2;
            fscanf ( pFile, "%lf %lf %lf %lf\n", &coefs[i4],&coefs[i4+1],&coefs[i4+2],&coefs[i4+3] );
            printf(         "%lf %lf %lf %lf\n",  coefs[i4], coefs[i4+1], coefs[i4+2], coefs[i4+3] );
        }
        fclose(pFile);
        return n;
    }

};

#endif



