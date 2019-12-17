
#ifndef Radiosity_h
#define Radiosity_h

#include "Vec3.h"
#include "Mat3.h"
#include "geom3D.h"
#include "raytrace.h"
#include "CMesh.h"

#include "TriangleRayTracer.h"

#include "Lingebra.h"

// See:
// https://en.wikipedia.org/wiki/Radiosity_(computer_graphics)

class Radiosity : public TriangleRayTracer, public LinSolver { public:
    double couplingTrashold  = 1e-8;

    double* M=0;

    // work arrays
    double* vals=0;    // a
    double* sources=0; //

    inline void allocateWork(int n){
        _realloc(vals   ,n);
        _realloc(sources,n);
    }

    inline void prepare(){
        int n = elements.size();
        allocateWork( n );
        setLinearProblem( n, vals, sources );
    }

    inline int makeCouplingMatrix(){
        int n = elements.size();
        _realloc(M,n*n);
        for(int i=0; i<n; i++) M[i*n+i] = 0.0;
        //for(int i=0; i<n; i++) M[i*n+i] = 1.0;
        int nvalid = 0;
        for(int i=0; i<n; i++){
            SurfElement& eli = elements[i];
            for(int j=0; j<i; j++){
                SurfElement& elj = elements[j];
                //if( eli.isurf == elj.isurf ){ M[i*n+j]=0.0; continue; }

                Vec3d  d = elj.pos - eli.pos;
                //double r = d.normalize();
                double r2 = d.normalize();
                double coupling = d.dot(eli.normal)*d.dot(elj.normal)/r2;
                if( fabs(coupling) < couplingTrashold ){ M[i*n+j]=0.0; continue; };
                coupling /= ( r2 + eli.area + elj.area );


                //coupling *=  eli.area * elj.area;
                coupling *= 2* eli.area/(4*M_PI);

                double occlusion = getOcclusion( eli.pos, d, sqrt(r2), eli.isurf, elj.isurf );

                coupling*=(1-occlusion);

                /*
                //glColor3f( occlusion, 0.0, 1-occlusion );
                //Draw::colorScale()
                //float c = sqrt(coupling)*0.3;
                float c = coupling*0.3;
                //glColor3f( c, c, c );
                glColor3f( c, 0.5, c );
                if(occlusion>0.5) glColor3f( 1.0, 0.0, 0.0 );
                Draw3D::drawLine( eli.pos, elj.pos );
                */

                coupling = fabs( coupling );
                if( occlusion <0.5 ){
                    Draw3D::drawLine( eli.pos, elj.pos );
                    //printf( " %i %i %g \n", i, j, coupling );
                }

                M[i*n+j] = coupling;
                M[j*n+i] = coupling;
                nvalid++;
            }
        }
        return nvalid;
    }

    inline void processTriangles( int ntri, Triangle3D* tris, double sz ){
        //nObstacles=ntri;
        //obstacles = tris;

        for(int i=0; i<ntri; i++){
            int nsamp = 5; // TODO :  automatially determine element size
            //trinagleToElements( {nsamp,nsamp}, tris[i], i );
            //trinagleToElements( nsamp, tris[i], i );
        }
        makeCouplingMatrix();
    }

    virtual void dotFunc( int n, double * x, double * Ax ) override {
        for(int i=0; i<n; i++){
            double Axi = 0.0;
            for(int j=0; j<n; j++){
                Axi += M[i*n+j] * x[j];
            }
            Ax[i]=Axi;
        }
    }

    inline void step_Direct(){
        for(int i=0; i<n; i++){
            double Axi = 0.0;
            for(int j=0; j<n; j++){
                Axi += M[i*n+j] * ( vals[j] + sources[j] );
            }
            vals[i] = Axi; // - sources[i];
        }
    }

};

#endif
