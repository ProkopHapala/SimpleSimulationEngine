
#ifndef Scatterer_h
#define Scatterer_h

#include "Vec3.h"
#include "Mat3.h"
#include "geom3D.h"
#include "raytrace.h"
#include "CMesh.h"

#include "TriangleRayTracer.h"

#include "Lingebra.h"

//#include "fastmath.h"

// See:
//
//
//   coupling probability
//  p_hit =  cross.dot(rd);
//  thick =  thick.dot(rd);
//  scatterOnSite( thick, freq, amp );  // store scattered rays into bins

class ScatterElement{ public:

    Vec3d pos;     //= Vec3dZero;

    int nfreq = 0; // number of frequeny(energy) spectrum subdivision
    int n     = 0; // number of subdivision per u,v axis
    int ndir  = 0;
    int ntot  = 0; // total number of rays stored
    //double* cross =0; // store projection of element crossection in different directions  ( to evaluate hit probability )
    //double* thick =0; // store projection of element thickness in different directions    ( to evaluate intercation probability )

    Vec2d dsamp;
    int ia=0,ib=0,iface=0;

    Mat3d cross;
    Mat3d thick;

    double* rays  =0; // store population of rays in particular direction

    // ==== Methods

    //Vec3d getDir(int idir){
    //    int iface = idir&7;
    //    int ip    = idir>>3;
    //};

    inline void getDir(int ia, int ib, int iface, Vec3d& h ){
        //Vec3d h;
        h.a = dsamp.a*ia;
        h.b = dsamp.b*ib;
        h.c = 1-h.a-h.b;
        h.octDir(iface);
        h.normalize_taylor3();
    };

    inline void getNextDir( Vec3d& h ){
        ia++;
        if(ia>=n){ia=0;ib++;}
        if(ib>=n){ib=0;iface++;}
        getDir(ia,ib,iface, h );
        //if(iface>=8){iface=0;iface++;}
        //idir++;
    };

    inline double scatterFunc(double cosa, double thick){
        // it is basically diffusion
        // exp(2*cos(x)-2)   ~=~   exp(-x^2)
        return exp( 2*(cosa-1)/thick );
    };

    void scatterCoefs_elastic( Vec3d dir, double thick, double Apre, double* coefs ){
        double sum = 0;
        ia=0;ib=0;iface=0;
        for(int i=0; i<ndir; i++){
            Vec3d   hi;
            getNextDir(hi);
            double cosa = dir.dot(hi);
            double amp  = scatterFunc(cosa,thick) * Apre;
            sum += amp;
            coefs[i]    += amp * Apre;
        }
        // ToDo : check sum == 1.0   !!!
    };

};

class Scattering : public TriangleRayTracer, public LinSolver { public:
    double couplingTrashold  = 1e-8;

    //int nObstacles=0;
    //Triangle3D* obstacles=0;

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

                coupling = fabs( coupling );
                if( occlusion <0.5 ){
                    Draw3D::drawLine( eli.pos, elj.pos );
                    //printf( " %i %i %g \n", i, j, coupling );
                }

                M[i*n+j] = coupling;
                M[j*n+i] = coupling;
            }
        }
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
