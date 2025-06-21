#ifndef Scatterer_h
#define Scatterer_h

/// @file
/// @brief Defines the `Scattering` class, a matrix-based solver for flux transport using radiosity-like approach.

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



/*

Theory:
=======

We work with bunches of particles (photons, neutrons) with some distribution in phase space defined by position and direction (p,d)
There are two fundamental processes and therefore two views of the problem:
1) through space flight ( where the spatial size of bunch spread due ot its finite spread in angular space ), this is problem of geometric optics
2) On-site scattering where the ray is scattered from one direction to an other.

1) flux disribution in ements element








*/


/// @struct HalfRay
/// @brief Represents a ray channel with angular spread and effective thickness.
struct HalfRay{
    double width; /// @brief Angular spread of the channel, how large the complementary elements seem from this point of view
    // area/(4*pi*distance^2)
};

/// @struct FluxRay
/// @brief Represents a flux ray with direction, flux value, and geometric coupling.
struct FluxRay{
    //int i,j; // scattering elements between which the flux is
    HalfRay left;
    HalfRay right;
    Vec3d  dir;
    double flux;
    double Kgeom; //  through space gemetric coupling |  ToDo : how it should be normalized? How to choose angular width ?
};

/// @struct ScatterElem
/// @brief Represents a scattering element with position, rotation, thickness, and area.
struct ScatterElem{

// the Scatterer scatter flux between rays attached to the node
// scheme 2 - assume there is external matrix "rays" which store flux along

    Vec3d pos;     /// @brief Position in space
    Mat3d rot;     /// @brief Rotation of main axes (elipsoide maping)
    Vec3d thicks;  /// @brief Thickness along main axes
    Vec3d areas;   /// @brief Crossection area along mains axes
    double beta;   /// @brief Decay exponent


    /// @brief Projects a direction onto the scatterer's main axes.
    //double project(const Vec3d& dir, const Vec3d& property){
    //    rot.dot
    //    return property.dot(V);
    //}

    /// @brief Computes the scattering probability between two directions.
    double scatter( Vec3d h0, Vec3d h1, double  ){
        // this scattering function has to be normalized so that sum of flux scattered oto all channels is = 1

        // to-do - this can be precomputed for each ray
        Vec3d cs0,cs1;
        rot.dot_to( h0, cs0); // decomposition of direction 1 in eigen direction of scatterer
        rot.dot_to( h1, cs1); // ---,,---      of direstion 2
        double thick0 = cs0.dot(thicks);
        double thick1 = cs1.dot(thicks);
        double area0  = cs0.dot(areas );
        double area1  = cs1.dot(areas );

        //double thick  = thick0 + thick1;
        double cosa     = h0.dot(h1);
        return exp( beta*cosa*(thick0 + thick1) ); // ToDo: scattering probabilities should depend on width of channels between which it scatters
    }

};

/// @class ScatterElement
/// @brief Represents a scattering element with position, frequency, and direction subdivisions.
class ScatterElement{ public:

    Vec3d pos;     /// @brief Position of scattering element
    int nfreq = 0; // number of frequeny(energy) spectrum subdivision
    int n     = 0; // number of subdivision per u,v axis
    int ndir  = 0;
    int ntot  = 0; // total number of rays stored
    //double* cross =0; // store projection of element crossection in different directions  ( to evaluate hit probability )
    //double* thick =0; // store projection of element thickness in different directions    ( to evaluate intercation probability )

    Vec2d dsamp;   /// @brief Direction sampling intervals
    int ia=0,ib=0,iface=0; /// @brief Current direction indices
    Mat3d cross;
    Mat3d thick;

    double* rays  =0; /// @brief Array storing ray populations in different directions

    /// @brief Gets the direction of a ray given its indices.
    //Vec3d getDir(int idir){
    //    int iface = idir&7;
    //    int ip    = idir>>3;
    //};

    /// @brief Gets the direction of a ray given its indices.
    inline void getDir(int ia, int ib, int iface, Vec3d& h ){
        //Vec3d h;
        h.a = dsamp.a*ia;
        h.b = dsamp.b*ib;
        h.c = 1-h.a-h.b;
        h.octDir(iface);
        h.normalize_taylor3();
    };

    /// @brief Gets the next direction in the iteration.
    inline void getNextDir( Vec3d& h ){
        ia++;
        if(ia>=n){ia=0;ib++;}
        if(ib>=n){ib=0;iface++;}
        getDir(ia,ib,iface, h );
        //if(iface>=8){iface=0;iface++;}
        //idir++;
    };

    /// @brief Computes the scattering function value.
    inline double scatterFunc(double cosa, double thick){
        // it is basically diffusion
        // exp(2*cos(x)-2)   ~=~   exp(-x^2)
        return exp( 2*(cosa-1)/thick );
    };

    /// @brief Computes the scattering coefficients for elastic scattering.
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

/// @class Scattering
/// @brief A matrix-based solver for flux transport using radiosity-like approach.
class Scattering : public TriangleRayTracer, public LinSolver { public:
    double couplingTrashold  = 1e-8; /// @brief Minimum coupling value to consider in matrix construction

    //int nObstacles=0;
    //Triangle3D* obstacles=0;

    double* coupling =0; /// @brief Geometric coupling matrix between elements (distance dependent)
    double* rays     =0; /// @brief Flux stored in each coupling channel

    // work arrays
    double* vals=0;    /// @brief Current flux values (solution vector)
    double* sources=0; /// @brief Source terms for each element (right-hand side)

    /// @brief Allocates work arrays for the solver.
    inline void allocateWork(int n){
        _realloc(vals   ,n);
        _realloc(sources,n);
    }

    /// @brief Prepares the solver for computation.
    inline void prepare(){
        int n = elements.size();
        allocateWork( n );
        setLinearProblem( n, vals, sources );
    }

    /// @brief Constructs the coupling matrix M[i][j] = geomCoupling(el_i,el_j) accounting for occlusion.
    inline int makeCouplingMatrix(){
        int n = elements.size();
        _realloc(M,n*n);
        int nvalid = 0;
        for(int i=0; i<n; i++) M[i*n+i] = 0.0;
        //for(int i=0; i<n; i++) M[i*n+i] = 1.0;
        for(int i=0; i<n; i++){
            SurfElement& eli = elements[i];
            for(int j=0; j<i; j++){
                SurfElement& elj = elements[j];
                //if( eli.isurf == elj.isurf ){ M[i*n+j]=0.0; continue; }

                double coupling = eli.geomCoupling( elj );
                /*
                Vec3d  d = elj.pos - eli.pos;
                //double r = d.normalize();
                double r2 = d.normalize();
                double coupling = d.dot(eli.normal)*d.dot(elj.normal)/r2;
                if( fabs(coupling) < couplingTrashold ){ M[i*n+j]=0.0; continue; };
                coupling /= ( r2 + eli.area + elj.area );
                */


                //coupling *=  eli.area * elj.area;
                coupling *= 2* eli.area/(4*M_PI);

                Vec3d d = elj.pos - eli.pos;
                double r = d.normalize();
                double occlusion = getOcclusion( eli.pos, d, r, eli.isurf, elj.isurf );

                coupling*=(1-occlusion);

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

    /// @brief Direct propagation step: vals[i] = sum_j M[i][j]*(vals[j] + sources[j])
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
