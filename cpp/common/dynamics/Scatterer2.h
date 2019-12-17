
#ifndef Scatterer2_h
#define Scatterer2_h

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


struct HalfChannel{
    double angWidth; // angular spread of the channel, how large the complementary elements seem from this point of view
    double thick;    // thickness of end scattere from this direction
    double fluxIn;
    double fluxOut;
    // area/(4*pi*distance^2)

    HalfChannel()=default;
    HalfChannel(double angWidth_, double thick_):angWidth(angWidth_),thick(thick_),fluxIn(0),fluxOut(0){};

    inline void influx  ( double val ){ fluxIn += val; }
    inline void transfer( double val ){ fluxOut = val*angWidth; }
};

struct Channel{
    //int i,j; // scattering elements between which the flux is
    Vec3d  dir;
    //double flux;
    //double Kgeom; //  through space gemetric coupling |  ToDo : how it should be normalized? How to choose angular width ?
    HalfChannel ends[2];

    inline void transfer(){
        ends[1].transfer(ends[0].fluxIn);
        ends[0].transfer(ends[1].fluxIn);
    }

    Channel()=default;
    Channel(const Vec3d& dir_, double w0, double w1, double th0, double th1 ):dir(dir_),ends{HalfChannel(w0,th0),HalfChannel(w1,th1)}{};

};


struct ScatterElem{

// the Scatterer scatter flux between rays attached to the node
// scheme 2 - assume there is external matrix "rays" which store flux along

    int isurf;

    Vec3d pos;     // position in space
    Mat3d rot;     // rotation of main axes (elipsoide maping)
    Vec3d thicks;  // thickness along main axes
    Vec3d areas;   // crossection area along mains axes
    double beta;   // decay exponent

    int ichan    ,nchan;
    int ibackChan,nbackChan;   // indexes in array of channels


    //double project(const Vec3d& dir, const Vec3d& property){
    //    rot.dot
    //    return property.dot(V);
    //}

    inline double scatterAmp( const Vec3d& h0, const  Vec3d h1 ) const {
        // this scattering function has to be normalized so that sum of flux scattered oto all channels is = 1

        // to-do - this can be precomputed for each ray
        Vec3d cs0,cs1;
        rot.dot_to( h0, cs0); // decomposition of direction 1 in eigen direction of scatterer
        rot.dot_to( h1, cs1); // ---,,---      of direstion 2
        double thick0 = cs0.dot(thicks); // ToDo: this can be precomputed in half-channel
        double thick1 = cs1.dot(thicks);
        double area0  = cs0.dot(areas );
        double area1  = cs1.dot(areas );

        //double thick  = thick0 + thick1;
        double cosa     = h0.dot(h1);
        return exp( beta*cosa*(thick0 + thick1) ); // ToDo: scattering probabilities should depend on width of channels between which it scatters
    }

    void scatterBiDir( Channel& chi, Channel& chj, int iend, int jend ){
        double amp = scatterAmp( chi.dir, chj.dir );
        chj.ends[jend].influx( amp * chi.ends[iend].fluxOut );
        chi.ends[iend].influx( amp * chj.ends[jend].fluxOut );
        //chj.ends[jend].fluxIn += amp * chi.ends[iend].fluxOut;
        //chi.ends[iend].fluxIn += amp * chj.ends[jend].fluxOut;
    }

};



class Scattering : public TriangleRayTracer, public LinSolver { public:

    std::vector<ScatterElem> elems;
    std::vector<Channel>     channels;
    int*       backChans = 0; // backward maping of channles

    virtual void dotFunc( int n, double * x, double * Ax ) override {}

    void makeChannles(){
        int n = elements.size();
        for(int i=0; i<n; i++){
            ScatterElem& eli = elems[i];
            for(int j=0; j<i; j++){
                ScatterElem& elj = elems[j];

                Vec3d   d = elj.pos - eli.pos;
                double  r = d.normalize();

                double occlusion = getOcclusion( eli.pos, d, r, eli.isurf, elj.isurf );

                if( occlusion <0.5 ){

                    Vec3d csi,csj;
                    eli.rot.dot_to( d, csi );
                    elj.rot.dot_to( d, csj );
                    double thi = csi.dot( eli.thicks );
                    double thj = csj.dot( elj.thicks );
                    double Si  = csi.dot( eli.areas );
                    double Sj  = csj.dot( elj.areas );
                    double invS = 1/(4*M_PI*r*r+Si+Sj);
                    channels.push_back( Channel( d, Si*invS, Si*invS, thi,thj ) );

                }

            }
        }
    }


    inline void step_Direct(){
        // ---- scattering phase
        for(int k=0; k<n; k++){ // node i
            ScatterElem& elem = elems[k];
            int jmin =      elem.ichan;
            int jmax = jmin+elem.nchan;
            int imin =      elem.ibackChan;
            int imax = imin+elem.nbackChan;
            for(int i=imin; i<imax; i++){ // all out-going channels
                Channel& chi = channels[backChans[i]];
                for(int j=jmin; j<jmax; j++){ // all incomming channels
                    elem.scatterBiDir( chi, channels[j], 0, 1 );
                }
            }
        }
        // ---- transfer phase
        for(Channel& ch: channels){ ch.transfer(); }
    }


    /*
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

                double coupling = eli.geomCoupling( elj );
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

    */

};

#endif
