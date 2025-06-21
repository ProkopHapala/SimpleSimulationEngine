
#ifndef Scatterer2_h
#define Scatterer2_h
/// @file @brief Defines the `Scattering2` class, a direct iterative solver for flux transport using a network of channels.

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

/// @brief Represents one end of a flux transport channel.
/// @brief It stores incoming and outgoing flux, along with geometric properties like angular width and thickness from this viewpoint, and the index of the associated scattering element.
struct HalfChannel{
    
    double angWidth; /// @brief Angular spread of the channel, how large the complementary element seems from this point of view.
    double thick;    /// @brief Thickness of the scattering element at this end of the channel, from this direction.
    double fluxIn;   /// @brief Flux accumulated into this half-channel from scattering events.
    double fluxOut;  /// @brief Flux propagated out from this half-channel after a transfer step.
    int elemIdx;     /// @brief Index of the ScatterElem2 this half-channel belongs to.
    // area/(4*pi*distance^2)

    HalfChannel()=default;
    HalfChannel(double angWidth_, double thick_, int elemIdx_ = -1):angWidth(angWidth_),thick(thick_),fluxIn(0),fluxOut(0),elemIdx(elemIdx_){};

    /// @brief Accumulates incoming flux from a scattering event.
    inline void influx  ( double val ){ fluxIn += val; }
    /// @brief Calculates the outgoing flux to be transferred along the channel.
    inline void transfer( double val ){ fluxOut = val*angWidth; }
};

/// @brief Represents a bidirectional connection for flux transport between two scattering elements.
struct Channel{
    //int i,j; // scattering elements between which the flux is
    Vec3d  dir;
    //double flux;
    //double Kgeom; //  through space gemetric coupling |  ToDo : how it should be normalized? How to choose angular width ?
    HalfChannel ends[2];

    /// @brief Performs the flux transfer step, propagating flux from one end to the other.
    inline void transfer(){
        ends[1].transfer(ends[0].fluxIn);
        ends[0].transfer(ends[1].fluxIn);
    }

    Channel()=default; // Default constructor
    Channel(const Vec3d& dir_, double w0, double w1, double th0, double th1, int idx0, int idx1 ):dir(dir_),ends{HalfChannel(w0,th0,idx0),HalfChannel(w1,th1,idx1)}{};

};

/// @brief Represents a discrete scattering element in the simulation.
/// @brief It has anisotropic properties and is connected to other elements via a set of channels.
struct ScatterElem2{

// the Scatterer scatter flux between rays attached to the node
// scheme 2 - assume there is external matrix "rays" which store flux along

    int isurf;  /// @brief Surface ID, used for occlusion checks.    
    Vec3d pos;     /// @brief Position in space.
    Mat3d rot;     /// @brief Rotation of main axes, defining the orientation of anisotropic properties.
    Vec3d thicks;  /// @brief Thickness along main axes.
    Vec3d areas;   /// @brief Cross-section area along main axes.
    double beta;   /// @brief Decay exponent for the scattering kernel.

    /// @brief Index and count for forward and backward channels associated with this element.
    int ichan    ,nchan;
    int ibackChan,nbackChan;   // indexes in array of channels

    double fluxIn;  /// @brief Total flux accumulated by this scattering element.
    double fluxOut; /// @brief Total flux emitted by this scattering element.


    //double project(const Vec3d& dir, const Vec3d& property){
    //    rot.dot
    //    return property.dot(V);
    //}

    /// @brief Calculates the scattering amplitude between two directions (h0, h1).
    /// @brief The kernel is an exponential function of the angle and effective thickness: A(h0,h1) = exp( beta*cosa*(thick0 + thick1) ).
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

    /// @brief Simulates a bidirectional scattering event between two channels connected to this element.
    void scatterBiDir( Channel& chi, Channel& chj, int iend, int jend ){
        double amp = scatterAmp( chi.dir, chj.dir );
        chj.ends[jend].influx( amp * chi.ends[iend].fluxOut );
        chi.ends[iend].influx( amp * chj.ends[jend].fluxOut );
        //chj.ends[jend].fluxIn += amp * chi.ends[iend].fluxOut;
        //chi.ends[iend].fluxIn += amp * chj.ends[jend].fluxOut;
    }

};

/// @brief A direct iterative solver for flux transport problems.
/// @brief This class simulates the explicit propagation and scattering of flux through a network of discrete channels connecting scattering elements. It avoids building a global matrix, making it potentially more memory-efficient and scalable for certain problems than a Radiosity-like approach.
class Scattering2 : public TriangleRayTracer, public LinSolver { public:

    std::vector<ScatterElem2> scattelems;
    std::vector<Channel>     channels;
    int*       backChans = 0; // backward maping of channles

    /// @brief Empty override, as this class does not use the `LinSolver`'s matrix-solving capabilities.
    virtual void dotFunc( int n, double * x, double * Ax ) override {}

    /// @brief Constructs the network of `Channel`s between all pairs of `ScatterElem`s.
    /// @brief For each pair, it checks for occlusion. If the path is clear, a channel is created with geometric properties derived from the elements' positions and areas.
    void makeChannles(){
        int n = scattelems.size();
        for(int i=0; i<n; i++){
            ScatterElem2& eli = scattelems[i];
            for(int j=0; j<i; j++){
                ScatterElem2& elj = scattelems[j];
                if( eli.isurf == elj.isurf ){ continue; }
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
                    channels.push_back( Channel( d, Si*invS, Si*invS, thi,thj, i, j ) ); // Pass element indices to HalfChannels
                }
            }
        }
    }

    /// @brief Populates the `scattelems` vector from the `elements` (SurfElement) vector.
    /// @brief Initializes the physical properties of each `ScatterElem2` with provided values.
    /// @param thicks The thickness vector for the scattering elements.
    /// @param areas The area vector for the scattering elements.
    /// @param beta The decay exponent for the scattering kernel.
    void setupScatterElems(const Vec3d& thicks, const Vec3d& areas, double beta) {
        scattelems.clear();
        for (const auto& surfElem : elements) { // `elements` is inherited from `TriangleRayTracer`
            ScatterElem2 elem;
            elem.pos = surfElem.pos;
            elem.isurf = surfElem.isurf;
            elem.rot = Mat3dIdentity; // Default orientation
            elem.thicks = thicks;
            elem.areas = areas;
            elem.beta = beta;
            elem.fluxIn = 0.0;  // Initialize flux
            elem.fluxOut = 0.0; // Initialize flux
            scattelems.push_back(elem);
        }
    }

    /// @brief Performs one full simulation step.
    /// @brief This involves two phases: 1) A scattering phase where flux is exchanged locally at each element between its connected channels, and 2) A transfer phase where the newly scattered flux is propagated along the channels to the other end.
    inline void step_Direct(){
        // ---- scattering phase
        for(int k=0; k<n; k++){ // node i
            ScatterElem2& elem = scattelems[k]; // Get the current scattering element
            elem.fluxOut = 0.0; // Reset outgoing flux for this step
            elem.fluxIn = 0.0; // Reset incoming flux for this step, it will be accumulated from channels
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
        for(Channel& ch: channels){ 
            ch.transfer(); 
            // Accumulate flux into ScatterElem2's fluxIn for visualization/next step
            scattelems[ch.ends[0].elemIdx].fluxIn += ch.ends[0].fluxOut;
            scattelems[ch.ends[1].elemIdx].fluxIn += ch.ends[1].fluxOut;
        }
    }

};

#endif
