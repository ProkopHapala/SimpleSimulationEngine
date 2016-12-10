
#ifndef  PointCloudComparator_h
#define  PointCloudComparator_h

#include <vector>
#include <set>
#include <unordered_map>
//#include "CubicRuler.h"

#include "fastmath.h"
#include "Vec3.h"


// https://en.wikipedia.org/wiki/Hungarian_algorithm

class Neigh{
    public:
    int i;     // assigned neighbor
    double r2; // distance

    friend bool operator<(const Neigh& l, const Neigh& r){ return l.r2<r.r2; }
};

class PointCloudComparator {
    public:
    int n;
    Vec3d  * points    = NULL;
    double * overlap_i = NULL;
    double * overlap_j = NULL;


    std::unordered_map<uint64_t,int> gridpoints;

    std::vector<int>   unassigned;
    std::vector<Neigh> neighs_ij;
    //std::vector neighs_ji;
    int   * ijoffs = NULL;
    //int   * jioffs = NULL;
    Neigh * assigned = NULL;


    double rmax, r2max, invRmax, invR2max;
    double rmin, r2min, invRmin, invR2min;
    Vec3d  pmin;
    Vec3d  pmax;
    double step;
    double invStep;

    // ========== functions

    inline void index2pos( const Vec3i& index, Vec3d& pos    ){ pos  .set( index.x*step+pmin.x,    index.y*step+pmin.y,    index.z*step+pmin.z    );  };
    inline void pos2index( const Vec3d& pos,   Vec3i& index  ){ index.set( (pos.x-pmin.x)*invStep, (pos.y-pmin.y)*invStep, (pos.z-pmin.z)*invStep );  };


    // rmin : minimal distance between different atoms
    // rmax : maximim distance between corresponding atoms,
    //  => rmax is often smaller than rmin !!! (that is not mistake)
    void setGrid(  double rmin_, double rmax_ ){
        rmax    = rmax_; r2max = rmax*rmax;  invR2max = 1/r2max; invRmax = 1/rmax;
        rmin    = rmin_; r2min = rmin*rmin;  invR2min = 1/r2min; invRmin = 1/rmin;
        step    = 0.95*rmin/( sqrt(3.0d) );
        invStep = 1/step;
        //posMin  = pos0_;
        printf( "rmin %f rmax %f step %f \n", rmin, rmax, step );
    }

    void allocate( int n_ ){
        if( overlap_i ) delete overlap_i;   overlap_i = new double[n_];
        if( overlap_j ) delete overlap_j;   overlap_j = new double[n_];
    }

    void setFindBounds(){
        pmin.set(+1.0e+300,+1.0e+300,+1.0e+300);
        pmax.set(-1.0e+300,-1.0e+300,-1.0e+300);
        printf("n %i\n",n);
        for(int i=0; i<n; i++ ){
            Vec3d& p = points[i];
            if(p.x<pmin.x){ pmin.x=p.x;} if(p.x>pmax.x){ pmax.x=p.x;}
            if(p.y<pmin.y){ pmin.y=p.y;} if(p.y>pmax.y){ pmax.y=p.y;}
            if(p.z<pmin.z){ pmin.z=p.z;} if(p.z>pmax.z){ pmax.z=p.z;}
            printf( "%i p: (%3.3f,%3.3f,%3.3f) \n", i, p.x, p.y, p.z );
        }
        printf( "pmin (%3.3f,%3.3f,%3.3f) pmax (%3.3f,%3.3f,%3.3f)\n", pmin.x, pmin.y, pmin.z, pmax.x, pmax.y, pmax.z );
        printf( "rmin %f rmax %f step %f \n", rmin, rmax, step );
        pmin.add(-2*rmax,-2*rmax,-2*rmax);
        pmax.add(+2*rmax,+2*rmax,+2*rmax);
        printf( "pmin (%3.3f,%3.3f,%3.3f) pmax (%3.3f,%3.3f,%3.3f)\n", pmin.x, pmin.y, pmin.z, pmax.x, pmax.y, pmax.z );
    }

    void setRefPoints( int n_, Vec3d * points_ ){
        n=n_; points = points_;
        setFindBounds();
        for(int i=0; i<n; i++){
            Vec3i vi;
            pos2index( points[i], vi );
            uint64_t ibox = scalar_id( vi );
            printf( "%i : %i %i %i %li \n", i, vi.x, vi.y, vi.z, ibox );
            //gridpoints[ibox] = i;
            auto p = gridpoints.insert( {ibox,i} );
            if( !p.second ){
                printf( "cannot insert point %i occupied by %i \n", i, *p.first );
            }
        }

    }

    int findNeighFor( int i ){
        // insert to best free position
        bool was_assigned = false;
        for( int ineigh=ijoffs[i]; ineigh<ijoffs[i+1]; ineigh++ ){
            Neigh& ng = neighs_ij[i];
            if( assigned[ng.i].i<0 ){
                assigned[ng.i].i  = i;
                assigned[ng.i].r2 = ng.r2;
                return -1;
            };
            if( assigned[ng.i].r2 > ng.r2 ){
                int displaced = assigned[ng.i].i;
                assigned[ng.i].i  = i;
                assigned[ng.i].r2 = ng.r2;
                return displaced;
            }
        }
    }

    // this search for assigment of correspoding atoms ( quite complex algorithm )
    // see https://en.wikipedia.org/wiki/Hungarian_algorithm
    double evalDistance_N2( int n2, Vec3d * points2 ){
        // build sparse distance matrix
        for( int i=0; i<n; i++ ){ assigned[i].i=-1; assigned[i].r2=1e+300; }
        neighs_ij.clear();
        std::set<Neigh> neighi;
        for( int i=0; i<n2; i++ ){
            int nfound = 0;
            Vec3d& pi = points2[i];
            neighi.clear();
            for( int j=0; j<n; j++ ){
                Vec3d dp     = pi - points[j];
                double r2 = dp.norm2();
                if( r2<r2max ){
                    neighi.insert( {j,r2} );
                    nfound++;
                };
            }
            if( nfound <= 0 ) return 1e+300;

            neighs_ij.reserve( neighs_ij.size()+neighi.size() );
            std::copy(neighi.begin(), neighi.end(), neighs_ij.begin());
            ijoffs[i] = neighs_ij.size();

        }
    }

    inline double dist_func( double r2 ){
        //return 1/(r2+1e-16);
        //return sq(1-(r2*invR2max))/(r2+1e-16);
        return (r2max-r2)/(r2+1e-16);
    }

    // compute overlap of atomic densities ... much simpler
    double overlap_N2( Vec3d * points2 ){
        double r2max = rmax*rmax;
        double dist = 0;
        for( int i=0; i<n; i++ ){
            overlap_i[i] = 0;
            overlap_j[i] = 0;
        }
        for( int i=0; i<n; i++ ){
            Vec3d& pi = points2[i];
            //printf( "%i (%3.3f,%3.3f,%3.3f)\n", i, pi.x, pi.y, pi.z );
            double sumi = 0;
            for( int j=0; j<n; j++ ){
                Vec3d dp     = pi - points[j];
                double r2 = dp.norm2();
                //printf( "%  i (%3.3f,%3.3f,%3.3f)  %f %f %f \n", j, points[j].x, points[j].y, points[j].z, r2, dist_func( r2 ), sumi );
                if (r2>r2max) continue;
                double overlap_ij = dist_func( r2 );
                //overlap_i[i] += overlap_ij;
                sumi         += overlap_ij;
                overlap_j[j] += overlap_ij;
            }
            if( sumi < 1e-16 ) return 1e+300;
            dist += 1/sumi;
        }
        for( int i=0; i<n; i++ ){
            double sumj = overlap_j[i];
            if( sumj < 1e-16 ) return 1e+300;
            dist += 1/sumj;
            //dist += 1/overlap_i[i] + 1/overlap_j[i];
        }
        return dist;
    }

};

#endif

