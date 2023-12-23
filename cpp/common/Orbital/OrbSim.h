
#ifndef  OrbSim_h
#define  OrbSim_h

#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
//#include "datatypes.h"  

Quat4f springForce( Vec3f d, Quat4f par ){
    float l  = d.norm();
    float dl = l - par.x;
    float k;
    if( dl > 0.0f ){
        k = -par.z;
    } else {
        k = par.y;
    }
    Quat4f fe; 
    fe.f = d*(k*dl/l);
    fe.e = k*dl*dl;
    return fe;
}

class OrbSim_f{ public:
    int nPoint=0, nNeighMax=0, nNeighTot=0;
    // cpu buffers
    Quat4f* points=0;
    Quat4f* forces=0;
    Quat4f* params=0;
    int*    neighs=0;

    void recalloc( int nPoint_, int nNeighMax_ ){
        nPoint = nPoint_; nNeighMax = nNeighMax_;
        nNeighTot = nPoint*nNeighMax;
        _realloc( points, nPoint    );
        _realloc( forces, nPoint    );
        _realloc( params, nNeighTot );
        _realloc( neighs, nNeighTot );
    }

    void evalTrussForce(){
        //#pragma omp paralel for 
        for(int iG=0; iG<nPoint; iG++){
            //const int iG = get_global_id(0);
            Quat4f p = points[iG];
            Quat4f f =Quat4f{0.0f,0.0f,0.0f,0.0f};
            //#pragma omp simd
            for(int ij=0; ij<nNeighTot; ij++){
                int j  = nNeighTot*iG + ij;
                int ja = neighs[j];
                if(j == -1) break;
                f.add( springForce( points[ja].f - p.f, params[j] ) );
            }
            forces[iG] = f; // we may need to do += in future
        } 
    }

    void printNeighs(int i){
        int j0 = i*nNeighMax;
        for(int jj=0;jj<nNeighMax;jj++){
            int j=j0+jj;
            int ing = neighs[j];
            if(ing<0) break;
            Quat4f par = params[j];
            printf( "ng[%i,%i|%i] l0,kP,kT,damp(%g,%g,%g,%g)\n", i, ing, jj, par.x,par.y,par.z,par.w );
        }
    }
    void printAllNeighs(){ printf("OrbSim_f::printAllNeighs(nPoint=%i,nNeighMax=%i)\n",nPoint,nNeighMax); for(int i=0;i<nPoint;i++){ printNeighs(i); }; };


};   // OrbSim_f

#endif
