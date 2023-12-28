
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

    int nBonds=0;
    int2*  bonds  =0;
    float* strain =0;
    float* l0s    =0; 
    Vec2f* maxStrain=0;

    void recalloc( int nPoint_, int nNeighMax_, int nBonds_=0 ){
        nPoint = nPoint_; nNeighMax = nNeighMax_;
        nNeighTot = nPoint*nNeighMax;
        _realloc( points, nPoint    );
        _realloc( forces, nPoint    );
        _realloc( params, nNeighTot );
        _realloc( neighs, nNeighTot );

        if(nBonds_>0){
            nBonds=nBonds_;
            _realloc( bonds,     nBonds );
            _realloc( strain,    nBonds );
            _realloc( l0s,       nBonds );
            _realloc( maxStrain, nBonds );
        }
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

    void evalBondTension(){
        for(int i=0;i<nBonds; i++ ){
            int2  b  = bonds[i];
            float l0 = l0s[i];
            float s  = ((points[b.y]-points[b.x]).norm() - l0)/l0;
            // ToDo: break the bond if strain > maxStrain;
        }
    }

    void applyCentrifugalForce( Vec3f p0, Vec3f ax, float omega ){
        double omega2 = omega*omega;
        for(int i=0;i<nPoint; i++ ){
            const Quat4f& p = points[i];
            //Vec3f f; f.set_cross(ax,p.f-p0);
            Vec3f f; f.set_sub(p.f,p0); f.mul( p.w*omega2 );
            f.makeOrthoU(ax);
            forces[i].f.add(f);
            //forces[i].f.add_mul( f, p.w*omega2 );
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

    void cleanForce(){ for (int i=0; i<nPoint; i++){ forces[i]=Quat4fZero; } };

    void move_GD(float dt){
        for(int i=0;i<nPoint; i++ ){
            Quat4f p = points[i];
            Quat4f f = forces[i];
            float dtm = dt/p.w;
            points[i].f.add_mul( f.f, dtm );
        }
    }

};   // OrbSim_f

#endif
