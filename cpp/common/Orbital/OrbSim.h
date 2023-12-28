
#ifndef  OrbSim_h
#define  OrbSim_h

#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
//#include "datatypes.h"  


float springForce( float l, float& f, Quat4f par ){
    float dl = l - par.x;
    float k;
    if( dl > 0.0f ){
        k = -par.z;
    } else {
        k = par.y;
    }
    //Quat4f fe; 
    //fe.f = d*(k*dl/l);
    //fe.e = k*dl*dl;
    f = k*dl;
    return k*dl*dl;
}

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
    Quat4f* points=0;  // position and mass
    Quat4f* forces=0;  // force and energy
    Quat4f* vel   =0;  // velocity

    Quat4f* params=0;  // neighbor parameters (l0,kP,kT,damp)
    int*    neighs=0;  // neighbor indices
    int*    neighBs=0; // neighbor bond indices

    int     nBonds =0; // number of bonds
    Quat4f* bparams=0; // bond parameters (l0,kP,kT,damp)
    int2*   bonds  =0; // indices of bonded points (i,j)
    float*  strain =0; // strain
    //float*  l0s    =0; // 
    Vec2f*  maxStrain=0;

    void recalloc( int nPoint_, int nNeighMax_, int nBonds_=0 ){
        nPoint = nPoint_; nNeighMax = nNeighMax_;
        nNeighTot = nPoint*nNeighMax;
        _realloc( points, nPoint    );
        _realloc( forces, nPoint    );
        _realloc( vel,    nPoint    );
        _realloc( params, nNeighTot );
        _realloc( neighs, nNeighTot );
        _realloc( neighBs,nNeighTot );

        if(nBonds_>0){
            nBonds=nBonds_;
            _realloc( bonds,     nBonds );
            _realloc( strain,    nBonds );
            //_realloc( l0s,       nBonds );
            _realloc( maxStrain, nBonds );
            _realloc( bparams,   nBonds );

        }
    }

    void evalTrussForce_neighs(){
        //#pragma omp paralel for 
        for(int iG=0; iG<nPoint; iG++){
            //const int iG = get_global_id(0);
            Quat4f p = points[iG];
            Quat4f f =Quat4f{0.0f,0.0f,0.0f,0.0f};
            //printf( "--- p[%i] \n", iG );
            //#pragma omp simd
            for(int ij=0; ij<nNeighMax; ij++){
                int j  = nNeighMax*iG + ij;
                int ja = neighs[j];
                if(ja == -1) break;
                //f.add( springForce( points[ja].f - p.f, params[j] ) );
                
                Vec3f d =  points[ja].f - p.f;
                float li = d.norm();
                /*
                float fi,ei = springForce( li, fi, params[j] );
                //f.add( Quat4f{ d*(fi/l), ei } );
                f.f.add_mul( d, fi/li );
                */
                float k = 1e+6;
                f.f.add_mul( d, (k*(li-params[j].x)/li) );

                //printf( "p[%i,ij=%i,j=%i] li=%7.3f dl=%8.5e fi=%8.5e e=%8.5e par(%7.3f,%8.5e,%8.5e,%8.5e) \n", iG,ij,ja, li, li-params[j].x, fi,ei, params[j].x,params[j].y,params[j].z,params[j].w );
            }
            forces[iG] = f; // we may need to do += in future
        } 
        //exit(0);
    }

    void evalTrussForce_neighs2(){
        //#pragma omp paralel for 
        for(int iG=0; iG<nPoint; iG++){
            //const int iG = get_global_id(0);
            Quat4f p = points[iG];
            Quat4f f =Quat4f{0.0f,0.0f,0.0f,0.0f};
            //printf( "--- p[%i] \n", iG );
            //#pragma omp simd
            for(int ij=0; ij<nNeighMax; ij++){
                int j  = nNeighMax*iG + ij;
                int ja = neighs[j];
                if(ja == -1) break;

                int ib = neighBs[j];
                Quat4f par = bparams[ib];
                //f.add( springForce( points[ja].f - p.f, params[j] ) );
                
                Vec3f d =  points[ja].f - p.f;
                float li = d.norm();
                /*
                float fi,ei = springForce( li, fi, params[j] );
                //f.add( Quat4f{ d*(fi/l), ei } );
                f.f.add_mul( d, fi/li );
                */
                float k = 1e+6;
                f.f.add_mul( d, (k*(li-par.x)/li) );

                //printf( "p[%i,ij=%i,j=%i] li=%7.3f dl=%8.5e fi=%8.5e e=%8.5e par(%7.3f,%8.5e,%8.5e,%8.5e) \n", iG,ij,ja, li, li-params[j].x, fi,ei, params[j].x,params[j].y,params[j].z,params[j].w );
            }
            forces[iG] = f; // we may need to do += in future
        } 
        //exit(0);
    }

    void evalTrussForce_bonds(){
        //#pragma omp paralel for 
        for(int i=0; i<nBonds; i++){
            int2  b = bonds[i];
            Vec3f d = points[b.y].f - points[b.x].f;
            float li = d.norm();
            //float fi,ei = springForce( li, fi, bparams[i] );
            float k = 1e+6;
            d.mul( (k*(li-bparams[i].x)/li) );
            forces[b.x].f.add(d);
            forces[b.y].f.sub(d);
        } 
        //exit(0);
    }

    void evalBondTension(){
        for(int i=0;i<nBonds; i++ ){
            int2  b  = bonds[i];
            float l0 = bparams[i].x;
            //float l0 = l0s[i];
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
    void cleanVel  (){ for (int i=0; i<nPoint; i++){ vel   [i]=Quat4fZero; } };

    void move_GD(float dt){
        for(int i=0;i<nPoint; i++ ){
            Quat4f p = points[i];
            Quat4f f = forces[i];
            p.f.add_mul( f.f, dt/p.w );
            //printf( "move_GD[%i] |d|=%g |f|=%g dt/m=%g m=%g \n", i, f.f.norm() * dt/p.w, f.f.norm(), dt/p.w, p.w );
            points[i]=p;
        }
    }

    void move_MD(float dt, float damp=0.0f ){
        float cdamp = 1.0f - damp;
        for(int i=0;i<nPoint; i++ ){
            Quat4f p = points[i];
            Quat4f f = forces[i];
            Quat4f v = vel   [i];
            v.f.mul( cdamp );
            v.f.add_mul( f.f, dt/p.w );
            p.f.add_mul( v.f, dt     );
            //printf( "move_GD[%i] |d|=%g |f|=%g dt/m=%g m=%g \n", i, f.f.norm() * dt/p.w, f.f.norm(), dt/p.w, p.w );
            vel   [i]=v;
            points[i]=p;
        }
    }

};   // OrbSim_f

#endif
