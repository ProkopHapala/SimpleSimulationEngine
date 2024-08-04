
#ifndef _h
#define NBFF_SR_h

#include <vector>

#include "fastmath.h"
//#include "Vec2.h"
#include "Vec3.h"
//#include "SMat3.h"
#include "quaternion.h"
//#include "Forces.h"
#include "Buckets3D.h"

void fitAABB( Vec6d& bb, int n, int* c2o, Vec3d* ps ){
    //Quat8d bb;
    //bb.lo = bb.lo = ps[c2o[0]];
    for(int i=0; i<n; i++){ 
        //printf( "fitAABB() i %i \n", i );
        int ip = c2o[i];
        //printf( "fitAABB() i=%i ip=%i \n", i, ip );
        Vec3d p = ps[ip];
        bb.lo.setIfLower  ( p );
        bb.hi.setIfGreater( p );
    }; 
    //return bb;
}


double getSR( const Vec3d& d, Vec3d& f, double Rcut, double R0, double E0, double K  ){
    double r2 = d.norm2();
    if(r2>(Rcut*Rcut)){
        f = Vec3dZero;
        return 0.0;
    }
    double r  = sqrt(r2)-R0;
    double dr = r-R0;
    double E  = 0.5*K*(dr*dr) - E0;
    f         = d*(-K*dr);
    return E;
}

class NBFF_SR{ public:
// non bonded forcefield
    int natom     = 0; // [natom] 
    Vec3d* REQs   = 0; // [natom] parameters {RvdW, EvdW, Q_charge}
    Vec3d* apos   = 0; // [natom] position
    Vec3d* avel   = 0; // [natom] velocity
    Vec3d* aforce = 0; // [natom] forces

    Vec3d* apos_old = 0; // [nAtom] atomic position before last step

    //int nCluster   = 0;
    //int nPerBox    = 16; // how many atoms per cluster

    // global parameters of bonded interactions
    double L_bond = 1.0; // bond length
    double K_bond = 1.0; // stiffnes of bonds (springs) within cluster

    // global parameters of non-bonded interactions
    double Rcol  = 1.0;   // collision radius
    double Rcut  = 5.0;   // cut-off radius for non-bonded interactions
    double R0_nb = 1.0;  // bond length
    double E0_nb = 1.0; // stiffnes of bonds (springs) within cluster
    double K_nb  = 1.0;  // stiffnes of bonds (springs) within cluster

    // clusters and their bounding boxes
    int      nBBs=0;
    Vec6d*   BBs=0;      // bounding boxes (can be either AABB, or cylinder, capsula) 
    Buckets  pointBBs;   // buckets for collision detection

    void realloc( int natom_, int nCluster ){
        //nPerBox  = nPerBox_;
        //nCluster = nCluster_;
        natom = natom_;
        pointBBs.realloc(nCluster, natom );
        _realloc( REQs   , natom );
        _realloc( apos   , natom );
        _realloc( avel   , natom );
        _realloc( aforce , natom );

        _realloc( apos_old, natom );
        
    }

    void clean_forces  ( Vec3d zero = Vec3dZero ){ for(int ia=0; ia<natom; ia++){ aforce[ia] = zero; } }
    void clean_velocity( Vec3d zero = Vec3dZero ){ for(int ia=0; ia<natom; ia++){ avel[ia] = zero;   } }

    void step( double dt ){
        clean_forces();
        //evalSortRange_BBs( Rcut, perCluster );
        evalSortRange_brute();   
        move_MD( dt, 0 );       
        solve_sphere_colision_brute( ); 
        solve_cluster_bonds  ( );  
        update_velocity( dt );     
    };

    void update_velocity(double dt){
        double inv_dt = 1.0/dt;
        for(int ia=0; ia<natom; ia++){
            avel    [ia] = ( apos[ia] = apos_old[ia] ) * inv_dt;
            apos_old[ia] = apos[ia];
        }
    }

    void move_MD( double dt, double damping ){
        double cdamp = 1-damping;
        for(int ia=0; ia<natom; ia++){
            Vec3d& p = apos[ia];
            Vec3d& v = avel[ia];
            Vec3d& f = aforce[ia];
            v.mul( cdamp );
            v.add_mul( f, dt);
            p.add_mul( v, dt);
        }
    }

    void solve_cluster_bonds  ( ){
        // Gauss-Seidel method
        for(int ib=0; ib<nBBs; ib++){
            Vec6d& bb = BBs[ib];
            const int  npi = pointBBs.cellNs[ib];
            const int  ip0 = pointBBs.cellI0s[ib];
            //const int* ips = pointBBs.cell2obj + ip0;
            int iao =  ip0 + npi - 1;
            for(int ia=ip0; ia<npi; ia++){
                Vec3d  d = apos[ia] - apos[iao];
                double l = d.norm();
                d.mul( (L_bond-l)/l );
                apos[ia ].add( d );
                apos[iao].sub( d );
            }
        }
    }

    void solve_sphere_colision_brute( ){
        // brute force
        double Rcol2 = Rcol*Rcol;
        for(int ia=0; ia<natom; ia++){
            Vec3d& pi = apos[ia];
            for(int ja=ia; ja<natom; ja++){
                if(ia==ja) continue;
                Vec3d& pj  = apos[ja];
                Vec3d  d   = pj-pi;
                double r2  = d.norm2();
                if(r2<Rcol2){
                    double r = sqrt(r2);
                    d.mul( (Rcol-r)/r );
                    pi.add( d );
                    pj.sub( d );
                }
            }
        }
    }

    double evalSortRange_brute( ){
        // brute force
        double Rcut2 = Rcut*Rcut;
        double E = 0;
        for(int ia=0; ia<natom; ia++){
            Vec3d& pi = apos[ia];
            for(int ja=ia; ja<natom; ja++){
                if(ia==ja) continue;
                Vec3d& pj  = apos[ja];
                Vec3d  d   = pj-pi;
                double r2  = d.norm2();
                if(r2<Rcut2){
                    double r = sqrt(r2);
                    Vec3d f; E += getSR( d, f, R0_nb, R0_nb, K_nb, Rcut );
                    aforce[ia].add( f );
                    aforce[ja].sub( f );
                }
            }
        }
        return E;
    }   

/*
    __attribute__((hot))  
    inline void updatePointBBs( bool bInit=true){
        const Buckets& buckets = pointBBs;
        //printf( "updatePointBBs() START \n" );
        for(int ib=0; ib<buckets.ncell; ib++){
            //printf( "updatePointBBs() ib %i \n", ib );
            if(bInit){ BBs[ib].lo = Vec3dmax; BBs[ib].hi = Vec3dmin; }
            int n = buckets.cellNs[ib];
            if(n>0){
                int i0 = buckets.cellI0s[ib];
                //printf( "updatePointBBs() ib %i n %i i0 %i \n", ib, n, i0 );
                fitAABB( BBs[ib], n, buckets.cell2obj+i0, vapos );
            }
        }
        //printf( "updatePointBBs() DONE \n" );
    }

    __attribute__((hot))  
    inline int selectInBox( const Vec6d& bb, const int ib, Vec3d* ps, Quat4d* paras, int* inds ){
        const int  npi = pointBBs.cellNs[ib];
        const int* ips = pointBBs.cell2obj +pointBBs.cellI0s[ib];
        int  n   = 0;
        for(int i=0; i<npi; i++){
            int ia = ips[i];
            const Vec3d& p = vapos[ ia ];
            if( (p.x>bb.lo.x)&&(p.x<bb.hi.x)&&
                (p.y>bb.lo.y)&&(p.y<bb.hi.y)&&
                (p.z>bb.lo.z)&&(p.z<bb.hi.z) 
            ){
                ps[i]    = p;
                paras[i] = REQs[ ips[i] ];
                inds[i]  = ia;
                n++;
            }
        }
        return n;
    }

    __attribute__((hot))  
    double evalSortRange_BBs( double Rcut, int ngmax ){
        //printf( "evalSortRange_BBs() START \n" );
        double E=0;
        const double R2cut = Rcut*Rcut;
        Quat4d REQis[ngmax];
        Quat4d REQjs[ngmax];
        Vec3d  pis  [ngmax];
        Vec3d  pjs  [ngmax];
        int   nis   [ngmax];
        int   njs   [ngmax];

        for(int ib=0; ib<nBBs; ib++){

            // --- Within the same bucket
            const int  ni =  selectInBox( BBs[ib], ib, pis, REQis, nis );  // this is maybe not optimal (we just select all atoms in the bucket)
            for(int ip=0; ip<ni; ip++){
                const Vec3d&  pi   = pis[ip];
                const Quat4d& REQi = REQis[ip];
                const int     ia   = nis[ip];
                for(int jp=ip+1; jp<ni; jp++){
                    const Vec3d& d  = pjs[jp] - pi;
                    int          ja = nis[jp];
                    Vec3d f; getSR( d, f, REQij.x-drSR, REQij.x, 1.0, 1.0 );
                    fapos[ia].sub(f);
                    fapos[ja].add(f);
                }
            }

            // --- between different buckets 
            for(int jb=0; jb<ib; jb++){
                Vec6d bb;
                bb = BBs[jb]; bb.lo.add(-Rcut); bb.hi.add(Rcut); const int ni =  selectInBox( bb, ib, pis, REQis, nis );  // atoms in bucket i overlapping with bucket j
                bb = BBs[ib]; bb.lo.add(-Rcut); bb.hi.add(Rcut); const int nj =  selectInBox( bb, jb, pjs, REQjs, njs );  // atoms in bucket j overlapping with bucket i
                for(int jp=0; jp<nj; jp++){
                    const Vec3d&  pj   = pjs[jp];
                    const Quat4d& REQj = REQjs[jp];
                    const int     ja   = nis[jp];
                    for(int ip=0; ip<ni; ip++){
                        const Vec3d& d = pis[jp] - pj;
                        const int ia   = nis[ip];
                        Quat4d REQij; combineREQ( REQis[ip], REQj, REQij );
                        Vec3d f; repulsion_R4( d, f, REQij.x-drSR, REQij.x, REQij.y*ampSR );
                        fapos[ia].sub(f);
                        fapos[ja].add(f);
                    }
                }
            }

        }
        //printf( "evalSortRange_BBs() DONE \n" );
        return E;

    }

*/

};

#endif
