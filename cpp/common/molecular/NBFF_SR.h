
#ifndef _h
#define NBFF_SR_h

#include <vector>

#include "fastmath.h"
//#include "Vec2.h"
#include "Vec3.h"
//#include "SMat3.h"
#include "quaternion.h"
#include "Forces.h"
#include "Buckets3D.h"

void fitAABB( Vec6d& bb, int n, int* c2o, Vec3d* ps ){
    //Quat8d bb;
    bb.lo = bb.hi = ps[c2o[0]];
    for(int i=1; i<n; i++){ 
        //printf( "fitAABB() i %i \n", i );
        int ip = c2o[i];
        Vec3d p = ps[ip];
        //printf( "fitAABB() i=%i ip=%i p(%g,%g,%g) \n", i, ip, p.x, p.y, p.z );
        bb.lo.setIfLower  ( p );
        bb.hi.setIfGreater( p );
    }; 
    //return bb;
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
    double Rcol  = 1.5;   // collision radius
    double Rcut  = 5.0;   // cut-off radius for non-bonded interactions
    double Rf    = 3.0;
    double R0_nb = 2.0;  // bond length
    double E0_nb =  0.1; // stiffnes of bonds (springs) within cluster
    double K_nb  =  0.1;  // stiffnes of bonds (springs) within cluster

    // clusters and their bounding boxes
    int      nBBs=0;
    Vec6d*   BBs=0;      // bounding boxes (can be either AABB, or cylinder, capsula) 
    Buckets  pointBBs;   // buckets for collision detection

    // picking and ray-casting
    int ipicked=-1;
    Vec3d pick_hray, pick_ray0; 
    double Kpick = -10.0;

    void realloc( int natom_, int nBBs_ ){
        //nPerBox  = nPerBox_;
        //nCluster = nCluster_;
        natom = natom_;
        nBBs = nBBs_;
        pointBBs.realloc(nBBs, natom );
        _realloc( BBs, nBBs );
        _realloc( REQs   , natom );
        _realloc( apos   , natom );
        _realloc( avel   , natom );
        _realloc( aforce , natom );

        _realloc( apos_old, natom );
        
    }

    void clean_forces  ( Vec3d zero = Vec3dZero ){ for(int ia=0; ia<natom; ia++){ aforce[ia] = zero; } }
    void clean_velocity( Vec3d zero = Vec3dZero ){ for(int ia=0; ia<natom; ia++){ avel[ia] = zero;   } }

    void step( double dt ){
        
        backup_positions();
        update_bboxes();

        //exit(0);
        clean_forces();
        //evalSortRange_BBs( Rcut, perCluster );
        evalSortRange_brute();   

        if(ipicked>=0){ 
            aforce[ipicked].add( getForceSpringRay( apos[ipicked], pick_hray, pick_ray0,  Kpick ) ); 
        }


        //backup_positions();
        move_MD( dt, 0.1 );      
        //move_GD( dt ); 
        //solve_sphere_colision_brute( );
        int nConstrIter = 10; 
        for(int i=0;i<nConstrIter;i++){
            solve_cluster_bonds_GSB_back();
            //solve_cluster_bonds_GSB( );  
            
        } 
        update_velocity( dt );     

        //update_velocity_imp(dt);
    };


    void update_bboxes(){
        //printf( "@cellI0s=%li @cellNs=%li @cell2obj=%li \n", (long)pointBBs.cellI0s, (long)pointBBs.cellNs, (long)pointBBs.cell2obj );
        for(int ib = 0; ib<nBBs; ib++){
            int i0 = pointBBs.cellI0s[ib];
            int ni = pointBBs.cellNs[ib];
            //printf( "ib=%i, i0=%i, ni=%i \n", ib, i0, ni );
            fitAABB( BBs[ib], ni, pointBBs.cell2obj + i0, apos );
            //printf( "BBs[%i] i0%i ni=%i min(%g,%g,%g) max(%g,%g,%g) \n", ib, i0, ni, BBs[ib].lo.x, BBs[ib].lo.y, BBs[ib].lo.z, BBs[ib].hi.x, BBs[ib].hi.y, BBs[ib].hi.z );
        }
    }

    void backup_positions(){
        for(int ia=0; ia<natom; ia++){ apos_old[ia] = apos[ia];}
    }

    void update_velocity(double dt){
        double inv_dt = 1.0/dt;
        for(int ia=0; ia<natom; ia++){
            avel    [ia] = ( apos[ia] - apos_old[ia] ) * inv_dt;
            //apos_old[ia] = apos[ia];
        }
    }

    void update_velocity_imp(double dt){
        double inv_dt = 1.0/dt;
        double inv_dt2 = inv_dt*inv_dt;
        for(int ia=0; ia<natom; ia++){
            avel    [ia].add_mul( apos[ia] - apos_old[ia], -0.5*inv_dt );
        }
    }

    void move_GD( double dt ){
        for(int ia=0; ia<natom; ia++){
            Vec3d& p = apos[ia];
            Vec3d& f = aforce[ia];
            p.add_mul( f, dt);
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

    void solve_cluster_bonds_GSB  ( ){
        //printf( "solve_cluster_bonds_GSB() \n" );
        // Gauss-Seidel method
        for(int ib=0; ib<nBBs; ib++){
            Vec6d& bb = BBs[ib];
            const int  npi = pointBBs.cellNs[ib];
            const int  ip0 = pointBBs.cellI0s[ib];
            //printf( "solve_cluster_bonds() [%i] npi=%i ip0=%i \n", ib, npi, ip0 );
            //const int* ips = pointBBs.cell2obj + ip0;
            int iao =  ip0 + npi - 1;
            for(int ia=ip0; ia<ip0+npi; ia++){
                //printf( "solve_cluster_bonds() [%i|%i,%i] dl=%g d(%g,%g,%g) \n", ib,ia,iao, l-L_bond, d.x, d.y, d.z );

                Vec3d  d = apos[ia] - apos[iao];
                double l = d.norm();
                d.mul( 0.5*(L_bond-l)/l );

                apos[ia ].add( d );
                apos[iao].sub( d );

                iao = ia;
                
            }
        }
        //exit(0);
    }

    void solve_cluster_bonds_GSB_back  ( ){
        //printf( "solve_cluster_bonds_GSB() \n" );
        // Gauss-Seidel method
        for(int ib=0; ib<nBBs; ib++){
            Vec6d& bb = BBs[ib];
            const int  npi = pointBBs.cellNs[ib];
            const int  ip0 = pointBBs.cellI0s[ib];
            //printf( "solve_cluster_bonds() [%i] npi=%i ip0=%i \n", ib, npi, ip0 );
            //const int* ips = pointBBs.cell2obj + ip0;
            int iao =  ip0;
            for(int ia=ip0+npi-1; ia>=ip0; ia--){
                //printf( "solve_cluster_bonds() [%i|%i,%i] dl=%g d(%g,%g,%g) \n", ib,ia,iao, l-L_bond, d.x, d.y, d.z );

                Vec3d  d = apos[ia] - apos[iao];
                double l = d.norm();
                d.mul( 0.5*(L_bond-l)/l );

                apos[ia ].add( d );
                apos[iao].sub( d );

                iao = ia;
                
            }
        }
        //exit(0);
    }

    void solve_cluster_bonds_GSP( ){
        printf( "solve_cluster_bonds_GSP() \n" );
        // Gauss-Seidel method
        for(int ib=0; ib<nBBs; ib++){
            Vec6d& bb = BBs[ib];
            const int  npi = pointBBs.cellNs[ib];
            const int  ip0 = pointBBs.cellI0s[ib];
            //printf( "solve_cluster_bonds() [%i] npi=%i ip0=%i \n", ib, npi, ip0 );
            //const int* ips = pointBBs.cell2obj + ip0;
            //int iao =  ip0 + npi - 1;
            for(int ia=ip0; ia<ip0+npi; ia++){
                int ia1 = ia-1; if(ia1<ip0     ) ia1 = ip0 + npi - 1;
                int ia2 = ia+1; if(ia2>=ip0+npi) ia2 = ip0;
                //printf( "solve_cluster_bonds() [%i|%i,%i] dl=%g d(%g,%g,%g) \n", ib,ia,iao, l-L_bond, d.x, d.y, d.z );

                Vec3d  d1 = apos[ia] - apos[ia1];
                double l1 = d1.norm();
                d1.mul( 0.5*(L_bond-l1)/l1 );

                Vec3d  d2 = apos[ia] - apos[ia2];
                double l2 = d2.norm();
                d2.mul( 0.5*(L_bond-l2)/l2 );

                apos[ia ].add( d1+d2 );
                //apos[ia1].sub( d1 );
                //apos[ia2].sub( d2 );
                
            }
        }
        //exit(0);
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
                    Vec3d f; 
                    // double Rcut, double R0, double E0, double K
                    //double Ei = getSR( d, f, Rcut, R0_nb, E0_nb, K_nb  );
                    double Ei = getSR2( d, f, Rcut, R0_nb, E0_nb, K_nb, Rf  );
                    //printf( "evalSortRange_brute[%i,%i] r=%g E=%g d(%g,%g,%g) f(%g,%g,%g) \n", ia, ja, r, Ei, d.x, d.y, d.z, f.x, f.y, f.z );
                    aforce[ia].sub( f );
                    aforce[ja].add( f );
                    E += Ei;
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
