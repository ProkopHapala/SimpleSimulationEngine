
#ifndef  TrussDynamics_d_cpp
#define  TrussDynamics_d_cpp

#include <stdio.h>
#include <string.h>

#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
//#include "datatypes.h"  
#include "VecN.h"

#include "Buckets.h"

#include "raytrace.h"
#include "geom3D.h"
#include "Interfaces.h"

#include "Cholesky.h"
#include "CG.h"
#include "arrayAlgs.h"

#include "SparseMatrix.h"
#include "SparseMatrix2.h"

#include "TrussDynamics_d.h"

#include "IO_utils.h"
#include "testUtils.h"


bool checkMatrixNaN(int n, int m, const double* A, const char* name, bool bExit=false){
    for(int i=0; i<n; i++){
        for(int j=0; j<m; j++){
            double Aij = A[i*m+j];
            if(isnan(Aij)||isinf(Aij)){ printf("ERROR in checkMatrixNaN() %s[%i,%i]=%g\n", name, i,j,Aij); if(bExit)exit(0); return true; };
        }
    }
    return false;
}

double checkDist(int n, const Vec3d* vec, const Vec3d* ref, int verb, double tol ){
    double err=0;
    for(int i=0; i<n; i++){ 
        double ei = (vec[i]-ref[i]).norm();
        err = fmax( err, ei );
        if(verb>1){
           if(ei>tol) printf( "[%i] vec(%20.10f,%20.10f,%20.10f) ref(%20.10f,%20.10f,%20.10f)\n", i, vec[i].x,vec[i].y,vec[i].z,  ref[i].x,ref[i].y,ref[i].z ); 
        }
    }; 
    if( verb>0 )printf( "Err| v - ref |=%g\n", err );
    return err;
}

void fitAABB( Quat8d& bb, int n, int* c2o, Quat4d * ps ){
    //Quat8d bb;
    //bb.lo = bb.lo = ps[c2o[0]];
    for(int i=0; i<n; i++){ 
        //printf( "fitAABB() i %i \n", i );
        int ip = c2o[i];
        //printf( "fitAABB() i=%i ip=%i \n", i, ip );
        Quat4d p = ps[ip];
        bb.lo.f.setIfLower  ( p.f );
        bb.hi.f.setIfGreater( p.f );
    }; 
    //return bb;
}

void fitAABB_edge( Quat8d& bb, int n, int* c2o, int2* edges, Quat4d * ps ){
    //Quat8d bb;
    //bb.lo = bb.lo = ps[c2o[0]];
    for(int i=0; i<n; i++){ 
        //printf( "fitAABB() i %i \n", i );
        int ie = c2o[i];
        int2 b = edges[ie];
        //printf( "fitAABB() i=%i ip=%i \n", i, ip );
        Quat4d pi = ps[b.x];
        bb.lo.f.setIfLower  ( pi.f );
        bb.hi.f.setIfGreater( pi.f );
        Quat4d pj = ps[b.y];
        bb.lo.f.setIfLower  ( pj.f );
        bb.hi.f.setIfGreater( pj.f );
    }; 
    //return bb;
}

void updatePointBBs(const Buckets& buckets, Quat8d* BBs, Quat4d* points, bool bInit){
    //printf( "updatePointBBs() START \n" );
    for(int ib=0; ib<buckets.ncell; ib++){
        //printf( "updatePointBBs() ib %i \n", ib );
        if(bInit){ BBs[ib].lo.f = Vec3dMax; BBs[ib].hi.f = Vec3dMin; }
        int n = buckets.cellNs[ib];
        if(n>0){
            int i0 = buckets.cellI0s[ib];
            //printf( "updatePointBBs() ib %i n %i i0 %i \n", ib, n, i0 );
            fitAABB( BBs[ib], n, buckets.cell2obj+i0, points );
        }
    }
    //printf( "updatePointBBs() DONE \n" );
}

void updateEdgeBBs(const Buckets& buckets, Quat8d* BBs, int2* edges, Quat4d* points, bool bInit){
    //printf( "updateEdgeBBs() START \n" );
    for(int ib=0; ib<buckets.ncell; ib++){
        if(bInit){ BBs[ib].lo.f = Vec3dMax; BBs[ib].hi.f = Vec3dMin; }
        int n = buckets.cellNs[ib];
        if(n>0){
            int i0 = buckets.cellI0s[ib];
            fitAABB_edge( BBs[ib], n, buckets.cell2obj+i0, edges, points );
        }
    }
    //printf( "updateEdgeBBs() DONE \n" );
}

void print_vector( int n, double * a, int pitch, int j0, int j1 ){
    for(int i=0; i<n; i++){
        double* ai = a + i*pitch;	
        for (int j=j0; j<j1; j++ ){	printf( "%f ", ai[j] );	}
    }    
    printf( "\n" );
}


// ============ Functions of TrussDynamics_d ============

    void TrussDynamics_d::recalloc( int nPoint_, int nNeighMax_, int nBonds_){
        nPoint = nPoint_; nNeighMax = nNeighMax_;
        nNeighTot = nPoint*nNeighMax;
        _realloc0( points, nPoint   , Quat4dZero );
        _realloc0( forces, nPoint   , Quat4dZero );
        _realloc0( vel,    nPoint   , Quat4dZero );
        _realloc0( vel0,    nPoint   , Quat4dZero );
        _realloc0( bvec,   nPoint   , Quat4dZero );
        _realloc0( bvec0,  nPoint   , Vec3dZero  );
        _realloc0( ps_0,   nPoint   , Vec3dZero  );
        _realloc0( ps_cor       ,nPoint, Vec3dZero );
        _realloc0( ps_pred      ,nPoint, Vec3dZero );
        _realloc0( linsolve_yy , nPoint, Vec3dZero ); // used also to store momentum in MomentumDiff solver
        //_realloc( impuls, nPoint  , Quat4dZero ); );
        _realloc0( params,  nNeighTot , Quat4dZero );
        _realloc0( neighs,  nNeighTot , 0 );
        _realloc0( neighBs, nNeighTot , {0,0} );
        _realloc0( neighB2s,nNeighTot , 0 );

        if(nBonds_>0){
            nBonds=nBonds_;
            _realloc0( bonds,     nBonds, {0,0} );
            _realloc0( strain,    nBonds, 0.0 );
            //_realloc( l0s,       nBonds );
            _realloc0( maxStrain, nBonds, {0.0,0.0} );
            _realloc0( bparams,   nBonds, Quat4dZero );
        }
    }

    void TrussDynamics_d::realloc_LinearSystem( bool bCG, bool bCholesky, int nNeighMaxLDLT_, bool bDens ){
        printf( "TrussDynamics_d_d::realloc_LinearSystem() nPoint=%i nNeighMaxLDLT=%i nNeighMax=%i\n", nPoint, nNeighMaxLDLT_, nNeighMax );
        nNeighMaxLDLT = nNeighMaxLDLT_;  if(nNeighMaxLDLT<nNeighMax){ printf("ERROR in TrussDynamics_d::prepare_Cholesky(): nNeighMaxLDLT(%i)<nNeighMax(%i) \n", nNeighMaxLDLT, nNeighMax); exit(0); }
        int n2 = nPoint*nPoint;
        // --- sparse Linear system Matrix and its Cholesky L*D*L^T decomposition
        if(bDens){ _realloc0( PDmat,  n2,     0.0 ); }
        _realloc0( linsolve_b ,nPoint, Vec3dZero );
        if(bCholesky){
            _realloc0( LDLT_L,      n2,     0.0             );
            _realloc0( LDLT_D,      nPoint, 0.0             );
            _realloc0( neighsLDLT, nPoint*nNeighMaxLDLT, -1 );
        }
        if(bCG){
            cgSolver.realloc(nPoint,3);
        }
    }

    void TrussDynamics_d::realloc_lin( bool bLinSolver){
        _realloc0( kFix, nPoint, 0.0 );
        if(bLinSolver){
            _realloc( kDirs, nBonds );
            //_realloc( f0s,   nBonds );
        }else{
            _realloc( hbs,   nBonds );
            _realloc( dpos,  nBonds );
            _realloc( fdpos, nBonds );
            _realloc( vdpos, nBonds );  
        }
    }

    void TrussDynamics_d::recallocBBs( int n, bool bPoint, bool bEdge, bool bFace, bool bClean ){
        nBBs = n;
        _realloc( BBs, nBBs );
        if(bPoint &&(nPoint>0) ){ pointBBs.realloc( nBBs, nPoint, true ); if(bClean)pointBBs.clean(); }
        if(bEdge  &&(nBonds>0) ){ edgeBBs .realloc( nBBs, nBonds, true ); if(bClean)pointBBs.clean(); }
        if(bFace  &&(nFaces>0) ){ faceBBs .realloc( nBBs, nFaces, true ); if(bClean)pointBBs.clean(); }
         //_realloc( , nBBs, nPoint );
    }


    double TrussDynamics_d::norm_butFixed( Vec3d* ps ){
        double r2 = 0.0;
        for (int i=0; i<nPoint; i++){ 
            if(kFix){ if(kFix[i]>1e-8) continue; }
            r2 += ps[i].norm2();
        }
        return sqrt(r2);;
    }
    double TrussDynamics_d::norm_butFixed( Quat4d* ps ){
        double r2 = 0.0;
        for (int i=0; i<nPoint; i++){ 
            if(kFix){ if(kFix[i]>1e-8) continue; }
            r2 += ps[i].f.norm2();
        }
        return sqrt(r2);;
    }

    // ================= Bounding Boxes

    void TrussDynamics_d::edgesToBBs(){
        // assign to cell with smaller number of points
        for(int i=0; i<nBonds; i++){
            int2 b = bonds[i];
            int ci = pointBBs.obj2cell[b.x];
            int cj = pointBBs.obj2cell[b.y];
            int ni = pointBBs.cellNs[ci]; 
            int nj = pointBBs.cellNs[cj];
            if(ni<nj){ edgeBBs.obj2cell[i]=ci; }else{ edgeBBs.obj2cell[i]=cj; }
        }
    }

    void TrussDynamics_d::printBBs(){
        for(int i=0; i<nBBs; i++){ 
            int np = (pointBBs.cellNs) ? pointBBs.cellNs[i] : 0;
            int ne = (edgeBBs .cellNs) ? edgeBBs .cellNs[i] : -1;
            int nf = (faceBBs .cellNs) ? faceBBs .cellNs[i] : -1;
            printf( "BBs[%i] (%g,%g,%g) (%g,%g,%g) np=%i ned=%i nfc=%i \n", i, BBs[i].lo.x, BBs[i].lo.y, BBs[i].lo.z, BBs[i].hi.x, BBs[i].hi.y, BBs[i].hi.z, np, ne, nf );
        }
    };

    // =================== Picking

    int TrussDynamics_d::pick_point_brute( const Vec3d& ray0, const Vec3d& hray, double Rmax, int i0, int i1 ) const {
        if(i0<0){ i0=0;      }
        if(i1<0){ i1=nPoint; }
        double r2min =  Rmax*Rmax;
        int imin    = -1;
        for(int i=i0; i<i1; i++){
            double t;
            double r2 = rayPointDistance2( ray0, hray, points[i].f, t );
            //printf( "pick_point_brute ipick %i r %g p(%g,%g,%g)\n", i, sqrt(r2), points[i].f.x,points[i].f.y,points[i].f.z );
            if(r2<r2min){ imin=i; r2min=r2; }
            //double ti = raySphere( ray0, hRay, R, ps[i] );
            //if(ti<tmin){ imin=i; tmin=ti; }
        }
        //printf( "pick_point_brute ipick %i r2min %g \n", ipick, r2min );
        return imin;
    }

    int TrussDynamics_d::pick_bond_brute( const Vec3d& ray0, const Vec3d& hRay, double Rmax, int i0, int i1 ) const {
        if(i0<0){ i0=0;      }
        if(i1<0){ i1=nBonds; }
        double dist_min =  Rmax;
        int    imin     = -1;
        for(int ib=i0; ib<i1; ib++){
            int2 b = bonds[ib];
            double t1,t2;
            Vec3d p0 = (Vec3d)points[b.x].f;
            Vec3d d  = (Vec3d)points[b.y].f - p0;
            double l = d.normalize();
            double dist = rayLine( ray0, hRay, p0, d, t1, t2 );
            if( (dist<dist_min) && (t2>0) && (t2<l) ){
                imin=ib; dist_min=dist;
            }
        }
        return imin;
    };

    int TrussDynamics_d::pick_BBox( const Vec3d& ray0, const Vec3d& hRay, double tmax, int i0, int i1 ){
        if(i0<0){ i0=0;     }
        if(i1<0){ i1=nBBs; }
        if(i0<0){ i0=0; i1=nBonds; }
        double tmin = tmax;
        int    imin = -1;
        for(int ib=i0; ib<i1; ib++){
            const Quat8d& bb = BBs[ib];
            double t = rayBox( ray0, hRay, (Vec3d)bb.lo.f, (Vec3d)bb.hi.f, hit_pos, hit_normal );
            if( t<tmin ){ imin=ib; tmin=t; }
        }
        return imin;
    };

    int TrussDynamics_d::pick_nearest(Vec3d ray0, Vec3d hray, int& ipick, int mask, double Rmax ) {
        if     (mask==1){ ipick=pick_point_brute((Vec3d)ray0,(Vec3d)hray,Rmax); return 1; }
        else if(mask==2){ ipick=pick_bond_brute ( ray0, hray, Rmax );           return 2; }
        return -1;
    };
    
    int TrussDynamics_d::pick_all(Vec3d ray0, Vec3d hray, int* out, int mask, double Rmax ) { return 0; };
    

    // =================== Solver Using Projective Dynamics and Cholesky Decomposition

    bool TrussDynamics_d::checkMasses(double mass_tolerance, bool bExit){
        int nUnconnected = 0;
        //printf("TrussDynamics_d::checkMasses() checking %i points with mass_tolerance=%g\n", nPoint, mass_tolerance);
        for(int i=0; i<nPoint; i++){
            if(points[i].w < mass_tolerance){
                printf("WARNING: TrussDynamics_d::checkMasses() point %i has mass %g (< %g), likely unconnected.\n", i, points[i].w, mass_tolerance);
                nUnconnected++;
            }
        }
        if(nUnconnected > 0){
            printf("ERROR in TrussDynamics_d::checkMasses: Found %i points with mass < %g. This would lead to a singular PD matrix. \n", nUnconnected, mass_tolerance);
            if(bExit) exit(0);
            return false;
        }
        return true;
    }

    __attribute__((hot)) 
    void TrussDynamics_d::make_PD_Matrix( double* A, double dt ) {
        printf( "TrussDynamics_d_d::make_PD_Matrix() dt=%g\n", dt );
        //for(int i=0; i<nPoint; i++){ printf("mass[%i] %g\n", points[i].w ); }; // debug
        //for(int i=0; i<nBonds; i++){ printf("ks[%i] %g\n",   params[i].y ); }; // debug
        int n = nPoint;
        int n2 = n*n;
        //for(int i=0; i<np2; i++){ PDmat[i]=0.0; }
        double idt2 = 1.0 / (dt * dt);
        for (int i = 0; i < n; ++i) { // point i
            double   Aii  = points[i].w * idt2; // Assuming points[i].w stores the mass
            if(kFix) Aii += kFix[i];            // fixed point
            //printf( "==i,i %i %g %g %g \n", i, Aii, points[i].w, idt2  );
            int2* ngi = neighBs + (i*nNeighMax);
            for( int ing=0; ing<nNeighMax; ing++ ){
                int2  ng = ngi[ing];
                int   ib = ng.y;
                if(ib<0) break;
                //if(ib<0){ printf("ERROR int TrussDynamics_d_d::make_PD_Matrix() invalid neighbor [i=%i,ing=%i] ng(%i,%i) neigh=%i\n", i, ing, ib, ng.x, neighs[i*nNeighMax+ing] ); exit(0); }  
                //printf( "make_PD_Matrix()[i=%i,ing=%i] ng(%i,%i)\n", i, ing, ib, ng.x );
                int2   b = bonds[ib];
                //int j  = neighs[b.x];
                int j    = b.y;
                double k = getLinearBondStiffness( ib );  

                //printf( "i,ib %i %i %g \n", i, ib+1, k  );

                Aii += k;
                // if (j > i) {
                //     A[i+j*n] = -k;
                //     A[j+i*n] = -k;
                // }

                if     ( b.y > i ){
                    A[i  *n+b.y] = -k;
                    A[b.y*n+i  ] = -k;
                }else if ( b.x > i ){
                    A[i  *n+b.x] = -k;
                    A[b.x*n+i  ] = -k;
                }
            }
            A[i+i*n] += Aii;
        }
    }

    __attribute__((hot)) 
    void TrussDynamics_d::make_PDmat_sparse( SparseMatrix<double>& A, double dt, bool bRealloc ) {
        printf( "TrussDynamics_d::make_PDmat_sparse() dt=%g nNeighMax=%i \n", dt, nNeighMax );
        double idt2 = 1.0 / (dt * dt);
        // --- count neighbors
        if( bRealloc ){ A.realloc( nPoint, nNeighMax+1 ); }
        
        for (int i=0; i<nPoint; i++) {
            double   Aii  = points[i].w * idt2; 
            if(kFix) Aii += kFix[i];  // fixed point
            int2* ngi = neighBs + (i*nNeighMax);
            for( int ing=0; ing<nNeighMax; ing++ ){
                int2  ng = ngi[ing];
                int   ib = ng.y;
                if(ib<0) break;
                int2   b = bonds[ib];
                int j    = b.y;
                double k = getLinearBondStiffness( ib );  
                if(isnan(k)||isinf(k)){ printf("ERROR in TrussDynamics_d::make_PDmat_sparse() [i=%i,ing=%i] ib=%i k=%g\n", i, ing, ib, k); exit(0); }
                //printf( "i,ib %i %i %g \n", i, ib+1, k  );
                Aii += k;
                if       ( b.y > i ){
                    A.set( i,   b.y, -k );
                    A.set( b.y, i  , -k );
                }else if ( b.x > i ){
                    A.set( i,   b.x, -k );
                    A.set( b.x, i  , -k );
                }
            }
            A.set(i,i, Aii );
        }
    }

    __attribute__((hot)) 
    void TrussDynamics_d::rhs_ProjectiveDynamics(const Vec3d* pnew, Vec3d* b) {
        double idt2 = 1.0 / (dt * dt);
        for (int i=0; i<nPoint; i++) {
            double   Ii  = points[i].w*idt2; 
            if(kFix) Ii += kFix[i];
            Vec3d bi; bi.set_mul(pnew[i], Ii );  // points[i].w is the mass
            int2* ngi = neighBs + (i*nNeighMax);
            int ni = 0;
            for( int ing=0; ing<nNeighMax; ing++ ){
                int2  ng = ngi[ing];
                int   ib = ng.y;
                if(ib<0) break;
                double k  = getLinearBondStiffness( ib );  
                double l0 = bparams[ib].x;
                int2   b  = bonds[ib];
                int j     = (b.x == i) ? b.y : b.x;
                //printf( "rhs[i=%2i,ib=%2i,j=%i] k=%g l0=%g\n", i, ib, j, k, l0 );
                Vec3d d = pnew[i] -  pnew[j];
                bi.add_mul(d, k * l0 / d.norm());  // params[ib].x is l0
                ni++;
            }
            //printf( "rhs[i=%2i] ni=%i bi(%g,%g,%g)\n", i, ni, bi.x,bi.y,bi.z );
            b[i] = bi;
        }
    }

    __attribute__((hot)) 
    Vec3d TrussDynamics_d::rhs_ProjectiveDynamics_i(int i, const Vec3d* pnew) {
        double idt2 = 1.0 / (dt * dt);
        double   Ii  = points[i].w*idt2; 
        if(kFix) Ii += kFix[i];
        Vec3d bi; bi.set_mul(pnew[i], Ii );  // points[i].w is the mass
        int2* ngi = neighBs + (i*nNeighMax);
        int ni = 0;
        for( int ing=0; ing<nNeighMax; ing++ ){
            int2  ng = ngi[ing];
            int   ib = ng.y;
            if(ib<0) break;
            double k  = getLinearBondStiffness( ib );  
            double l0 = bparams[ib].x;
            int2   b  = bonds[ib];
            int j     = (b.x == i) ? b.y : b.x;
            //printf( "rhs[i=%2i,ib=%2i,j=%i] k=%g l0=%g\n", i, ib, j, k, l0 );
            Vec3d d = pnew[i] -  pnew[j];
            bi.add_mul(d, k * l0 / d.norm());  // params[ib].x is l0
            ni++;
        }
        //printf( "rhs[i=%2i] ni=%i bi(%g,%g,%g)\n", i, ni, bi.x,bi.y,bi.z );
        return bi;
    }

    __attribute__((hot)) 
    void TrussDynamics_d::rhs_ProjectiveDynamics_(const Vec3d* pnew, Vec3d* b) {
        double idt2 = 1.0 / (dt * dt);
        for (int i=0; i<nPoint; i++) {
            b[i] = rhs_ProjectiveDynamics_i(i,pnew);
        }
    }

    void TrussDynamics_d::updatePD_RHS( const Vec3d* pnew, Quat4d* bvec ){  // RHS for projective dynamics
        for(int i=0; i<nPoint; i++){
            const Vec3d    pi = pnew[i];
            const double idt2 = 1.0 / (dt * dt);
            double         Ii = points[i].w*idt2; 
            if(kFix)       Ii += kFix[i];
            Vec3d bi; bi.set_mul(pi, Ii );   // b_i    =  M_i/dt^2 p'_i   +  \sum_j ( K_{ij} d_{ij} )         
            double Aii = Ii;                 // A_{ii} =  M_i/dt^2        +  \sum_j   K_{ij} 
            const int j0 = i * nNeighMax;
            for(int jj=0; jj<nNeighMax; jj++){
                const int j  = j0 + jj;
                const int jG = neighs[j];
                if(jG == -1) break;
                const Quat4d par = params[j];
                const Vec3d  pj  = pnew[jG];
                const Vec3d  dij = pi - pj;
                const double l   = dij.norm();
                const double k   = par.z;
                bi   .add_mul( dij, k*par.x/l );   // RHS      bi  = \sum_j K_{ij} d_{ij} (+ inertial term M_i/dt^2 p'_i )
                Aii  += k;                         // diagonal Aii =  \sum_j k_{ij}       (+ inertial term M_i/dt^2      )
            }
            bvec[i].f = bi;  // use forces array to store RHS vector
            bvec[i].w = Aii;
        }
    }

    // evaluate f = A p
    void TrussDynamics_d::dotPD( const Vec3d* p, Vec3d* f ){  // RHS for projective dynamics
        for(int i=0; i<nPoint; i++){
            const Vec3d    pi = p[i];
            const double idt2 = 1.0 / (dt * dt);
            double         Ii = points[i].w*idt2; 
            if(kFix)       Ii += kFix[i];
            Vec3d fi           = Vec3dZero;   // sum_j  =  \sum_j ( K_{ij} p_j} )   + Aii * pi        
            double Aii  = Ii;                 // A_{ii} =  \sum_j   K_{ij}          + M_i/dt^2 
            const int j0 = i * nNeighMax;
            for(int jj=0; jj<nNeighMax; jj++){
                const int j  = j0 + jj;
                const int jG = neighs[j];
                if(jG == -1) break;
                const Quat4d par = params[j];
                const Vec3d  pj  = p[jG];
                const double k   = par.z;
                fi.add_mul( pj, -k ); 
                Aii  += k;                  
            }
            fi .add_mul( pi, Aii );
            f[i] = fi;
        }
    }

    void TrussDynamics_d::updatePD_dRHS( const Vec3d* pnew, Quat4d* bvec ){  // RHS for projective dynamics
        // Calculates db = b - A p0 ( i.e. combines dotPD and updatePD_RHS into a single loop) 
        for(int i=0; i<nPoint; i++){
            const Vec3d    pi = pnew[i];
            const double idt2 = 1.0 / (dt * dt);
            double         Ii = points[i].w*idt2; 
            if(kFix)       Ii += kFix[i];
            Vec3d bi = Vec3dZero; 
            //bi.set_mul(pi,  Ii );   // b_i    =  M_i/dt^2 p'_i   +  \sum_j ( K_{ij} d_{ij} )   
            //bi.add_mul(pi, -Ii ); 
            double Aii = Ii;                 // A_{ii} =  M_i/dt^2        +  \sum_j   K_{ij} 
            const int j0 = i * nNeighMax;
            for(int jj=0; jj<nNeighMax; jj++){
                const int j  = j0 + jj;
                const int jG = neighs[j];
                if(jG == -1) break;
                const Quat4d par = params[j];
                const Vec3d  pj  = pnew[jG];
                const Vec3d  dij = pi - pj;
                const double l   = dij.norm();
                const double k   = par.z;
                //bi   .add_mul( pi, -k ); 
                //bi   .add_mul( pj,  k ); 
                //bi   .add_mul( dij,-k );
                //bi   .add_mul( dij, k*par.x/l );   // RHS      bi  = \sum_j K_{ij} d_{ij} (+ inertial term M_i/dt^2 p'_i )
                bi   .add_mul( dij, k*(par.x/l - 1) ); 
                Aii  += k;                         // diagonal Aii =  \sum_j k_{ij}       (+ inertial term M_i/dt^2      )
            }
            //bi .add_mul( pi, -Aii );
            bvec[i].f = bi;  // use forces array to store RHS vector
            bvec[i].w = Aii;
        }
    }

    void TrussDynamics_d::updateJacobi_lin( Vec3d* ps_in, Vec3d* ps_out, Quat4d* bvec, Vec3d* rs ){
        for(int i=0; i<nPoint; i++){
            const int j0 = i * nNeighMax;
            Vec3d sum_j  = Vec3dZero;
            for(int jj=0; jj<nNeighMax; jj++){
                const int j  = j0 + jj;
                const int jG = neighs[j];
                if(jG == -1) break;
                sum_j.add_mul( ps_in[jG],  params[j].z );   // kPull
            }
            const Quat4d bi     = bvec[i];  // forces array stores RHS vector
            if(ps_out ) ps_out[i] = (bi.f + sum_j       )*(1/bi.w);
            if(rs     ) rs    [i] = (bi.f + sum_j - ps_in[i]*bi.w );
        }
    }

    void TrussDynamics_d::updateGaussSeidel_lin( Vec3d* ps, Quat4d* bvec ){
        for(int i=0; i<nPoint; i++){
            const int j0 = i * nNeighMax;
            Vec3d sum_j  = Vec3dZero;
            for(int jj=0; jj<nNeighMax; jj++){
                const int j  = j0 + jj;
                const int jG = neighs[j];
                if(jG == -1) break;
                sum_j.add_mul( ps[jG],  params[j].z );   // kPull
            }
            const Quat4d bi     = bvec[i];  // forces array stores RHS vector
            ps[i]  = (bi.f + sum_j)*(1/bi.w);
        }
    }

    void TrussDynamics_d::evalTrussForce( Vec3d* ps_in, Vec3d* force ){
        // This solver calculates bi and Aii on-the fly, preventing updatePD_RHS call and possibly improving numerical stability (?)
        // NOTE: this is goes beyond linar solution of Ap=b because we update b every iteration, and use non-linear operations like d_{ij}=(pi-pj)/|pi-pj| 
        for(int i=0; i<nPoint; i++){
            const Vec3d    pi = ps_in[i];
            Vec3d          fi = Vec3dZero;
            const int j0 = i * nNeighMax;
            for(int jj=0; jj<nNeighMax; jj++){
                const int j  = j0 + jj;
                const int jG = neighs[j];
                if(jG == -1) break;
                
                const Quat4d par = params[j];
                const Vec3d  pj  = ps_in[jG];
                const Vec3d  dij = pi - pj;
                const double l   = dij.norm();
                const double dl  = l - par.x; 
                const double k   = par.z;
                fi.add_mul( dij, k*dl/l ); 
            }
            force[i]  = fi;
        }
    }

    void TrussDynamics_d::updateJacobi_fly( Vec3d* ps_in, Vec3d* ps_out ){
        // This solver calculates bi and Aii on-the fly, preventing updatePD_RHS call and possibly improving numerical stability (?)
        // NOTE: this is goes beyond linar solution of Ap=b because we update b every iteration, and use non-linear operations like d_{ij}=(pi-pj)/|pi-pj| 
        for(int i=0; i<nPoint; i++){
            const Vec3d    pi = ps_in[i];
            const double idt2 = 1.0 / (dt * dt);
            double         Ii = points[i].w*idt2; 
            if(kFix)       Ii += kFix[i];
            Vec3d sum_j; sum_j.set_mul(pi, Ii );   // b_i    =  M_i/dt^2 p'_i   +  \sum_j ( K_{ij} d_{ij} )         
            double Aii = Ii;                       // A_{ii} =  M_i/dt^2        +  \sum_j   K_{ij} 
            const int j0 = i * nNeighMax;
            for(int jj=0; jj<nNeighMax; jj++){
                const int j  = j0 + jj;
                const int jG = neighs[j];
                if(jG == -1) break;
                
                //   update RHS bi
                const Quat4d par = params[j];
                const Vec3d  pj  = ps_in[jG];
                const Vec3d  dij = pi - pj;
                const double l   = dij.norm();
                const double k   = par.z;
                sum_j.add_mul( dij, k*par.x/l );   //  b_i  +=  \sum_j ( K_{ij} d_{ij} )   
                sum_j.add_mul( pj , k         );   //  s_j  +=  \sum_j ( K_{ij} p_j    )
                Aii  += k;

            }
            sum_j.mul( 1/Aii );
            ps_out[i]  = sum_j;
        }
    }

    void TrussDynamics_d::updateGaussSeidel_fly( Vec3d* ps ){
        // This solver calculates bi and Aii on-the fly, preventing updatePD_RHS call and possibly improving numerical stability (?)
        // NOTE: this is goes beyond linar solution of Ap=b because we update b every iteration, and use non-linear operations like d_{ij}=(pi-pj)/|pi-pj| 
        for(int i=0; i<nPoint; i++){
            const Vec3d    pi = ps[i];
            const double idt2 = 1.0 / (dt * dt);
            double         Ii = points[i].w*idt2; 
            if(kFix)       Ii += kFix[i];
            Vec3d sum_j; sum_j.set_mul(pi, Ii );   // b_i    =  M_i/dt^2 p'_i   +  \sum_j ( K_{ij} d_{ij} )         
            double Aii = Ii;                       // A_{ii} =  M_i/dt^2        +  \sum_j   K_{ij} 
            const int j0 = i * nNeighMax;
            for(int jj=0; jj<nNeighMax; jj++){
                const int j  = j0 + jj;
                const int jG = neighs[j];
                if(jG == -1) break;
                
                //   update RHS bi
                const Quat4d par = params[j];
                const Vec3d  pj  = ps[jG];
                const Vec3d  dij = pi - pj;
                const double l   = dij.norm();
                const double k   = par.z;
                sum_j.add_mul( dij, k*par.x/l );   //  b_i  +=  \sum_j ( K_{ij} d_{ij} )   
                sum_j.add_mul( pj , k         );   //  s_j  +=  \sum_j ( K_{ij} p_j    )
                Aii  += k;

            }
            sum_j.mul( 1/Aii );
            ps[i]  = sum_j;
        }
    }

    void TrussDynamics_d::updateIterativeMomentum( Vec3d* psa, Vec3d* psb ){
        //Vec3d* psa = ps_in;
        //Vec3d* psb = ps_out;
        //if(LinSolveMethod::JacobiMomentum == (LinSolveMethod)linSolveMethod  ){
        updatePD_RHS(psa, bvec ); // TODO: we need to calculate this only whem we use linear solvers, not if we use _fly solvers
        //}
        for (int i = 0; i < nSolverIters; i++) {  
            switch( (LinSolveMethod)linSolveMethod ){
                case LinSolveMethod::JacobiMomentum:    { updateJacobi_lin( psa, psb, bvec ); } break;
                case LinSolveMethod::JacobiFlyMomentum: { updateJacobi_fly( psa, psb );       } break;
                case LinSolveMethod::GSMomentum:     { 
                    for (int j=0; j<nPoint; j++){ psb[j]=psa[j]; }
                    updateGaussSeidel_lin( psb, bvec ); 
                } break; 
                case LinSolveMethod::GSFlyMomentum:     { 
                    for (int j=0; j<nPoint; j++){ psb[j]=psa[j]; }
                    updateGaussSeidel_fly( psb ); 
                } break; 
            }
            //updateJacobi_lin( psa, psb, bvec );      //    p'_{k+1}    = Solve(A,p_k,b)
            double bmix = mixer.get_bmix( i );       //    bmix_k = bmix(k)
            if( (i==0) || (i>=(nSolverIters-1) ) ) bmix = 0.0; // make sure we stop momentum mixing for the last iteration      
            for(int j=0; j<nPoint; j++){ 
                Vec3d p = psb[j]; 
                p.add_mul( linsolve_yy[j], bmix );   //    p_{k+1} = p'_k + bmix d_k
                linsolve_yy[j] = p - psa[j];         //    d_{k+1} = p_{k+1} - p_k
                psa[j]         = p;
            }   // we use linsolve_yy to store momentum
            //for (int i=0; i<nPoint; i++){ ps_pred[i]=ps_cor[i]; }
        }

        if( bApplyResudualForce ){
            double dt2 = dt * dt;
            evalTrussForce( psa, linsolve_b );
            //updateJacobi_lin( psa, 0, bvec, linsolve_b );
            // for (int i=0; i<nPoint; i++){ 
            //     if(kFix){ if(kFix[i]>1e-8) continue; }
            //     psa[i].add_mul( linsolve_b[i], residualForceFactor*dt2/points[i].w ); 
            // }
        }

        for (int i=0; i<nPoint; i++){ psb[i]=psa[i]; }
    }

    void TrussDynamics_d::updateIterativeJacobi( Vec3d* psa, Vec3d* psb ){
        updatePD_RHS(psa, bvec );
        for (int i = 0; i < nSolverIters; i++) {  
            updateJacobi_lin( psa, psb, bvec ); 
            for (int i=0; i<nPoint; i++){ psa[i]=psb[i]; }
        }
    }

    // Here we replace Ap=b   by    A(p-p0) = b-Ap0
    void TrussDynamics_d::updateIterativeJacobiDiff( Vec3d* psa, Vec3d* psb ){
        //printf( "updateIterativeJacobiDiff() \n" );
        const bool bUse_dRHS = true; 
        //const bool bUse_dRHS = false; 
        if(bUse_dRHS){
            updatePD_dRHS( psa, bvec  );    // b = Kd - Ap' = \sum_j{ k_{ij} ( l0_{ij}/l_{ij} - 1 ) ( p_i - p_j )
        }else{
            dotPD       ( psa, bvec0 );     // b0 = Ap0
            updatePD_RHS( psa, bvec  );     // b = Kd + Ip'
            for (int i=0; i<nPoint; i++){ bvec[i].f.sub(bvec0[i]);  }       //  db = b - b0 
        }
    //    updatePD_RHS( psa, bvec  );    // b = Kd + Ip'
    //    dotPD       ( psa, bvec0 );     // b0 = Ap0
    //    for (int i=0; i<nPoint; i++){ bvec0[i]= bvec[i].f - bvec0[i];  } 
    //    updatePD_dRHS( psa, bvec  ); 
    //    for (int i=0; i<nPoint; i++){ printf( "bvec %i |dr|: %12.6e b(%12.6e,%12.6e,%12.6e) b0(%12.6e,%12.6e,%12.6e)\n", i, (bvec[i].f-bvec0[i]).norm(), bvec[i].x, bvec[i].y, bvec[i].z, bvec0[i].x, bvec0[i].y, bvec0[i].z ); } 
    //    exit(0);
        for (int i=0; i<nPoint; i++){ ps_0[i]=psa[i]; psa[i]=Vec3dZero; }   //  dp = p - p0      
        for (int i = 0; i < nSolverIters; i++) {  
            updateJacobi_lin( psa, psb, bvec );
            for (int i=0; i<nPoint; i++){ psa[i]=psb[i]; }
        }
        for (int i=0; i<nPoint; i++){ psb[i].add(ps_0[i]); }   //  p = p0 + dp 
        // Debug:
        //    double db_norm  = norm_butFixed( bvec );
        //    double dp_norm  = norm_butFixed( psa );
        //    printf( "updateIterativeJacobiDiff() p0_norm: %g b0_norm: %g b_norm: %g db_norm: %g dp_norm: %g \n", p0_norm, b0_norm, b_norm, db_norm, dp_norm );
    }

    void TrussDynamics_d::updateIterativeMomentumDiff( Vec3d* psa, Vec3d* psb ){
        const bool bUse_dRHS = true; 
        updatePD_dRHS( psa, bvec  );    // b = Kd - Ap' = \sum_j{ k_{ij} ( l0_{ij}/l_{ij} - 1 ) ( p_i - p_j )
        for (int i=0; i<nPoint; i++){ ps_0[i]=psa[i]; psa[i]=Vec3dZero; }   //  dp = p - p0 , start with zero      
        for (int i = 0; i < nSolverIters; i++) {  
            updateJacobi_lin( psa, psb, bvec );
            double bmix = mixer.get_bmix( i );       //    bmix_k = bmix(k)
            if( (i==0) || (i>=(nSolverIters-1) ) ) bmix = 0.0; // make sure we stop momentum mixing for the last iteration      
            for(int j=0; j<nPoint; j++){ 
                Vec3d p = psb[j]; 
                p.add_mul( linsolve_yy[j], bmix );   //    p_{k+1} = p'_k + bmix d_k
                linsolve_yy[j] = p - psa[j];         //    d_{k+1} = p_{k+1} - p_k
                psa[j]         = p;
            }
        }
        for (int i=0; i<nPoint; i++){ psb[i].add(ps_0[i]); }   //  p = p0 + dp 
    }

    // Here we replace Ap=b   by    A(p-p0) = b-Ap0
    void TrussDynamics_d::updateIterativeExternDiff( Vec3d* psa, Vec3d* psb ){
        dotPD       ( psa, bvec0 );     // b0 = Ap0
        updatePD_RHS( psa, bvec  );     // b = Kd + Ip'
        for (int i=0; i<nPoint; i++){ Quat4d b=bvec[i]; b.f.sub(bvec0[i]); extern_b[i]=(Quat4f)b; }   //  db = b - b0 
        for (int i=0; i<nPoint; i++){ extern_x[i] = Quat4fZero;           }       //  dp = p - p0   
        extern_solve();
        for (int i=0; i<nPoint; i++){ psb[i] = (Vec3d)extern_x[i].f + psa[i]; }   //  p = p0 + dp 
    }

    /*    double run_updateJacobi_smart( int psa, int psb, int itr ) {
        //run_updateJacobi_neighs(psa, psb);
        // I want some function of form     A/( B + itr ) which interpolates between  bmix_end and bmix_start when itr transitions from nitr_end to nitr_start
        if(itr<nitr_start){ run_updateJacobi_neighs(psa, psb);  }
        else              { 
            if ( itr > nitr_end ) { bmix.y = bmix_end; }
            else                  { bmix.y = bmix_start + (bmix_end - bmix_start ) * (itr - nitr_start) / (nitr_end - nitr_start); }
            run_updateJacobi_mix(psa, psb); 
        }
        
       return bmix.y;
    }
    */

    void TrussDynamics_d::prepare_LinearSystem( bool bRealloc, bool bCG, bool bCholesky, int nNeighMaxLDLT_, bool bDens ){ 
        double dt = this->dt;
        printf( "TrussDynamics_d::prepare_LinearSystem() nPoint=%i dt=%g nNeighMaxLDLT_=%i\n", nPoint, dt, nNeighMaxLDLT_ );
        //nNeighMaxLDLT=nNeighMaxLDLT_;
        if(bRealloc)realloc_LinearSystem( bCG, bCholesky, nNeighMaxLDLT_, bDens );

        this->dt = dt;
        int n2 = nPoint*nPoint;

        //long t0;
        //t0=getCPUticks(); make_PD_Matrix(    PDmat,    dt );         printf("@time make_PD_Matrix()    t= %g [MTicks] \n", (getCPUticks()-t0)*1.0e-6 );
        //t0=getCPUticks(); make_PDmat_sparse( PDsparse, dt, true );   printf("@time make_PD_Matrix()    t= %g [MTicks] \n", (getCPUticks()-t0)*1.0e-6 );

        timeit( "TIME make_PDmat_sparse() t= %g [MTicks]\n", 1e-6, [&](){ make_PDmat_sparse( PDsparse, dt, true ); });
        mat2file<int>   ( "PDsparse_inds.log",  nPoint, nNeighMax+1, (int*)   PDsparse.inds, "%5i " );
        mat2file<double>( "PDsparse_vals.log",  nPoint, nNeighMax+1, (double*)PDsparse.vals );

        if(PDsparse.checkNaN(true)){ printf("ERROR in TrussDynamics_d::prepare_LinearSystem() PDsparse contains NaN => exit \n"); exit(0); }

        if(bDens){
            timeit( "TIME make_PD_Matrix()    t= %g [MTicks]\n", 1e-6, [&](){ make_PD_Matrix(    PDmat,    dt );       });
            if( PDsparse.checkDens( PDmat, 1e-16, true            ) ){ printf("ERROR in TrussDynamics_d::prepare_LinearSystem() PDsparse does not match (dense)PDmat => exit \n"); exit(0); }
            if( checkMatrixNaN(nPoint, nPoint,PDmat, "PDmat", true) ){ printf("ERROR in TrussDynamics_d::prepare_LinearSystem() (dense)PDmat contains NaN => exit \n"); exit(0); }
        }

        // { // Check M.dot(vec) sparse
        //     for(int i=0; i<nPoint; i++){ ps_pred[i] = points[i].f *1.1; }; 
        //     dotM_ax                   ( nPoint,3, PDmat, (double*)ps_pred, (double*)linsolve_b );
        //     PDsparse.dot_dens_vector_m(   3,             (double*)ps_pred, (double*)ps_cor     );
        //     for(int i=0; i<nPoint; i++){ printf( "[%i] PDmat*v(%20.10f,%20.10f,%20.10f) PDsparse(%20.10f,%20.10f,%20.10f)\n", i, linsolve_b[i].x,linsolve_b[i].y,linsolve_b[i].z,  ps_cor[i].x,ps_cor[i].y,ps_cor[i].z ); }; 
        //     exit(0);
        // }

        if( bCG ){
            // setup Conjugate-Gradient solver
            cgSolver.setLinearProblem(  nPoint, 3, (double*)ps_cor, (double*)linsolve_b, false );
            cgSolver.initDiagPrecond( PDmat );
            //cgSolver.dotFunc = [&](int n, double* x, double* y){  dotM_ax( n,3, PDmat, x, y ); };
            //cgSolver.dotFunc = [&](int n, double* x, double* y){  PDsparse.dot_dens_vector_m( 3, x, y ); };

            cgSolver.dotFunc = [&](int n, double* x, double* y){  
                long t0 = getCPUticks();
                PDsparse.dot_dens_vector_m( 3, x, y ); 
                time_cg_dot += (getCPUticks() - t0)*1e-6; 
            };

            // {   printf( "# ---------- Test Conjugate-Gradient ");
            //     cgSolver.dotFunc = [&](int n, double* x, double* y){   dotM_ax( n,1, PDmat, x, y );  };
            //     for(int i=0; i<nPoint; i++){ ps_pred[i] = points[i].f + vel[i].f*dt; ps_cor[i]=ps_pred[i];  }; 
            //     rhs_ProjectiveDynamics( ps_pred, linsolve_b );
            //     for(int i=0; i<nPoint; i++){ 
            //         ((double*)linsolve_b)[i] = linsolve_b[i].x;
            //         ((double*)ps_cor    )[i] = ps_cor    [i].x;
            //     }
            //     cgSolver.solve_m1(1e-6,2);
            //     //exit(0);
            // }
        }
        
        if( bCholesky ){
            for(int i=0; i<nPoint; i++){  for(int j=0; j<nNeighMaxLDLT; j++){ if(j>=nNeighMax){neighsLDLT[i*nNeighMaxLDLT+j]=-1;}else{ neighsLDLT[i*nNeighMaxLDLT+j]=neighs[i*nNeighMax+j]; };  } }
            
            bool bSparse = false;
            //bool bSparse = true;
            if(bSparse){
                Lingebra::CholeskyDecomp_LDLT_sparse( PDmat, LDLT_L, LDLT_D, neighsLDLT, nPoint, nNeighMaxLDLT );
                sortNeighs( nPoint, nNeighMaxLDLT, neighsLDLT);
            }else{
                //Lingebra::CholeskyDecomp_LDLT( PDmat, LDLT_L, LDLT_D, nPoint, true );
                //LsparseT.fromDense( nPoint, LDLT_L, 1.e-16, true );
                //for(int i=0; i<n2; i++){ LDLT_L[i]=0; }
                Lingebra::CholeskyDecomp_LDLT( PDmat, LDLT_L, LDLT_D, nPoint );

                if( checkMatrixNaN(nPoint, nPoint, LDLT_L, "LDLT_L") || checkMatrixNaN(nPoint, 1, LDLT_D, "LDLT_D") ){ exit(0); }


                //Lsparse.fromDense( nPoint, LDLT_L, 1.e-16 );
                Lsparse.fromDense   ( nPoint, LDLT_L, 1e-300 );
                LsparseT.fromFwdSubT( nPoint, LDLT_L);
                //mat2file<double>( "PDmat.log",  nPoint,nPoint, PDmat  );
                //mat2file<double>( "LDLT_L.log", nPoint,nPoint, LDLT_L );

                

            }
           

            // {   printf( "# ---------- Test Cholesky Lsparse \n");
            //     for(int i=0; i<nPoint; i++){ ps_pred[i] = points[i].f*1.1; }; 
            //     rhs_ProjectiveDynamics( ps_pred, linsolve_b );
            //     Lingebra::forward_substitution_m( LDLT_L, (double*)linsolve_b,  (double*)linsolve_yy, nPoint, 3 );
            //     Lsparse.fwd_subs_m( 3,  (double*)linsolve_b,  (double*)ps_cor );
            //     checkDist( nPoint, ps_cor, linsolve_yy, 2 );
            // }

            //mat2file<double>( "LDLT_L.log", nPoint,nPoint, LDLT_L );
            Lsparse.fprint_inds("Lsparse_inds.log");
            //Lsparse.fprint_vals("Lsparse_vals.log");

            LsparseT.fprint_inds("LsparseT_inds.log");
            //LsparseT.fprint_vals("LsparseT_vals.log");
            //exit(0);

            double* L_reconstructed  = Lsparse .reconstruct_dense( );
            double* LT_reconstructed = LsparseT.reconstruct_dense( );
            writeMatrix("LDLT_L_dense.txt",         nPoint, nPoint,  LDLT_L, false);
            writeMatrix("LDLT_LT_reconstructed.txt", nPoint, nPoint, LT_reconstructed, false);
            writeMatrix("LDLT_L_reconstructed.txt", nPoint, nPoint,  L_reconstructed, false);
            delete[] L_reconstructed;
            delete[] LT_reconstructed;
        }


    }

    void TrussDynamics_d::dampPoints( double Vdamping ){
        for (int i : damped_points){
            Vec3d vi  = vel [i].f;
            Vec3d vi0 = vel0[i].f;
            Vec3d acc = (vi0-vi) * Vdamping;
            vel0[i].f.add_mul( acc, -dt );
            forces[i].f.add_mul( acc, points[i].w );
        }
    }

    __attribute__((hot)) 
    void TrussDynamics_d::run_LinSolve(int niter) {
        printf( "TrussDynamics_d::run_LinSolve()  linSolveMethod=%i nSolverIters=%i ps_cor=%p points=%p \n", linSolveMethod, nSolverIters, ps_cor, points  );
        const int m=3;
        memcpy(ps_cor, points, nPoint * sizeof(Vec3d));
        double dt2 = dt * dt;
        double inv_dt = 1/dt;
        double cdamp = fmax( 1-damping*dt, 0. );
        cleanForce();

        for (int iter = 0; iter < niter; iter++) {
            // Evaluate forces (assuming you have a method for this)
            
            // ---- Update Forces
            cleanForce();
            //evalForces();
            if (user_update){ user_update(dt);}
            evalEdgeVerts();
            dampPoints( 2.5 );   // Damping - this is not rigorous - it does not perfectly preserve momentum / angular momentum, resp. points have like double the mass

            //  Predictor step
            for (int i=0;i<nPoint;i++){ 
                forces[i].f.add( getPointForce( i ) );
                //forces[i].f.add_mul( Gravity, points[i].w ); // Grafivity
                //forces[i].f.add_mul( vel[i].f, Cdrag ); // Drag force
                ps_pred[i] = points[i].f + vel[i].f*dt + forces[i].f*(dt2/points[i].w); 
                //printf( "ps_pred[%3i](%10.6f,%10.6f,%10.6f) v(%10.6f,%10.6f,%10.6f) p(%10.6f,%10.6f,%10.6f) dt=%g \n", i, ps_pred[i].x,ps_pred[i].y,ps_pred[i].z, vel[i].x,vel[i].y,vel[i].z, points[i].x,points[i].y,points[i].z, dt );
            }

            // Apply fixed constraints
            if(kFix) for (int i = 0; i < nPoint; i++) { if (kFix[i] > 1e-8 ) { ps_pred[i] = points[i].f; } }
            // Compute right-hand side
            
            long t0 = getCPUticks();
            switch( (LinSolveMethod)linSolveMethod ){
                case LinSolveMethod::Jacobi:{
                    //printf( "Jacobi nSolverIters=%i \n", nSolverIters  );
                    updateIterativeJacobi( ps_pred, ps_cor  );
                } break;
                case LinSolveMethod::JacobiDiff:{
                    updateIterativeJacobiDiff( ps_pred, ps_cor );
                } break;
                case LinSolveMethod::MomentumDiff:{
                    updateIterativeMomentumDiff( ps_pred, ps_cor );
                }break;
                case LinSolveMethod::ExternDiff:{
                    updateIterativeExternDiff( ps_pred, ps_cor );
                } break;
                case LinSolveMethod::GaussSeidel:{
                    //printf( "Jacobi nSolverIters=%i \n", nSolverIters  );
                    for (int i=0; i<nPoint; i++){ ps_cor[i]=ps_pred[i]; }
                    updatePD_RHS(ps_cor, bvec );
                    for (int i = 0; i < nSolverIters; i++) {  
                        updateGaussSeidel_lin( ps_cor, bvec ); 
                    }
                } break;
                case LinSolveMethod::JacobiMomentum:
                case LinSolveMethod::GSMomentum: 
                case LinSolveMethod::JacobiFlyMomentum:
                case LinSolveMethod::GSFlyMomentum:
                    updateIterativeMomentum( ps_pred, ps_cor );
                    break;
                case LinSolveMethod::CholeskySparse:{
                    //printf("TrussDynamics_d::run_LinSolve()  LinSolveMethod::CholeskySparse \n");
                    // Solve using LDLT decomposition (assuming you have this method)
                    //solve_LDLT_sparse(b, ps_cor);
                    rhs_ProjectiveDynamics(ps_pred, linsolve_b );
                    Lingebra::forward_substitution_sparse           ( nPoint,m,  LDLT_L, (double*)linsolve_b,  (double*)linsolve_yy, neighsLDLT, nNeighMaxLDLT );
                    for (int i=0; i<nPoint; i++){ linsolve_yy[i].mul(1/LDLT_D[i]); } // Diagonal 
                    Lingebra::forward_substitution_transposed_sparse( nPoint,m,  LDLT_L, (double*)linsolve_yy, (double*)ps_cor,      neighsLDLT, nNeighMaxLDLT );
                } break;
                case LinSolveMethod::Cholesky:{
                    //printf("TrussDynamics_d::run_LinSolve()  LinSolveMethod::Cholesky \n");
                    //Lingebra::forward_substitution_m( LDLT_L, (double*)linsolve_b,  (double*)linsolve_yy, nPoint,m );
                    //Lsparse.fwd_subs_m( m,  (double*)linsolve_b,  (double*)ps_cor );
                    //if( checkDist( nPoint, ps_cor, linsolve_yy, 1 ) ){ printf("ERROR run_LinSolve.checkDist() => exit()"); exit(0); };
                    rhs_ProjectiveDynamics(ps_pred, linsolve_b );
                    Lsparse.fwd_subs_m( m,  (double*)linsolve_b,  (double*)linsolve_yy );
                    for (int i=0; i<nPoint; i++){ linsolve_yy[i].mul(1/LDLT_D[i]); } // Diagonal 
                    LsparseT.fwd_subs_T_m( m,  (double*)linsolve_yy,  (double*)ps_cor );
                    //Lingebra::forward_substitution_T_m( LDLT_L, (double*)linsolve_yy, (double*)ps_cor,      nPoint,m );
                    //if( checkDist( nPoint, ps_pred, ps_cor, 2 ) ){ printf("ERROR run_LinSolve.checkDist() => exit()"); exit(0); };

                } break;
                case LinSolveMethod::CG:{
                    //printf("TrussDynamics_d_d::run_LinSolve()  LinSolveMethod::CG \n");
                    rhs_ProjectiveDynamics(ps_pred, linsolve_b );
                    for(int i=0; i<nPoint; i++){ ps_cor[i]=ps_pred[i]; };
                    cg_tol = 1e-2;
                    cgSolver_niterdone += cgSolver.solve( cg_tol );
                    
                } break;
                case LinSolveMethod::CGsparse:{} break;
            }
            time_LinSolver += (getCPUticks()-t0)*1e-6;
            //mat2file<double>( "points.log",  nPoint,4, (double*)points      );
            //mat2file<double>( "vel.log",     nPoint,4, (double*)vel         );
            //mat2file<double>( "ps_pred.log", nPoint,3, (double*)ps_pred     );
            //mat2file<double>( "b.log",       nPoint,3, (double*)linsolve_b  );
            //mat2file<double>( "yy.log",      nPoint,3, (double*)linsolve_yy );
            //mat2file<double>( "ps_cor.log",  nPoint,3, (double*)ps_cor      );

            //exit(0);

            // Compute residual
            // double res = 0.0;
            // for (int i=0;i<nPoint;i++) {
            //     Vec3d  d = ps_cor[i] - points[i].f;
            //     double l = d.norm();
            //     if (l > res) res = l;
            // }
            // printf("residual[%d] : %f\n", iter, res);
            // Update velocity and points
            
            //   Corrector step
            double l2sum = 0.0;
            for (int i=0;i<nPoint;i++) {
                //Vec3d dv = ps_cor[i] - ps_pred[i];
                //dv.mul(1.0 / dt);
                //vel   [i].f.add( dv );

                // To-Do : We need more rigorous approach how to update velocity

                // position-based velocity update
                //double vr2 = vel[i].norm2();
                Vec3d v    = ps_cor[i] - points[i].f;
                double l2 = v.norm2();
                l2sum += l2;
                v.mul(inv_dt);
                v.mul( cdamp );
                //v.mul( sqrt(  vel[i].norm2()/( v.norm2() + 1e-8 ) ) ); // This is unphysicsl
                if(kFix){ if( kFix[i] > 1e-8){ v=Vec3dZero; } }
                vel[i].f = v;
                // update positions
                points[i].f = ps_cor[i];

                //printf( "STEP: %6i time: %16.8f p.y: %1

                if( isnan(l2)){ printf("ERROR in TrussDynamics_d::run_LinSolve() iter=%i #points=%i  l2=%g ps_cor(%g,%g,%g) points(%g,%g,%g)\n", iter, i, l2, ps_cor[i].x, ps_cor[i].y, ps_cor[i].z, points[i].f.x, points[i].f.y, points[i].f.z); exit(0); }
                points[i].f = ps_cor[i];
            }            
            time += dt;
            //printf( "STEP: %6i time: %16.8f l2sum: %20.8e\n", iter, time, l2sum );
        } // for iter ( time stepping )
        //exit(0);
    }


    __attribute__((hot)) 
    void TrussDynamics_d::run_Cholesky_omp_simd(int niter) {
        //printf( "TrussDynamics_d::run_Cholesky_omp_simd() \n" );
        const int m=3;
        memcpy(ps_cor, points, nPoint * sizeof(Vec3d));
        double dt2 = dt * dt;
        double inv_dt = 1/dt;
        cleanForce();
        int iter=0;
        for (iter = 0; iter < niter; iter++) {
            #pragma omp simd
            for (int i=0;i<nPoint;i++){ ps_pred[i] = points[i].f + vel[i].f*dt + forces[i].f*dt2;  }
            #pragma omp simd
            for (int i=0;i<nPoint;i++){ linsolve_b[i] = rhs_ProjectiveDynamics_i(i,ps_pred); }
            Lsparse.fwd_subs_m( m,  (double*)linsolve_b,  (double*)linsolve_yy );            
            #pragma omp simd
            for(int i=0; i<nPoint; i++){ linsolve_yy[i].mul(1/LDLT_D[i]); } // Diagonal 
            LsparseT.fwd_subs_T_m( m,  (double*)linsolve_yy,  (double*)ps_cor );
            #pragma omp simd
            for (int i=0;i<nPoint;i++) {
                Vec3d    v   = (ps_cor[i] - points[i].f) * inv_dt;
                if(kFix){ if( kFix[i] > 1e-8){ v=Vec3dZero; } }
                vel[i].f    = v;
                points[i].f = ps_cor[i];
            }
        }
    }

    __attribute__((hot)) 
    void TrussDynamics_d::run_Cholesky_omp(int niter) {
        //printf( "TrussDynamics_d::run_Cholesky_omp() \n" );
        const int m=3;
        memcpy(ps_cor, points, nPoint * sizeof(Vec3d));
        double dt2 = dt * dt;
        double inv_dt = 1/dt;
        cleanForce();
        const int ns=3;
        double sum[ns];
        double* bb  = (double*)linsolve_b; 
        double* yy  = (double*)linsolve_yy;
        double* xx  = (double*)ps_cor ;
        int iter=0;
        #pragma omp parallel shared(sum) private(iter)
        for (iter = 0; iter < niter; iter++) {
            #pragma omp for
            for (int i=0;i<nPoint;i++){ ps_pred[i] = points[i].f + vel[i].f*dt + forces[i].f*dt2;  }
            #pragma omp for 
            for (int i=0;i<nPoint;i++){ linsolve_b[i] = rhs_ProjectiveDynamics_i(i,ps_pred); }
            #pragma omp single
            {
                Lsparse.fwd_subs_m( m,  (double*)linsolve_b,  (double*)linsolve_yy );
                //for(int i=0; i<nPoint; i++){ Lsparse.fwd_subs_mi( m,  (double*)linsolve_b,  (double*)linsolve_yy ); }
            }
            // for (int i=0; i<Lsparse.n; i++){
            //     for(int s=0;s<ns;s++){ sum[s]=0.0; }
            //     const int i0 = Lsparse.i0s [i];
            //     const int ni = Lsparse.nngs[i];
            //     #pragma omp for reduction(+:sum)
            //     for (int k=0; k<ni; k++){
            //         const int ik = i0+k;
            //         const int j  = Lsparse.inds[ik]*ns;
            //         const double l  = Lsparse.vals[ik];  
            //         for(int s=0;s<ns;s++){ sum[s]+=l*yy[j+s]; }
            //     }
            //     for(int s=0;s<ns;s++){ 
            //         const int ii=i*ns+s; yy[ii]=bb[ii]-sum[s]; 
            //     }
            //     //#pragma omp barrier
            // }
            
            #pragma omp for
            for(int i=0; i<nPoint; i++){ linsolve_yy[i].mul(1/LDLT_D[i]); } // Diagonal 
            #pragma omp single
            {
                LsparseT.fwd_subs_T_m( m,  (double*)linsolve_yy,  (double*)ps_cor );
                //for(int i=0; i<nPoint; i++){ LsparseT.fwd_subs_T_mi( m,  (double*)linsolve_yy,  (double*)ps_cor ); }
            }
            // for (int i=Lsparse.n-1; i>=0; i--){
            //     for(int s=0;s<ns;s++){ sum[s]=0.0; }
            //     const int i0 = LsparseT.i0s [i];
            //     const int ni = LsparseT.nngs[i];
            //     #pragma omp for reduction(+:sum)
            //     for (int k=0; k<ni; k++){
            //         const int ik = i0+k;
            //         const int j  = LsparseT.inds[ik]*ns;
            //         const double   l  = LsparseT.vals[ik];  
            //         for(int s=0;s<ns;s++){ sum[s]+=l*xx[j+s]; }
            //     }
            //     for(int s=0;s<ns;s++){ 
            //         const int ii=i*ns+s; xx[ii]=yy[ii]-sum[s]; 
            //     }
            // }

            #pragma omp for
            for (int i=0;i<nPoint;i++) {
                vel[i].f    = (ps_cor[i] - points[i].f) * inv_dt;
                points[i].f = ps_cor[i];
            }
        }
    }

    // =================== Linearized Elasticity Truss Simulation

    void TrussDynamics_d::prepareLinearizedTruss(){
        for(int ib=0; ib<nBonds; ib++){ 
            int2    b = bonds[ib];
            Vec3d   d = points[b.y].f-points[b.x].f;
            double  l = d.normalize();
            hbs [ib].f = (Vec3f)d;
            hbs [ib].e = l-bparams[ib].x; 
            dpos[ib]   = Quat4fZero; 
        }
    }

    void TrussDynamics_d::evalTrussForcesLinearized(){ 
        const double k = kGlobal;
        for(int ib=0;ib<nBonds;ib++){
            const int2   b  = bonds[ib];
            Quat4f h        = hbs[ib];
            float        dl = h.w + h.f.dot( dpos[b.y].f-dpos[b.x].f );
            h.f.mul(dl*k);
            fdpos[b.x].f.sub(h.f);
            fdpos[b.y].f.add(h.f);
        }
    }

    float TrussDynamics_d::evalTrussForceLinearized_neighs2(int iG){
        const Quat4f dp = dpos[iG];
        Quat4f f = Quat4f{0.0f,0.0f,0.0f,0.0f};
        const double k = kGlobal;
        for(int ij=0; ij<nNeighMax; ij++){
            const int j  = nNeighMax*iG + ij;
            const int2 b = neighBs[j];
            if(b.x == -1) break;
            //const Quat4d par = bparams[b.y];
            const Quat4f h  = hbs[b.y];
            float  dl = h.w + h.f.dot( dpos[b.x].f-dp.f );
            f.f.add_mul( h.f, k*dl );
        }
        fdpos[iG] = f; // we may need to do += in future
        return f.f.norm2();
    }
    float TrussDynamics_d::evalTrussForcesLinearized_neighs2(){
        float F2=0;
        for(int iG=0; iG<nPoint; iG++){ F2+=evalTrussForceLinearized_neighs2(iG); } 
        return F2;
    }

    void TrussDynamics_d::prepareLinearizedTruss_ling( double* bvec, bool bAddForce ){
        //printf( "prepareLinearizedTruss_ling nBonds %i \n", nBonds );
        Vec3d* f0s = (Vec3d*)bvec;
        //double k = kGlobal;
        double k = 50.0;
        for(int i=0; i<nPoint; i++){ f0s[i]=Vec3dZero; }
        for(int ib=0; ib<nBonds; ib++){ 
            int2    b = bonds[ib];
            Vec3d   d = points[b.y].f-points[b.x].f;
            double  l = d.normalize();
            kDirs[ib] = d; 
            Vec3d f = d*(l-bparams[ib].x)*k;
            f0s[b.x].add(f);
            f0s[b.y].sub(f);
            //printf( "prepare_lin[%i|%i,%i] f0=%g l=%g hdir(%g,%g,%g) \n", ib,b.x,b.y, d.dot(f), l, kDirs[ib].x, kDirs[ib].y, kDirs[ib].z );
            //dpos[ib]  = Quat4fZero; 
        }
        //printf("fdl: "); print_vector(nPoint, (double*)f0s, 3,0,3 );
        if(bAddForce){
            for(int iG=0; iG<nPoint; iG++){ f0s[iG].add( forces[iG].f ); }
        }
        //printf("f:  "); print_vector(nPoint, (double*)f0s, 3,0,3 );
    }

    void TrussDynamics_d::dot_Linearized_bonds(int n, double* x, double * Ax){
        //printf( "TrussDynamics_d_d::dot_Linearized_neighs2(n=%i) @x=%li @Ax=%li\n", n, (long*)x, (long*)Ax );
        Vec3d* dpos  = (Vec3d*)x;
        Vec3d* fdpos = (Vec3d*)Ax;
        int nG = n/3;
        for(int iG=0; iG<nG; iG++){
            double kp = kFix[iG] + kLinRegularize;
            //kp = 0;
            //printf( "Kdp[%i] k=%g dp(%g,%g,%g)\n", iG, kp, dpos[iG].x, dpos[iG].y, dpos[iG].z );
            fdpos[iG] = dpos[iG]*kp;
        }
        //printf("fdpos:"); print_vector(nG, (double*)fdpos, 3,0,3 );
        //double k = kGlobal;
        double k = 50.0;
        for(int ib=0; ib<nBonds; ib++){
            const int2  b = bonds[ib];
            const Vec3d h = kDirs[ib];
            // double k = bparams.x;
            double dl = h.dot( dpos[b.y]-dpos[b.x] );
            //h.mul( -dl*k );
            Vec3d f = h*(-dl*k);
            //printf( "dot_bond[%i|%i,%i] k=%g dl=%g h(%g,%g,%g) f(%g,%g,%g)\n", ib,b.x,b.y, k, dl,  h.x,h.y,h.z,   f.x,f.y,f.z );
            fdpos[b.x].add( f );
            fdpos[b.y].sub( f );
        }
    }

    void TrussDynamics_d::dot_Linearized_neighs2(int n, double* x, double * Ax){
        //printf( "TrussDynamics_d_d::dot_Linearized_neighs2(n=%i) @x=%li @Ax=%li\n", n, (long*)x, (long*)Ax );
        Vec3d * dpos  = (Vec3d*)x;
        Vec3d * fdpos = (Vec3d*)Ax;
        int nG = n/3;
        for(int iG=0; iG<nG; iG++){
            const Vec3d dp = dpos[iG];
            //double kOnSide = kLinRegularize;
            double kp = kFix[iG] + kLinRegularize;
            //kp = 0;
            //kp.add(kLinRegularize);
            Vec3d f = dp*kp;   // regularize
            //printf( "dot_Linearized_neighs2[%i] dp(%g,%g,%g) f(%g,%g,%g) kp=%g \n", iG, dp.x, dp.y, dp.z, f.x, f.y, f.z, kp );
            //double k = kGlobal;
            double k = 50.0;
            for(int ij=0; ij<nNeighMax; ij++){
                const int j  = nNeighMax*iG + ij;
                const int2 b = neighBs[j];
                if(b.x == -1) break;
                //const Quat4d par = bparams[b.y];
                const Vec3d h  = kDirs[b.y];
                double dl = h.dot( dpos[b.x]-dp );
                f.add_mul( h, -dl*k );
            }
            fdpos[iG] = f; // we may need to do += in future
            //fdpos[iG] = Vec3d{iG,1,2};
        }
    }

    void TrussDynamics_d::findGuassSeidelPivotingPriority( double* priority ){
        // The lightest points should be solved first
        // also the points with the most neighbors should be solved first
        // also the points which are lighter than their neighbors should be solved first

        // Perhaps better strategy would be to solve first the points which move the most in previous iteration

        for(int iG=0; iG<nPoint; iG++){
            double m = points[iG].w;
            double neigh_mass = 0;
            for(int ij=0; ij<nNeighMax; ij++){
                const int j  = nNeighMax*iG + ij;
                const int2 b = neighBs[j];
                if(b.x == -1) break;
                neigh_mass += points[b.x].w;
            }
            priority[iG] = (neigh_mass + 1e-3)/m;
        }
    }

    void TrussDynamics_d::move_dpos( double dt ){
        for(int iG=0; iG<nPoint; iG++){ 
            Quat4d p = points[iG];
            Quat4d v = vel   [iG];
            const Quat4d f = forces[iG];
            points_bak[iG] = p;
            v.f.add_mul( f.f, dt);
            //Quat4d dp; dp.f.set_mul(v.f,dt); dp.w=p.w;
            //p.f.add( dp.f );
            p.f.add_mul( v.f, dt);
            points[iG] = p;
            vel   [iG] = v;
            //dpos  [iG] = (Quat4f)dp;
        }
    }

    void TrussDynamics_d::apply_dpos( double sc ){
        for(int iG=0; iG<nPoint; iG++){ 
            points[iG].f.add_mul( (Vec3d)dpos[iG].f, sc );
        }
    }

    void TrussDynamics_d::update_velocity( double dt ){
        double invdt = 1.0/dt;
        for(int iG=0; iG<nPoint; iG++){ 
            //Vec3d dp = (Vec3d)dpos[iG].f;
            //points[iG].f.add( dp );
            // warrning  points are already updated by p_k1 =  p_k + ( v_k + f_k *dt )*dt
            // but to update v we need    v_k1 = ( p_k1 - p_k )/dt, therefore we need to store p_k1, or add ( v_k + f_k *dt )*dt to dp
            //   so dp = ( v_k + f_k *dt )*dt + dp_constrain
            Vec3d dp = points[iG].f - points_bak[iG].f;
            vel[iG].f.set_mul( dp, invdt );
        }
    }

    double TrussDynamics_d::constr_jacobi_neighs2_absolute(){
        // This is constrain-solver
        //   we know how long should be each bond
        //   we need to solve for dpos (i.e. how much we should shift each point to make the bond correct length)
        //   we shift more the points which are lighter (i.e. we use inverse mass as a weight)
        //   it should correspond to to position based dynamics (PBD) developed by Muller et al. 2007 (https://matthias-research.github.io/pages/) https://doi.org/10.1016/j.jvcir.2007.01.005  https://matthias-research.github.io/pages/publications/posBasedDyn.pdf
        double dlmax = 0;
        for(int iG=0; iG<nPoint; iG++){
            const Quat4d pi = points[iG];
            //const double k = 50.0;
            //double im = 1.0/pi.w;
            Vec3d dp = Vec3dZero;
            double wsum = 0;
            for(int ij=0; ij<nNeighMax; ij++){
                const int j  = nNeighMax*iG + ij;
                const int2 b = neighBs[j];
                if(b.x == -1) break;
                Quat4d pj  = points[b.x];
                Vec3d d; d.set_sub( pj.f, pi.f);
                double l  = d.normalize();
                double dl = l - bparams[b.y].x; // l0
                //double w = im/(im+1./pj.w);
                double w = pj.w/(pi.w+pj.w);  // (1/mi)/(1/mi+1/mj) = mi*mj/( mi*(mi + mj)) = mj/(mi+mj)
                dp.add_mul( d, dl*w*w );
                wsum += w;
                {
                    dl    = (dl>0    ) ? dl : -dl;
                    dlmax = (dl>dlmax) ? dl : dlmax;
                }
            }
            // Note about weighting: Sum_i{ d*w_i^2 }/Sum_i{ w_i } 
            // if we have just one neighbor then d*w_i^2/w_i = d*w_i
            // if we have two same neighbors then 2*(d*w_i^2)/2*w_i = d*w_i
            dp.mul( 1.0/wsum );
            //points[iG].f.add( dp ); // This would be gauss-seidel, but we need to do it in parallel, and we do-not want have dependency on order of points (pivottinh)
            //dpos[iG].f.add( (Vec3f)dp ); // we may need to do += in future
            dpos[iG].f = (Vec3f)dp; // we may need to do += in future
            //fdpos[iG] = Vec3d{iG,1,2};
        }
        return dlmax;
    }

    double TrussDynamics_d::constr_jacobi_neighs2_diff(){
        // This is constrain-solver
        //   we know how much longer/shorter should be each bond (dls)
        //   we need to solve for dpos (i.e. how much we should shift each point to make the bond correct length)
        double dlmax = 0;
        for(int iG=0; iG<nPoint; iG++){
            const Quat4f dpi = dpos[iG];
            Vec3f dp         = Vec3fZero;
            double wsum      = 0;
            for(int ij=0; ij<nNeighMax; ij++){
                const int j  = nNeighMax*iG + ij;
                const int2 b = neighBs[j];
                if(b.x == -1) break;
                const Vec3f  h  = (Vec3f)kDirs[b.y];
                const Quat4f dpj = dpos[b.x];
                float dl = h.dot( dpj.f - dpi.f ) + dls[ b.y ];
                float mj = dpj.x;
                float w  = dpi.w/(dpi.w+dpj.w);  // (1/mi)/(1/mi+1/mj) = mi*mj/( mi*(mi + mj)) = mj/(mi+mj)
                dp.add_mul( h, dl*w*w );
                wsum += w;
                {
                    dl    = (dl>0    ) ? dl : -dl;
                    dlmax = (dl>dlmax) ? dl : dlmax;
                }
            }
            dp.mul( 1.0/wsum );
            dpos[iG].f = (Vec3f)dp;
            //dpos[iG].f.add( dp ); // we may need to do += in future
            //fdpos[iG] = Vec3d{iG,1,2};
        }
        return dlmax;
    }


    void TrussDynamics_d::run_constr_dynamics(int nitr, double dt ){
        int nsolver = 3;
        double dlconv = 1e-3;
        //dt = 0.001;
        dt = 1.0;
        for(int i=0; i<nitr; i++){
            move_dpos( dt );
            for(int isolver=0; isolver<nsolver; isolver++){
                double dlmax = constr_jacobi_neighs2_absolute();
                //printf( "run_constr_dynamics[%i,%i] dlmax=%g \n", i,isolver,  dlmax );
                apply_dpos( 1.5 );
                //constr_jacobi_neighs2_diff();
            }
            update_velocity( dt );
        }
    }


    void TrussDynamics_d::solveLinearizedConjugateGradient(){
        // modified from    genLinSolve_CG()      /cpp/common/math/Lingebra.h
        // We can just re-prase conjugate gradient to dynamics by putting
        // x_k     = pos_k
        // p_k     = vel_k
        // r_k     = force_k
        // alpha_k = dt_k 


        // the idea is that we chose alpha_k to minimize the residual force and make it orthogonal to the previous residual force
        // dot( r_k1, r_k ) = 0
        // dot( r_k1, r_k ) = dot( r_k - alpha_k * A*p_k, r_k ) = dot( r_k, r_k ) - alpha_k * dot( A*p_k, r_k ) = 0
        // alpha_k = dot( r_k, r_k ) / dot( A*p_k, r_k )
        // Alternatively
        // dot( r_k1, r_k ) = dot( b - A*(x_k+d_k) , b - A*x_k ) = dot(b,b) - dot(b,A*d_k) - 2*dot(b,A*x_k) + dot(A*d_k,A*x_k)  = bb + dot(A*d_k, A*x_k - b  ) - 2*dot(b,A*x_k) = 0
        // dot( A*d_k, A*x_k - b  )  = 2*dot(b,A*x_k) - bb
        // dot( A*d_k,   r_k      )  = 2*dot(b,A*x_k) - bb
        // d_k = ak * p_k        
        // ak  = dot(r_k,r_k) / dot( r_k , A * p_k )
        // 1   = dot(r_k,r_k) / dot( r_k , A * d_k )

        // https://en.wikipedia.org/wiki/Conjugate_gradient_method
        // l = dot(r_k,r_k) / (p_k * Ap_k)
        // x_k1 = x_k + l_k * p_k
        // r_k1 = b - A*x_k1               = r_k          - l_k * Ap_k     (in exact arithmetic it is equivalent, because r_k = b - A*x_k )
        //      = b - A*(x_k + l_k * p_k)  =  (b - A*x_k) - l_k * Ap_k    
        // fac = dot(r_k1,r_k1) / dot(r_k,r_k)
        // p_k1 = r_k1 + fac * p_k                  // it is like velocity damping assuming that r_k1 is the force and p_k is the velocity      


        //   l = dot(r_k1,r_k1) / (p_k1 * Ap_k1)

          
        // The problem:
        //  * to calculate l_k = dot(r_k,r_k) / (p_k * Ap_k)     we need to calculate  A * p_k
        //  * however to calculate   r_k1 = b   -       A*x_k1   we neeed to calculate A * x_k1
        //  * we may also calculater r_k1 = r_k - l_k * A*p_k    but this may be numerically unstable 
        //  *  if there is some way how to calculate ak using A*x_k without calculating A*p_k then we can avoid this problem



    };
    

    // =================== Truss Simulation

    void TrussDynamics_d::updateInveriants(bool bPrint){
        mass=0;
        cog =Vec3dZero;
        vcog=Vec3dZero;
        for(int i=0; i<nPoint; i++){ double mi=points[i].w; mass+=mi;  cog.add_mul( points[i].f, mi ); vcog.add_mul( vel[i].f, mi ); }
        cog .mul( 1.0/mass );
        vcog.mul( 1.0/mass );
        I=Mat3dZero;
        L=Vec3dZero;
        torq=Vec3dZero;
        for(int i=0; i<nPoint; i++){ 
            double mi=points[i].w; 
            Vec3d d; 
            d   .set_sub   ( points[i].f, cog    );
            L   .add_crossw( d, vel[i]   .f, mi  );
            torq.add_cross ( d, forces[i].f      );
            I   .add_outer ( d, d, mi            );
        }
        if(bPrint){
            printf( "TrussDynamics_d::updateInveriants() mass %g cog(%g,%g,%g) vcog(%g,%g,%g)  L(%g,%g,%g) torq(%g,%g,%g) \n", mass, cog.x,cog.y,cog.z, vcog.x,vcog.y,vcog.z, L.x,L.y,L.z, torq.x,torq.y,torq.z );
            printf( "TrussDynamics_d::updateInveriants() I \n" ); I.print();
        }
    }

    void TrussDynamics_d::evalTrussForce_neighs(int iG){
        Quat4d p = points[iG];
        Quat4d f =Quat4d{0.0f,0.0f,0.0f,0.0f};
        //printf( "--- p[%i] \n", iG );
        //#pragma omp simd
        for(int ij=0; ij<nNeighMax; ij++){
            int j  = nNeighMax*iG + ij;
            int ja = neighs[j];
            if(ja == -1) break;
            //f.add( springForce( points[ja].f - p.f, params[j] ) );
            
            Vec3d d =  points[ja].f - p.f;
            double li = d.norm();
            /*
            double fi,ei = springForce( li, fi, params[j] );
            //f.add( Quat4d{ d*(fi/l), ei } );
            f.f.add_mul( d, fi/li );
            */
            double k = kGlobal;
            f.f.add_mul( d, (k*(li-params[j].x)/li) );

            //printf( "p[%i,ij=%i,j=%i] li=%7.3f dl=%8.5e fi=%8.5e e=%8.5e par(%7.3f,%8.5e,%8.5e,%8.5e) \n", iG,ij,ja, li, li-params[j].x, fi,ei, params[j].x,params[j].y,params[j].z,params[j].w );
        }
        forces[iG] = f; // we may need to do += in future
    }

    void TrussDynamics_d::evalTrussForce_neighs2(int iG){
        //const double Adamp = collision_damping*0.5/dt;
        //const int iG = get_global_id(0);
        const Quat4d p = points[iG];
        const Quat4d v = vel   [iG];
        Quat4d f = Quat4d{0.0f,0.0f,0.0f,0.0f};
        //printf( "--- p[%i] \n", iG );
        //#pragma omp simd
        for(int ij=0; ij<nNeighMax; ij++){
            const int j  = nNeighMax*iG + ij;
            const int2 b = neighBs[j];
            if(b.x == -1) break;
            const Quat4d par = bparams[b.y];
            //f.add( springForce( points[ja].f - p.f, params[j] ) );
            Quat4d d = points[b.x];
            d.f.sub( p.f );
            const double l  = d.f.norm();
            //double k        = kGlobal;    // This works
            const double k    = par.z;    // This crashes
            const double fl   = k*(l-par.x);
            f.f.add_mul( d.f, fl/l );

            // const double dv  = d.f.dot( vel[b.x].f - v.f );
            // if(dv<0){ fl*=1-dv; }

            // // collision damping
            // const double invL = 1/l;
            // const double dv  = d.f.dot( vel[b.x].f - v.f );
            // double imp = Adamp * p.w*d.w*dv/(p.w+d.w);
            // //double imp = 0.1* 0.5 * p.w*d.w*dv/(p.w+d.w);
            // //const double imp = 0;
            // //imp/=dt;
            //f.f.add_mul( d.f, ( imp + fl )*invL );
            //printf( "p[%i,ij=%i,j=%i] fl: %10.2e dl: %10.2e l: %10.2e l0: %10.2e k: %10.2e \n", iG,ij,j, fl, l-par.x, l, par.x, k );
        }
        forces[iG] = f; // we may need to do += in future
    }

    void TrussDynamics_d::evalTrussForces_neighs2(){
        for(int iG=0; iG<nPoint; iG++){ evalTrussForce_neighs2(iG); } 
    }

    void TrussDynamics_d::evalTrussForces_neighs(){
        for(int iG=0; iG<nPoint; iG++){ evalTrussForce_neighs(iG); } 
    }

    

    void TrussDynamics_d::evalTrussForces_bonds(){
        for(int i=0; i<nBonds; i++){
            int2  b = bonds[i];
            const Quat4d& pi = points[b.x];
            const Quat4d& pj = points[b.y];
            Vec3d d = pj.f - pi.f;
            double l = d.norm();
            //double fi,ei = springForce( li, fi, bparams[i] );
            double k = kGlobal;
            double f = k*(l-bparams[i].x);

            // Limit acceleration force to improve stability   -- it does not seem to help
            // double mmin = (pi.w < pj.w) ? pi.w : pj.w;
            // if( (fabs(f)/mmin) > maxAcc ){   // here it would be more efficient to use momentum rather than force 
            //     f = maxAcc*mmin;
            // }

            // collision damping
            //double vi   = d.dot( vel[b.x].f );
            //double vj   = d.dot( vel[b.y].f );
            //double dv   = vj   - vi;
            double invL = 1./l;
            double dv   = d.dot( vel[b.y].f - vel[b.x].f )*invL;
            double mcog = pj.w + pi.w;
            double imp  = collision_damping * pi.w*pi.w*dv/mcog;

            d.mul( (imp + f)*invL );          
            forces[b.x].f.add(d);
            forces[b.y].f.sub(d);
        } 
    }

    void TrussDynamics_d::evalEdgeVerts(){
        //for(int i=0; i<nEdgeVertBonds; i++){ evalEdgeVert( edgeVertBonds[i].verts, edgeVertBonds[i].c, edgeVertBonds[i].K ); }
        for(int i=0; i<nEdgeVertBonds; i++){ evalEdgeVert( i ); }
    }

    void TrussDynamics_d::evalTrussCollisionImpulses_bonds( double rate ){
        // Collision forces are based on momentum-conserving impulses, it does not need to know anything about the interaction potential (e.g. spring constants)
        for(int i=0; i<nBonds; i++){
            int2  b = bonds[i];
            const Quat4d& pi = points[b.x];
            const Quat4d& pj = points[b.y];
            Vec3d d  = pi.f - pj.f;
            double l  = d.normalize();
            double vi = d.dot( vel[b.x].f );
            double vj = d.dot( vel[b.y].f );
            double mcog = pj.w + pi.w;
            double dv   = vj   - vi;
            //double vcog = (pi.w*vi + pj.w*vj)/(pi.w+pj.w);
            // ---- Inelastic collision
            //double dvi  = vcog - vi;
            //double imp1 = dvi*pi.w;
            // Analytical Derivation:
            // dvi = ( pi.w*vi + pj.w*vj                        )/(pi.w+pj.w) - (pi.w+pj.w)*vi/(pi.w+pj.w)
            //     = ( pi.w*vi + pj.w*vj -  pi.w*vi - pj.w*vi   )/(pi.w+pj.w)
            //     = (           pj.w*vj            - pj.w*vi   )/(pi.w+pj.w)
            //     = (           pj.w*(vj-vi)                   )/(pi.w+pj.w)
            // imp1 = dvi*pi.w = pi.w*pj.w*(vj-vi)/(pi.w+pj.w)
            double imp = pi.w*pi.w*dv/mcog;
            //if( i==146 ){
            //    double dvj  = vcog - vj;
            //    double imp2 = dvj*pj.w; // shoould be the same as imp1 
            //    printf( "evalTrussCollisionForces_bonds[%i] imp1 %g imp2 %g \n", i, imp1, imp2  );
            //}
            // apply force to points
            d.mul( imp*rate );
            forces[b.x].f.add( d );
            forces[b.y].f.sub( d );
        } 
    }

    double TrussDynamics_d::evalBondTension(){
        double E = 0.0;
        for(int i=0;i<nBonds; i++ ){
            int2   b  = bonds[i];
            double l0 = bparams[i].x;
            //double l0 = l0s[i];
            double l   = (points[b.y].f-points[b.x].f).norm();
            double dl  = l - l0;
            double s   = dl/l0;
            //if( fabs(s)>0.5 ){ printf( "evalBondTension[%i] strain=%g l=%g l0=%g\n", i, s, l, l0 ); }
            //if(i==6272){ printf( "evalBondTension[%i](%i,%i) strain=%g l=%g l0=%g\n", i, b.x,b.y, s, l, l0 ); }
            strain[i] = s;
            E+= 0.5*dl*dl*bparams[i].z;
            // ToDo: break the bond if strain > maxStrain;
        }
        //exit(0);
        return E;
    }

    void TrussDynamics_d::applyForceRotatingFrame( Vec3d p0, Vec3d ax, double omega ){
        double omega2 = omega*omega;
        Vec3d omega_ax = ax*omega*2.0;
        for(int i=0;i<nPoint; i++ ){
            const Quat4d& p = points[i];
            const Quat4d& v = vel   [i];
            //Vec3d f; f.set_cross(ax,p.f-p0);
            Vec3d d,f;
            d.set_sub(p.f,p0);
            d.makeOrthoU(ax);
            f.set_mul( d, omega2 );     // centrifugal force  = r*m*omega^2
            f.add_cross(omega_ax,v.f);  // Coriolis force     = 2*m*omega*v
            forces[i].f.add_mul(f, p.w );
            //forces[i].f.add_mul( f, p.w*omega2 );
        }
    }

    void TrussDynamics_d::applyForceCentrifug( Vec3d p0, Vec3d ax, double omega ){
        for(int i=0;i<nPoint; i++ ){ applyForceCentrifug_i( i, p0, ax,omega ); }
    }

    void TrussDynamics_d::printNeighs(int i){
        int j0 = i*nNeighMax;
        for(int jj=0;jj<nNeighMax;jj++){
            int j=j0+jj;
            int ing = neighs[j];
            if(ing<0) break;
            Quat4d par = params[j];
            printf( "ng[%i,%i|%i] l0,kP,kT,damp(%g,%g,%g,%g)\n", i, ing, jj, par.x,par.y,par.z,par.w );
        }
    }
    void TrussDynamics_d::printAllNeighs(){ printf("TrussDynamics_d::printAllNeighs(nPoint=%i,nNeighMax=%i)\n",nPoint,nNeighMax); for(int i=0;i<nPoint;i++){ printNeighs(i); }; };

    double TrussDynamics_d::getFmax(){ 
        double fmax=0;
        for(int i=0; i<nPoint; i++){ double f=forces[i].norm(); fmax=fmax>f?fmax:f; }   
        //printf( "|fmax|=%g\n", fmax );
        return fmax;
    }

    void TrussDynamics_d::setFixPoints( int n, int* fixPoints, double Kfix, bool bRealloc ){
        if(bRealloc) reallocFixed();
        for(int i=0;i<n;i++){
            int ip = fixPoints[i]; 
            if(ip>nPoint){ printf("ERROR in TrussDynamics_d::setFixPoints fixing point %i > sim.nPoint(%i) \n", ip, nPoint); exit(0); }
            printf("TrussDynamics_d_d::setFixPoints()[%3i]: %6i \n", i, ip ); 
            kFix[ip] = Kfix ; 
        }
    }

    void TrussDynamics_d::addAngularVelocity( Vec3d p0, Vec3d ax ){
        for(int i=0; i<nPoint; i++){
            Vec3d dp = points[i].f - p0;
            Vec3d v; 
            v.set_cross(ax,dp);
            //v.set_cross(dp,ax);
            //cross(axis, p)
            vel[i].f.add(v);
        }
    };

    void TrussDynamics_d::addAngularVelocity2( Vec3d p0, Vec3d ax ){
        for(int i=0; i<nPoint; i++){
            Vec3d dp = points[i].f - p0;
            Vec3d v; 
            v.set_cross(ax,dp);
            v.mul( dp.norm() );
            vel[i].f.add(v);
        }
    };

    double TrussDynamics_d::move_GD(double dt){
        double ff=0;
        for(int i=0;i<nPoint; i++ ){
            Quat4d p = points[i];
            Quat4d f = forces[i];
            ff += f.f.norm2();
            p.f.add_mul( f.f, dt/p.w );
            //printf( "move_GD[%i] |d|=%g |f|=%g dt/m=%g m=%g \n", i, f.f.norm() * dt/p.w, f.f.norm(), dt/p.w, p.w );
            points[i]=p;
        }
        return ff;
    }

    inline double TrussDynamics_d::move_i_MD(int i, double dt, double cdamp ){
        Quat4d f = forces[i];
        Quat4d p = points[i];
        Quat4d v = vel   [i];

        //double vf = v.f.dot(f.f);
        //if( vf>0.0 ){ f.f.mul( 0.99 ); }

        v.f.mul( cdamp );
        v.f.add_mul( f.f, dt/p.w );
        p.f.add_mul( v.f, dt     );
        //printf( "move_GD[%i] |d|=%g |f|=%g dt/m=%g m=%g \n", i, f.f.norm() * dt/p.w, f.f.norm(), dt/p.w, p.w );
        vel   [i]=v;
        points[i]=p;
        return f.f.norm2();
    }

    double TrussDynamics_d::move_MD(double dt, double damp ){
        double cdamp = 1.0f - damp;
        double ff = 0.0; 
        for(int i=0;i<nPoint; i++ ){ move_i_MD(i, dt, cdamp ); }
        return ff;
    }

int TrussDynamics_d::run( int niter, double dt, double damp  ){
    double f2 = -1;
    niter=2;
    for(int itr=0; itr<niter; itr++){
        cleanForce();   
        //evalTrussForces_neighs();
        evalTrussForces_neighs2();
        //evalTrussForces_bonds();
        //applyCentrifugalForce( {0.,0.,0.}, {0.0,0.0,1.0}, 1e-2 );
        applyForceRotatingFrame( rot0.f, omega.f, omega.w );
        //move_GD( 0.00001 );
        //move_MD( 1e-3, 1e-5 );
        //move_GD( 1e-7 );
        f2 = move_MD( dt, damp );
        printf( "TrussDynamics_d_f::run[%i] |F|=%g\n", itr, sqrt(f2) );
    }
    exit(0);
    return niter;
}

inline void TrussDynamics_d::setOpt( double dt_, double damp_ ){
    dt      = dt_max   = dt_;  dt_min=0.1*dt_max;
    damping = damp_max = damp_;
    cleanForce( );
    cleanVel  ( );
}

void TrussDynamics_d::FIRE_update( double& vf, double& vv, double& ff, double& cv, double& cf ){
    double cs = vf/sqrt(vv*ff);
    //if( vf < 0.0 ){
    if( vv>1600.0 )damping  = damp_max;
    if( cs<0.02 ){
		dt       = fmax( dt * fdec, dt_min );
	  	damping  = damp_max;
		lastNeg  = 0;
        cv=0.0; cf=0.0;
        //printf( "FIRE<0 cv,cf(%g,%g) cs=%g   vf,vv,ff(%g,%g,%g) \n", cv,cf, cs, vf,vv,ff );
	}else{
		cf   =     damping * sqrt(vv/(ff+ff_safety));
		cv   = 1 - damping;
		if( lastNeg > minLastNeg ){
			dt      = fmin( dt * finc, dt_max );
			damping = damping  * falpha;
		}
		lastNeg++;
        
	}
    //printf( "FIRE>0 cv,cf(%g,%g) cs=%g dt=%g damp=%g/%g lastNeg=%i vf,vv,ff(%g,%g,%g) \n", cv,cf, cs,  dt, damping,damp_max, lastNeg,  vf,vv,ff );
}


int TrussDynamics_d::run_omp( int niter_max, bool bDynamic, double dt_, double damp_ ){
    //printf( "run_omp() niter_max %i dt %g Fconv %g Flim %g timeLimit %g outE %li outF %li \n", niter_max, dt, Fconv, Flim, timeLimit, (long)outE, (long)outF );
    double cdamp = 1.0f - damp_;
    double E,F2,ff,vv,vf;
    //long T0 = getCPUticks();
    int itr=0,niter=niter_max;
    #pragma omp parallel shared(niter,itr,cdamp)
    while(itr<niter){
        if(itr<niter){
        //#pragma omp barrier
        #pragma omp single
        {E=0;ff=0;vv=0;vf=0;}
        //------ eval forces
        //#pragma omp barrier
        //#pragma omp for reduction(+:E)
        #pragma omp for 
        for(int iG=0; iG<nPoint; iG++){ 
            forces[iG] = Quat4dZero;
            //evalTrussForce_neighs(iG);
            evalTrussForce_neighs2(iG);
            
            //applyForceCentrifug_i    ( iG, rot0.f, omega.f, omega.w );
            //if(bDynamic){ applyForceRotatingFrame_i( iG, rot0.f, omega.f, omega.w ); }
            //else        { applyForceCentrifug_i    ( iG, rot0.f, omega.f, omega.w ); }

            // ToDo: Through-space interactions (e.g. magnetic damping of ropes by induced currents)   
            //  1) Magneto-static F = (mu_0/2pi)(L/r)dot(I_i,I_j)    // see http://hyperphysics.phy-astr.gsu.edu/hbase/magnetic/wirfor.html
            //  2) Volatage induced in wire moving in magnetic field: https://en.wikipedia.org/wiki/Moving_magnet_and_conductor_problem#Conductor_frame
            //  3) Magnetic Breaking: 

        }
        // ---- assemble (we need to wait when all atoms are evaluated)
        //#pragma omp barrier
        /*
        if(!bDynamic){    // FIRE if not dynamic
            #pragma omp for reduction(+:vv,vf,ff)
            for(int i=0;i<nPoint; i++ ){
                Quat4d p = points[i];
                Quat4d f = forces[i];
                Quat4d v = vel   [i];
                vv += v.f.norm2();
                vf += v.f.dot(f.f);
                ff += f.f.norm2();
            }
            #pragma omp single
            { 
                FIRE_update( vf, vv, ff, cv, cf );
                //printf( "FIRE cv,cf(%g,%g)   vf,vv,ff(%g,%g,%g) \n", cv,cf, vf,vv,ff );
            }
        }
        */
        #pragma omp single
        { 
            if(user_update) user_update(dt);  // call e.g. spacecraft control (move the wheels etc. )
            evalEdgeVerts();
        }
        #pragma omp for reduction(+:ff)
        for(int i=0;i<nPoint; i++ ){
            //move_i_MD( i, dt, 1.0 );
            ff += move_i_MD( i, dt, cdamp );
            /*
            if(bDynamic){ 
                move_i_MD( i, dt, cdamp );
            }else{
                vel[i].f = vel[i].f*cv  + forces[i].f*cf;  // FIRE
                move_i_MD( i, dt, 1.0 );
            }
            */
        }
        //#pragma omp barrier
        #pragma omp single
        { 
            F_residual = sqrt(ff);
            //printf( "TrussDynamics_d::run_omp() itr %i/%i E=%g |F|=%g time= %g [ms]( %g [us/%i iter]) \n", itr,niter_max, E, F_residual, time*1e+3, time*1e+6/itr, itr );
            time+=dt;
            itr++; 
        }
        } // if(itr<niter){
    }
    //{
    //double t = (getCPUticks() - T0)*tick2second;
    //if(itr>=niter_max)if(verbosity>0)printf( "run_omp() NOT CONVERGED in %i/%i E=%g |F|=%g time= %g [ms]( %g [us/%i iter]) \n", itr,niter_max, E, sqrt(F2), t*1e+3, t*1e+6/itr, itr );
    //}
    return itr;
}

//void TrussDynamics_d::print_points(){}

#endif
