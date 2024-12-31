
#ifndef  OrbSim_f_h
#define  OrbSim_f_h

#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
//#include "datatypes.h"  

#include "Buckets.h"

#include "raytrace.h"
#include "geom3D.h"
#include "Interfaces.h"

#include "Cholesky.h"
#include "SparseMatrix.h"
#include "SparseMatrix2.h"

float checkDist(int n, const Vec3f* vec, const Vec3f* ref, int verb=1, float tol=1e-12 ){
    float err=0;
    for(int i=0; i<n; i++){ 
        float ei = (vec[i]-ref[i]).norm();
        err = fmax( err, ei );
        if(verb>1){
           if(ei>tol) printf( "[%i] vec(%20.10f,%20.10f,%20.10f) ref(%20.10f,%20.10f,%20.10f)\n", i, vec[i].x,vec[i].y,vec[i].z,  ref[i].x,ref[i].y,ref[i].z ); 
        }
    }; 
    if( verb>0 )printf( "Err| v - ref |=%g\n", err );
    return err;
}

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
struct EdgeVertBond_f{ 
    Vec3i verts; 
    float c; 
    float K; 
    Vec3f f=Vec3fZero; 
};

void fitAABB( Quat8f& bb, int n, int* c2o, Quat4f * ps ){
    //Quat8f bb;
    //bb.lo = bb.lo = ps[c2o[0]];
    for(int i=0; i<n; i++){ 
        //printf( "fitAABB() i %i \n", i );
        int ip = c2o[i];
        //printf( "fitAABB() i=%i ip=%i \n", i, ip );
        Quat4f p = ps[ip];
        bb.lo.f.setIfLower  ( p.f );
        bb.hi.f.setIfGreater( p.f );
    }; 
    //return bb;
}

void fitAABB_edge( Quat8f& bb, int n, int* c2o, int2* edges, Quat4f * ps ){
    //Quat8f bb;
    //bb.lo = bb.lo = ps[c2o[0]];
    for(int i=0; i<n; i++){ 
        //printf( "fitAABB() i %i \n", i );
        int ie = c2o[i];
        int2 b = edges[ie];
        //printf( "fitAABB() i=%i ip=%i \n", i, ip );
        Quat4f pi = ps[b.x];
        bb.lo.f.setIfLower  ( pi.f );
        bb.hi.f.setIfGreater( pi.f );
        Quat4f pj = ps[b.y];
        bb.lo.f.setIfLower  ( pj.f );
        bb.hi.f.setIfGreater( pj.f );
    }; 
    //return bb;
}

inline void updatePointBBs(const Buckets& buckets, Quat8f* BBs, Quat4f* points, bool bInit=true){
    //printf( "updatePointBBs() START \n" );
    for(int ib=0; ib<buckets.ncell; ib++){
        //printf( "updatePointBBs() ib %i \n", ib );
        if(bInit){ BBs[ib].lo.f = Vec3fMax; BBs[ib].hi.f = Vec3fMin; }
        int n = buckets.cellNs[ib];
        if(n>0){
            int i0 = buckets.cellI0s[ib];
            //printf( "updatePointBBs() ib %i n %i i0 %i \n", ib, n, i0 );
            fitAABB( BBs[ib], n, buckets.cell2obj+i0, points );
        }
    }
    //printf( "updatePointBBs() DONE \n" );
}

inline void updateEdgeBBs(const Buckets& buckets, Quat8f* BBs, int2* edges, Quat4f* points, bool bInit=true){
    //printf( "updateEdgeBBs() START \n" );
    for(int ib=0; ib<buckets.ncell; ib++){
        if(bInit){ BBs[ib].lo.f = Vec3fMax; BBs[ib].hi.f = Vec3fMin; }
        int n = buckets.cellNs[ib];
        if(n>0){
            int i0 = buckets.cellI0s[ib];
            fitAABB_edge( BBs[ib], n, buckets.cell2obj+i0, edges, points );
        }
    }
    //printf( "updateEdgeBBs() DONE \n" );
}

class OrbSim_f : public Picker { public:
    float time=0;
    int nPoint=0, nNeighMax=0, nNeighTot=0;
    // cpu buffers
    Quat4f* points=0;  // position and mass
    Quat4f* forces=0;  // force and energy
    Quat4f* vel   =0;  // velocity
    Quat4f* impuls=0;  // accumulated impulse from corrector

    float* kFix=0;   // force constant for fixed points
    float  kLinRegularize = 1.0;

    Quat4f* params=0;  // neighbor parameters ( l0,kP,kT,damp )
    float * kngs  =0;  // neighbor stiffness  ( for linear solver )
    int*    neighs=0;  // neighbor indices
    int2*   neighBs=0; // neighbor bond indices
    //int*    neighBs=0; // neighbor bond indices
    int*    neighB2s=0;  // neighbor indices


    // Cholesky / Projective Dynamics 
    //float* LDLT_L=0;
    //float* LDLT_D=0; 
    //int* neighsLDLT=0;

    // Cholesky / Projective Dynamics 
    float* PDmat=0;
    float* LDLT_L=0;
    float* LDLT_D=0; 
    int* neighsLDLT=0;
    int  nNeighMaxLDLT=0;
    SparseMatrix2<float>  Lsparse; 
    SparseMatrix2<float>  LsparseT; 
    SparseMatrix<float> PDsparse;
    //CGsolver cgSolver;
    double time_LinSolver     = 0;

    int linSolveMethod = 0;
    enum class LinSolveMethod{ CG,CGsparse,Cholesky,CholeskySparse };

    Vec3f*  ps_cor      =0; // new Vec3d[nPoint];
    Vec3f*  ps_pred     =0; // new Vec3d[nPoint];
    Vec3f*  linsolve_b  =0; // new Vec3d[nPoint];
    Vec3f*  linsolve_yy =0; // new Vec3d[nPoint];



    int     nBonds =0; // number of bonds
    Quat4f* bparams=0; // bond parameters (l0,kP,kT,damp)
    int2*   bonds  =0; // indices of bonded points (i,j)
    float*  strain =0; // strain
    //float*  l0s    =0; // 
    Vec2f*  maxStrain=0;

    // ====== Invairiants

    float mass = 0;
    Vec3f cog  = Vec3fZero;
    Vec3f vcog = Vec3fZero;
    Mat3f I    = Mat3fZero;
    Vec3f L    = Vec3fZero;
    Vec3f torq = Vec3fZero;

    float F_residual = 0.0;

    // ====== Collision
    // ToDo: this should be moved to a separate class ?
    int       nBBs=0;
    Quat8f*  BBs=0; // bounding boxes (can be either AABB, or cylinder, capsula) 
    Buckets  pointBBs;    // buckets for collision detection
    Buckets  edgeBBs;
    Buckets  faceBBs;
    Buckets  pointChunks;  // chunks for parallelization, these points are copied to local memory when solving one edgeBBs chunk of bonds  

    // Faces are used just for ray-tracing, collision detection, etc
    int   nFaces=0;
    int4* faces=0; // indices of points tringles or quads, the last index is -1 is the face is a triangle, ot it can be used also to store face type if just triangles are used

    int nEdgeVertBonds=0;
    EdgeVertBond_f* edgeVertBonds=0; // indices of bonded points (i,j)

    // callback function pointer what to do in between iterations
    void (*user_update)(float dt);

    // Rotating frame
    //Vec3f p0{0.,0.,0.};
    //Vec3f ax{0.0,0.0,1.0};
    //float omega = 0.05;
    Quat4f accel{0.0f,0.0f,0.0f,0.0f};    // acceleration
    Quat4f rot0 {0.0f,0.0f,0.0f,0.0f};    // center of rotation
    //Quat4f omega{0.0f,0.0f,1.0f,0.05f}; // angular velocity, (xyz=axisxyz,w=magnitude)
    Quat4f omega{0.0f,0.0f,1.0f,0.05f};


    Vec3d hit_pos, hit_normal;


    //float maxAcc = 1e+6;
    float maxAcc = 1.0;
    float collision_damping = 0.002;
    //float collision_damping = 1.0;
    //float collision_damping = 1.1;   // if collision_damping > 1.0 then it is like successive over-relaxation (SOR) method ? https://en.wikipedia.org/wiki/Successive_over-relaxation

    //float kGlobal = 1e+6;
    
    float  kGlobal = 1e+7;
    float  dt      = 2e-3;    
    float  inv_dt2 = 1.0f / (dt * dt);

    //float damping = 1e-4;
    float damping  = 0.05;
    int    lastNeg = 0;
    // FIRE
    int    minLastNeg   = 5;
    float finc         = 1.1;
    float fdec         = 0.5;
    float falpha       = 0.98;
    float dt_max       = dt;
    float dt_min       = 0.1 * dt;
    float damp_max     = damping;
    float ff_safety    = 1e-16;
    float cv,cf;

    // ============ Functions

    inline void set_time_step( float dt_ ){
        dt      = dt_;
        inv_dt2 = 1.0f / (dt * dt);
    }

    void recalloc( int nPoint_, int nNeighMax_, int nBonds_=0){
        nPoint = nPoint_; nNeighMax = nNeighMax_;
        nNeighTot = nPoint*nNeighMax;
        _realloc( points, nPoint    );
        _realloc( forces, nPoint    );
        _realloc( vel,    nPoint    );
        _realloc( impuls, nPoint    );
        _realloc( params, nNeighTot );
        _realloc( kngs,   nNeighTot );
        _realloc( neighs, nNeighTot );
        _realloc( neighBs,nNeighTot );
        _realloc( neighB2s,nNeighTot );

        if(nBonds_>0){
            nBonds=nBonds_;
            _realloc( bonds,     nBonds );
            _realloc( strain,    nBonds );
            //_realloc( l0s,       nBonds );
            _realloc( maxStrain, nBonds );
            _realloc( bparams,   nBonds );
        }
    }

    void recallocBBs( int n, bool bPoint=true, bool bEdge=true, bool bFace=true, bool bClean=true ){
        nBBs = n;
        _realloc( BBs, nBBs );
        if(bPoint &&(nPoint>0) ){ pointBBs.realloc( nBBs, nPoint, true ); if(bClean)pointBBs.clean(); }
        if(bEdge  &&(nBonds>0) ){ edgeBBs .realloc( nBBs, nBonds, true ); if(bClean)pointBBs.clean(); }
        if(bFace  &&(nFaces>0) ){ faceBBs .realloc( nBBs, nFaces, true ); if(bClean)pointBBs.clean(); }
         //_realloc( , nBBs, nPoint );
    }

    // ================= Bounding Boxes

    void edgesToBBs(){
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

    void printBBs(){
        for(int i=0; i<nBBs; i++){ 
            int np = (pointBBs.cellNs) ? pointBBs.cellNs[i] : 0;
            int ne = (edgeBBs .cellNs) ? edgeBBs .cellNs[i] : -1;
            int nf = (faceBBs .cellNs) ? faceBBs .cellNs[i] : -1;
            printf( "BBs[%i] (%g,%g,%g) (%g,%g,%g) np=%i ned=%i nfc=%i \n", i, BBs[i].lo.x, BBs[i].lo.y, BBs[i].lo.z, BBs[i].hi.x, BBs[i].hi.y, BBs[i].hi.z, np, ne, nf );
        }
    };

    // =================== Picking

    int pick_point_brute( const Vec3f& ray0, const Vec3f& hray, float Rmax, int i0=-1, int i1=-1 ){
        if(i0<0){ i0=0;      }
        if(i1<0){ i1=nPoint; }
        float r2min =  Rmax*Rmax;
        int imin    = -1;
        for(int i=i0; i<i1; i++){
            float t;
            float r2 = rayPointDistance2( ray0, hray, points[i].f, t );
            //printf( "pick_point_brute ipick %i r %g p(%g,%g,%g)\n", i, sqrt(r2), points[i].f.x,points[i].f.y,points[i].f.z );
            if(r2<r2min){ imin=i; r2min=r2; }
            //float ti = raySphere( ray0, hRay, R, ps[i] );
            //if(ti<tmin){ imin=i; tmin=ti; }
        }
        //printf( "pick_point_brute ipick %i r2min %g \n", ipick, r2min );
        return imin;
    }

    int pick_bond_brute( const Vec3d& ray0, const Vec3d& hRay, double Rmax, int i0=-1, int i1=-1 ) const {
        if(i0<0){ i0=0;      }
        if(i1<0){ i1=nBonds; }
        double dist_min =  Rmax*Rmax;
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

    int pick_BBox( const Vec3d& ray0, const Vec3d& hRay, double tmax, int i0=-1, int i1=-1 ){
        if(i0<0){ i0=0;     }
        if(i1<0){ i1=nBBs; }
        if(i0<0){ i0=0; i1=nBonds; }
        double  tmin = tmax;
        int    imin = -1;
        for(int ib=i0; ib<i1; ib++){
            const Quat8f& bb = BBs[ib];
            double t = rayBox( ray0, hRay, (Vec3d)bb.lo.f, (Vec3d)bb.hi.f, hit_pos, hit_normal );
            if( t<tmin ){ imin=ib; tmin=t; }
        }
        return imin;
    };

    virtual int pick_nearest(Vec3d ray0, Vec3d hray, int& ipick, int mask, double Rmax ) override {
        if     (mask==1){ ipick=pick_point_brute((Vec3f)ray0,(Vec3f)hray,Rmax); return 1; }
        else if(mask==2){ ipick=pick_bond_brute ( ray0, hray, Rmax );           return 2; }
        return -1;
    };
    
    virtual int pick_all(Vec3d ray0, Vec3d hray, int* out, int mask, double Rmax ) override { return 0; };
    
    virtual void* getPickedObject(int picked, int mask) override { 
        if     (mask==1){ return (void*)&points[picked]; }
        else if(mask==2){ return (void*)&bonds [picked]; }
        return 0; 
    };

    // =================== Solver Using Projective Dynamics and Cholesky Decomposition

    void rhs_ProjectiveDynamics(Vec3f* pnew, Vec3f* b) {
        float idt2 = 1.0 / (dt * dt);
        for (int i = 0; i < nPoint; i++) {
            Vec3f bi;
            bi.set_mul(pnew[i], points[i].w * idt2);  // points[i].w is the mass
            int neighB_start = (i == 0) ? 0 : neighB2s[i-1];
            int neighB_end = neighB2s[i];
            for (int nb = neighB_start; nb < neighB_end; nb++) {
                int ib = neighs[nb];
                float k = params[ib].y;  // assuming params.y is the spring constant
                int i_ = bonds[ib].x;
                int j_ = bonds[ib].y;
                int j = (i_ == i) ? j_ : i_;
                Vec3f d = points[i].f -  points[j].f;
                bi.add_mul(d, k * params[ib].x / d.norm());  // params[ib].x is l0
            }
            b[i] = bi;
        }
    }

    __attribute__((hot)) 
    Vec3f rhs_ProjectiveDynamics_i(int i, const Vec3f* pnew) {
        float idt2 = 1.0 / (dt * dt);
        Vec3f bi; bi.set_mul(pnew[i], points[i].w*idt2);  // points[i].w is the mass
        int2* ngi = neighBs + (i*nNeighMax);
        int ni = 0;
        for( int ing=0; ing<nNeighMax; ing++ ){
            int2  ng = ngi[ing];
            int   ib = ng.y;
            if(ib<0) break;
            float k  = bparams[ib].y;  // assuming params.y is the spring constant
            float l0 = bparams[ib].x;
            int2   b  = bonds[ib];
            int j     = (b.x == i) ? b.y : b.x;
            //printf( "rhs[i=%2i,ib=%2i,j=%i] k=%g l0=%g\n", i, ib, j, k, l0 );
            Vec3f d = pnew[i] -  pnew[j];
            bi.add_mul(d, k * l0 / d.norm());  // params[ib].x is l0
            ni++;
        }
        //printf( "rhs[i=%2i] ni=%i bi(%g,%g,%g)\n", i, ni, bi.x,bi.y,bi.z );
        return bi;
    }

    __attribute__((hot)) 
    void rhs_ProjectiveDynamics_(const Vec3f* pnew, Vec3f* b) {
        float idt2 = 1.0 / (dt * dt);
        for (int i=0; i<nPoint; i++) {
            b[i] = rhs_ProjectiveDynamics_i(i,pnew);
        }
    }

    void run_Cholesky(int niter) {
        Vec3f*  ps_cor  = new Vec3f[nPoint];
        Vec3f*  ps_pred = new Vec3f[nPoint];
        Vec3f*  b       = new Vec3f[nPoint];
        Vec3f*  yy      = new Vec3f[nPoint];

        const int m=3;
        memcpy(ps_cor, points, nPoint * sizeof(Vec3f));
        float idt = 1.0 / dt;
        float dt2 = idt * idt;
        for (int iter = 0; iter < niter; iter++) {
            // Evaluate forces (assuming you have a method for this)
            //evalForces();
            for (int i = 0; i < nPoint; i++) { ps_pred[i] = points[i].f + vel[i].f*dt + forces[i].f*dt2; }  // Predict step
            for (int i = 0; i < nPoint; i++) { if( kFix[i]>0 ) { ps_pred[i] = points[i].f; } }              // Apply fixed constraints
            rhs_ProjectiveDynamics(ps_pred, b);             // Compute right-hand side
            // Solve using LDLT decomposition (assuming you have this method)
            //solve_LDLT_sparse(b, ps_cor);
            Lingebra::forward_substitution_sparse           ( nPoint,m,  LDLT_L, (float*)b, (float*)yy, neighsLDLT, nNeighMax );
            for (int i=0; i<nPoint; i++){ yy[i].mul(1/LDLT_D[i]); } // Diagonal 
            Lingebra::forward_substitution_transposed_sparse( nPoint,m, LDLT_L, (float*)yy, (float*)ps_cor, neighsLDLT, nNeighMax );

            // Compute residual
            float res = 0.0;
            for (int i=0;i<nPoint;i++) {
                Vec3f  d = ps_cor[i] - points[i].f;
                float l = d.norm();
                if (l > res) res = l;
            }
            //printf("residual[%d] : %f\n", iter, res);

            // Update velocity and points
            for (int i=0;i<nPoint;i++) {
                Vec3f dv = ps_cor[i] - ps_pred[i];
                dv.mul(1.0 / dt);
                vel   [i].f.add( dv );
                points[i].f = ps_cor[i];
            }

            // Call user update function if set
            if (user_update) {
                user_update(dt);
            }
        }

        delete[] ps_cor;
        delete[] ps_pred;
        delete[] b;
    }


    void realloc_LinearSystem( bool bCG=true, bool bCholesky=true, int nNeighMaxLDLT_=32, bool bDens=true ){
        printf( "OrbSim_f::realloc_LinearSystem() nPoint=%i nNeighMaxLDLT=%i nNeighMax=%i\n", nPoint, nNeighMaxLDLT_, nNeighMax );
        nNeighMaxLDLT = nNeighMaxLDLT_;  if(nNeighMaxLDLT<nNeighMax){ printf("ERROR in OrbSim::prepare_Cholesky(): nNeighMaxLDLT(%i)<nNeighMax(%i) \n", nNeighMaxLDLT, nNeighMax); exit(0); }
        int n2 = nPoint*nPoint;
        // --- sparse Linear system Matrix and its Cholesky L*D*L^T decomposition
        if(bDens){ _realloc0( PDmat      ,n2    , 0.0f      ); }
        _realloc0( ps_cor     ,nPoint, Vec3fZero );
        _realloc0( ps_pred    ,nPoint, Vec3fZero );
        _realloc0( linsolve_b ,nPoint, Vec3fZero );
        if(bCholesky){
            _realloc0(  linsolve_yy,nPoint, Vec3fZero);
            _realloc0( LDLT_L, n2,     0.0f );
            _realloc0( LDLT_D, nPoint, 0.0f );
            _realloc0( neighsLDLT, nPoint*nNeighMaxLDLT, -1 );
        }
        //if(bCG){  cgSolver.realloc(nPoint,3);}
    }

    __attribute__((hot)) 
    void make_PD_Matrix( float* A, float dt ) {
        printf( "OrbSim_f::make_PD_Matrix() dt=%g\n", dt );
        //for(int i=0; i<nPoint; i++){ printf("mass[%i] %g\n", points[i].w ); }; // debug
        //for(int i=0; i<nBonds; i++){ printf("ks[%i] %g\n",   params[i].y ); }; // debug
        int n = nPoint;
        int n2 = n*n;
        //for(int i=0; i<np2; i++){ PDmat[i]=0.0; }
        float idt2 = 1.0f / (dt * dt);
        for (int i = 0; i < n; ++i) { // point i
            float Aii = points[i].w * idt2; // Assuming points[i].w stores the mass
            //printf( "==i,i %i %g %g %g \n", i, Aii, points[i].w, idt2  );
            int2* ngi = neighBs + (i*nNeighMax);
            for( int ing=0; ing<nNeighMax; ing++ ){
                int2  ng = ngi[ing];
                int   ib = ng.y;
                if(ib<0) break;
                //if(ib<0){ printf("ERROR int OrbSim_d::make_PD_Matrix() invalid neighbor [i=%i,ing=%i] ng(%i,%i) neigh=%i\n", i, ing, ib, ng.x, neighs[i*nNeighMax+ing] ); exit(0); }  
                //printf( "make_PD_Matrix()[i=%i,ing=%i] ng(%i,%i)\n", i, ing, ib, ng.x );
                int2   b = bonds[ib];
                //int j  = neighs[b.x];
                int j    = b.y;
                double k = bparams[ib].y; // Assuming params[ib].y stores the spring constant

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
    void make_PDmat_sparse( SparseMatrix<float>& A, float dt, bool bRealloc ) {
        printf( "OrbSim_f::make_PDmat_sparse() dt=%g nNeighMax=%i \n", dt, nNeighMax );
        float idt2 = 1.0f / (dt * dt);
        // --- count neighbors
        if( bRealloc ){ A.realloc( nPoint, nNeighMax+1 ); }
        for (int i=0; i<nPoint; i++) {
            float Aii = points[i].w * idt2; 
            int2* ngi = neighBs + (i*nNeighMax);
            for( int ing=0; ing<nNeighMax; ing++ ){
                int2  ng = ngi[ing];
                int   ib = ng.y;
                if(ib<0) break;
                int2   b = bonds[ib];
                int j    = b.y;
                double k = bparams[ib].y; // Assuming params[ib].y stores the spring constant
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

    void prepare_LinearSystem( float dt, bool bRealloc=true, bool bCG=true, bool bCholesky=true, int nNeighMaxLDLT_=32, bool bDens=true ){ 
        printf( "OrbSim_d::prepare_LinearSystem() dt=%g nNeighMaxLDLT_=%i\n", dt, nNeighMaxLDLT_ );
        //nNeighMaxLDLT=nNeighMaxLDLT_;
        if(bRealloc)realloc_LinearSystem( bCG, bCholesky, nNeighMaxLDLT_, bDens );
        this->dt = dt;
        int n2 = nPoint*nPoint;

        timeit( "TIME make_PDmat_sparse() t= %g [MTicks]\n", 1e-6, [&](){ make_PDmat_sparse( PDsparse, dt, true ); });
        mat2file<int>  ( "PDsparse_inds.log",  nPoint, nNeighMax+1, (int*)  PDsparse.inds, "%5i " );
        mat2file<float>( "PDsparse_vals.log",  nPoint, nNeighMax+1, (float*)PDsparse.vals );

        if(bDens){
            timeit( "TIME make_PD_Matrix()    t= %g [MTicks]\n", 1e-6, [&](){ make_PD_Matrix(    PDmat,    dt );       });
            if( PDsparse.checkDens( PDmat, 1e-16, true ) ){ printf("ERROR in OrbSim_d::prepare_LinearSystem() PDsparse does not match PDmat => exit \n"); exit(0); }
        }
        
        // if( bCG ){
        //     // setup Conjugate-Gradient solver
        //     cgSolver.setLinearProblem(  nPoint, 3, (float*)ps_cor, (float*)linsolve_b, false );
        //     cgSolver.initDiagPrecond( PDmat );
        //     cgSolver.dotFunc = [&](int n, float* x, float* y){  
        //         long t0 = getCPUticks();
        //         PDsparse.dot_dens_vector_m( 3, x, y ); 
        //         time_cg_dot += (getCPUticks() - t0)*1e-6; 
        //     };
        // }
        
        if( bCholesky ){
            for(int i=0; i<nPoint; i++){  for(int j=0; j<nNeighMaxLDLT; j++){ if(j>=nNeighMax){neighsLDLT[i*nNeighMaxLDLT+j]=-1;}else{ neighsLDLT[i*nNeighMaxLDLT+j]=neighs[i*nNeighMax+j]; };  } }
            
            bool bSparse = false;
            //bool bSparse = true;
            if(bSparse){
                Lingebra::CholeskyDecomp_LDLT_sparse( PDmat, LDLT_L, LDLT_D, neighsLDLT, nPoint, nNeighMaxLDLT );
                sortNeighs( nPoint, nNeighMaxLDLT, neighsLDLT);
            }else{
                Lingebra::CholeskyDecomp_LDLT( PDmat, LDLT_L, LDLT_D, nPoint );
                Lsparse.fromDense( nPoint, LDLT_L, 1.e-16 );
                LsparseT.fromFwdSubT( nPoint, LDLT_L);
            }
            {   printf( "# ---------- Test Cholesky Lsparse \n");
                for(int i=0; i<nPoint; i++){ ps_pred[i] = points[i].f*1.1; }; 
                rhs_ProjectiveDynamics( ps_pred, linsolve_b );
                Lingebra::forward_substitution_m( LDLT_L, (float*)linsolve_b,  (float*)linsolve_yy, nPoint, 3 );
                Lsparse.fwd_subs_m( 3,  (float*)linsolve_b,  (float*)ps_cor );
                checkDist( nPoint, ps_cor, linsolve_yy, 2 );
            }
            //mat2file<float>( "LDLT_L.log", nPoint,nPoint, LDLT_L );
            Lsparse.fprint_inds("Lsparse_inds.log");
            //Lsparse.fprint_vals("Lsparse_vals.log");
            LsparseT.fprint_inds("LsparseT_inds.log");
            //LsparseT.fprint_vals("LsparseT_vals.log");
            //exit(0);
        }

    }

    __attribute__((hot)) 
    void run_LinSolve(int niter) {
        //printf( "OrbSim::run_LinSolve() \n" );
        const int m=3;
        memcpy(ps_cor, points, nPoint * sizeof(Vec3f));
        float dt2 = dt * dt;
        float inv_dt = 1/dt;
        cleanForce();
        for (int iter = 0; iter < niter; iter++) {
            // Evaluate forces (assuming you have a method for this)
            //evalForces();
            // Predict step
            for (int i=0;i<nPoint;i++){ 
                ps_pred[i] = points[i].f + vel[i].f*dt + forces[i].f*dt2; 
                //printf( "ps_pred[%3i](%10.6f,%10.6f,%10.6f) v(%10.6f,%10.6f,%10.6f) p(%10.6f,%10.6f,%10.6f) dt=%g \n", i, ps_pred[i].x,ps_pred[i].y,ps_pred[i].z, vel[i].x,vel[i].y,vel[i].z, points[i].x,points[i].y,points[i].z, dt );
            }

            // Apply fixed constraints
            //for (int i = 0; i < nPoint; i++) {    if (kFix[i] > 0) { ps_pred[i] = points[i].f; } }
            // Compute right-hand side
            rhs_ProjectiveDynamics(ps_pred, linsolve_b );

            long t0 = getCPUticks();
            switch( (LinSolveMethod)linSolveMethod ){
                case LinSolveMethod::CholeskySparse:{
                    // Solve using LDLT decomposition (assuming you have this method)
                    //solve_LDLT_sparse(b, ps_cor);
                    Lingebra::forward_substitution_sparse           ( nPoint,m,  LDLT_L, (float*)linsolve_b,  (float*)linsolve_yy, neighsLDLT, nNeighMaxLDLT );
                    for (int i=0; i<nPoint; i++){ linsolve_yy[i].mul(1/LDLT_D[i]); } // Diagonal 
                    Lingebra::forward_substitution_transposed_sparse( nPoint,m,  LDLT_L, (float*)linsolve_yy, (float*)ps_cor,      neighsLDLT, nNeighMaxLDLT );
                } break;
                case LinSolveMethod::Cholesky:{
                    //Lingebra::forward_substitution_m( LDLT_L, (float*)linsolve_b,  (float*)linsolve_yy, nPoint,m );
                    //Lsparse.fwd_subs_m( m,  (float*)linsolve_b,  (float*)ps_cor );
                    //if( checkDist( nPoint, ps_cor, linsolve_yy, 1 ) ){ printf("ERROR run_LinSolve.checkDist() => exit()"); exit(0); };
                    Lsparse.fwd_subs_m( m,  (float*)linsolve_b,  (float*)linsolve_yy );
                    for (int i=0; i<nPoint; i++){ linsolve_yy[i].mul(1/LDLT_D[i]); } // Diagonal 
                    LsparseT.fwd_subs_T_m( m,  (float*)linsolve_yy,  (float*)ps_cor );
                    //Lingebra::forward_substitution_T_m( LDLT_L, (float*)linsolve_yy, (float*)ps_cor,      nPoint,m );
                    //if( checkDist( nPoint, ps_pred, ps_cor, 2 ) ){ printf("ERROR run_LinSolve.checkDist() => exit()"); exit(0); };

                } break;
            }
            time_LinSolver += (getCPUticks()-t0)*1e-6;            
            for (int i=0;i<nPoint;i++) {
                // To-Do : We need more rigorous approach how to update velocity
                // position-based velocity update
                float vr2 = vel[i].norm2();
                Vec3f v    = (ps_cor[i] - points[i].f);
                v.mul(inv_dt);
                //v.mul( sqrt(  vel[i].norm2()/( v.norm2() + 1e-8 ) ) ); // This is unphysicsl
                vel[i].f = v;
                // update positions
                points[i].f = ps_cor[i];
            }
            // Call user update function if set
            //if (user_update){ user_update(dt);}
        }
    }

    __attribute__((hot)) 
    void run_Cholesky_omp_simd(int niter) {
        //printf( "OrbSim::run_LinSolve() \n" );
        const int m=3;
        memcpy(ps_cor, points, nPoint * sizeof(Vec3f));
        float dt2 = dt * dt;
        float inv_dt = 1/dt;
        cleanForce();
        int iter=0;
        for (iter = 0; iter < niter; iter++) {
            #pragma omp simd
            for (int i=0;i<nPoint;i++){ ps_pred[i] = points[i].f + vel[i].f*dt + forces[i].f*dt2;  }
            #pragma omp simd
            for (int i=0;i<nPoint;i++){ linsolve_b[i] = rhs_ProjectiveDynamics_i(i,ps_pred); }
            Lsparse.fwd_subs_m( m,  (float*)linsolve_b,  (float*)linsolve_yy );            
            #pragma omp simd
            for(int i=0; i<nPoint; i++){ linsolve_yy[i].mul(1/LDLT_D[i]); } // Diagonal 
            LsparseT.fwd_subs_T_m( m,  (float*)linsolve_yy,  (float*)ps_cor );
            #pragma omp simd
            for (int i=0;i<nPoint;i++) {
                vel[i].f    = (ps_cor[i] - points[i].f) * inv_dt;
                points[i].f = ps_cor[i];
            }
        }
    }

    // =================== Truss Simulation

    void updateInveriants(bool bPrint=false){
        mass=0;
        cog =Vec3fZero;
        vcog=Vec3fZero;
        for(int i=0; i<nPoint; i++){ float mi=points[i].w; mass+=mi;  cog.add_mul( points[i].f, mi ); vcog.add_mul( vel[i].f, mi ); }
        cog .mul( 1.0/mass );
        vcog.mul( 1.0/mass );
        I=Mat3fZero;
        L=Vec3fZero;
        torq=Vec3fZero;
        for(int i=0; i<nPoint; i++){ 
            float mi=points[i].w; 
            Vec3f d; 
            d   .set_sub   ( points[i].f, cog    );
            L   .add_crossw( d, vel[i]   .f, mi  );
            torq.add_cross ( d, forces[i].f      );
            I   .add_outer ( d, d, mi            );
        }
        if(bPrint){
            printf( "OrbSim::updateInveriants() mass %g cog(%g,%g,%g) vcog(%g,%g,%g)  L(%g,%g,%g) torq(%g,%g,%g) \n", mass, cog.x,cog.y,cog.z, vcog.x,vcog.y,vcog.z, L.x,L.y,L.z, torq.x,torq.y,torq.z );
            printf( "OrbSim::updateInveriants() I \n" ); I.print();
        }
    }

    void evalTrussForce_neighs(int iG){
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
            float k = kGlobal;
            f.f.add_mul( d, (k*(li-params[j].x)/li) );

            //printf( "p[%i,ij=%i,j=%i] li=%7.3f dl=%8.5e fi=%8.5e e=%8.5e par(%7.3f,%8.5e,%8.5e,%8.5e) \n", iG,ij,ja, li, li-params[j].x, fi,ei, params[j].x,params[j].y,params[j].z,params[j].w );
        }
        forces[iG] = f; // we may need to do += in future
    }

    /*
    void evalTrussForces_neighs(){
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
                
                // float fi,ei = springForce( li, fi, params[j] );
                // //f.add( Quat4f{ d*(fi/l), ei } );
                // f.f.add_mul( d, fi/li );
                
                float k = kGlobal;
                f.f.add_mul( d, (k*(li-params[j].x)/li) );

                //printf( "p[%i,ij=%i,j=%i] li=%7.3f dl=%8.5e fi=%8.5e e=%8.5e par(%7.3f,%8.5e,%8.5e,%8.5e) \n", iG,ij,ja, li, li-params[j].x, fi,ei, params[j].x,params[j].y,params[j].z,params[j].w );
            }
            forces[iG] = f; // we may need to do += in future
        } 
        //exit(0);
    }
    */

    inline void evalTrussForce_neighs2(int iG){
        const float Adamp = collision_damping*0.5/dt;
        //const int iG = get_global_id(0);
        const Quat4f p = points[iG];
        const Quat4f v = vel   [iG];
        Quat4f f = Quat4f{0.0f,0.0f,0.0f,0.0f};
        //printf( "--- p[%i] \n", iG );
        //#pragma omp simd
        for(int ij=0; ij<nNeighMax; ij++){
            const int j  = nNeighMax*iG + ij;
            const int2 b = neighBs[j];
            if(b.x == -1) break;
            const Quat4f par = bparams[b.y];
            //f.add( springForce( points[ja].f - p.f, params[j] ) );
            Quat4f d = points[b.x];
            d.f.sub( p.f );
            const float l  = d.f.norm();
            float k        = kGlobal;
            float fl       = k*(l-par.x);
            const float invL = 1/l;


            // const float dv  = d.f.dot( vel[b.x].f - v.f );
            // if(dv<0){ fl*=1-dv; }


            // // collision damping
            
            // const float dv  = d.f.dot( vel[b.x].f - v.f );
            // float imp = Adamp * p.w*d.w*dv/(p.w+d.w);
            // //float imp = 0.1* 0.5 * p.w*d.w*dv/(p.w+d.w);
            // //const float imp = 0;
            // //imp/=dt;
            float imp = 0;

            f.f.add_mul( d.f, ( imp + fl )*invL );
            //printf( "p[%i,ij=%i,j=%i] li=%7.3f dl=%8.5e fi=%8.5e e=%8.5e par(%7.3f,%8.5e,%8.5e,%8.5e) \n", iG,ij,ja, li, li-params[j].x, fi,ei, params[j].x,params[j].y,params[j].z,params[j].w );
        }
        forces[iG] = f; // we may need to do += in future
    }

    void evalTrussForces_neighs2(){
        for(int iG=0; iG<nPoint; iG++){ evalTrussForce_neighs2(iG); } 
    }

    void evalTrussForces_neighs(){
        for(int iG=0; iG<nPoint; iG++){ evalTrussForce_neighs(iG); } 
    }

    

    void evalTrussForces_bonds(){
        for(int i=0; i<nBonds; i++){
            int2  b = bonds[i];
            const Quat4f& pi = points[b.x];
            const Quat4f& pj = points[b.y];
            Vec3f d = pj.f - pi.f;
            float l = d.norm();
            //float fi,ei = springForce( li, fi, bparams[i] );
            float k = kGlobal;
            float f = k*(l-bparams[i].x);

            // Limit acceleration force to improve stability   -- it does not seem to help
            // float mmin = (pi.w < pj.w) ? pi.w : pj.w;
            // if( (fabs(f)/mmin) > maxAcc ){   // here it would be more efficient to use momentum rather than force 
            //     f = maxAcc*mmin;
            // }

            // collision damping
            //float vi   = d.dot( vel[b.x].f );
            //float vj   = d.dot( vel[b.y].f );
            //float dv   = vj   - vi;
            float invL = 1./l;
            float dv   = d.dot( vel[b.y].f - vel[b.x].f )*invL;
            float mcog = pj.w + pi.w;
            float imp  = collision_damping * pi.w*pi.w*dv/mcog;

            d.mul( (imp + f)*invL );          
            forces[b.x].f.add(d);
            forces[b.y].f.sub(d);
        } 
    }

    Vec3f getEdgeVertPos( int i ){
        Vec3i b = edgeVertBonds[i].verts;
        const Quat4f& pa = points[b.x];
        const Quat4f& pb = points[b.y];
        const Quat4f& pc = points[b.z];
        float c = edgeVertBonds[i].c;
        float mc = 1-c;
        return pa.f*mc + pb.f*c;
    }

    //void evalEdgeVert( Vec3i b, float c, float K ){
    void evalEdgeVert( int i ){
        Vec3i b = edgeVertBonds[i].verts;
        float c = edgeVertBonds[i].c;
        float K = edgeVertBonds[i].K;
        // ToDo: perhaps we should interpolate it by B-spline to make the path more smooth
        // ToDo: Apply Force in the direction of the edge, constrain perpendicular to the edge
        // ToDo: Damping ( collision damping )
        
        // fit vert to edge
        const Quat4f& pa = points[b.x];
        const Quat4f& pb = points[b.y];
        const Quat4f& pc = points[b.z];
        float mc = 1-c;
        Vec3f d = pc.f - (pa.f*mc + pb.f*c);
        //glColor3f(1.0,0.0,0.0); Draw3D::drawVecInPos( d, pc.f*1.0 );
        
        float invL = 1./d.norm();
        float dv   = d.dot( vel[b.z].f - vel[b.x].f*mc - vel[b.y].f*c )*invL;
        float mab  = pa.w*mc + pb.w*c;
        float imp  = collision_damping * pc.w*mab*dv/( pc.w  + mab );
        //imp/=dt; 

        d.mul( K + imp*invL ); 
        edgeVertBonds[i].f = d;
        //printf( "evalEdgeVert[%i,%i,%i] d(%g,%g,%g) c=%g \n", b.x,b.y,b.z, d.x,d.y,d.z, c );   
        forces[b.x].f.add_mul(d,mc);
        forces[b.y].f.add_mul(d, c);
        forces[b.z].f.sub(d);

        // Force
        // const Quat4f& pa = points[b.x];
        // const Quat4f& pc = points[b.z];
        // Vec3f d = pc.f - pa.f;          
        // // Damping
        // //float dv   = d.dot( vel[b.x].f + vel[b.y].f - vel[b.z].f );
        // d.mul( K );
        // forces[b.x].f.add(d);
        // //forces[b.y].f.add_mul(d, c);
        // forces[b.z].f.sub(d);
    }

    void evalEdgeVerts(){
        //for(int i=0; i<nEdgeVertBonds; i++){ evalEdgeVert( edgeVertBonds[i].verts, edgeVertBonds[i].c, edgeVertBonds[i].K ); }
        for(int i=0; i<nEdgeVertBonds; i++){ evalEdgeVert( i ); }
    }

    void evalTrussCollisionImpulses_bonds( float rate=1.0 ){
        // Collision forces are based on momentum-conserving impulses, it does not need to know anything about the interaction potential (e.g. spring constants)
        for(int i=0; i<nBonds; i++){
            int2  b = bonds[i];
            const Quat4f& pi = points[b.x];
            const Quat4f& pj = points[b.y];
            Vec3f d  = pi.f - pj.f;
            float l  = d.normalize();
            float vi = d.dot( vel[b.x].f );
            float vj = d.dot( vel[b.y].f );
            float mcog = pj.w + pi.w;
            float dv   = vj   - vi;
            //float vcog = (pi.w*vi + pj.w*vj)/(pi.w+pj.w);
            // ---- Inelastic collision
            //float dvi  = vcog - vi;
            //float imp1 = dvi*pi.w;
            // Analytical Derivation:
            // dvi = ( pi.w*vi + pj.w*vj                        )/(pi.w+pj.w) - (pi.w+pj.w)*vi/(pi.w+pj.w)
            //     = ( pi.w*vi + pj.w*vj -  pi.w*vi - pj.w*vi   )/(pi.w+pj.w)
            //     = (           pj.w*vj            - pj.w*vi   )/(pi.w+pj.w)
            //     = (           pj.w*(vj-vi)                   )/(pi.w+pj.w)
            // imp1 = dvi*pi.w = pi.w*pj.w*(vj-vi)/(pi.w+pj.w)
            float imp = pi.w*pi.w*dv/mcog;
            //if( i==146 ){ // Debug
            //    float dvj  = vcog - vj;
            //    float imp2 = dvj*pj.w; // shoould be the same as imp1 
            //    printf( "evalTrussCollisionForces_bonds[%i] imp1 %g imp2 %g \n", i, imp1, imp2  );
            //}
            // apply force to points
            d.mul( imp*rate );
            forces[b.x].f.add( d );
            forces[b.y].f.sub( d );
        } 
    }

    void evalBondTension(){
        for(int i=0;i<nBonds; i++ ){
            int2  b  = bonds[i];
            float l0 = bparams[i].x;
            //float l0 = l0s[i];
            float l   = (points[b.y].f-points[b.x].f).norm();
            float s   = (l-l0)/l0;
            //if( fabs(s)>0.5 ){ printf( "evalBondTension[%i] strain=%g l=%g l0=%g\n", i, s, l, l0 ); }
            //if(i==6272){ printf( "evalBondTension[%i](%i,%i) strain=%g l=%g l0=%g\n", i, b.x,b.y, s, l, l0 ); }
            strain[i] = s;
            // ToDo: break the bond if strain > maxStrain;
        }
        //exit(0);
    }

    void applyForceRotatingFrame_i( int i, Vec3f p0, Vec3f ax, float omega ){
        const Quat4f& p = points[i];
        const Quat4f& v = vel   [i];
        Vec3f d,f;
        // Coriolis force     = 2*m*omega*v
        f.set_cross(ax,v.f);        
        f.mul( 2.0*omega );
        // centrifugal force  = r*m*omega^2
        d.set_sub(p.f,p0);
        d.makeOrthoU(ax);
        f.add_mul( d, omega*omega );     
        // apply force
        forces[i].f.add_mul(f, p.w );
    }

    void applyForceCentrifug_i( int i, Vec3f p0, Vec3f ax, float omega ){
        const Quat4f& p = points[i];
        Vec3f d;
        d.set_sub(p.f,p0);
        d.makeOrthoU(ax);   
        forces[i].f.add_mul(d, p.w*omega*omega );
    }

    void applyForceRotatingFrame( Vec3f p0, Vec3f ax, float omega ){
        float omega2 = omega*omega;
        Vec3f omega_ax = ax*omega*2.0;
        for(int i=0;i<nPoint; i++ ){
            const Quat4f& p = points[i];
            const Quat4f& v = vel   [i];
            //Vec3f f; f.set_cross(ax,p.f-p0);
            Vec3f d,f;
            d.set_sub(p.f,p0);
            d.makeOrthoU(ax);
            f.set_mul( d, omega2 );     // centrifugal force  = r*m*omega^2
            f.add_cross(omega_ax,v.f);  // Coriolis force     = 2*m*omega*v
            forces[i].f.add_mul(f, p.w );
            //forces[i].f.add_mul( f, p.w*omega2 );
        }
    }

    void applyForceCentrifug( Vec3f p0, Vec3f ax, float omega ){
        for(int i=0;i<nPoint; i++ ){ applyForceCentrifug_i( i, p0, ax,omega ); }
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

    float getFmax(){ 
        float fmax=0;
        for(int i=0; i<nPoint; i++){ float f=forces[i].norm(); fmax=fmax>f?fmax:f; }   
        //printf( "|fmax|=%g\n", fmax );
        return fmax;
    }

    void cleanForce (){ for (int i=0; i<nPoint; i++){ forces[i]=Quat4fZero;   } };
    void cleanVel   (){ for (int i=0; i<nPoint; i++){ vel   [i]=Quat4fZero;   } };
    void cleanImpuls(){ for (int i=0; i<nPoint; i++){ impuls[i]=Quat4fZero;   } };
    void setKngs    (){ for (int i=0; i<nNeighTot; i++){ kngs[i]=params[i].z; } };   // kPull

    float move_GD(float dt){
        float ff=0;
        for(int i=0;i<nPoint; i++ ){
            Quat4f p = points[i];
            Quat4f f = forces[i];
            ff += f.f.norm2();
            p.f.add_mul( f.f, dt/p.w );
            //printf( "move_GD[%i] |d|=%g |f|=%g dt/m=%g m=%g \n", i, f.f.norm() * dt/p.w, f.f.norm(), dt/p.w, p.w );
            points[i]=p;
        }
        return ff;
    }

    inline float move_i_MD(int i, float dt, float cdamp ){
        Quat4f f = forces[i];
        Quat4f p = points[i];
        Quat4f v = vel   [i];

        //float vf = v.f.dot(f.f);
        //if( vf>0.0 ){ f.f.mul( 0.99 ); }

        v.f.mul( cdamp );
        v.f.add_mul( f.f, dt/p.w );
        p.f.add_mul( v.f, dt     );
        //printf( "move_GD[%i] |d|=%g |f|=%g dt/m=%g m=%g \n", i, f.f.norm() * dt/p.w, f.f.norm(), dt/p.w, p.w );
        vel   [i]=v;
        points[i]=p;
        return f.f.norm2();
    }

    float move_MD(float dt, float damp=0.0f ){
        float cdamp = 1.0f - damp;
        float ff = 0.0; 
        for(int i=0;i<nPoint; i++ ){ move_i_MD(i, dt, cdamp ); }
        return ff;
    }

int run( int niter, float dt, float damp  ){
    float f2 = -1;
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
        //printf( "OrbSim_f::run[%i] |F|=%g\n", itr, sqrt(f2) );
    }
    return niter;
}

inline void setOpt( float dt_, float damp_ ){
    dt      = dt_max   = dt_;  dt_min=0.1*dt_max;
    damping = damp_max = damp_;
    cleanForce( );
    cleanVel  ( );
}

void FIRE_update( float& vf, float& vv, float& ff, float& cv, float& cf ){
    float cs = vf/sqrt(vv*ff);
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


int run_omp( int niter_max, bool bDynamic, float dt_, float damp_ ){
    //printf( "run_omp() niter_max %i dt %g Fconv %g Flim %g timeLimit %g outE %li outF %li \n", niter_max, dt, Fconv, Flim, timeLimit, (long)outE, (long)outF );
    float cdamp = 1.0f - damp_;
    float E,F2,ff,vv,vf;
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
            forces[iG] = Quat4fZero;
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
                Quat4f p = points[i];
                Quat4f f = forces[i];
                Quat4f v = vel   [i];
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
            //printf( "OrbSim::run_omp() itr %i/%i E=%g |F|=%g time= %g [ms]( %g [us/%i iter]) \n", itr,niter_max, E, F_residual, time*1e+3, time*1e+6/itr, itr );
            time+=dt;
            itr++; 
        }
        } // if(itr<niter){
    }
    //{
    //float t = (getCPUticks() - T0)*tick2second;
    //if(itr>=niter_max)if(verbosity>0)printf( "run_omp() NOT CONVERGED in %i/%i E=%g |F|=%g time= %g [ms]( %g [us/%i iter]) \n", itr,niter_max, E, sqrt(F2), t*1e+3, t*1e+6/itr, itr );
    //}
    return itr;
}

};   // OrbSim_f

#endif
