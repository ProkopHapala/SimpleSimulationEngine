
#ifndef  TrussDynamics_f_h
#define  TrussDynamics_f_h

#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "datatypes.h"  

#include "Buckets.h"

#include "raytrace.h"
#include "geom3D.h"
#include "Interfaces.h"

#include "Cholesky.h"
#include "SparseMatrix.h"
#include "SparseMatrix2.h"

// ========== Free Functions ==========

inline float springForce( float l, float& f, Quat4f par ){
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

inline Quat4f springForce( Vec3f d, Quat4f par ){
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

float checkDist(int n, const Vec3f* vec, const Vec3f* ref, int verb=1, float tol=1e-12 );
void fitAABB( Quat8f& bb, int n, int* c2o, Quat4f * ps );
void fitAABB_edge( Quat8f& bb, int n, int* c2o, int2* edges, Quat4f * ps );
void updatePointBBs(const Buckets& buckets, Quat8f* BBs, Quat4f* points, bool bInit=true);
void updateEdgeBBs(const Buckets& buckets, Quat8f* BBs, int2* edges, Quat4f* points, bool bInit=true);

// ========== Helper classes ==========

struct EdgeVertBond_f{ 
    Vec3i verts; 
    float c; 
    float K; 
    Vec3f f=Vec3fZero; 
};

// ========== main class ==========

class TrussDynamics_f : public Picker { public:
    float time=0;
    int nPoint=0, nNeighMax=0, nNeighTot=0;
    // cpu buffers
    Quat4f* points=0;  // position and mass
    Quat4f* forces=0;  // force and energy
    Quat4f* vel   =0;  // velocity
    Quat4f* impuls=0;  // accumulated impulse from corrector
    Quat4f* bvec  =0;

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
    SparseMatrix<float>   PDsparse;
    //CGsolver cgSolver;
    double time_LinSolver     = 0;

    int nSolverIters   = 10;
    int linSolveMethod = 2; // 0=CG, 1=CGsparse, 2=Cholesky
    enum class LinSolveMethod{ CG,CGsparse,Cholesky,CholeskySparse };

    Quat4f* ps_cor      =0; // new Vec3d[nPoint];
    Quat4f* ps_pred     =0; // new Vec3d[nPoint];
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

    // ============ inline Functions

    inline void set_time_step( float dt_ ){
        dt      = dt_;
        inv_dt2 = 1.0f / (dt * dt);
    }

    inline float getLinearBondStiffness( int ib ){
        // (l0,kPress,kTens,damp)
        //return  bparams[ib].y;  //k = kPress;
        return  bparams[ib].z;  //k = kTens;  
    }

    inline void reallocFixed(){ _realloc0( kFix, nPoint, 0.0f ); }
    inline void cleanForce (){ for (int i=0; i<nPoint; i++){ forces[i]=Quat4fZero;   } };
    inline void cleanVel   (){ for (int i=0; i<nPoint; i++){ vel   [i]=Quat4fZero;   } };
    inline void cleanImpuls(){ for (int i=0; i<nPoint; i++){ impuls[i]=Quat4fZero;   } };
    inline void setKngs    (){ for (int i=0; i<nNeighTot; i++){ kngs[i]=params[i].z; } };   // kPull

    inline void applyForceRotatingFrame_i( int i, Vec3f p0, Vec3f ax, float omega ){
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

    inline void applyForceCentrifug_i( int i, Vec3f p0, Vec3f ax, float omega ){
        const Quat4f& p = points[i];
        Vec3f d;
        d.set_sub(p.f,p0);
        d.makeOrthoU(ax);   
        forces[i].f.add_mul(d, p.w*omega*omega );
    }
    
    inline Vec3f getEdgeVertPos( int i ){
        Vec3i b = edgeVertBonds[i].verts;
        const Quat4f& pa = points[b.x];
        const Quat4f& pb = points[b.y];
        const Quat4f& pc = points[b.z];
        float c = edgeVertBonds[i].c;
        float mc = 1-c;
        return pa.f*mc + pb.f*c;
    }

    //void evalEdgeVert( Vec3i b, float c, float K ){
    inline void evalEdgeVert( int i ){
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

    // ============ virtual Functions
    
    inline virtual int pick_nearest(Vec3d ray0, Vec3d hray, int& ipick, int mask, double Rmax ) override {
        if     (mask==1){ ipick=pick_point_brute((Vec3f)ray0,(Vec3f)hray,Rmax); return 1; }
        else if(mask==2){ ipick=pick_bond_brute ( ray0, hray, Rmax );           return 2; }
        return -1;
    };
    
    inline virtual int pick_all(Vec3d ray0, Vec3d hray, int* out, int mask, double Rmax ) override { return 0; };
    
    inline virtual void* getPickedObject(int picked, int mask) override { 
        if     (mask==1){ return (void*)&points[picked]; }
        else if(mask==2){ return (void*)&bonds [picked]; }
        return 0; 
    };


    // ============ Functions

    float springForce(float l, float& f, Quat4f par );
    Quat4f springForce(Vec3f d, Quat4f par );
    void recalloc(int nPoint_, int nNeighMax_, int nBonds_=0);
    void recallocBBs(int n, bool bPoint=true, bool bEdge=true, bool bFace=true, bool bClean=true );
    void edgesToBBs();
    void printBBs();
    int pick_point_brute(const Vec3f& ray0, const Vec3f& hray, float Rmax, int i0=-1, int i1=-1 ) const;
    int pick_bond_brute( const Vec3d& ray0, const Vec3d& hRay, double Rmax, int i0=-1, int i1=-1 ) const;
    int pick_BBox(const Vec3d& ray0, const Vec3d& hRay, double tmax, int i0=-1, int i1=-1 );
    void updatePD_RHS(const Vec3f* pnew, Quat4f* bvec );
    void dotPD(const Vec3f* p, Vec3f* f );
    void updateJacobi_lin(Quat4f* ps_in, Quat4f* ps_out, Quat4f* bvec, Quat4f* rs=0 );
    void make_PD_Matrix(float* A, float dt );
    void make_PDmat_sparse(SparseMatrix<float>& A, float dt, bool bRealloc );
    void rhs_ProjectiveDynamics(const Quat4f* pnew, Quat4f* b);
    Quat4f rhs_ProjectiveDynamics_i(int i, const Quat4f* pnew);
    void rhs_ProjectiveDynamics_(const Quat4f* pnew, Quat4f* b);
    void realloc_LinearSystem(bool bCG=true, bool bCholesky=true, int nNeighMaxLDLT_=32, bool bDens=true );
    void prepare_LinearSystem(bool bRealloc=true, bool bCG=true, bool bCholesky=true, int nNeighMaxLDLT_=32, bool bDens=true );
    void run_LinSolve(int niter);
    void run_Cholesky_omp_simd(int niter);
    void updateInveriants(bool bPrint=false);
    void evalTrussForce_neighs(int iG);
    void evalTrussForce_neighs2(int iG);
    void evalTrussForces_neighs2();
    void evalTrussForces_neighs();
    void evalTrussForces_bonds();
    void evalEdgeVerts();
    void evalTrussCollisionImpulses_bonds(float rate=1.0 );
    float evalBondTension();
    void applyForceRotatingFrame(Vec3f p0, Vec3f ax, float omega );
    void applyForceCentrifug(Vec3f p0, Vec3f ax, float omega );
    void printNeighs(int i);
    void printAllNeighs();
    float getFmax();
    void setFixPoints(int n, int* fixPoints, double Kfix=1e12, bool bRealloc=true );
    void addAngularVelocity(Vec3f p0, Vec3f ax );
    float move_GD(float dt);
    float move_i_MD(int i, float dt, float cdamp );
    float move_MD(float dt, float damp=0.0f );
    int run(int niter, float dt, float damp  );
    void setOpt(float dt_, float damp_ );
    void FIRE_update(float& vf, float& vv, float& ff, float& cv, float& cf );
    int run_omp(int niter_max, bool bDynamic, float dt_, float damp_ );

};   // TrussDynamics_f

#endif
