
#ifndef  TrussDynamics_d_h
#define  TrussDynamics_d_h

#include "datatypes.h"
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

//#include "Cholesky.h"
#include "CG.h"
#include "arrayAlgs.h"

#include "SparseMatrix.h"
#include "SparseMatrix2.h"
#include <unordered_set>

// ========== Free Functions ==========

inline double springForce( double l, double& f, Quat4d par ){
    double dl = l - par.x;
    double k;
    if( dl > 0.0f ){
        k = -par.z;
    } else {
        k = par.y;
    }
    //Quat4d fe; 
    //fe.f = d*(k*dl/l);
    //fe.e = k*dl*dl;
    f = k*dl;
    return k*dl*dl;
}

inline Quat4d springForce( Vec3d d, Quat4d par ){
    double l  = d.norm();
    double dl = l - par.x;
    double k;
    if( dl > 0.0f ){
        k = -par.z;
    } else {
        k = par.y;
    }
    Quat4d fe; 
    fe.f = d*(k*dl/l);
    fe.e = k*dl*dl;
    return fe;
}

double checkDist(int n, const Vec3d* vec, const Vec3d* ref, int verb=1, double tol=1e-12 );
void fitAABB(Quat8d& bb, int n, int* c2o, Quat4d * ps );
void fitAABB_edge(Quat8d& bb, int n, int* c2o, int2* edges, Quat4d * ps );
void updatePointBBs(const Buckets& buckets, Quat8d* BBs, Quat4d* points, bool bInit=true);
void updateEdgeBBs(const Buckets& buckets, Quat8d* BBs, int2* edges, Quat4d* points, bool bInit=true);
void print_vector(int n, double * a, int pitch, int j0, int j1 );

// ========== Helper classes ==========

struct EdgeVertBond{ 
    Vec3i verts; 
    double c; 
    double K; 
    Vec3d f=Vec3dZero; 
};

struct SmartMixer{
    float b_end   = 0.75f;
    float b_start = 0.55f;
    int   istart = 3;
    int   iend   = 10; 

    // float get_bmix(int itr){
    //     if      ( itr < istart ) { return 0; }
    //     else if ( itr > iend   ) { return b_end; }
    //     else                     { return b_start + (b_end - b_start ) * (itr - istart) / (iend - istart); }
    // }

    inline float get_bmix(int itr){
        if ( itr < istart ) { return 0; }
        else                { return b_end; }
    }
};

// ========== main class ==========

class TrussDynamics_d: public Picker { public:
    double time=0;
    double Cdrag   = 0.0;
    //Vec3d Gravity = {0,-9.81,0};
    // Rotating frame
    Vec3d pos0{0.,0.,0.};
    Vec3d ax{0.0,0.0,1.0};
    //double omega = 0.05;
    //Quat4d accel{ 0.0,-9.81 , 0.0 , 0.0 };    // acceleration
    Quat4d accel{ 0.0, 0.0 , 0.0 , 0.0  };    // acceleration
    Quat4d rot0 { 0.0, 0.0 , 0.0 , 0.0  };    // center of rotation
    Quat4d omega{ 0.0, 0.0 , 0.0 , 0.05 };    // angular velocity

    //double dt      = 2e-3; //double kGlobal = 1e+6;
    //double dt      = 2e-3;    double kGlobal = 1e+7;
    double dt      = 0.5e-3;  double kGlobal = 1e+8;
    //double damping = 1e-4;
    double damping  = 0.0;
    int nSolverIters = 10;

    std::vector<int> damped_bonds;
    std::unordered_set<int> damped_points;


    int nPoint=0, nNeighMax=0, nNeighTot=0;
    // cpu buffers
    Quat4d* points=0;  // position and mass
    Quat4d* forces=0;  // force and energy
    Quat4d* vel   =0;  // velocity
    Quat4d* vel0 = 0;
    //Quat4d* impuls=0;  // accumulated impulse from corrector
    Quat4d* bvec  =0;  // right hand side of linear system Ap=b ( it has meaning of internal force due to projective dybamics matrix A = K + D + I ) 
    Vec3d * bvec0 =0;

    //Solver<float>* solver_f = 0;
    
    // Binding to external linear solver of system Ap=b ( where A is matrix of projective dynamics )
    Quat4f* extern_b = 0;
    Quat4f* extern_x = 0;
    void (*extern_solve)();   // function pointer extern_solve

    double* kFix=0;   // force constant for fixed points
    double  kLinRegularize = 1.0;

    Quat4d* params=0;  // neighbor parameters (l0,kP,kT,damp)
    int*    neighs=0;  // neighbor indices
    int2*   neighBs=0; // neighbor bond indices
    //int*    neighBs=0; // neighbor bond indices
    int*    neighB2s=0;  // neighbor indices

    // Cholesky / Projective Dynamics 
    double* PDmat=0;
    double* LDLT_L=0;
    double* LDLT_D=0; 
    int* neighsLDLT=0;
    int  nNeighMaxLDLT=0;
    SparseMatrix2<double>  Lsparse; 
    SparseMatrix2<double>  LsparseT; 
    SparseMatrix<double> PDsparse;
    
    SmartMixer mixer;

    // variables realted to Conjugate Gradient
    CGsolver cgSolver;
    double cg_tol             = 1e-3;
    int    cgSolver_niterdone = 0;
    double time_LinSolver     = 0;
    double time_cg_dot        = 0;

    // choice of linear solver method
    int linSolveMethod = 2; // 0=CG, 1=CGsparse, 2=Cholesky
    enum class LinSolveMethod{ CG=0, CGsparse=1, Cholesky=2, CholeskySparse=3, Jacobi=4, GaussSeidel=5, JacobiMomentum=6, GSMomentum=7, JacobiFlyMomentum=8, GSFlyMomentum=9, Force=10, JacobiDiff=11, MomentumDiff=12, ExternDiff=13 };
    bool   bApplyResudualForce = true;
    double residualForceFactor = 1.0;

    Vec3d*  ps_cor      =0; // new Vec3d[nPoint];
    Vec3d*  ps_pred     =0; // new Vec3d[nPoint]; 
    Vec3d*  linsolve_b  =0; // new Vec3d[nPoint];
    Vec3d*  linsolve_yy =0; // new Vec3d[nPoint];
    Vec3d*  ps_0        =0;

    int     nBonds  =0;    // number of bonds
    Quat4d* bparams =0;    // bond parameters (l0,kPress,kTens,damp)
    int2*   bonds   =0;    // indices of bonded points (i,j)
    double*  strain =0;    // strain
    //double*  l0s  =0;    // 
    Vec2d*  maxStrain=0;

    // ====== Linearized Truss

    Quat4f* hbs  =0;  // normalized direction of the stick, and initial distortion of the length from neutral length
    Quat4f* dpos =0;  // distortion of point from neutral postion
    Quat4f* fdpos=0;  // force on distortion of point from neutral postion
    Quat4f* vdpos=0;  // velocity of distortion of point from neutral postion

    Quat4d* points_bak = 0; // backup of points for position based dynamics

    float* dls =0;  // distortion of point from neutral postion
    Vec3d* kDirs=0;  // d/|d| where d = p1-p0 normalized direction of the stick, and initial distortion of the length from neutral length
    //double* f0s  =0;  // (|d|-l0) *k  force on distortion of point from neutral postion


    // ====== Invairiants

    double mass = 0;
    Vec3d cog  = Vec3dZero;
    Vec3d vcog = Vec3dZero;
    Mat3d I    = Mat3dZero;
    Vec3d L    = Vec3dZero;
    Vec3d torq = Vec3dZero;

    double F_residual = 0.0;

    // ====== Collision
    // ToDo: this should be moved to a separate class ?
    int       nBBs=0;
    Quat8d*  BBs=0; // bounding boxes (can be either AABB, or cylinder, capsula) 
    Buckets  pointBBs;    // buckets for collision detection
    Buckets  edgeBBs;
    Buckets  faceBBs;
    Buckets  pointChunks;  // chunks for parallelization, these points are copied to local memory when solving one edgeBBs chunk of bonds  

    // Faces are used just for ray-tracing, collision detection, etc
    int   nFaces=0;
    int4* faces=0; // indices of points tringles or quads, the last index is -1 is the face is a triangle, ot it can be used also to store face type if just triangles are used

    int nEdgeVertBonds=0;
    EdgeVertBond* edgeVertBonds=0; // indices of bonded points (i,j)

    // callback function pointer what to do in between iterations
    void (*user_update)(double dt);

    Vec3d hit_pos, hit_normal;

    //double maxAcc = 1e+6;
    double maxAcc = 1.0;
    double collision_damping = 0.002;
    //double collision_damping = 1.0;
    //double collision_damping = 1.1;   // if collision_damping > 1.0 then it is like successive over-relaxation (SOR) method ? https://en.wikipedia.org/wiki/Successive_over-relaxation

    int    lastNeg = 0;
    // FIRE
    int    minLastNeg   = 5;
    double finc         = 1.1;
    double fdec         = 0.5;
    double falpha       = 0.98;
    double dt_max       = dt;
    double dt_min       = 0.1 * dt;
    double damp_max     = damping;
    double ff_safety    = 1e-16;
    double cv,cf;

    // ============ inline  Functions

    inline void set_time_step( float dt_ ){
        dt      = dt_;
        //inv_dt2 = 1.0f / (dt * dt);
    }

    inline Vec3d getPointForce( int i ){
        return (vel[i].f*Cdrag) + (accel.f*points[i].w);
    }

    inline void reallocFixed(){ _realloc0( kFix, nPoint, 0.0 ); }
    inline void cleanForce(){ for (int i=0; i<nPoint; i++){ forces[i]=Quat4dZero; } };
    inline void cleanVel  ( Quat4d v0=Quat4dZero ){ for (int i=0; i<nPoint; i++){ vel[i]=v0; } };

    inline double getLinearBondStiffness( int ib ){
        // (l0,kPress,kTens,damp)
        //return  bparams[ib].y;  //k = kPress;
        return  bparams[ib].z;  //k = kTens;  
    }

    inline void applyForceRotatingFrame_i( int i, Vec3d p0, Vec3d ax, double omega ){
        const Quat4d& p = points[i];
        const Quat4d& v = vel   [i];
        Vec3d d,f;
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

    inline void applyForceCentrifug_i( int i, Vec3d p0, Vec3d ax, double omega ){
        const Quat4d& p = points[i];
        Vec3d d;
        d.set_sub(p.f,p0);
        d.makeOrthoU(ax);   
        forces[i].f.add_mul(d, p.w*omega*omega );
    }

    inline Vec3d getEdgeVertPos( int i ){
        Vec3i b = edgeVertBonds[i].verts;
        const Quat4d& pa = points[b.x];
        const Quat4d& pb = points[b.y];
        const Quat4d& pc = points[b.z];
        double c = edgeVertBonds[i].c;
        double mc = 1-c;
        return pa.f*mc + pb.f*c;
    }

    //void evalEdgeVert( Vec3i b, double c, double K ){
    inline void evalEdgeVert( int i ){
        Vec3i  b = edgeVertBonds[i].verts; // vert indexes (edge.a,edge.b,ivert);
        double c = edgeVertBonds[i].c;     //  interpolation parameter
        double K = edgeVertBonds[i].K;     //  stiffness constant 
        // ToDo: perhaps we should interpolate it by B-spline to make the path more smooth
        // ToDo: Apply Force in the direction of the edge, constrain perpendicular to the edge
        // ToDo: Damping ( collision damping )
        
        // fit vert to edge
        const Quat4d& pa = points[b.x];
        const Quat4d& pb = points[b.y];
        const Quat4d& pc = points[b.z];
        double mc = 1-c;
        Vec3d d = pc.f - (pa.f*mc + pb.f*c);
        //glColor3f(1.0,0.0,0.0); Draw3D::drawVecInPos( d, pc.f*1.0 );
        
        double invL = 1./d.norm();
        double dv   = d.dot( vel[b.z].f - vel[b.x].f*mc - vel[b.y].f*c )*invL;
        double mab  = pa.w*mc + pb.w*c;
        double imp  = collision_damping * pc.w*mab*dv/( pc.w  + mab );
        //imp/=dt; 

        d.mul( K + imp*invL ); 
        edgeVertBonds[i].f = d;
        //printf( "evalEdgeVert[%i,%i,%i] d(%g,%g,%g) c=%g \n", b.x,b.y,b.z, d.x,d.y,d.z, c );   
        forces[b.x].f.add_mul(d,mc);  // edge vert 1
        forces[b.y].f.add_mul(d, c);  // edge vert 2
        forces[b.z].f.sub(d);         // vert

        // Force
        // const Quat4d& pa = points[b.x];
        // const Quat4d& pc = points[b.z];
        // Vec3d d = pc.f - pa.f;          
        // // Damping
        // //double dv   = d.dot( vel[b.x].f + vel[b.y].f - vel[b.z].f );
        // d.mul( K );
        // forces[b.x].f.add(d);
        // //forces[b.y].f.add_mul(d, c);
        // forces[b.z].f.sub(d);
    }

    // ============ virtual Functions

    virtual int pick_nearest(Vec3d ray0, Vec3d hray, int& ipick, int mask, double Rmax ) override;
    virtual int pick_all(Vec3d ray0, Vec3d hray, int* out, int mask, double Rmax ) override;
    //virtual void* getPickedObject(int picked, int mask) override;

    virtual inline void* getPickedObject(int picked, int mask) { 
        if     (mask==1){ return (void*)&points[picked]; }
        else if(mask==2){ return (void*)&bonds [picked]; }
        return 0; 
    };

    // ============ Functions
 
void recalloc(int nPoint_, int nNeighMax_, int nBonds_=0);
void realloc_LinearSystem(bool bCG=true, bool bCholesky=true, int nNeighMaxLDLT_=32, bool bDens=true );
void realloc_lin(bool bLinSolver=true);
void recallocBBs(int n, bool bPoint=true, bool bEdge=true, bool bFace=true, bool bClean=true );
double norm_butFixed(Vec3d* ps );
double norm_butFixed(Quat4d* ps );
void edgesToBBs();
void printBBs();
int pick_point_brute(const Vec3d& ray0, const Vec3d& hray, double Rmax, int i0=-1, int i1=-1 ) const;
int pick_bond_brute( const Vec3d& ray0, const Vec3d& hRay, double Rmax, int i0=-1, int i1=-1 ) const ;
int pick_BBox(const Vec3d& ray0, const Vec3d& hRay, double tmax, int i0=-1, int i1=-1 );
void make_PD_Matrix(double* A, double dt );
void make_PDmat_sparse(SparseMatrix<double>& A, double dt, bool bRealloc );
void rhs_ProjectiveDynamics(const Vec3d* pnew, Vec3d* b);
Vec3d rhs_ProjectiveDynamics_i(int i, const Vec3d* pnew);
void rhs_ProjectiveDynamics_(const Vec3d* pnew, Vec3d* b);
void updatePD_RHS(const Vec3d* pnew, Quat4d* bvec );
void dotPD(const Vec3d* p, Vec3d* f );
void updatePD_dRHS(const Vec3d* pnew, Quat4d* bvec );
void updateJacobi_lin(Vec3d* ps_in, Vec3d* ps_out, Quat4d* bvec, Vec3d* rs=0 );
void updateGaussSeidel_lin(Vec3d* ps, Quat4d* bvec );
void evalTrussForce(Vec3d* ps_in, Vec3d* force );
void updateJacobi_fly(Vec3d* ps_in, Vec3d* ps_out );
void updateGaussSeidel_fly(Vec3d* ps );
void updateIterativeMomentum(Vec3d* psa, Vec3d* psb );
void updateIterativeJacobi(Vec3d* psa, Vec3d* psb );
void updateIterativeJacobiDiff(Vec3d* psa, Vec3d* psb );
void updateIterativeMomentumDiff(Vec3d* psa, Vec3d* psb );
void updateIterativeExternDiff(Vec3d* psa, Vec3d* psb );
void prepare_LinearSystem(bool bRealloc=true, bool bCG=true, bool bCholesky=true, int nNeighMaxLDLT_=32, bool bDens=true );
void dampPoints(double Vdamping );
void linsolve( Vec3d* ps_pred, Vec3d* ps_cor );
void run_LinSolve(int niter);
void run_Cholesky_omp_simd(int niter);
void run_Cholesky_omp(int niter);
void prepareLinearizedTruss();
void evalTrussForcesLinearized();
float evalTrussForceLinearized_neighs2(int iG);
float evalTrussForcesLinearized_neighs2();
void prepareLinearizedTruss_ling(double* bvec, bool bAddForce=true );
void dot_Linearized_bonds(int n, double* x, double * Ax);
void dot_Linearized_neighs2(int n, double* x, double * Ax);
void findGuassSeidelPivotingPriority(double* priority );
void move_dpos(double dt=1.0 );
void apply_dpos(double sc=1.0 );
void update_velocity(double dt=1.0 );
double constr_jacobi_neighs2_absolute();
double constr_jacobi_neighs2_diff();
void run_constr_dynamics(int nitr, double dt );
void solveLinearizedConjugateGradient();
void updateInveriants(bool bPrint=false);
void evalTrussForce_neighs(int iG);
void evalTrussForce_neighs2(int iG);
void evalTrussForces_neighs2();
void evalTrussForces_neighs();
void evalTrussForces_bonds();
//Vec3d getEdgeVertPos(int i );
//void evalEdgeVert(int i );
void evalEdgeVerts();
void evalTrussCollisionImpulses_bonds(double rate=1.0 );
double evalBondTension();
//void applyForceRotatingFrame_i(int i, Vec3d p0, Vec3d ax, double omega );
//void applyForceCentrifug_i(int i, Vec3d p0, Vec3d ax, double omega );
void applyForceRotatingFrame(Vec3d p0, Vec3d ax, double omega );
void applyForceCentrifug(Vec3d p0, Vec3d ax, double omega );
void printNeighs(int i);
void printAllNeighs();
double getFmax();
void setFixPoints(int n, int* fixPoints, double Kfix=1e12, bool bRealloc=true );
// void reallocFixed();
// void cleanForce();
// void cleanVel(Quat4d v0=Quat4dZero );
void addAngularVelocity(Vec3d p0, Vec3d ax );
void addAngularVelocity2(Vec3d p0, Vec3d ax );
double move_GD(double dt);
double move_i_MD(int i, double dt, double cdamp );
double move_MD(double dt, double damp=0.0f );
int run(int niter, double dt, double damp  );
void setOpt(double dt_, double damp_ );
void FIRE_update(double& vf, double& vv, double& ff, double& cv, double& cf );
int run_omp(int niter_max, bool bDynamic, double dt_, double damp_ );

bool checkMasses(double mass_tolerance=1e-9, bool bExit=true);

};   // TrussDynamics_d

// --- standalone functions using TrussDynamics_d

void loadSimFromFile( const char* fname, TrussDynamics_d& sim );

#endif
