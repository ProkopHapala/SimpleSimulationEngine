
#ifndef  OrbSim_d_h
#define  OrbSim_d_h

#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
//#include "datatypes.h"  

#include "Buckets.h"

#include "raytrace.h"
#include "geom3D.h"
#include "Interfaces.h"

double springForce( double l, double& f, Quat4d par ){
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

Quat4d springForce( Vec3d d, Quat4d par ){
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
struct EdgeVertBond{ 
    Vec3i verts; 
    double c; 
    double K; 
    Vec3d f=Vec3dZero; 
};

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

inline void updatePointBBs(const Buckets& buckets, Quat8d* BBs, Quat4d* points, bool bInit=true){
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

inline void updateEdgeBBs(const Buckets& buckets, Quat8d* BBs, int2* edges, Quat4d* points, bool bInit=true){
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

class OrbSim: public Picker { public:
    double time=0;
    int nPoint=0, nNeighMax=0, nNeighTot=0;
    // cpu buffers
    Quat4d* points=0;  // position and mass
    Quat4d* forces=0;  // force and energy
    Quat4d* vel   =0;  // velocity

    Quat4d* params=0;  // neighbor parameters (l0,kP,kT,damp)
    int*    neighs=0;  // neighbor indices
    int2*   neighBs=0; // neighbor bond indices
    //int*    neighBs=0; // neighbor bond indices
    int*    neighB2s=0;  // neighbor indices

    int     nBonds =0; // number of bonds
    Quat4d* bparams=0; // bond parameters (l0,kP,kT,damp)
    int2*   bonds  =0; // indices of bonded points (i,j)
    double*  strain =0; // strain
    //double*  l0s    =0; // 
    Vec2d*  maxStrain=0;

    // ====== Linearized Truss

    Quat4f* hbs  =0;  // normalized direction of the stick, and initial distortion of the length from neutral length
    Quat4f* dpos =0;  // distortion of point from neutral postion
    Quat4f* fdpos=0;  // force on distortion of point from neutral postion
    Quat4f* vdpos=0;  // velocity of distortion of point from neutral postion

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

    // Rotating frame
    //Vec3d p0{0.,0.,0.};
    //Vec3d ax{0.0,0.0,1.0};
    //double omega = 0.05;
    Quat4d accel{0.0f,0.0f,0.0f,0.0f};    // acceleration
    Quat4d rot0 {0.0f,0.0f,0.0f,0.0f};    // center of rotation
    //Quat4d omega{0.0f,0.0f,1.0f,0.05f}; // angular velocity, (xyz=axisxyz,w=magnitude)
    Quat4d omega{0.0f,0.0f,1.0f,0.05f};


    Vec3d hit_pos, hit_normal;


    //double maxAcc = 1e+6;
    double maxAcc = 1.0;
    double collision_damping = 0.002;
    //double collision_damping = 1.0;
    //double collision_damping = 1.1;   // if collision_damping > 1.0 then it is like successive over-relaxation (SOR) method ? https://en.wikipedia.org/wiki/Successive_over-relaxation

    //double kGlobal = 1e+6;
    
    double dt      = 2e-3;    double kGlobal = 1e+7;
    //double dt      = 0.5e-3;  double kGlobal = 1e+8;

    //double damping = 1e-4;
    double damping  = 0.05;
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

    void recalloc( int nPoint_, int nNeighMax_, int nBonds_=0){
        nPoint = nPoint_; nNeighMax = nNeighMax_;
        nNeighTot = nPoint*nNeighMax;
        _realloc( points, nPoint    );
        _realloc( forces, nPoint    );
        _realloc( vel,    nPoint    );
        _realloc( params, nNeighTot );
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

    int pick_point_brute( const Vec3d& ray0, const Vec3d& hray, double Rmax, int i0=-1, int i1=-1 ){
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

    int pick_bond_brute( const Vec3d& ray0, const Vec3d& hRay, double Rmax, int i0=-1, int i1=-1 ) const {
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

    int pick_BBox( const Vec3d& ray0, const Vec3d& hRay, double tmax, int i0=-1, int i1=-1 ){
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

    virtual int pick_nearest(Vec3d ray0, Vec3d hray, int& ipick, int mask, double Rmax ) override {
        if     (mask==1){ ipick=pick_point_brute((Vec3d)ray0,(Vec3d)hray,Rmax); return 1; }
        else if(mask==2){ ipick=pick_bond_brute ( ray0, hray, Rmax );           return 2; }
        return -1;
    };
    
    virtual int pick_all(Vec3d ray0, Vec3d hray, int* out, int mask, double Rmax ) override { return 0; };
    
    virtual void* getPickedObject(int picked, int mask) override { 
        if     (mask==1){ return (void*)&points[picked]; }
        else if(mask==2){ return (void*)&bonds [picked]; }
        return 0; 
    };

    // =================== Linearized Elasticity Truss Simulation

    void prepareLinearizedTruss(){
        for(int ib=0; ib<nBonds; ib++){ 
            int2    b = bonds[ib];
            Vec3d   d = points[b.y].f-points[b.x].f;
            double  l = d.normalize();
            hbs [ib].f = (Vec3f)d;
            hbs [ib].e = l-bparams[ib].x; 
            dpos[ib]   = Quat4fZero; 
        }
    }

    void evalTrussForcesLinearized(){ 
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

    inline float evalTrussForceLinearized_neighs2(int iG){
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
    float evalTrussForcesLinearized_neighs2(){
        float F2=0;
        for(int iG=0; iG<nPoint; iG++){ F2+=evalTrussForceLinearized_neighs2(iG); } 
        return F2;
    }

    void solveLinearizedConjugateGradient(){
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

    void updateInveriants(bool bPrint=false){
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
            printf( "OrbSim::updateInveriants() mass %g cog(%g,%g,%g) vcog(%g,%g,%g)  L(%g,%g,%g) torq(%g,%g,%g) \n", mass, cog.x,cog.y,cog.z, vcog.x,vcog.y,vcog.z, L.x,L.y,L.z, torq.x,torq.y,torq.z );
            printf( "OrbSim::updateInveriants() I \n" ); I.print();
        }
    }

    void evalTrussForce_neighs(int iG){
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

    inline void evalTrussForce_neighs2(int iG){
        const double Adamp = collision_damping*0.5/dt;
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
            double k        = kGlobal;
            double fl       = k*(l-par.x);
            const double invL = 1/l;


            // const double dv  = d.f.dot( vel[b.x].f - v.f );
            // if(dv<0){ fl*=1-dv; }


            // // collision damping
            
            // const double dv  = d.f.dot( vel[b.x].f - v.f );
            // double imp = Adamp * p.w*d.w*dv/(p.w+d.w);
            // //double imp = 0.1* 0.5 * p.w*d.w*dv/(p.w+d.w);
            // //const double imp = 0;
            // //imp/=dt;
            double imp = 0;

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

    Vec3d getEdgeVertPos( int i ){
        Vec3i b = edgeVertBonds[i].verts;
        const Quat4d& pa = points[b.x];
        const Quat4d& pb = points[b.y];
        const Quat4d& pc = points[b.z];
        double c = edgeVertBonds[i].c;
        double mc = 1-c;
        return pa.f*mc + pb.f*c;
    }

    //void evalEdgeVert( Vec3i b, double c, double K ){
    void evalEdgeVert( int i ){
        Vec3i b = edgeVertBonds[i].verts;
        double c = edgeVertBonds[i].c;
        double K = edgeVertBonds[i].K;
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
        forces[b.x].f.add_mul(d,mc);
        forces[b.y].f.add_mul(d, c);
        forces[b.z].f.sub(d);

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

    void evalEdgeVerts(){
        //for(int i=0; i<nEdgeVertBonds; i++){ evalEdgeVert( edgeVertBonds[i].verts, edgeVertBonds[i].c, edgeVertBonds[i].K ); }
        for(int i=0; i<nEdgeVertBonds; i++){ evalEdgeVert( i ); }
    }

    void evalTrussCollisionImpulses_bonds( double rate=1.0 ){
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
            //if( i==146 ){ // Debug
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

    void evalBondTension(){
        for(int i=0;i<nBonds; i++ ){
            int2  b  = bonds[i];
            double l0 = bparams[i].x;
            //double l0 = l0s[i];
            double l   = (points[b.y].f-points[b.x].f).norm();
            double s   = (l-l0)/l0;
            //if( fabs(s)>0.5 ){ printf( "evalBondTension[%i] strain=%g l=%g l0=%g\n", i, s, l, l0 ); }
            //if(i==6272){ printf( "evalBondTension[%i](%i,%i) strain=%g l=%g l0=%g\n", i, b.x,b.y, s, l, l0 ); }
            strain[i] = s;
            // ToDo: break the bond if strain > maxStrain;
        }
        //exit(0);
    }

    void applyForceRotatingFrame_i( int i, Vec3d p0, Vec3d ax, double omega ){
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

    void applyForceCentrifug_i( int i, Vec3d p0, Vec3d ax, double omega ){
        const Quat4d& p = points[i];
        Vec3d d;
        d.set_sub(p.f,p0);
        d.makeOrthoU(ax);   
        forces[i].f.add_mul(d, p.w*omega*omega );
    }

    void applyForceRotatingFrame( Vec3d p0, Vec3d ax, double omega ){
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

    void applyForceCentrifug( Vec3d p0, Vec3d ax, double omega ){
        for(int i=0;i<nPoint; i++ ){ applyForceCentrifug_i( i, p0, ax,omega ); }
    }

    void printNeighs(int i){
        int j0 = i*nNeighMax;
        for(int jj=0;jj<nNeighMax;jj++){
            int j=j0+jj;
            int ing = neighs[j];
            if(ing<0) break;
            Quat4d par = params[j];
            printf( "ng[%i,%i|%i] l0,kP,kT,damp(%g,%g,%g,%g)\n", i, ing, jj, par.x,par.y,par.z,par.w );
        }
    }
    void printAllNeighs(){ printf("OrbSim::printAllNeighs(nPoint=%i,nNeighMax=%i)\n",nPoint,nNeighMax); for(int i=0;i<nPoint;i++){ printNeighs(i); }; };

    double getFmax(){ 
        double fmax=0;
        for(int i=0; i<nPoint; i++){ double f=forces[i].norm(); fmax=fmax>f?fmax:f; }   
        //printf( "|fmax|=%g\n", fmax );
        return fmax;
    }

    void cleanForce(){ for (int i=0; i<nPoint; i++){ forces[i]=Quat4dZero; } };
    void cleanVel  (){ for (int i=0; i<nPoint; i++){ vel   [i]=Quat4dZero; } };

    double move_GD(double dt){
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

    inline double move_i_MD(int i, double dt, double cdamp ){
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

    double move_MD(double dt, double damp=0.0f ){
        double cdamp = 1.0f - damp;
        double ff = 0.0; 
        for(int i=0;i<nPoint; i++ ){ move_i_MD(i, dt, cdamp ); }
        return ff;
    }

int run( int niter, double dt, double damp  ){
    double f2 = -1;
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

inline void setOpt( double dt_, double damp_ ){
    dt      = dt_max   = dt_;  dt_min=0.1*dt_max;
    damping = damp_max = damp_;
    cleanForce( );
    cleanVel  ( );
}

void FIRE_update( double& vf, double& vv, double& ff, double& cv, double& cf ){
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


int run_omp( int niter_max, bool bDynamic, double dt_, double damp_ ){
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
            //printf( "OrbSim::run_omp() itr %i/%i E=%g |F|=%g time= %g [ms]( %g [us/%i iter]) \n", itr,niter_max, E, F_residual, time*1e+3, time*1e+6/itr, itr );
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

};   // OrbSim

#endif
