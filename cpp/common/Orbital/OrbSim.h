
#ifndef  OrbSim_h
#define  OrbSim_h

#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
//#include "datatypes.h"  

#include "raytrace.h"
#include "geom3D.h"
#include "Interfaces.h"

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

class OrbSim_f : public Picker { public:
    int nPoint=0, nNeighMax=0, nNeighTot=0;
    // cpu buffers
    Quat4f* points=0;  // position and mass
    Quat4f* forces=0;  // force and energy
    Quat4f* vel   =0;  // velocity

    Quat4f* params=0;  // neighbor parameters (l0,kP,kT,damp)
    int*    neighs=0;  // neighbor indices
    int2*   neighBs=0; // neighbor bond indices
    //int*    neighBs=0; // neighbor bond indices

    int     nBonds =0; // number of bonds
    Quat4f* bparams=0; // bond parameters (l0,kP,kT,damp)
    int2*   bonds  =0; // indices of bonded points (i,j)
    float*  strain =0; // strain
    //float*  l0s    =0; // 
    Vec2f*  maxStrain=0;

    // Rotating frame
    Vec3f p0{0.,0.,0.};
    Vec3f ax{0.0,0.0,1.0};
    float omega = 0.05;

    float dt      = 2e-3;
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
    double ff_safety    = 1e-32;


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


    // =================== Picking

    int pick_point_brute( const Vec3f& ray0, const Vec3f& hray, float Rmax ){
        float r2min =  Rmax*Rmax;
        int imin    = -1;
        for(int i=0; i<nPoint; i++){
            float t;
            float r2 = rayPointDistance2( ray0, hray, points[i].f, t );
            //printf( "pick_point_brute ipick %i r %g p(%g,%g,%g)\n", i, sqrt(r2), points[i].f.x,points[i].f.y,points[i].f.z );
            if(r2<r2min){ imin=i; r2min=r2; }
            //double ti = raySphere( ray0, hRay, R, ps[i] );
            //if(ti<tmin){ imin=i; tmin=ti; }
        }
        //printf( "pick_point_brute ipick %i r2min %g \n", ipick, r2min );
        return imin;
    }

    int pick_bond_brute( const Vec3d& ray0, const Vec3d& hRay, double Rmax ) const {
        double dist_min =  Rmax;
        int    imin     = -1;
        for(int ib=0; ib<nBonds; ib++){
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

    // =================== Truss Simulation

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

    inline void evalTrussForce_neighs2(int iG){
        //const int iG = get_global_id(0);
        Quat4f p = points[iG];
        Quat4f f =Quat4f{0.0f,0.0f,0.0f,0.0f};
        //printf( "--- p[%i] \n", iG );
        //#pragma omp simd
        for(int ij=0; ij<nNeighMax; ij++){
            int j  = nNeighMax*iG + ij;
            int2 b = neighBs[j];
            if(b.x == -1) break;
            Quat4f par = bparams[b.y];
            //f.add( springForce( points[ja].f - p.f, params[j] ) );
            Vec3f d =  points[b.x].f - p.f;
            float li = d.norm();
            float k = 1e+6;
            f.f.add_mul( d, (k*(li-par.x)/li) );
            //printf( "p[%i,ij=%i,j=%i] li=%7.3f dl=%8.5e fi=%8.5e e=%8.5e par(%7.3f,%8.5e,%8.5e,%8.5e) \n", iG,ij,ja, li, li-params[j].x, fi,ei, params[j].x,params[j].y,params[j].z,params[j].w );
        }
        forces[iG] = f; // we may need to do += in future
    }

    void evalTrussForces_neighs2(){
        for(int iG=0; iG<nPoint; iG++){
            evalTrussForce_neighs2(iG);
        } 
    }

    void evalTrussForces_bonds(){
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
    }

    void evalBondTension(){
        for(int i=0;i<nBonds; i++ ){
            int2  b  = bonds[i];
            float l0 = bparams[i].x;
            //float l0 = l0s[i];
            float l   = (points[b.y]-points[b.x]).norm();
            float s   = (l-l0)/l0;
            //if( fabs(s)>0.5 ){ printf( "evalBondTension[%i] strain=%g l=%g l0=%g\n", i, s, l, l0 ); }
            strain[i] = s;
            // ToDo: break the bond if strain > maxStrain;
        }
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
        double omega2 = omega*omega;
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

    double getFmax(){ 
        double fmax=0;
        for(int i=0; i<nPoint; i++){ float f=forces[i].norm(); fmax=fmax>f?fmax:f; }   
        //printf( "|fmax|=%g\n", fmax );
        return fmax;
    }

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

    inline void move_i_MD(int i, float dt, float cdamp ){
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

    double move_MD(float dt, float damp=0.0f ){
        float cdamp = 1.0f - damp;
        double ff = 0.0; 
        for(int i=0;i<nPoint; i++ ){
            Quat4f p = points[i];
            Quat4f f = forces[i];
            Quat4f v = vel   [i];
            ff += f.f.norm2();
            v.f.mul( cdamp );
            v.f.add_mul( f.f, dt/p.w );
            p.f.add_mul( v.f, dt     );
            //printf( "move_GD[%i] |d|=%g |f|=%g dt/m=%g m=%g \n", i, f.f.norm() * dt/p.w, f.f.norm(), dt/p.w, p.w );
            vel   [i]=v;
            points[i]=p;
        }
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
        applyForceRotatingFrame( p0, ax, omega );
        //move_GD( 0.00001 );
        //move_MD( 1e-3, 1e-5 );
        //move_GD( 1e-7 );
        f2 = move_MD( dt, damp );
        printf( "OrbSim_f::run[%i] |F|=%g\n", itr, sqrt(f2) );
    }
    return niter;
}

inline void setOpt( double dt_, double damp_ ){
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
    float E,F2,ff,vv,vf, cf,cv;
    //long T0 = getCPUticks();
    int itr=0,niter=niter_max;
    #pragma omp parallel shared(niter,itr,cdamp)
    while(itr<niter){
        if(itr<niter){
        //#pragma omp barrier
        #pragma omp single
        {E=0;F2=0;ff=0;vv=0;vf=0;}
        //------ eval forces
        //#pragma omp barrier
        //#pragma omp for reduction(+:E)
        #pragma omp for 
        for(int iG=0; iG<nPoint; iG++){ 
            forces[iG] = Quat4fZero;
            evalTrussForce_neighs2(iG);
            if(bDynamic){ applyForceRotatingFrame_i( iG, p0, ax, omega ); }
            else        { applyForceCentrifug_i    ( iG, p0, ax, omega ); }
        }
        // ---- assemble (we need to wait when all atoms are evaluated)
        //#pragma omp barrier
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
        #pragma omp for
        for(int i=0;i<nPoint; i++ ){
            if(bDynamic){ 
                move_i_MD( i, dt, cdamp );
            }else{
                vel[i].f = vel[i].f*cv  + forces[i].f*cf;  // FIRE
                move_i_MD( i, dt, 1.0 );
            }
        }
        //#pragma omp barrier
        #pragma omp single
        { 
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

};   // OrbSim_f

#endif
