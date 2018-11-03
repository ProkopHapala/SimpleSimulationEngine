
#ifndef libFlight_h
#define libFlight_h

#include  <vector>

#include "arrayAlgs.h"

#include "AeroSurf.h"
#include "AeroCraft.h"
#include "geom3D.h"

#include "ShooterCommon.h"
#include "Projectile3D.h"

struct ControlState{
    double x;
    double dxdtmax, xmax, xmin;

    void setLimits( double dxdtmax_, double xmax_, double xmin_ ){
        dxdtmax=dxdtmax_; xmax=xmax_; xmin=xmin_;
    }

    double setDiff( double dxdt, double dt ){
        dxdt  = _clamp( dxdt, -dxdtmax, dxdtmax );
        x     = _clamp( x+dxdt*dt,  xmax, xmin );
        return x;
    }

};

enum { iPitch=0, iYaw=1, iRoll=2, iThrottle=3, iTrigger=4 };

class FlightWorld{ public:

    double time    = 0;
    double dt      = 0.01;
    Vec3d gravity  = (Vec3d){0.0,-9.81,0.0};
    Vec3d vwind    = (Vec3d){0.0, 0.00,0.0};
    AeroCraft craft1;
    AeroCraft craft_bak;

    int       nTargets=0;
    Sphere3d* targets =0;

    constexpr static int  nControls = 4;
    ControlState  controlsState[nControls];

    double         airDensity    = 1.22;
    int            nShotPerBurst = 16;
    double         shotLifeTime  = 10.0; // [s]
    double         maxBurstTime  = 0.5;  // [s]
    ProjectileType shotType;
    std::vector<Burst3d> bursts;

    Vec2d * targetHits = 0;

    int burstTargetCollision( int iTarget, Burst3d& burst, double dt ){
        Sphere3d& obj = targets[iTarget];
        int nhit=0;
        if( sq( obj.r)> burst.bbox.dist2_Cilinder( obj.p ) ){
            //printf( " bboxHit %i target %i \n", burst.id, iTarget );
            for( int j=0; j<burst.shots.size(); j++ ){
            //for( Particle3d& p : burst.shots ){
                Particle3d& p  = burst.shots[j];
                if(p.age > shotLifeTime) continue;
                //p.getOldPos(dt,hdir);
                Vec3d hdir = p.vel*-dt; // toward old
                double l   = hdir.normalize();
                if( sq(obj.r) > linePointDistance2( p.pos, hdir, obj.p, l ) ){
                    Vec2d& th = targetHits[iTarget];
                    th.x = fmin(p.age,th.x);
                    th.y = fmax(p.age,th.y);
                    //printf( "hit target %i box %i prj %i %g (%g,%g)\n", iTarget, burst.id, j, p.age, th.x, th.y );
                    burst.hit(j);
                    nhit++;
                }
            }
        }
        return nhit;
    }

    void shoot( const AeroCraft& craft1 ){
        //printf( "shoot \n" );
        double vMuzzle = 800.0; // [m/s]
        if(
            (bursts.size()==0)
            || ( bursts.back().shots.size() > nShotPerBurst )
            || ( bursts.back().age > maxBurstTime )
        ){
            //printf("New Burst %i\n", bursts.size() );
            bursts.push_back( Burst3d( &shotType, bursts.size() ) );
        }
        bursts.back().addShot( craft1.pos, craft1.vel + craft1.rotMat.c*vMuzzle );
    }

    void clearProjectiles(){
        /*
        double ageMax = 0.0;
        for(int ib=0; ib<bursts.size(); ib++){
            Burst3d& b = bursts[ib];
            printf( "burst %i np %i \n", ib, b.shots.size() );
            for( Particle3d& p : b.shots ){
                ageMax = fmax(ageMax, p.age);
            }
        }
        printf(  "ageMax %g \n", ageMax );
        */
        for(Burst3d& b: bursts){
            int n = prune( b.shots.size(), &b.shots[0], [&](const Particle3d& p){ return p.age > shotLifeTime; } );
            /*
            if( (ageMax>shotLifeTime) && (b.shots.size()>0)  ){
                printf( "nlf %i nsz %i ", n, b.shots.size() );
                for(int i=0; i<b.shots.size(); i++){
                    if(n==i)printf( "|" );
                    printf( " %5.5e", b.shots[i].age );
                }
                printf("\n");
            }
            */
            if( n<b.shots.size() ) b.shots.resize(n);
        }
        int nbleft = prune( bursts.size(), &bursts[0], [&](const Burst3d& b){ return b.shots.size()<=0; } );
        if( nbleft<bursts.size()){
            //printf( "nbleft %i bursts.size %i \n", nbleft, bursts.size() );
            bursts.resize(nbleft);
        }
    }

    void setWing( int i, double angle ){
        // TODO: we should use better system for rotation
        craft1.panels[i].lrot = craft_bak.panels[i].lrot;
        craft1.panels[i].lrot.rotate(angle,craft1.panels[i].lrot.a);
    }

    /*
    void setControlDiff( int i, double dangle, double dt, double vmax, double amin, double amax ){
        // TODO: we should use better system for rotation
        craft1.panels[i].lrot = craft_bak.panels[i].lrot;
        craft1.panels[i].lrot.rotate(angle,craft1.panels[i].lrot.a);
    }
    */

    void fly( int n, int nsub, double dt, Vec3d* pos, Vec3d* vel, Mat3d* rot ){
        //printf( "n %i nsub %i ntot \n", n, nsub, n*nsub );
        for(int i=0; i<n; i++){
            //printf( "fly %i\n", i );
            for(int j=0; j<nsub; j++){
                //printf( "fly i %i j %i \n ", i, j );
                craft1.clean_temp();
                craft1.applyAeroForces( vwind );
                craft1.force.add( gravity*craft1.mass );
                craft1.move(dt);
                //printf( "fly i %i j %i pos (%g,%g,%g)  \n ", i, j, craft1.pos.x, craft1.pos.y, craft1.pos.z );
                time+=dt;
            }
            pos[i] = craft1.pos;
            vel[i] = craft1.vel;
            rot[i] = craft1.rotMat;
        }
    }

    void flyAndShootTargets( int nsub, double dt_, double* controlBuff, double* stateBuff, double* targetHits_ ){
        dt=dt_;
        double ndt = nsub * dt;
        setWing( 2,   controlsState[iPitch].setDiff( controlBuff[iPitch], ndt ) );
        setWing( 3,   controlsState[iYaw  ].setDiff( controlBuff[iYaw  ], ndt ) );
        double roll = controlsState[iRoll ].setDiff( controlBuff[iRoll ], ndt );
        setWing( 0, +roll );
        setWing( 1, -roll );

        double enginePowerMax = 1e+6; // W
        craft1.propelers[0].power = enginePowerMax * controlsState[iThrottle].setDiff( controlBuff[iThrottle], ndt );
        if( controlBuff[iTrigger] > 0.5 ){ shoot( craft1 ); };

        targetHits = (Vec2d*)targetHits_;
        for( int i=0; i<nTargets; i++ ){
            targetHits[i].x = +1e+300;
            targetHits[i].y = -1e+300;
        }

        for(int j=0; j<nsub; j++){
            for(Burst3d& b: bursts){
                //if( b.age > shotLifeTime ) continue;
                b.move( dt, gravity, airDensity );
                for(int i=0; i<nTargets; i++){
                    burstTargetCollision( i, b, dt );
                }
                for(Particle3d& p : b.shots){ if(p.pos.y<0.0) p.age=1e+300; };
            }
            clearProjectiles();

            craft1.clean_temp();
            craft1.applyAeroForces( Vec3dZero );
            craft1.force.add( gravity*craft1.mass );
            craft1.move(dt);

            time+=dt;
        }
        Vec3d * v3stateBuff = (Vec3d*)stateBuff;
        v3stateBuff[0] = craft1.pos;      // 0
        v3stateBuff[1] = craft1.vel;      // 3
        v3stateBuff[2] = craft1.rotMat.c; // 6
        v3stateBuff[3] = craft1.rotMat.b; // 9
        v3stateBuff[4] = craft1.omega;    // 12
        v3stateBuff[5] = craft1.force;    // 15
        v3stateBuff[6] = craft1.torq;     // 18
        for(int i=0; i<nControls; i++){ stateBuff[21+i] = controlsState[i].x; } // 22
        // enum { iPitch=0, iYaw=1, iRoll=2, iThrottle=3, iTrigger=4 };
    }

};


#endif
