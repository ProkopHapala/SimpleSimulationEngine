
#ifndef libFlight_h
#define libFlight_h

#include  <vector>

#include "arrayAlgs.h"

#include "AeroSurf.h"
#include "AeroCraft.h"
#include "geom3D.h"

#include "ShooterCommon.h"
#include "Projectile3D.h"

class FlightWorld{ public:

    double time    = 0;
    Vec3d gravity  = (Vec3d){0.0,-9.81,0.0};
    Vec3d vwind    = (Vec3d){0.0, 0.00,0.0};
    AeroCraft craft1;
    AeroCraft craft_bak;

    int       nTargets=0;
    Sphere3d* targets =0;

    double         airDensity    = 1.22;
    int            nShotPerBurst = 16;
    double         shotLifeTime  = 10.0; // [s]
    ProjectileType shotType;
    std::vector<Burst3d> bursts;

    int burstTargetCollision( Sphere3d& obj, Burst3d& burst, double dt ){
        int nhit=0;
        if( sq( obj.r)> burst.bbox.dist2_Cilinder( obj.p ) ){
            for( int j=0; j<burst.shots.size(); j++ ){
            //for( Particle3d& p : burst.shots ){
                Particle3d& p  = burst.shots[j];
                if(p.age > shotLifeTime  ) continue;
                //p.getOldPos(dt,hdir);
                Vec3d hdir = p.vel*-dt; // toward old
                double l   = hdir.normalize();
                if( sq(obj.r) > linePointDistance2( p.pos, hdir, obj.p, l ) ){
                    burst.hit(j);
                    nhit++;
                }
            }
        }
        return nhit;
    }

    void shoot( const AeroCraft& craft1 ){
        double vMuzzle = 800.0; // [m/s]
        if( (bursts.size()==0) || ( bursts.back().shots.size() > nShotPerBurst ) ){
            bursts.push_back( Burst3d( &shotType, bursts.size() ) );
        }
        bursts.back().addShot( craft1.pos, craft1.vel + craft1.rotMat.c*vMuzzle );
    }

    void clearProjectiles(){
        auto pruneShots = [&](Burst3d& b) {
            auto pruneShots = [&](Particle3d& p){ return p.age > shotLifeTime; };
            int n = prune( b.shots.size(), &b.shots[0], pruneShots );
            b.shots.resize(n);
            return b.shots.size()==0;
        };
        //std::vector<Particle3d> particleTemp;
        int n = prune( bursts.size(), &bursts[0], pruneShots );
        bursts.resize(n);
    }

    void setWing( int i, double angle ){
        // TODO: we should use better system for rotation
        craft1.panels[i].lrot = craft_bak.panels[i].lrot;
        craft1.panels[i].lrot.rotate(angle,craft1.panels[i].lrot.a);
    }


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

    void flyAndShootTargets( int nsub, double dt, double* controlBuff, double* stateBuff, int* targetShot ){
        enum { iPitch=0, iYaw=1, iRoll=2, iThrottle=3, iTrigger=4 };
        setWing( 0,  controlBuff[iPitch] );
        setWing( 1,  controlBuff[iYaw  ] );
        setWing( 2, +controlBuff[iRoll ] );
        setWing( 3, -controlBuff[iRoll ] );
        craft1.propelers[0].power = controlBuff[iThrottle];
        if( controlBuff[iTrigger] > 0.5 ){ shoot( craft1 ); };
        for(int j=0; j<nsub; j++){
            for(Burst3d& b: bursts){
                if( b.age > shotLifeTime ) continue;
                b.move( dt, gravity, airDensity );
                for(int i=0; i<nTargets; i++){
                    targetShot += burstTargetCollision( targets[i], b, dt );
                }
            }
            craft1.clean_temp();
            craft1.applyAeroForces( Vec3dZero );
            craft1.force.add( {0.0,-9.81,0.0} );
            craft1.move(dt);

            time+=dt;
        }
        Vec3d * v3stateBuff = (Vec3d*)stateBuff;
        v3stateBuff[0] = craft1.pos;
        v3stateBuff[1] = craft1.vel;
        v3stateBuff[2] = craft1.rotMat.c;
        v3stateBuff[3] = craft1.rotMat.b;
    }

};


#endif
