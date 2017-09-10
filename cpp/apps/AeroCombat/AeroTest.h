
#ifndef AeroTester_h
#define AeroTester_h

#include "testUtils.h"

class AeroTester{ public:
    AeroCraftControler* autoPilot1 = NULL;
    AeroCraft*          myCraft    = NULL;
    double             gravityG    = 9.81;

    int ntrj=0,ntrjMax=0;
    Vec3d  *trjPos=NULL,*trjVel=NULL,*trjForce=NULL,*trjFw=NULL,*trjUp=NULL;
    double *trjT=NULL;

    void reallocateTrj(int n);
    void evalAircraftTrajectory( int n, int nsub, int msub, double dt );
    //void doStaticTesting();
    void doStaticTesting( int ntrj_, double vmin, double vmax, double dv );

    //AeroTester( AeroCraftControler* autoPilot1, AeroCraft* myCraft, int n ){}

};

void AeroTester::reallocateTrj(int n){
  ntrjMax=n;
  if(trjPos  ) delete trjPos;   trjPos   = new Vec3d[n];
  if(trjVel  ) delete trjVel;   trjVel   = new Vec3d[n];
  if(trjForce) delete trjForce; trjForce = new Vec3d[n];
  if(trjFw   ) delete trjFw;    trjFw    = new Vec3d[n];
  if(trjUp   ) delete trjUp;    trjUp    = new Vec3d[n];
  if(trjT    ) delete trjT;     trjT     = new double[n];
}

void AeroTester::evalAircraftTrajectory( int n, int nsub, int msub, double dt ){
    ntrj=n;
    double t=0;
    long ticks1 = getCPUticks();
    for(int i=0; i<n; i++){
        trjPos  [i] = myCraft->pos;
        trjVel  [i] = myCraft->vel;
        trjForce[i] = myCraft->force;
        trjFw   [i] = myCraft->rotMat.c;
        trjUp   [i] = myCraft->rotMat.b;
        trjT    [i] = t;
        for(int j=0; j<nsub; j++){
            //resetSteer();
            autoPilot1->control(dt*msub);
            for(int itr=0; itr<msub; itr++){
                myCraft->clean_temp();
                myCraft->force.set      ( { 0, gravityG*myCraft->mass, 0 } );
                myCraft->applyAeroForces( {0,0,0} );
                myCraft->move(dt);
                t+=dt;
            }
		}
    }
    double  dticks        = getCPUticks() - ticks1;
    double  ticksPerIter  = dticks/(n*nsub*msub);
    double  tickPerSec    = dticks/(t-trjT[0]);
    printf( "Ticks %g /iter %g /sec %g\n", dticks, ticksPerIter, tickPerSec );
}

void AeroTester::doStaticTesting( int ntrj_, double vmin, double vmax, double dv ){
    //--- propeller characterisic
    //double vmin=0.01;
    //double vmax=300.0;
    //double dv  =5.0;
    for(double v=vmin; v<vmax; v+=dv){
        double thrust = myCraft->propelers[0].getThrust(v);
        printf(" v=%f [m/s] thrust=%f [N] \n",  v, thrust );
    }
    //int ntrj_=500;
    reallocateTrj(ntrj_);
    evalAircraftTrajectory(ntrj_, 100, 10, 0.001);  // correspond to interactive simulation
    //evalAircraftTrajectory(ntrj_, 25, 4, 0.01);
    //exit(0);
}


#endif
