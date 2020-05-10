
#ifndef MechEuler2D_h
#define MechEuler2D_h

#include "fastmath.h"
#include "Vec2.h"
//#include "Vec3.h"





class MechEuler2D{ public:
    double cohesive_pressure = 1.0; // if pressure lower than cohesive pressure, material tend to stick toghether
    constexpr static double const_Rgas    = 8.31446261815324; // J/(K*mol)

    double step    = 1.0;
    double invStep = 1.0;
    Vec2i nc;
    int nctot=0;

    int nmat = 0;
    //double*  Ns_buff=0;
    double** matNs=0; // molar amount of different materials
    double*  Ntot=0;    // total molar amount
    double*  Ntmp=0;
    double*  Vs=0;
    double*  Ts=0;      // Temperature
    double*  ps=0;      // pressure
    Vec2d*   vs=0;      // drift velocities (or momentum?)

    double sumN;

    void realloc( Vec2i nc_, int nmat_ ){
        nc    = nc_;
        nctot = nc.x*nc.y;
        _realloc( Ntot, nctot );
        _realloc( Vs,   nctot );
        _realloc( Ts,   nctot );
        _realloc( ps,   nctot );
        _realloc( vs,   nctot );

        _realloc( Ntmp, nctot );

        nmat = nmat_;
        if(matNs){
            for(int imat=0; imat<nmat; imat++){ delete [] matNs[imat]; };
            delete [] matNs;
            matNs = 0;
        }
        matNs = new double*[nmat];
        for(int imat=0; imat<nmat; imat++){
            matNs[imat] = new double[nctot];
        }
    }

    void setBoundary( double val, double* arr ){
        for(int i=0; i<nc.x  ;i++){ arr[i]=val; };               arr+=nc.x;
        for(int i=1; i<nc.y-1;i++){ arr[0]=val; arr[nc.x-1]=val; arr+=nc.x; }
        for(int i=0; i<nc.x  ;i++){ arr[i]=val; };
    }

    void clearVelocitis(){ for(int i=0; i<nctot;i++){ vs[i]=Vec2dZero; }; }
    void clearMater    (){
        for(int imat=0; imat<nmat; imat++){
            double* Ns = matNs[imat];
            for(int i=0; i<nctot;i++){
                Ns[i]=0.0;
            }
        }
    }

    void material2cell(){
        for(int i=0; i<nctot; i++){ Ntot[i]=0; }
        for(int imat=0; imat<nmat; imat++){
            double* Ns = matNs[imat];
            for(int i=0; i<nctot;i++){
                Ntot[i] += Ns[i];
            }
        }
    }

    void updateCellTermodynamics(){
        for(int i=0; i<nctot;i++){
            // ToDo : Include Cell Volume for non-flat grid symmetry (spherical, cylindical)
            double oldT = Ts[i];
            double oldN = ps[i]/( oldT * const_Rgas);
            double N    = Ntot[i];
            double T    = pow( N/oldN, 2./3.);
            double p    = N*T*const_Rgas/Vs[i];
            Ts[i]       = T;
            ps[i]       = p;
        }
    }

    void updateCellPressure(){
        for(int i=0; i<nctot;i++){
            ps[i] = Ntot[i]*Ts[i]*const_Rgas/Vs[i]; //  p = nRT/V
            //if(i==(nc.x*100+100)){printf( "updateCellPressure[100,100]  p %g N %g T %g V %g \n", ps[i], Ntot[i],Ts[i],Vs[i] );}
        }
    }

    /*
    void accelerate( double dt ){
        int i=nc.x;
        for(int iy=1;iy<nc.y; iy++){
            i++;
            for(int ix=1;ix<nc.x; ix++){
                // ToDo : we may blure it here  use 8 neighbors instead 4 neighbor cells
                double fx = (ps[i] - ps[i-1   ])*invStep;
                double fy = (ps[i] - ps[i-nc.x])*invStep;
                if((ix==120)&&(iy==120))printf( "p %g f (%g,%g)\n", ps[i], fx,fy );
                vs[i].add_mul( {fx,fy}, dt );
                i++;
            }
        }
    }
    */

    void accelerate2( double dt ){
        int i=nc.x;
        for(int iy=1;iy<nc.y-1; iy++){
            i++;
            for(int ix=1;ix<nc.x-1; ix++){
                // ToDo : we may blure it here  use 8 neighbors instead 4 neighbor cells
                double fx = (ps[i-1   ]-ps[i+1   ])*invStep;
                double fy = (ps[i-nc.x]-ps[i+nc.x])*invStep;
                //if((ix==120)&&(iy==120))printf( "p %g f (%g,%g)\n", ps[i], fx,fy );
                vs[i].add_mul( {fx,fy}, dt );
                i++;
            }
            i++;
        }
    }

    inline double evalCohesion( double& Ni, double& Nj, double pi, double pj, double dt ){
        //double dp = (pi-pj-cohesive_pressure);
        //pi-=cohesive_pressure;
        //pj-=cohesive_pressure;
        double dp = pi-pj; // in reverse so we save some minus signs
        double dN = Ni-Nj;
        double Nmove;
        if(dN>0){
            dp=cohesive_pressure-dp;
            //dN*=(dp-cohesive_pressure);

        }else{  // dN<0
            dp=cohesive_pressure+dp;
            //dN*=(dp+cohesive_pressure);
        }
        //printf( " dp %g \n", dp );
        dN*=dp*dt;
        //dN*=dt;
        if(dN>0){ if(dN>Nj)dN=Nj; }else{ if(-dN>Ni)dN=-Ni; }   // make sure we dont move more than there is
        //return dN;
        Ni+=dN; Nj-=dN; // update
        return dN;
    }

    void cohesion(double dt){
        // This function tries to prevent bluring of density blobs
        // ToDo : impulses caused by moving material should update velocity !!!!!!!
        for(int imat=0; imat<nmat; imat++){
            double* Ns = matNs[imat];
            setBoundary( 0, Ntmp );
            for(int iy=1;iy<nc.y; iy++){
                for(int ix=1;ix<nc.x; ix++){
                    int i = iy*nc.x + ix;
                    double p00,p01,p10;
                    double N00,N01,N10;
                    const int i01=i-1,i10=i-nc.x;
                    N00=Ns[i  ]; p00=ps[i  ];
                    N01=Ns[i01]; p01=ps[i01];
                    N10=Ns[i10]; p10=ps[i10];

                    double dNx = evalCohesion( N00, N01, p00, p01, dt );
                    double dNy = evalCohesion( N00, N10, p00, p10, dt );
                    Ns[i  ]=N00;
                    Ns[i01]=N01;
                    Ns[i10]=N10;
                    //double dNx = evalCohesion( N00, N01, p00, p01, dt );
                    //double dNy = evalCohesion( N00, N10, p00, p10, dt );
                    //Ns[i  ]=N00+dNx+dNy;
                    //Ns[i01]=N01-dNx;
                    //Ns[i10]=N10-dNy;

                    //if((ix==8)&&(iy==8)){ printf("p00 %g \n", p00); };
                }
            }
        }
    }

    void cohesion2(double dt){
        // This function tries to prevent bluring of density blobs
        // ToDo : impulses caused by moving material should update velocity !!!!!!!
        for(int imat=0; imat<nmat; imat++){
            double* Ns = matNs[imat];
            setBoundary( 0, Ntmp );
            for(int iy=1;iy<nc.y; iy++){
                for(int ix=1;ix<nc.x; ix++){
                    int i = iy*nc.x + ix;
                    double p00,p01,p10,p11;
                    double N00,N01,N10,N11;
                    const int i01=i-1,i10=i-nc.x,i11=i10-1;
                    N00=Ns[i  ]; p00=ps[i  ];
                    N01=Ns[i01]; p01=ps[i01];
                    N10=Ns[i10]; p10=ps[i10];
                    N11=Ns[i11]; p11=ps[i11];
                    double dNx = evalCohesion( N00, N01, p00, p01, dt           );
                    double dNy = evalCohesion( N00, N10, p00, p10, dt           );
                    double dNa = evalCohesion( N00, N11, p00, p11, dt*M_SQRT1_2 );
                    double dNb = evalCohesion( N01, N10, p01, p10, dt*M_SQRT1_2 );
                    Ns[i  ]=N00;
                    Ns[i01]=N01;
                    Ns[i10]=N10;
                    Ns[i11]=N11;
                    //Ns[i  ]=N00+dNx+dNy+dNa;
                    //Ns[i01]=N01-dNx+dNb;
                    //Ns[i10]=N10-dNy-dNb;
                    //Ns[i11]=N11-dNa;
                    //if((ix==8)&&(iy==8)){ printf("p00 %g \n", p00); };
                }
            }
        }
    }

    void viscosity(double dt){ // smooth out velocities
    }

    /*
    double advec(double dt){
        //printf( "nc(%i,%i) \n", nc.x, nc.y );
        sumN = 0;
        for(int imat=0; imat<nmat; imat++){
            double* Ns = matNs[imat];
            setBoundary( 0, Ntmp );
            for(int iy=1;iy<nc.y-1; iy++){
                for(int ix=1;ix<nc.x-1; ix++){
                    int i = iy*nc.x + ix;
                    //printf( "(%i,%i)\n", ix, iy );
                    //if( (i<0)||(i>=nctot) ){ printf( "ERROR i=%i (%i,%i) |  \n", i, ix, iy, nctot ); exit(0); }
                    Vec2d v = vs[i];
                    //Vec2d v = (vs[i] + vs[i+1] + vs[i+nc.x] + vs[i+nc.x+1])*0.25;
                    v.mul(-dt);
                    int dix=1,diy=nc.x;
                    if(v.x<0){dix=-dix;v.x=-v.x;};
                    if(v.y<0){diy=-diy;v.y=-v.y;};
                    if(v.x>1)v.x=1;
                    if(v.y>1)v.y=1;
                    double mx = 1-v.x;
                    double my = 1-v.y;
                    //Ns[i]   = Ns[i-1];
                    //if( ((i        )<0)||(i>=nctot) ){ printf( "ERROR i=%i (%i,%i) |  \n", i, ix, iy, nctot ); exit(0); }
                    //if( ((i+diy    )<0)||(i>=nctot) ){ printf( "ERROR i=%i (%i,%i) |  \n", i, ix, iy, nctot ); exit(0); }
                    //if( ((i    +dix)<0)||(i>=nctot) ){ printf( "ERROR i=%i (%i,%i) |  \n", i, ix, iy, nctot ); exit(0); }
                    //if( ((i+diy+dix)<0)||(i>=nctot) ){ printf( "ERROR i=%i (%i,%i) |  \n", i, ix, iy, nctot ); exit(0); }
                    Ntmp[i] = mx *my *Ns[i        ]
                            + mx *v.y*Ns[i+diy    ]
                            + v.x*my *Ns[i    +dix]
                            + v.x*v.y*Ns[i+diy+dix];
                    //if((ix==8)&&(iy==7)){
                    //    //printf( "%i vs %i \n", i, iy*nc.x + iy );
                    //    printf( "m[%i](%g,%g) (%g,%g) (%i,%i) | %g %g %g %g -> %g \n", imat,  v.x, v.y,  mx,my, dix, diy,  Ns[i], Ns[i+diy], Ns[i+dix], Ns[i+diy+dix], Ntmp[i]  );
                    //}
                }
            }
            for(int i=0; i<nctot; i++){
                double N = Ntmp[i];
                //printf( "N[%i] %g \n", i, N );
                //if((N<0)||(N>1)) printf("ERROR [i] %g \n", i, N);
                sumN+=N; Ns[i]=N;
            }
        }
        return sumN;
    }
    */


    double advec2(double dt){
        //printf( "nc(%i,%i) \n", nc.x, nc.y );
        sumN = 0;
        for(int imat=0; imat<nmat; imat++){
            double* Ns = matNs[imat];
            setBoundary( 0, Ntmp );
            for(int iy=1;iy<nc.y-1; iy++){
                for(int ix=1;ix<nc.x-1; ix++){
                    int i = iy*nc.x + ix;
                    //printf( "(%i,%i)\n", ix, iy );
                    //if( (i<0)||(i>=nctot) ){ printf( "ERROR i=%i (%i,%i) |  \n", i, ix, iy, nctot ); exit(0); }
                    Vec2d v = vs[i];
                    //Vec2d v = (vs[i] + vs[i+1] + vs[i+nc.x] + vs[i+nc.x+1])*0.25;
                    v.mul(dt);
                    int dix=1,diy=nc.x;
                    if(v.x<0){dix=-dix;v.x=-v.x;};
                    if(v.y<0){diy=-diy;v.y=-v.y;};
                    if(v.x>1)v.x=1;
                    if(v.y>1)v.y=1;
                    double mx = 1-v.x;
                    double my = 1-v.y;
                    //Ns[i]   = Ns[i-1];
                    //if( ((i        )<0)||(i>=nctot) ){ printf( "ERROR i=%i (%i,%i) |  \n", i, ix, iy, nctot ); exit(0); }
                    //if( ((i+diy    )<0)||(i>=nctot) ){ printf( "ERROR i=%i (%i,%i) |  \n", i, ix, iy, nctot ); exit(0); }
                    //if( ((i    +dix)<0)||(i>=nctot) ){ printf( "ERROR i=%i (%i,%i) |  \n", i, ix, iy, nctot ); exit(0); }
                    //if( ((i+diy+dix)<0)||(i>=nctot) ){ printf( "ERROR i=%i (%i,%i) |  \n", i, ix, iy, nctot ); exit(0); }
                    double& N00 = Ns[i        ];
                    double& N01 = Ns[i    +dix];
                    double& N10 = Ns[i+diy    ];
                    double& N11 = Ns[i+diy+dix];
                    double dN00 = mx *my *N00;
                    double dN01 = v.x*my *N01;
                    double dN10 = mx *v.y*N10;
                    double dN11 = v.x*v.y*N11;
                    N00  = dN00 + dN01 + dN10 + dN11;
                    N01 -= dN01;
                    N10 -= dN10;
                    N11 -= dN11;
                    //if((ix==8)&&(iy==7)){
                    //    //printf( "%i vs %i \n", i, iy*nc.x + iy );
                    //    printf( "m[%i](%g,%g) (%g,%g) (%i,%i) | %g %g %g %g -> %g \n", imat,  v.x, v.y,  mx,my, dix, diy,  Ns[i], Ns[i+diy], Ns[i+dix], Ns[i+diy+dix], Ntmp[i]  );
                    //}
                }
            }
        }
        return sumN;
    }


    void update(double dt){
        material2cell();
        updateCellPressure();
        accelerate2  (dt);
        //viscosity    (dt);
        //cohesion     (dt);
        //cohesion2    (dt);
        //advec        (dt);
        advec2        (dt);
    }

};

#endif

