
#ifndef  TrajectoryVariation_h
#define  TrajectoryVariation_h

#include "Vec3.h"
#include "Mat3.h"
#include "CubicBSpline.h"

// Gauss-Legendre 6 // Gauss-Legendre converge fastest on test functions, for polynoms of degree <=6 6 points gives exact integral
//static const double IntegrationNodes  [6] ={ 0.96623475710157603d,  0.83060469323313224d, 0.61930959304159838d, 0.38069040695840151d, 0.1693953067668677d,  0.033765242898423975d };
static const double IntegrationNodes  [6] ={   0.033765242898423975d, 0.1693953067668677d, 0.38069040695840151d, 0.61930959304159838d, 0.83060469323313224d, 0.96623475710157603d };
static const double IntegrationWeights[6] ={ 0.085662246189585275d, 0.18038078652406944d, 0.23395696728634541d, 0.23395696728634555d, 0.18038078652406908d, 0.085662246189584831d };

/*

 Tx   = cx*ddx - Gx(x,y)
 Ty   = cy*ddy - Gy(x,y)
 x_cx = bx

 ( Tx^2 + Ty^2 )_cx1 = Tx*( ddx - Gx_x * x_cx    )  + Ty*( 0 - Gy_x * x_cx  )
                     = Tx*( ddx - Gx_x * bx      )  + Ty*( 0 - Gy_x * bx    )

*/

bool DEBUG_PLOT = true;

inline void addGravity( const Vec3d p, Vec3d& T, Vec3d& T2_b ){
    // T is thrust vector; T2_b is variation of |T|^2 with basis function
    constexpr double GM = -1.0d;  // TODO read mass of body
    Vec3d  dp  = p;               // TODO substract position of body

    double ir2 = 1.0d/dp.norm2();
    double ir  = sqrt( ir2 );
    double ir3 = ir2*ir;
    Vec3d  G   = p*(ir3*GM);

    // this can be prhaps optimized
    Mat3d  DG;
    DG.set_outer(dp,G);
    //DG.set_outer(G,dp);
    DG.mul      (3*ir2);
    DG.diag_add (  ir3);

    if(DEBUG_PLOT){
        glColor3f(0.0f,0.0f,0.0f); Draw3D::drawPoint   (       p);
        glColor3f(1.0f,0.0f,1.0f); Draw3D::drawVecInPos(T*-1.0,p); // Centrifugal force
        glColor3f(0.0f,0.0f,1.0f); Draw3D::drawVecInPos(G     ,p); // Gravity
        glColor3f(1.0f,0.0f,0.0f); Draw3D::drawVecInPos(T-G   ,p); // Thrust
    }

    T   .sub(G);
    DG.dot_to( T, T2_b );
    //Vec3d T2_b_; DG.dot_to( T, T2_b_ );
    //T2_b.add(T2_b_);
}


double getTrajectoryVariation(int ncoef, Vec3d * CPs, Vec3d * dCPs, double u_t ){
//double getDerivs( double [][] Cs, double [][] dCs){

    //double T      = PI*0.25;
    //double u_t    = nseg/T;
    //double u_t    = T/nseg;
    //double t_u    = 1.0/u_t;
    double u_t2   = u_t*u_t;

    for (int i=0;i<ncoef;i++){ dCPs[i].set( 0.0d );}
    double sumT2  = 0;
    for (int i=0;i<ncoef-3;i++){
        Vec3d cp0 = CPs[i  ];
        Vec3d cp1 = CPs[i+1];
        Vec3d cp2 = CPs[i+2];
        Vec3d cp3 = CPs[i+3];
        Vec3d f0,f1,f2,f3; f0.set(0.0d); f1.set(0.0d); f2.set(0.0d); f3.set(0.0d);
        for (int iu=0; iu<6; iu++){
            double      u =  IntegrationNodes  [iu];
            double      w =  IntegrationWeights[iu];

            // --- evauleate acceleration
            double ddb0,ddb1,ddb2,ddb3;
            CubicBSpline::ddbasis(u,ddb0,ddb1,ddb2,ddb3);
            Vec3d  T = cp0*ddb0 + cp1*ddb1 + cp2*ddb2 + cp3*ddb3;

            // --- Evaluate Gravity
            double b0,b1,b2,b3;
            CubicBSpline::basis( u, b0, b1, b2, b3 );
            Vec3d  p = cp0*b0 + cp1*b1 + cp2*b2 + cp3*b3;
            Vec3d T2_b;
            addGravity( p, T, T2_b );

            // correction of boundary conditions
            if(i<5){
                if(i==4){
                    b1   += 0.5*  b0;
                    ddb1 += 0.5*ddb0;
                }else{
                    b2   += 0.5*  b1  -   b0;
                    ddb2 += 0.5*ddb1  - ddb0;
                }
            }else if( i>ncoef-3 ){
                if(i==ncoef-2){
                    b2   += 0.5*  b3;
                    ddb2 += 0.5*ddb3;
                }else{
                    b1   += 0.5*  b2  +   b3;
                    ddb1 += 0.5*ddb2  + ddb3;
                }
            }

            // --- acum variational derivs
            sumT2      += T.norm2()*w;

            w *=-2.0;
            f0.add_mul( T*ddb0 + T2_b*b0, w );
            f1.add_mul( T*ddb1 + T2_b*b1, w );
            f2.add_mul( T*ddb2 + T2_b*b2, w );
            f3.add_mul( T*ddb3 + T2_b*b3, w );

            /*
            sumT2+= w*( T.x*T.x + T.y*T.y ); // length of thrust vector
            double T2_ddbx =  2*w*  Tx                 ;
            double T2_bx   = -2*w*( Tx*Gx_x + Ty*Gy_x );
            fx0+= T2_ddbx * ddb0 + T2_bx * b0;
            fx1+= T2_ddbx * ddb1 + T2_bx * b1;
            fx2+= T2_ddbx * ddb2 + T2_bx * b2;
            fx3+= T2_ddbx * ddb3 + T2_bx * b3;
            double T2_ddby =  2*w*  Ty                 ;
            double T2_by   = -2*w*( Tx*Gx_y + Ty*Gy_y );
            fy0+= T2_ddby * ddb0 + T2_by * b0;
            fy1+= T2_ddby * ddb1 + T2_by * b1;
            fy2+= T2_ddby * ddb2 + T2_by * b2;
            fy3+= T2_ddby * ddb3 + T2_by * b3;
            */

        } // ncoef

        dCPs[i  ].add( f0 );
        dCPs[i+1].add( f1 );
        dCPs[i+2].add( f2 );
        dCPs[i+3].add( f3 );
    }
    return sumT2;
}

#endif
