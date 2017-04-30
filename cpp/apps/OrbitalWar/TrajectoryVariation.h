
#ifndef  TrajectoryVariation_h
#define  TrajectoryVariation_h

// Gauss-Legendre 6 // Gauss-Legendre converge fastest on test functions, for polynoms of degree <=6 6 points gives exact integral
static const double IntegrationNodes  [6] ={ 0.96623475710157603d,  0.83060469323313224d, 0.61930959304159838d, 0.38069040695840151d, 0.1693953067668677d,  0.033765242898423975d };
static const double IntegrationWeights[6] ={ 0.085662246189585275d, 0.18038078652406944d, 0.23395696728634541d, 0.23395696728634555d, 0.18038078652406908d, 0.085662246189584831d };

/*

 Tx   = cx*ddx - Gx(x,y)
 Ty   = cy*ddy - Gy(x,y)
 x_cx = bx

 ( Tx^2 + Ty^2 )_cx1 = Tx*( ddx - Gx_x * x_cx    )  + Ty*( 0 - Gy_x * x_cx  )
                     = Tx*( ddx - Gx_x * bx      )  + Ty*( 0 - Gy_x * bx    )

*/




double getDerivs(int ncoef, Vec3d * CPs, Vec3d * dCPs, double u_t ){
//double getDerivs( double [][] Cs, double [][] dCs){

    //double T      = PI*0.25;
    //double u_t    = nseg/T;
    //double u_t    = T/nseg;
    //double t_u    = 1.0/u_t;
    double u_t2   = u_t*u_t;

    for (int i=0;i<ncoef;i++){ dCs[i].set( 0.0d );}
    double sumT2  = 0;
    for (int i=3;i<ncoef;i++){
        double i6 = 1.0d/6.0d; // TODO shoudl be removed

        Vec3d p0 = CPs[i  ];
        Vec3d p1 = CPs[i+1];
        Vec3d p2 = CPs[i+2];
        Vec3d p3 = CPs[i+3];

        Vec3d f0,f1,f2,f3; f0.set(0.0d); f1.set(0.0d); f2.set(0.0d); f3.set(0.0d);

        for (int iu=0; iu<6; iu++){
            double      u =  IntegrationNodes  [iu];
            double      w =  IntegrationWeights[iu];

            // Evaluate Gravity
            // CubicBSpline::basis(u,ddb0,ddb1,ddb2,ddb2);
            double u2  = u*u; double u3 = u2*u;
            double b0  = i6*(                    u3 );
            double b1  = i6*( 1 + 3*u + 3*u2 - 3*u3 );
            double b2  = i6*( 4       - 6*u2 + 3*u3 );
            double b3  = i6*( 1 - 3*u + 3*u2 -   u3 );

            // getGravity( p, ); ... put in function
            double x   = cx0*b0 + cx1*b1 + cx2*b2 + cx3*b3;
            double y   = cy0*b0 + cy1*b1 + cy2*b2 + cy3*b3;
            double ir2 = 1.0d/( x*x + y*y );
            double ir  = Math.sqrt( ir2 );
            double ir3 =  ir2*ir;
            double Gx  =  x*ir3;
            double Gy  =  y*ir3;

            double Gx_x  =  3*x*Gx*ir2 - ir3;
            double Gx_y  =  3*y*Gx*ir2;
            double Gy_y  =  3*y*Gy*ir2 - ir3;
            double Gy_x  =  3*x*Gy*ir2;

            // evauleate acceleration
            // CubicBSpline::ddbasis(u,ddb0,ddb1,ddb2,ddb2);
            double ddb0 = i6*(       6*u)*u_t2;
            double ddb1 = i6*(  6 - 18*u)*u_t2;
            double ddb2 = i6*(-12 + 18*u)*u_t2;
            double ddb3 = i6*(  6 -  6*u)*u_t2;
            double ddx  = cx0*ddb0 + cx1*ddb1 + cx2*ddb2 + cx3*ddb3;
            double ddy  = cy0*ddb0 + cy1*ddb1 + cy2*ddb2 + cy3*ddb3;

            // put i all together  ... Thrust vector = inertial acceleration + gravity
            double Tx   = ddx + Gx;
            double Ty   = ddy + Gy;

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

            // F0.add( T2_ddb * ddb0 + T2_bx * b0);
            // F1.add( T2_ddb * ddb0 + T2_bx * b0);
            // F2.add( T2_ddb * ddb0 + T2_bx * b0);
            // F3.add( T2_ddb * ddb0 + T2_bx * b0);

            sumT2+= w*( Tx*Tx + Ty*Ty ); // length of thrust vector
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

        } // ncoef


        dCs[0][i  ] -= fx0;
        dCs[0][i-1] -= fx1;
        dCs[0][i-2] -= fx2;
        dCs[0][i-3] -= fx3;
        dCs[1][i  ] -= fy0;
        dCs[1][i-1] -= fy1;
        dCs[1][i-2] -= fy2;
        dCs[1][i-3] -= fy3;

    }
    return sumT2;
}

#endif
