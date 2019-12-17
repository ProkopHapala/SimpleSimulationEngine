
#ifndef  grids2d_h
#define  grids2d_h

//#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"


//  numlock assigment
//    7 8 9
//   4 5 6
//    2 3

inline void subLoopTile(
    double p7, double p8, double p9, double p4, double p5, double p6, double p2, double p3,
    double& o5, double& o6, double& o8, double& o9
){
    constexpr double c_para = 0.375d;
    constexpr double c_perp = 0.125d;
    constexpr double c_on  = 0.625d;
    constexpr double c_off = 0.0625d;
    o5 = p5*c_on + ( p7 + p8 + p4 + p6 + p2 + p3 )*c_off;
    o6 = (p5+p6)*c_para + (p3+p8)*c_perp;
    o8 = (p5+p8)*c_para + (p6+p7)*c_perp;
    o9 = (p6+p8)*c_para + (p5+p9)*c_perp;
}

double * subdivideLoopGrid( int nx, int ny, double * ps, double * ps_ ){
    int nx_=nx<<1;
    int ny_=ny<<1;
    printf( "ny %i nx %i \n", ny, nx );
    //double * ps_ = new double[nx_*ny_];
    // currently works only for PBC
    for (int iy=0; iy<ny; iy++){
        int iyp=iy+1; if(iyp>=ny) iyp=iyp-ny;
        int iym=iy-1; if(iym<0  ) iym=ny+iym;
        printf(  " >> iy %i %i %i \n", iy,iym,iyp );
        iyp*=nx; iym*=nx;
        for (int ix=0; ix<nx; ix++){
            int ixp=ix+1; if(ixp>=nx) ixp=ixp-nx;
            int ixm=ix-1; if(ixm<0  ) ixm=nx+ixm;
            int j0 = (iy<<1)*nx_+(ix<<1);
            printf(  " > ix %i %i %i \n", ix,ixm,ixp );
            printf( "nxy  %i 7:%i 8:%i 9:%i  4:%i 5:%i 6:%i  2:%i 3:%i \n", nx*ny, iyp+ixm,iyp+ix,iyp+ixp,  iy+ixm,iy+ix,iy+ixp,  iym, iym+ixp );
            printf( "nxy_ %i : %i   %i %i %i \n", nx_*ny_, j0,j0+1,j0+nx_,j0+nx_+1 );
            //printf( " %i %i   %i %i   %i %i \n", iy,ix,  iym,iyp,   ixm,ixp  );
            subLoopTile( ps[iyp+ixm],ps[iyp+ix], ps[iyp+ixp],   ps[iy+ixm],ps[iy+ix],ps[iy+ixp],  ps[iym],ps[iym+ixp],   ps_[j0],ps_[j0+1],ps_[j0+nx_],ps_[j0+nx_+1] );
            //int i0 = iy*nx+ix;
            //subLoopTile( ps[i0+nx-1],ps[i0+nx], ps[i0+nx+1],   ps[i0-1],ps[i0],ps[i0-1],  ps[i0-nx],ps[i0-nx+1],   ps_[j0],ps_[j0+1],ps_[j0+nx_],ps_[j0+nx_+1] );
        }
    }
    /*
    for (int iy=1; iy<ny-1; iy++){
            int i0 = iy*nx;
            int j0 = (iy<<1)*nx_;
            // left
            subLoopTile( ps[i0+nx+nx-1],ps[i0+nx], ps[i0+nx+1],   ps[i0+nx-1],ps[i0],ps[i0-1],  ps[i0-nx],ps[i0-nx+1],   ps_[j0],ps_[j0+1],ps_[j0+nx_],ps_[j0+nx_+1] );
            // right
            subLoopTile( ps[i0+nx-1],ps[i0+nx], ps[i0+nx+1],   ps[i0-1],ps[i0],ps[i0-1],  ps[i0-nx],ps[i0-nx+1],   ps_[j0],ps_[j0+1],ps_[j0+nx_],ps_[j0+nx_+1] );
    }
    for (int ix=1; ix<nx-1; ix++){
            // top
            subLoopTile( ps[i0+nx-1],ps[i0+nx], ps[i0+nx+1],   ps[i0-1],ps[i0],ps[i0-1],  ps[i0-nx],ps[i0-nx+1],   ps_[j0],ps_[j0+1],ps_[j0+nx_],ps_[j0+nx_+1] );
            // bottom
            subLoopTile( ps[i0+nx-1],ps[i0+nx], ps[i0+nx+1],   ps[i0-1],ps[i0],ps[i0-1],  ps[i0-nx],ps[i0-nx+1],   ps_[j0],ps_[j0+1],ps_[j0+nx_],ps_[j0+nx_+1] );
    }
    */
    return ps_;
}

#endif

