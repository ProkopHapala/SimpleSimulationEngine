
#ifndef GridUtils_h
#define GridUtils_h

#include "Grid.h"

double coulombGrid_brute( Vec3i ns1, Vec3i ns2, Vec3d pos1, Vec3d pos2, Mat3d rot1, Mat3d rot2, double* rho1, double* rho2, double offset, double damp ){
    int nyz1=ns1.y*ns1.z;
    int nyz2=ns2.y*ns2.z;
    double sum = 0;
    for(int ix=0; ix<ns1.x; ix++){
        //printf("ix %i \n", ix );
        for(int iy=0; iy<ns1.y; iy++){
            for(int iz=0; iz<ns1.z; iz++){
                Vec3d pi;  rot1.dot_to_T({ix+offset,iy+offset,iz+offset},pi);    //= rot.a*ix + rot.c*iy + rot.c*iz ;
                pi.add(pos1);
                double qi = rho1[ix*nyz1 + iy*ns1.z + iz];
                for(int jx=0; jx<ns2.x; jx++){
                    for(int jy=0; jy<ns2.y; jy++){
                        for(int jz=0; jz<ns2.z; jz++){
                            Vec3d pj;  rot2.dot_to_T({jx+offset,jy+offset,jz+offset},pj);
                            pj.add(pos2);
                            Vec3d p; p.set_sub(pi,pj);
                            //double r  = p.norm();
                            double qj = rho2[jx*nyz2 + jy*ns2.z + jz];
                            sum += qi*qj/sqrt( p.norm2() + damp );
                        }
                    }
                }
            }
        }
    }
    return sum;
}

#endif









