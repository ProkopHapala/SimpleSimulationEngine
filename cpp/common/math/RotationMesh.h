
#ifndef  RotationMesh_h
#define  RotationMesh_h

#include <vector>
#include <unordered_map>
#include <cstring>
#include <string>

// implementation
//#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <string>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
//#include "raytrace.h"
//#include "geom3D.h"

class RotationMesh{
    public:

    int nMax=0;
    int n   =0;           // number of rotations
    Quat4d * rots = NULL;
    Mat3d  * mrots = NULL;   // just DEBUG

    int nMaxNeigh  =0;
    int*  nneighs  = NULL;
    int** neighs   = NULL;

    void fromDirsNroll( int ndirs, int nroll, Vec3d * dirs, const Vec3d& pole ){
        n=0;
        for(int i=0; i<ndirs; i++){
            Mat3d mrot;
            mrot.c=dirs[i]; mrot.c.normalize();
            mrot.fromDirUp( mrot.c, pole );
            Quat4d qrot; qrot.fromMatrix(mrot);
            double a  = M_PI/nroll;
            double ca = cos(a);
            double sa = sin(a);
            for(int j=0; j<nroll; j++){
                //printf( "%i %i %i %i\n", i, j, n, nMax );
                rots[n]=qrot;
                qrot.roll2( ca, sa );
                //qrot.roll2( 2*M_PI/nroll );
                mrots[n]=mrot;
                double ca_ = ca*ca-sa*sa;
                double sa_ = 2*ca*sa;
                mrot.a.rotate_csa(ca_,sa_, mrot.c);
                mrot.b.rotate_csa(ca_,sa_, mrot.c);
                n++;
            }
        }
    }

    int findNeighs(double maxDist){
        int nntot=0;
        for(int i=0;i<n;i++){
            Quat4d roti = rots[i];
            int ni=0;
            for(int j=0;j<n;j++){
                if(i==j) continue;
                double d = roti.dist_cos(rots[j]);
                //printf( "%i %f %f (%3.3f,%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f,%3.3f)\n", j, d, roti.dot(rots[j]), rots[j].x, rots[j].y, rots[j].z, rots[j].w,  roti.x,roti.y,roti.z,roti.w );
                if(d<maxDist){
                    if(ni>=nMaxNeigh){ printf("ERROR: findNeighs : ni>=nMaxNeigh\n"); exit(0); }
                    neighs[i][ni] = j;
                    ni++;
                }
            }
            //we may also order neighbors by distance ?
            nneighs[i]=ni;
            nntot+=ni;
            //exit(0);
        }
        return nntot;
    }

    void allocate(int nMax_,int nMaxNeigh_){
        nMax=nMax_; nMaxNeigh=nMaxNeigh_;
        rots   = new Quat4d[nMax];
        mrots  = new Mat3d[nMax];
        neighs = new int*  [nMax];
        nneighs= new int   [nMax];
        for(int i=0; i<nMax; i++){
            nneighs[i]= 0;
            neighs[i]=new int[nMaxNeigh];
        }
    }

};

#endif



