
#ifndef BondAdaptedMesh_h
#define BondAdaptedMesh_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mesh.h"



class AtomGrid{ public:
    OMesh mesh;
    std::vector<Tetrahedron> tetrahedrons;

    void fromNeighbors( int ia0, int nneigh, int* neighs, Vec3d * poss, double Rmax ){
        //Mesh mesh;
        //Plane3D* planes = new Plane3D[nneigh*2];
        int nplanes = 2*nneigh;
        std::vector<Plane3D> planes{nplanes+6}; // we need space for additional planes
        Vec3d p0=poss[ia0];
        for(int ja=0; ja<nneigh; ja++){
            int jng=neighs[ja];
            int ip=ja*2;
            planes[ip  ].fromVec( p0-poss[jng] );
            planes[ip  ].iso*=-0.5;
            planes[ip+1].normal=planes[ip].normal*-1.00;
            planes[ip+1].iso = -Rmax;
        }
        // TODO: check if some direction is missing:
        //  - take planes[0].normal and find second far from colinear
        //  -
        // TODO: check if two neighbors are close to colinear
        mesh.fromPlanes( nplanes, &planes[0] );

        for(int i=0;i<nplanes; i++){
            int np  = mesh.polygons[i]->ipoints.size();
            if(np<3){
                printf( "ERROR in AtomGrid.fromNeighbors() degenerated polygon[%i]", i  );
                exit(0);
                continue;
            }
            if(np==3){
                int ip1 = mesh.polygons[i]->ipoints[0];
                int ip2 = mesh.polygons[i]->ipoints[1];
                int ip3 = mesh.polygons[i]->ipoints[2];
                tetrahedrons.push_back( {p0, mesh.points[ip1], mesh.points[ip2], mesh.points[ip3]} );
                continue;
            }
            int ipo = mesh.polygons[i]->ipoints[np-1];
            //Vec3d pb = (p0+poss[i])*0.5;
            Vec3d pb = p0 + planes[i].normal * planes[i].iso;

            for(int j=0; j<np; j++){
                int ip = mesh.polygons[i]->ipoints[j];
                tetrahedrons.push_back( {p0, pb, mesh.points[ipo], mesh.points[ip]} );
                ipo=ip;
            }
        }
        //delete planes;
    }

};

#endif
