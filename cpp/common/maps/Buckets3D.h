#ifndef  Buckets3D_h
#define  Buckets3D_h

#include "integerOps.h"
#include "fastmath.h"
#include "Vec3.h"
#include "grids3D.h"
#include "Buckets.h"

class GridRulerInterface{public:};

class Buckets3D : public Buckets, public CubeGridRuler { public:

    void pointsToCells( int np, Vec3d* ps ){
        resizeObjs ( np, true );
        for( int i=0; i<np; i++ ){
           obj2cell_bind[i] = icell( ps[i] );
        }
        updateCells( np, obj2cell_bind );
    }

    int getNeighbors( Vec3i ip, int* neighs ){
        int nfound = 0;
        nfound += getInCell( ixyz2i( {ip.x+1,ip.y+1,ip.z+1} ), neighs+nfound );  // 1
        nfound += getInCell( ixyz2i( {ip.x  ,ip.y+1,ip.z+1} ), neighs+nfound );  // 2
        nfound += getInCell( ixyz2i( {ip.x-1,ip.y+1,ip.z+1} ), neighs+nfound );  // 3
        nfound += getInCell( ixyz2i( {ip.x+1,ip.y  ,ip.z+1} ), neighs+nfound );  // 4
        nfound += getInCell( ixyz2i( {ip.x  ,ip.y  ,ip.z+1} ), neighs+nfound );  // 5
        nfound += getInCell( ixyz2i( {ip.x-1,ip.y  ,ip.z+1} ), neighs+nfound );  // 6
        nfound += getInCell( ixyz2i( {ip.x+1,ip.y-1,ip.z+1} ), neighs+nfound );  // 7
        nfound += getInCell( ixyz2i( {ip.x  ,ip.y-1,ip.z+1} ), neighs+nfound );  // 8
        nfound += getInCell( ixyz2i( {ip.x-1,ip.y-1,ip.z+1} ), neighs+nfound );  // 9

        nfound += getInCell( ixyz2i( {ip.x+1,ip.y+1,ip.z  } ), neighs+nfound );  // 10
        nfound += getInCell( ixyz2i( {ip.x  ,ip.y+1,ip.z  } ), neighs+nfound );  // 11
        nfound += getInCell( ixyz2i( {ip.x-1,ip.y+1,ip.z  } ), neighs+nfound );  // 12
        nfound += getInCell( ixyz2i( {ip.x+1,ip.y  ,ip.z  } ), neighs+nfound );  // 13
        //nfound += getInCell( ixyz2i( {ip.x  ,ip.y  ,ip.z  } ), neighs+nfound );  // 14 Allready inside
        nfound += getInCell( ixyz2i( {ip.x-1,ip.y  ,ip.z  } ), neighs+nfound ); // 15
        nfound += getInCell( ixyz2i( {ip.x+1,ip.y-1,ip.z  } ), neighs+nfound ); // 16
        nfound += getInCell( ixyz2i( {ip.x  ,ip.y-1,ip.z  } ), neighs+nfound ); // 17
        nfound += getInCell( ixyz2i( {ip.x-1,ip.y-1,ip.z  } ), neighs+nfound ); // 18

        nfound += getInCell( ixyz2i( {ip.x+1,ip.y+1,ip.z-1} ), neighs+nfound ); // 19
        nfound += getInCell( ixyz2i( {ip.x  ,ip.y+1,ip.z-1} ), neighs+nfound ); // 20
        nfound += getInCell( ixyz2i( {ip.x-1,ip.y+1,ip.z-1} ), neighs+nfound ); // 21
        nfound += getInCell( ixyz2i( {ip.x+1,ip.y  ,ip.z-1} ), neighs+nfound ); // 22
        nfound += getInCell( ixyz2i( {ip.x  ,ip.y  ,ip.z-1} ), neighs+nfound ); // 23
        nfound += getInCell( ixyz2i( {ip.x-1,ip.y  ,ip.z-1} ), neighs+nfound ); // 24
        nfound += getInCell( ixyz2i( {ip.x+1,ip.y-1,ip.z-1} ), neighs+nfound ); // 25
        nfound += getInCell( ixyz2i( {ip.x  ,ip.y-1,ip.z-1} ), neighs+nfound ); // 26
        nfound += getInCell( ixyz2i( {ip.x-1,ip.y-1,ip.z-1} ), neighs+nfound ); // 27
        return nfound;
    }


    int getForwardNeighbors( Vec3i ip, int* neighs ){
        int nfound = 0;
      //nfound += getInCell( ixyz2i( {ip.x+1,ip.y+1,ip.z+1} ), neighs+nfound );  // 1
        nfound += getInCell( ixyz2i( {ip.x-1,ip.y-1,ip.z-1} ), neighs+nfound ); // 27
      //nfound += getInCell( ixyz2i( {ip.x  ,ip.y+1,ip.z+1} ), neighs+nfound );  // 2
        nfound += getInCell( ixyz2i( {ip.x  ,ip.y-1,ip.z-1} ), neighs+nfound ); // 26
      //nfound += getInCell( ixyz2i( {ip.x-1,ip.y+1,ip.z+1} ), neighs+nfound );  // 3
        nfound += getInCell( ixyz2i( {ip.x+1,ip.y-1,ip.z-1} ), neighs+nfound ); // 25

      //nfound += getInCell( ixyz2i( {ip.x+1,ip.y  ,ip.z+1} ), neighs+nfound );  // 4
        nfound += getInCell( ixyz2i( {ip.x-1,ip.y  ,ip.z-1} ), neighs+nfound ); // 24
      //nfound += getInCell( ixyz2i( {ip.x  ,ip.y  ,ip.z+1} ), neighs+nfound );  // 5
        nfound += getInCell( ixyz2i( {ip.x  ,ip.y  ,ip.z-1} ), neighs+nfound ); // 23
      //nfound += getInCell( ixyz2i( {ip.x-1,ip.y  ,ip.z+1} ), neighs+nfound );  // 6
        nfound += getInCell( ixyz2i( {ip.x+1,ip.y  ,ip.z-1} ), neighs+nfound ); // 22

      //nfound += getInCell( ixyz2i( {ip.x+1,ip.y-1,ip.z+1} ), neighs+nfound );  // 7
        nfound += getInCell( ixyz2i( {ip.x-1,ip.y+1,ip.z-1} ), neighs+nfound ); // 21
      //nfound += getInCell( ixyz2i( {ip.x  ,ip.y-1,ip.z+1} ), neighs+nfound );  // 8
        nfound += getInCell( ixyz2i( {ip.x  ,ip.y+1,ip.z-1} ), neighs+nfound ); // 20
      //nfound += getInCell( ixyz2i( {ip.x-1,ip.y-1,ip.z+1} ), neighs+nfound );  // 9
        nfound += getInCell( ixyz2i( {ip.x+1,ip.y+1,ip.z-1} ), neighs+nfound ); // 19

      //nfound += getInCell( ixyz2i( {ip.x+1,ip.y+1,ip.z  } ), neighs+nfound );  // 10
        nfound += getInCell( ixyz2i( {ip.x-1,ip.y-1,ip.z  } ), neighs+nfound ); // 18
      //nfound += getInCell( ixyz2i( {ip.x  ,ip.y+1,ip.z  } ), neighs+nfound );  // 11
        nfound += getInCell( ixyz2i( {ip.x  ,ip.y-1,ip.z  } ), neighs+nfound ); // 17
      //nfound += getInCell( ixyz2i( {ip.x-1,ip.y+1,ip.z  } ), neighs+nfound );  // 12
        nfound += getInCell( ixyz2i( {ip.x+1,ip.y-1,ip.z  } ), neighs+nfound ); // 16

      //nfound += getInCell( ixyz2i( {ip.x+1,ip.y  ,ip.z  } ), neighs+nfound );  // 13
        nfound += getInCell( ixyz2i( {ip.x-1,ip.y  ,ip.z  } ), neighs+nfound ); // 15

        //nfound += getInCell( ixyz2i( {ip.x  ,ip.y  ,ip.z  } ), neighs+nfound );  // 14 Allready inside
        return nfound;
    }

};


#endif
