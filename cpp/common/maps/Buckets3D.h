#ifndef  Buckets3D_h
#define  Buckets3D_h

#include "integerOps.h"
#include "fastmath.h"
#include "Vec3.h"
#include "grids3D.h"
#include "Buckets.h"

class Buckets3D : public Buckets, public CubeGridRuler { public:

    bool bResizing=false;
    int* neighs_in=0;
    int  neighs_in_size=0; 

    void setup_Buckets3D( Vec3d pmin, Vec3d pmax_, double step ){
        CubeGridRuler::setup( pmin, pmax_, step  );
        Buckets::resizeCells( ntot );
    }

    void updateNeighsBufferSize(){
        int new_N = maxInBucket*14;
        if(neighs_in_size<new_N ){  printf("Buckets3D::updateNeighsBufferSize() new_N %i \n", new_N);  neighs_in_size=new_N;  _realloc(neighs_in,neighs_in_size); }
    }

    void pointsToCells( int np, Vec3d* ps, bool* ignore=0 ){
        //printf( "DEBUG pointsToCells() 1 \n"  );
        if(bResizing)resizeObjs ( np, true );
        if(obj2cell_bind==0){ Buckets::resizeObjs( np, true ); }
        //printf( "DEBUG pointsToCells() 2 obj2cell_bind %li \n", (long)obj2cell_bind  );
        if  (ignore){ for( int i=0; i<np; i++ ){ if(!ignore[i])obj2cell_bind[i] = icell( ps[i] );} }
        else        { for( int i=0; i<np; i++ ){               obj2cell_bind[i] = icell( ps[i] );} }
        //printf( "DEBUG pointsToCells() 3 np %i \n", np  );
        updateCells( np, obj2cell_bind );
        //printf( "DEBUG pointsToCells() 4 \n"  );
    }

    int getNeighbors( Vec3i ip, int* neighs ){
        int nfound = 0;
        nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y+1,ip.z+1} ), neighs+nfound );  // 1
        nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y+1,ip.z+1} ), neighs+nfound );  // 2
        nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y+1,ip.z+1} ), neighs+nfound );  // 3
        nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y  ,ip.z+1} ), neighs+nfound );  // 4
        nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y  ,ip.z+1} ), neighs+nfound );  // 5
        nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y  ,ip.z+1} ), neighs+nfound );  // 6
        nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y-1,ip.z+1} ), neighs+nfound );  // 7
        nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y-1,ip.z+1} ), neighs+nfound );  // 8
        nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y-1,ip.z+1} ), neighs+nfound );  // 9

        nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y+1,ip.z  } ), neighs+nfound );  // 10
        nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y+1,ip.z  } ), neighs+nfound );  // 11
        nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y+1,ip.z  } ), neighs+nfound );  // 12
        nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y  ,ip.z  } ), neighs+nfound );  // 13
        //nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y  ,ip.z  } ), neighs+nfound );  // 14 Allready inside
        nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y  ,ip.z  } ), neighs+nfound ); // 15
        nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y-1,ip.z  } ), neighs+nfound ); // 16
        nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y-1,ip.z  } ), neighs+nfound ); // 17
        nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y-1,ip.z  } ), neighs+nfound ); // 18

        nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y+1,ip.z-1} ), neighs+nfound ); // 19
        nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y+1,ip.z-1} ), neighs+nfound ); // 20
        nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y+1,ip.z-1} ), neighs+nfound ); // 21
        nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y  ,ip.z-1} ), neighs+nfound ); // 22
        nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y  ,ip.z-1} ), neighs+nfound ); // 23
        nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y  ,ip.z-1} ), neighs+nfound ); // 24
        nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y-1,ip.z-1} ), neighs+nfound ); // 25
        nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y-1,ip.z-1} ), neighs+nfound ); // 26
        nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y-1,ip.z-1} ), neighs+nfound ); // 27
        return nfound;
    }


    int getForwardNeighbors( Vec3i ip, int* neighs ){
        //printf("DEBUG getForwardNeighbors() neighs %li \n", (long)neighs );
        int nfound = 0;
      //nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y+1,ip.z+1} ), neighs+nfound );  // 1
        nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y-1,ip.z-1} ), neighs+nfound ); // 27
      //nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y+1,ip.z+1} ), neighs+nfound );  // 2
        nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y-1,ip.z-1} ), neighs+nfound ); // 26
      //nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y+1,ip.z+1} ), neighs+nfound );  // 3
        nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y-1,ip.z-1} ), neighs+nfound ); // 25
        //printf("DEBUG getForwardNeighbors() 3 \n");
      //nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y  ,ip.z+1} ), neighs+nfound );  // 4
        nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y  ,ip.z-1} ), neighs+nfound ); // 24
      //nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y  ,ip.z+1} ), neighs+nfound );  // 5
        nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y  ,ip.z-1} ), neighs+nfound ); // 23
      //nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y  ,ip.z+1} ), neighs+nfound );  // 6
        nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y  ,ip.z-1} ), neighs+nfound ); // 22
        //printf("DEBUG getForwardNeighbors() 6 \n");
      //nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y-1,ip.z+1} ), neighs+nfound );  // 7
        nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y+1,ip.z-1} ), neighs+nfound ); // 21
      //nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y-1,ip.z+1} ), neighs+nfound );  // 8
        nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y+1,ip.z-1} ), neighs+nfound ); // 20
      //nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y-1,ip.z+1} ), neighs+nfound );  // 9
        nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y+1,ip.z-1} ), neighs+nfound ); // 19
        //printf("DEBUG getForwardNeighbors() 9 \n");
      //nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y+1,ip.z  } ), neighs+nfound );  // 10
        nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y-1,ip.z  } ), neighs+nfound ); // 18
      //nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y+1,ip.z  } ), neighs+nfound );  // 11
        nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y-1,ip.z  } ), neighs+nfound ); // 17
      //nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y+1,ip.z  } ), neighs+nfound );  // 12
        nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y-1,ip.z  } ), neighs+nfound ); // 16
        //printf("DEBUG getForwardNeighbors() 12 \n");
      //nfound += getInCell( ixyz2i_wrap( {ip.x+1,ip.y  ,ip.z  } ), neighs+nfound );  // 13
        nfound += getInCell( ixyz2i_wrap( {ip.x-1,ip.y  ,ip.z  } ), neighs+nfound ); // 15
        //printf("DEBUG getForwardNeighbors() 13 \n");
        //nfound += getInCell( ixyz2i_wrap( {ip.x  ,ip.y  ,ip.z  } ), neighs+nfound );  // 14 Allready inside
        return nfound;
    }

};


#endif
