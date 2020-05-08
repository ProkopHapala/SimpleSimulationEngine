#ifndef Grid2DAlgs_h
#define Grid2DAlgs_h

#include <cstring> // memcpy

#include "Vec2.h"
#include "GridIndex2D.h"
//#include "Grid2DAlgs.cpp"


//    5   3   7
//    0       1
//    4   2   6
//                                           0     1    2     3     4        5          6           7
//static const int    SquareNeighs       []={-1,0, 1,0, 0,-1, 0,1, -1,-1,    -1,+1,     +1,-1,     +1,+1     };
//static const double SquareNeighsDist   []={ 1.0, 1.0,  1.0, 1.0, M_SQRT2,   M_SQRT2,   M_SQRT2,   M_SQRT2  };
//static const double SquareNeighsDistInv[]={ 1.0, 1.0,  1.0, 1.0, M_SQRT1_2, M_SQRT1_2, M_SQRT1_2, M_SQRT1_2};

//    6   5   4
//    7       3
//    0   1   2
//                                             0         1        2          3     4          5       6          7
static const int    SquareNeighs       []={ -1,-1,      0,-1,   1,-1,      1,0,   1,1,       0,+1,  -1,+1,     -1,0 };
static const double SquareNeighsDist   []={  M_SQRT2,   1.0,   M_SQRT2,    1.0,   M_SQRT2,   1.0,   M_SQRT2  ,  1.0 };
static const double SquareNeighsDistInv[]={  M_SQRT1_2, 1.0,   M_SQRT1_2,  1.0,   M_SQRT1_2, 1.0,   M_SQRT1_2,  1.0 };


//void bisecNoise( Vec2i n, double * hs, double frndMin, double frndMax );
void bisecNoise      ( int npow, double * hs, double frndMin, double frndMax );
void bisecPaternNoise( int npow, double * hs, double frndMin, double frndMax );

class Grid2DAlg : public GridIndex2D{ public:
    static const int nNeighMax = 8;
    int nneigh    = 0;
    Vec2i  neighs        [nNeighMax];
    double neigh_dist    [nNeighMax];
    double neigh_invDist [nNeighMax];

    void set( const Grid2DAlg& rf ){
        n   =rf.n;
        ntot=rf.ntot;
        nneigh  = rf.nneigh; printf( "nneigh %i %i \n" , nneigh, rf.nneigh );
        for( int i=0; i<nneigh; i++ ){
            neighs       [i] = rf.neighs        [i];
            neigh_dist   [i] = rf.neigh_dist    [i];
            neigh_invDist[i] = rf.neigh_invDist [i];
        }
    }

    void initNeighsSquareN( int n ){ // typically 4 or 8
        nneigh = n;
        memcpy(neighs,       SquareNeighs       ,2*nneigh*sizeof(int)   );
        memcpy(neigh_dist,   SquareNeighsDist   ,  nneigh*sizeof(double));
        memcpy(neigh_invDist,SquareNeighsDistInv,  nneigh*sizeof(double));
    }

    void initNeighs_6(bool switched){
        nneigh = 6;
        //int  sela[]={0,1,2,3,4,7};
        //int  selb[]={0,1,2,3,5,6};
        int  sela[]={0,1,3,4,5,7};
        int  selb[]={1,2,3,5,6,7};
        int* sel;
        if(switched){sel=sela;}else{sel=selb;}
        for(int i=0; i<nneigh; i++){
            int j=sel[i];
            neighs       [i]=((Vec2i*)SquareNeighs)[j];
            neigh_dist   [i]=1.0;
            neigh_invDist[i]=1.0;
        }
    }
    void initNeighsSquareMask( uint8_t mask ){
        uint8_t bit =1;
        nneigh=0;
        for(int i=0; i<8; i++){
            if(mask&bit){
                neighs       [nneigh]=((Vec2i*)SquareNeighs)[i];
                neigh_dist   [nneigh]=SquareNeighsDist[i];
                neigh_invDist[nneigh]=SquareNeighsDistInv[i];
                nneigh++;
            }
            bit<<=1;
        }
    }

};

/*
class PathFindingGrid2D :public Grid2DAlg{ public:

};
*/

#endif
