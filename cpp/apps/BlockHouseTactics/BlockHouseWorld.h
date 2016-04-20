
#ifndef BlockHouseWorld_h
#define BlockHouseWorld_h

#include <unordered_set>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

//#include "SoftBody.h" // dynamics

class WallType{
    public:
	int shape;
	double mass;
	double stiffness;
	double strength;
};

class Block{
    public:
    char ix,iy,iz;
    char content;
	char sides[6];
/*
    Uint8  roof;
    Uint8  roof;
    Uint8  floor;
    Uint8  wall_px;
    Uint8  wall_mx;
    Uint8  wall_py;
    Uint8  wall_my;
*/

    inline bool isEmpty(){
        return 0 == sides[0] == sides[1] == sides[2] == sides[3] == sides[4] == sides[5] == content;
    }

    inline void setEmpty(){
        sides[0] =0; sides[1]=0; sides[2]=0; sides[3]=0; sides[4]=0; sides[5]=0; content=0;
    }
};

class BlockHouseWorld{
	public:
    const static int nMaxBlocks = 1024;
    int nBlocks = 0;
    int iBlock  = 0;
	Block blocks[nMaxBlocks];

    const static int nMaxTypes = 256;
	int nTypes = 0;
	int iTypes = 0;
    WallType wallTypes[nMaxTypes];

    Vec3i nMax;
    Vec3d scaling;
    Vec3d invScaling;
    Vec3d pos0;
    Mat3d rotations[6];

    // =========== function implementation ( should be moved to .cpp )

    inline void index2pos( const Vec3i& index, Vec3d& pos    ){ pos  .set( index.x*scaling.x+pos0.x,    index.y*scaling.y+pos0.y,    index.z*scaling.z+pos0.z    );  };
    inline void pos2index( const Vec3d& pos,   Vec3i& index  ){ index.set( (pos.x-pos0.x)*invScaling.x, (pos.y-pos0.y)*invScaling.y, (pos.z-pos0.z)*invScaling.z );  };

/*
    void buildRotations( const Mat3d& rot ){
        rotations[0].set( rot.a, rot.b, rot.c );
        rotations[1].set( rot.a, rot.c, rot.b );
        rotations[2].set( rot.b, rot.a, rot.c );
        rotations[3].set( rot.b, rot.c, rot.a );
        rotations[4].set( rot.c, rot.a, rot.b );
        rotations[5].set( rot.c, rot.b, rot.a );
    }
*/

    int init(){
        /*
        rotations[0].set( {0.0,0.0,+1.0}, {+1.0,0.0,0.0}, {0.0,+1.0,0.0} ); // roof
        rotations[1].set( {0.0,0.0,-1.0}, {+1.0,0.0,0.0}, {0.0,+1.0,0.0} ); // floor
        rotations[2].set( {+1.0,0.0,0.0}, {0.0,0.0,+1.0}, {0.0,+1.0,0.0} ); // wall
        rotations[3].set( {-1.0,0.0,0.0}, {0.0,0.0,+1.0}, {0.0,+1.0,0.0} );
        rotations[4].set( {0.0,+1.0,0.0}, {0.0,0.0,+1.0}, {+1.0,0.0,0.0} );
        rotations[5].set( {0.0,-1.0,0.0}, {0.0,0.0,+1.0}, {+1.0,0.0,0.0} );
        */

        rotations[0].set(  {+1.0,0.0,0.0}, {0.0,+1.0,0.0}, {0.0,0.0,+1.0} ); // roof
        rotations[1].set(  {+1.0,0.0,0.0}, {0.0,+1.0,0.0}, {0.0,0.0,-1.0} ); // floor
        rotations[2].set(  {0.0,0.0,+1.0}, {0.0,+1.0,0.0}, {+1.0,0.0,0.0} ); // wall
        rotations[3].set(  {0.0,0.0,+1.0}, {0.0,+1.0,0.0}, {-1.0,0.0,0.0} );
        rotations[4].set(  {0.0,0.0,+1.0}, {+1.0,0.0,0.0}, {0.0,+1.0,0.0} );
        rotations[5].set(  {0.0,0.0,+1.0}, {+1.0,0.0,0.0}, {0.0,-1.0,0.0} );

        nMax   .set(256,256,256);
        pos0   .set( {-127.0d,-127.0d,-127.0d} );
        scaling.set( {1.0,1.0,1.0} );
    }

    int findBlock ( char ix, char iy, char iz ){
        for( int i=0; i<nBlocks; i++ ){
            if( ( blocks[i].ix == ix ) && ( blocks[i].iy == iy ) && ( blocks[i].iz == iz ) ){
                return i;
            };
        }
        return -1;
    };

    int changeBlock( char ix, char iy, char iz, char iSide, char type ){
        int i = findBlock ( ix, iy, iz );
        if( i < 0 ){
            if( nBlocks >= nMaxBlocks ) return -1;
            i=nBlocks;
            nBlocks++;
            blocks[i].ix = ix;
            blocks[i].iy = iy;
            blocks[i].iz = iz;
            blocks[i].setEmpty();
        };
        blocks[i].sides[iSide] = type;
        return i;
    };

    int eraseBlock( char ix, char iy, char iz ){
        int i = findBlock ( ix, iy, iz );
        if( i > 0 ){
            blocks[i].setEmpty();
        };
        return i;
    };

    int defragmentBlocks ( ){
        int nFound = 0;
        int nEmpty = 0;
        int i      = nBlocks;
        int j      = 0;
        while( i>j ){
            if( blocks[i].isEmpty() ){ nEmpty++; i--; continue; }
            int k;
            for( k = j; k<i; k++ ){
                if( blocks[k].isEmpty() ){
                    blocks[k] = blocks[i];
                    i--;
                    nEmpty++;
                    break;
                }
                nFound++;
            }
            nFound++;
            j = k+1;
        }
        nBlocks = nFound;
        return nEmpty;
    }


};

#endif
