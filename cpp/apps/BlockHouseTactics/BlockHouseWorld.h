
#ifndef BlockHouseWorld_h
#define BlockHouseWorld_h

#include <unordered_map>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

#include "SoftBody.h" // dynamics


const static int nMaxTypes = 255;

class Node{
    int    id;
    double mass;
};

class WallType{
    public:
	int shape;
	double mass;
	double stiffness;
	double strength;
};

class Block{
    public:
    uint8_t ix,iy,iz;
    uint8_t content;
	uint8_t sides[6];
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
        return 255 == sides[0] == sides[1] == sides[2] == sides[3] == sides[4] == sides[5] == content;
    }

    inline void setEmpty(){
        sides[0] =255; sides[1]=255; sides[2]=255; sides[3]=255; sides[4]=255; sides[5]=255; content=255;
    }
};

const static uint8_t wall_nodes [6][4][3] = {
    {{0,0,0},{1,0,0},{0,1,0},{1,1,0}},
    {{0,0,1},{1,0,1},{0,1,1},{1,1,1}},

    {{0,0,0},{0,0,1},{0,1,0},{0,1,1}},
    {{1,0,0},{1,0,1},{1,1,0},{1,1,1}},

    {{0,0,0},{1,0,0},{0,0,1},{1,0,1}},
    {{0,1,0},{1,1,0},{0,1,1},{1,1,1}},
};

inline uint32_t xyz2i( uint8_t ix, uint8_t iy, uint8_t iz ){
    return ix | ( iy << 8 ) | ( iz << 16 );
}

inline void i2xyz( uint32_t i,  uint8_t& ix, uint8_t& iy, uint8_t& iz ){
    iz=( i & 0xF00 ) >> 16;
    iy=( i & 0x0F0 ) >> 8;
    ix=( i & 0x00F );
}

class BlockHouseWorld{
	public:
    const static int nMaxBlocks = 1024;
    int nBlocks = 0;
    int iBlock  = 0;
	Block blocks[nMaxBlocks];

	int nTypes = 0;
	int iTypes = 0;
    WallType wallTypes[nMaxTypes];

    Vec3i nMax;
    Vec3d scaling;
    Vec3d invScaling;
    Vec3d pos0;
    Mat3d rotations[6];

    SoftBody truss;
    int nodeCount = 0,bondCount=0;
    std::unordered_map<uint32_t,int> nodes;
    std::unordered_map<uint32_t,int> bonds;

    // =========== function implementation ( should be moved to .cpp )

    inline void index2pos( const Vec3i& index, Vec3d& pos    ){ pos  .set( index.x*scaling.x+pos0.x,    index.y*scaling.y+pos0.y,    index.z*scaling.z+pos0.z    );  };
    inline void pos2index( const Vec3d& pos,   Vec3i& index  ){ index.set( (pos.x-pos0.x)*invScaling.x, (pos.y-pos0.y)*invScaling.y, (pos.z-pos0.z)*invScaling.z );  };

    int getValidNode( uint32_t key ){
        auto it = nodes.find( key );
        if ( it == nodes.end() ){
            int inode = nodeCount;
            nodes.insert( {key, nodeCount} );
            nodeCount++;
            return inode;
        }else{
            return it->second;
        };
    }

    int getValidBond(  int i, int j ){
        uint32_t key = i<<16 + j;
        auto it = bonds.find( key );
        if ( it == nodes.end() ){
            int ibond = bondCount;
            bonds.insert( {key, ibond} );
            bondCount++;
            return ibond;
        }else{
            return it->second;
        };
    }

    int insertNode( const Block& block, const uint8_t * corner ){
        int inod = getValidNode( xyz2i( block.ix + corner[0], block.iy + corner[1], block.iz + corner[2] ) );
        // nodes[inod].mass +=  ???
        // points nodes[ ] +=
        return inod;
    }

    int insertBond( int i, int j, const WallType& wtype ){
        int ibond = getValidBond( i, j );
        //insertBond( ib, i, j, wtype->kTens_diag, wtype->kPres_diag, wtype->l0_diag );
        return ibond;
    }

    void block2truss( const Block& block ){
        for(int iSide=0; iSide<6; iSide++){
            int type = block.sides[iSide];
            if( type < nMaxTypes ){
                int nod00 = insertNode( block, wall_nodes[iSide][0] );
                int nod01 = insertNode( block, wall_nodes[iSide][1] );
                int nod10 = insertNode( block, wall_nodes[iSide][2] );
                int nod11 = insertNode( block, wall_nodes[iSide][3] );

                insertBond( nod00, nod11, wallTypes[type] );
                insertBond( nod10, nod01, wallTypes[type] );

                insertBond( nod00, nod01, wallTypes[type] );
                insertBond( nod00, nod10, wallTypes[type] );
                insertBond( nod11, nod01, wallTypes[type] );
                insertBond( nod11, nod10, wallTypes[type] );

            }
        }
    }


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
        //pos0   .set( {0.0d,0.0d,-0.0d} );
        scaling.set( {1.0,1.0,1.0} );
    }

    int findBlock ( uint8_t ix, uint8_t iy, uint8_t iz ){
        for( int i=0; i<nBlocks; i++ ){
            if( ( blocks[i].ix == ix ) && ( blocks[i].iy == iy ) && ( blocks[i].iz == iz ) ){
                return i;
            };
        }
        return -1;
    };


    int changeBlock( uint8_t ix, uint8_t iy, uint8_t iz, uint8_t iSide, uint8_t type ){
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

    int eraseBlock( uint8_t ix, uint8_t iy, uint8_t iz ){
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
