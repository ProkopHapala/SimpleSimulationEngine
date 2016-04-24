
#ifndef BlockHouseWorld_h
#define BlockHouseWorld_h

#include <unordered_map>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

#include "SoftBody.h" // dynamics

class Node{
    public:
    uint16_t  id;    // unique indentifier
    //uint8_t   ix,iy,iz;  // grid point index
    Vec3d     pos;   //
    //double    mass;
};

class WallType{
    public:
	int shape;
	double mass;
    BondType edge,diag;
};

class Block{
    public:
    uint8_t ix,iy,iz;
    uint8_t content;
	uint8_t sides[6];

    inline bool isEmpty(){
        return 255 == sides[0] == sides[1] == sides[2] == sides[3] == sides[4] == sides[5] == content;
    }

    inline void setEmpty(){
        sides[0] =255; sides[1]=255; sides[2]=255; sides[3]=255; sides[4]=255; sides[5]=255; content=255;
    }
};


// FIXME : I'm not sure why the points are inverted
const static uint8_t wall_nodes [6][4][3] = {
    {{1,1,1},{0,1,1},{1,0,1},{0,0,1}},
    {{1,1,0},{0,1,0},{1,0,0},{0,0,0}},

    {{1,1,1},{1,1,0},{1,0,1},{1,0,0}},
    {{0,1,1},{0,1,0},{0,0,1},{0,0,0}},

    {{1,1,1},{0,1,1},{1,1,0},{0,1,0}},
    {{1,0,1},{0,0,1},{1,0,0},{0,0,0}},
};

inline uint32_t xyz2i( uint8_t ix, uint8_t iy, uint8_t iz ){
    return ix | ( iy << 8 ) | ( iz << 16 );
}

inline void i2xyz( uint32_t i,  uint8_t& ix, uint8_t& iy, uint8_t& iz ){
    iz=( i & 0xF00 ) >> 16;
    iy=( i & 0x0F0 ) >> 8;
    ix=( i & 0x00F );
}

// =====================
//    BlockHouseWorld
// =====================

class BlockHouseWorld{
	public:
    const static int nMaxBlocks = 1024;
    int nBlocks = 0;
    int iBlock  = 0;
	Block blocks[nMaxBlocks];

	const static int nMaxTypes = 255;
	int nTypes = 0;
	int iTypes = 0;
    WallType wallTypes[nMaxTypes];

    Vec3i nMax;
    Vec3d scaling;
    Vec3d invScaling;
    Vec3d pos0;
    Mat3d rotations[6];

    SoftBody truss;
    //int nodeCount = 0,bondCount=0;
    std::unordered_map<uint32_t,Node> nodes;
    std::unordered_map<uint32_t,Bond> bonds;


    BondType default_BondType_diag = {
        0,          // id
        7.8,        // density
        1e+5,1e+5,  // stiffness
        1e+9,1e+9  // strength
    };

    BondType default_BondType_edge = {
        0,          // id
        7.8,        // density
        1e+5,1e+5,  // stiffness
        1e+9,1e+9  // strength
    };

    // =========== function declaration

    Node& insertNode( const Block& block, const uint8_t * corner );
    Bond& insertBond( uint16_t i, uint16_t j, double l0, const BondType& type );
    void  block2truss( const Block& block );
    void  blocks2truss( );

    void setDefaultWallType();
    int  init();

    int findBlock        ( uint8_t ix, uint8_t iy, uint8_t iz );
    int changeBlock      ( uint8_t ix, uint8_t iy, uint8_t iz, uint8_t iSide, uint8_t type );
    int eraseBlock       ( uint8_t ix, uint8_t iy, uint8_t iz );
    int defragmentBlocks ( );

    // =========== inline functions

    inline void index2pos( const Vec3i& index, Vec3d& pos    ){ pos  .set( index.x*scaling.x+pos0.x,    index.y*scaling.y+pos0.y,    index.z*scaling.z+pos0.z    );  };
    inline void pos2index( const Vec3d& pos,   Vec3i& index  ){ index.set( (pos.x-pos0.x)*invScaling.x, (pos.y-pos0.y)*invScaling.y, (pos.z-pos0.z)*invScaling.z );  };

};

#endif
