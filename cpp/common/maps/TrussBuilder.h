
#ifndef  TrussBuilder_h
#define  TrussBuilder_h

#include <vector>
#include <unordered_map>

#include "fastmath.h"
//#include "Vec2.h"
#include "Vec3.h"
//#include "Mat3.h"
#include "quaternion.h"

#include "SoftBody.h"

class GridNode{
    public:
    int_fast16_t  ix,iy,iz;
    bool exist,fixed;
    //int_fast32_t  id;    // unique indentifier
    //uint8_t   ix,iy,iz;  // grid point index
    //uint8_t
    //Vec3d     pos;
    //double    mass;
};

// ===== point index conversion


inline int_fast64_t xyz2id( int_fast16_t ix, int_fast16_t iy, int_fast16_t iz ){
    return ix | ( ((int_fast64_t)iy) << 16 ) | ( ((int_fast64_t)iz) << 32 );
}

inline void id2xyz( int_fast64_t id,  int_fast16_t& ix, int_fast16_t& iy, int_fast16_t& iz ){
    iz=( id & 0xFF0000 ) >> 32;
    iy=( id & 0x00FF00 ) >> 16;
    ix=( id & 0x0000FF );
}

// ===== point index conversion

inline int_fast64_t xy2id( int_fast32_t ix, int_fast32_t iy ){
    return (int_fast64_t)ix | ( ((int_fast64_t)iy) << 32 );
    //return ((int_fast64_t)iy) << 32;
    //return iy << 32;
}

inline int_fast64_t id2xy( int_fast64_t id, int_fast32_t ix, int_fast32_t iy ){
    iy=( id & 0xFFFF0000 ) >> 32;
    ix=( id & 0x0000FFFF );
}


// ====================
// ==== TrussBuilder
// ====================

class TrussBuilder{
    public:

    //static Vec3i nMax;
    constexpr static int          nMax = 1<<16;
    constexpr static int_fast16_t ioff = 1<<15;

    Vec3d  scaling;
    Vec3d  invScaling;
    Vec3d  pos0;

    int nfixed = 0;

    std::vector<BondType> bondTypes;

    std::unordered_map<int_fast64_t,Bond>         bonds;
    std::unordered_map<int_fast64_t,int_fast32_t> nodeIs;
    std::vector<GridNode>                         nodes;

    // ======== function declarations

    //GridNode& insertNode( int_fast16_t ix, int_fast16_t iy, int_fast16_t iz );
    int_fast32_t insertNode  ( int_fast16_t ix, int_fast16_t iy, int_fast16_t iz );
    int_fast32_t getNodeIndex( int_fast16_t ix, int_fast16_t iy, int_fast16_t iz );
    bool         removeNode  ( int_fast16_t ix, int_fast16_t iy, int_fast16_t iz );

    //Bond& insertBond( int_fast32_t i, int_fast32_t j, double l0, const BondType& type );
    //Bond& insertBond( int_fast16_t ix0, int_fast16_t iy0, int_fast16_t iz0, int_fast16_t ix1, int_fast16_t iy1, int_fast16_t iz1, double l0, const BondType& type );
    //Bond& insertBond( int_fast16_t ix0, int_fast16_t iy0, int_fast16_t iz0, int_fast16_t ix1, int_fast16_t iy1, int_fast16_t iz1,            const BondType& type );

    Bond& insertBond( int_fast32_t i, int_fast32_t j, double l0, BondType* type );
    Bond& insertBond( int_fast16_t ix0, int_fast16_t iy0, int_fast16_t iz0, int_fast16_t ix1, int_fast16_t iy1, int_fast16_t iz1, double l0, BondType* type );
    Bond& insertBond( int_fast16_t ix0, int_fast16_t iy0, int_fast16_t iz0, int_fast16_t ix1, int_fast16_t iy1, int_fast16_t iz1,            BondType* type );

    void insertBox( int_fast16_t ix0, int_fast16_t iy0, int_fast16_t iz0, uint8_t mask, BondType* type );

    int_fast64_t getBondKey( int_fast16_t ix0, int_fast16_t iy0, int_fast16_t iz0, int_fast16_t ix1, int_fast16_t iy1, int_fast16_t iz1 );
    bool         removeBond( int_fast16_t ix0, int_fast16_t iy0, int_fast16_t iz0, int_fast16_t ix1, int_fast16_t iy1, int_fast16_t iz1 );
    bool         removeBond( int_fast32_t i,   int_fast32_t j );

    bool removeNodesWithoutBond( );

    void  init( int nNodesGuess, int nBondsGuess, int nBondsTypesGuess );
    void  toSoftBody( SoftBody& truss );

    void toFile  ( char * fname );
    void fromFile( char * fname );

    // =========== inline functions

    inline void index2pos( const Vec3i& index, Vec3d& pos    ){ pos  .set( index.x*scaling.x+pos0.x,    index.y*scaling.y+pos0.y,    index.z*scaling.z+pos0.z    );  };
    inline void pos2index( const Vec3d& pos,   Vec3i& index  ){ index.set( (pos.x-pos0.x)*invScaling.x, (pos.y-pos0.y)*invScaling.y, (pos.z-pos0.z)*invScaling.z );  };

    inline double dist2( const Vec3i& ip1, const Vec3i& ip2  ){
        Vec3d p1,p2,dp;
        index2pos( ip1, p1 );
        index2pos( ip2, p2 );
        dp.set_sub(p2,p1);
        return dp.norm2();
    };

    void  setScaling( const Vec3d& scaling_ ){
        scaling = scaling_; invScaling.set( 1/scaling.x, 1/scaling.y, 1/scaling.z );
    };

};

#endif

