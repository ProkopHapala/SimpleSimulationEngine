
#ifndef  TrussBuilder_h
#define  TrussBuilder_h

#include <unordered_map>

#include "fastmath.h"
//#include "Vec2.h"
#include "Vec3.h"
//#include "Mat3.h"
#include "quaternion.h"

class NodeIndex{
    public:
    union{
        struct{ int_fast16_t x,y,z,w; };
        //struct{ int_fast16_t a,b,c,d;     };
        int_fast16_t array[4];
        int_fast64_t i;
    };
};

class BondIndex{
    public:
    union{
        struct{ int_fast32_t a,b; };
        int_fast32_t array[2];
        int_fast64_t i;
    };
};

class Node{
    public:
    int_fast32_t  id;    // unique indentifier
    //uint8_t   ix,iy,iz;  // grid point index
    //uint8_t
    Vec3d     pos;
    //double    mass;
};

// ===== point index conversion

/*
inline uint32_t xyz2i( int_fast16_t ix, int_fast16_t iy, int_fast16_t iz ){
    return ix | ( iy << 8 ) | ( iz << 16 );
}

inline void i2xyz( int_fast64_t i,  int_fast16_t& ix, int_fast16_t& iy, int_fast16_t& iz ){
    iz=( i & 0xF00 ) >> 32;
    iy=( i & 0x0F0 ) >> 16;
    ix=( i & 0x00F );
}

// ===== point index conversion

inline uint32_t ij2ib( int_fast16_t ix, int_fast32_t iy, int_fast32_t iz ){
    return ix | ( iy << 8 ) | ( iz << 16 );
}

inline uint32_t ib2ij( int_fast16_t ix, int_fast16_t iy, int_fast16_t iz ){
    return ix | ( iy << 8 ) | ( iz << 16 );
}
*/

// ====================
// ==== TrussBuilder
// ====================

class TrussBuilder{

    Vec3i nMax;
    Vec3d scaling;
    Vec3d invScaling;
    Vec3d pos0;

    std::unordered_map<int_fast64_t,Node> nodes;
    std::unordered_map<int_fast64_t,Bond> bonds;

    Node& insertNode( int_fast16_t ix, int_fast16_t iy, int_fast16_t iz );
    Bond& insertBond( int_fast32_t i, int_fast32_t j, double l0, const BondType& type );
    void  toSoftBody( SoftBody& truss );

    // =========== inline functions

    inline void index2pos( const Vec3i& index, Vec3d& pos    ){ pos  .set( index.x*scaling.x+pos0.x,    index.y*scaling.y+pos0.y,    index.z*scaling.z+pos0.z    );  };
    inline void pos2index( const Vec3d& pos,   Vec3i& index  ){ index.set( (pos.x-pos0.x)*invScaling.x, (pos.y-pos0.y)*invScaling.y, (pos.z-pos0.z)*invScaling.z );  };

};

#endif

