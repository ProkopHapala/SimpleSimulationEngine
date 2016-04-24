
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

/*
const static uint8_t wall_nodes [6][4][3] = {
    {{0,0,0},{1,0,0},{0,1,0},{1,1,0}},
    {{0,0,1},{1,0,1},{0,1,1},{1,1,1}},

    {{0,0,0},{0,0,1},{0,1,0},{0,1,1}},
    {{1,0,0},{1,0,1},{1,1,0},{1,1,1}},

    {{0,0,0},{1,0,0},{0,0,1},{1,0,1}},
    {{0,1,0},{1,1,0},{0,1,1},{1,1,1}},
};
*/

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
        1e+9,1e+9,  // strength
    };

    BondType default_BondType_edge = {
        0,          // id
        7.8,        // density
        1e+5,1e+5,  // stiffness
        1e+9,1e+9,  // strength
    };

    // =========== function implementation ( should be moved to .cpp )

    inline void index2pos( const Vec3i& index, Vec3d& pos    ){ pos  .set( index.x*scaling.x+pos0.x,    index.y*scaling.y+pos0.y,    index.z*scaling.z+pos0.z    );  };
    inline void pos2index( const Vec3d& pos,   Vec3i& index  ){ index.set( (pos.x-pos0.x)*invScaling.x, (pos.y-pos0.y)*invScaling.y, (pos.z-pos0.z)*invScaling.z );  };


    Node& insertNode( const Block& block, const uint8_t * corner ){
        uint8_t ix = block.ix + corner[0];
        uint8_t iy = block.iy + corner[1];
        uint8_t iz = block.iz + corner[2];
        int key = xyz2i( ix, iy, iz );
        int sz0 = nodes.size();
        Node& node = nodes[key]; // get valid pointer ( if new alocate, if aold take it )
        if( nodes.size() > sz0 ){ // new element
            //node.ix = ix;
            //node.iy = iy;
            //node.iz = iz;
            index2pos( {ix,iy,iz}, node.pos ); node.pos.add_mul( scaling, -0.5 );
            node.id = sz0;
            //nodeCount++;
            printf( "insert node %i (%i,%i,%i) (%3.3f,%3.3f,%3.3f) \n", node.id, ix,iy, iz, node.pos.x, node.pos.y, node.pos.z );
        }else{
            printf( "found  node %i (%i,%i,%i) (%3.3f,%3.3f,%3.3f) \n", node.id, ix,iy, iz, node.pos.x, node.pos.y, node.pos.z );
        }
        return node;
    }

    Bond& insertBond( uint16_t i, uint16_t j, double l0, const BondType& type ){
        uint32_t key = (i<<16) + j;
        int sz0 = bonds.size();
        Bond& bond = bonds[key];   // get valid pointer ( if new alocate, if aold take it )
        if( bonds.size() > sz0 ){  // new element
            bond.i    = i;
            bond.j    = j;
            bond.id   = sz0;
            bond.type = type;
            bond.l0   = l0;
            printf( "insert bond (%i,%i) %i %i \n", i, j, key, bond.id );
        }else{
            if( type.sPress > bond.type.sPress ){
                bond.type = type;
            }
            printf( "found  bond (%i,%i) %i %i \n", i, j, key, bond.id );
        }
        return bond;
    }

    void block2truss( const Block& block ){
        for(int iSide=0; iSide<6; iSide++){
            int type = block.sides[iSide];
            if( type < nMaxTypes ){
                Node& nod00 = insertNode( block, wall_nodes[iSide][0] );
                Node& nod01 = insertNode( block, wall_nodes[iSide][1] );
                Node& nod10 = insertNode( block, wall_nodes[iSide][2] );
                Node& nod11 = insertNode( block, wall_nodes[iSide][3] );
                printf( " node.ids %i %i %i %i \n", nod00.id, nod01.id, nod10.id, nod11.id );

                double ldiag = 1.41421356237;
                double ledge = 1.0d;
                insertBond( nod00.id, nod11.id, ldiag, wallTypes[type].diag );
                insertBond( nod10.id, nod01.id, ldiag, wallTypes[type].diag );

                insertBond( nod00.id, nod01.id, ledge, wallTypes[type].edge );
                insertBond( nod00.id, nod10.id, ledge, wallTypes[type].edge );
                insertBond( nod11.id, nod01.id, ledge, wallTypes[type].edge );
                insertBond( nod11.id, nod10.id, ledge, wallTypes[type].edge );

            }
        }
    }

    void blocks2truss( ){
        printf( " blocks2truss 0 \n" );
        nodes.clear();  nodes.reserve ( 1024    );
        bonds.clear();  bonds.reserve ( 16*1024 );
        for( int i=0; i<nBlocks; i++ ){
            block2truss( blocks[i] );
        }
        printf( " blocks2truss 1 \n" );
        truss.allocate( nodes.size(), bonds.size(), 3, NULL, NULL, NULL, NULL );
        printf( " blocks2truss 2 \n" );
        for( auto it : nodes ){
            Node& node = it.second;
            truss.points[ node.id ] = node.pos; // do we need this at all ?
        }
        for( auto it : bonds ){
            Bond& bond = it.second;
            truss.bonds[ bond.id ] = bond;
        }
        truss.prepareBonds ( false            );
        truss.preparePoints( true, -1.0, 0.0 );

        for( int i=0; i<truss.nfix; i++){
            truss.fix[i] = i;   // FIXME ... this is just some stupid DEBUG default
        }
    }


    void setDefaultWallType(){
        for( int i=0;  i<nMaxTypes; i++ ){
            wallTypes[i].shape = 0;
            wallTypes[i].diag  = default_BondType_diag;
            wallTypes[i].edge  = default_BondType_edge;
            wallTypes[i].mass  = 1.0;
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

        setDefaultWallType();
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
