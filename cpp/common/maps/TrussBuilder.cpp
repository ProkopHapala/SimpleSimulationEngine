

#include "SoftBody.h"
#include "TrussBuilder.h" // THE HEADER

Node& TrussBuilder::insertNode( int_fast16_t ix, int_fast16_t iy, int_fast16_t iz ){
    NodeIndex inod; inod.x=ix; inod.y=iy; inod.z=iz;
    int sz0 = nodes.size();
    Node& node = nodes[inod.i]; // get valid pointer ( if new alocate, if aold take it )
    if( nodes.size() > sz0 ){ // new element
        //node.ix = ix;
        //node.iy = iy;
        //node.iz = iz;
        index2pos( {ix,iy,iz}, node.pos ); node.pos.add_mul( scaling, -0.5 );
        node.id = sz0;
        //nodeCount++;
        //printf( "insert node %i (%i,%i,%i) (%3.3f,%3.3f,%3.3f) \n", node.id, ix,iy, iz, node.pos.x, node.pos.y, node.pos.z );
    }else{
        //printf( "found  node %i (%i,%i,%i) (%3.3f,%3.3f,%3.3f) \n", node.id, ix,iy, iz, node.pos.x, node.pos.y, node.pos.z );
    }
    return node;
}

Bond& TrussBuilder::insertBond( int_fast32_t i, int_fast32_t j, double l0, const BondType& type ){
    BondIndex ib;
    int sz0 = bonds.size();
    Bond& bond = bonds[ib.i];   // get valid pointer ( if new alocate, if aold take it )
    if( bonds.size() > sz0 ){  // new element
        bond.i    = i;   // we actually don't need this
        bond.j    = j;
        bond.id   = sz0;
        bond.type = type;
        bond.l0   = l0;
        //printf( "insert bond (%i,%i) %i %i \n", i, j, key, bond.id );
    }else{
        if( type.sPress > bond.type.sPress ){
            bond.type = type;
        }
        //printf( "found  bond (%i,%i) %i %i \n", i, j, key, bond.id );
    }
    return bond;
}

void TrussBuilder::toSoftBody( SoftBody& truss ){
    truss.allocate( nodes.size(), bonds.size(), 3, NULL, NULL, NULL, NULL );
    for( auto it : nodes ){
        Node& node = it.second;
        truss.points[ node.id ] = node.pos; // do we need this at all ?
    }
    for( auto it : bonds ){
        Bond& bond = it.second;
        truss.bonds[ bond.id ] = bond;
    }
}

