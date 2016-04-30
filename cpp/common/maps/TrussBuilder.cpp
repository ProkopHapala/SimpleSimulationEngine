
#include <iostream>
#include <fstream>

#include "TrussBuilder.h" // THE HEADER

// ============= editor functions

int_fast32_t TrussBuilder::insertNode( int_fast16_t ix, int_fast16_t iy, int_fast16_t iz ){
    //ID64 index; index.x=ix; index.y=iy; index.z=iz;
    int_fast64_t key = xyz2id( ix, iy, iz );
    int sz0 = nodeIs.size();
    int_fast32_t& index = nodeIs[key]; // get valid pointer ( if new alocate, if aold take it )
    if( nodeIs.size() > sz0 ){ // new element
        //GridNode& node = nodes[index];
        //node.ix = ix;
        //node.iy = iy;
        //node.iz = iz;
        //nodes.push_back(GridNode());
        index = sz0;
        nodes.push_back({ix,iy,iz,true,false});
        //index2pos( {ix,iy,iz}, nodes[id].pos ); node.pos.add_mul( scaling, -0.5 );
        //node.id = sz0;
        //nodeCount++;
        printf( "insert node %li %i (%i,%i,%i) %i (%i,%i,%i)\n", key, index, ix,iy, iz,  nodes.size(), nodes[nodes.size()-1].ix, nodes[nodes.size()-1].iy, nodes[nodes.size()-1].iz );
        //printf( "insert node %i %i (%i,%i,%i) \n", id, index.id, index.x, index.y, index.z );
    }else{
        printf( "found  node %li %i (%i,%i,%i) \n", key, index, ix,iy, iz );
        //printf( "found  node %i %i (%i,%i,%i) \n", id, index.id, index.x, index.y, index.z );
    }
    return index;
}

int_fast32_t TrussBuilder::getNodeIndex( int_fast16_t ix, int_fast16_t iy, int_fast16_t iz ){
    int_fast64_t key = xyz2id( ix, iy, iz );
    auto it = nodeIs.find( key );
    if( it == nodeIs.end() ) return -1;
    return it->second;
}

bool TrussBuilder::removeNode( int_fast16_t ix, int_fast16_t iy, int_fast16_t iz ){
    int_fast64_t key = xyz2id( ix, iy, iz );
    auto it = nodeIs.find( key );
    if( it == nodeIs.end() )return false;
    nodes[it->second].exist = false;
    nodeIs.erase(it);
    return true;
}

Bond& TrussBuilder::insertBond( int_fast32_t i, int_fast32_t j, double l0, const BondType& type ){
    //ID64 ib;  ib.a=i; ib.b=j;
    if( i > j ) { int_fast32_t sw = i; i=j; j=sw; };  // permutations symmetry
    int_fast64_t key = xy2id( i, j );
    int sz0 = bonds.size();
    Bond& bond = bonds[key];   // get valid pointer ( if new alocate, if aold take it )
    if( bonds.size() > sz0 ){  // new element
        bond.i    = i;   // we actually don't need this
        bond.j    = j;
        bond.id   = sz0;
        bond.type = type;
        bond.l0   = l0;
        printf( "insert bond (%i,%i) %li %i \n", i, j, key, bond.id );
    }else{
        if( type.sPress > bond.type.sPress ){
            bond.type = type;
        }
        printf( "found  bond (%i,%i) %li %i \n", i, j, key, bond.id );
    }
    return bond;
}

Bond& TrussBuilder::insertBond( int_fast16_t ix0, int_fast16_t iy0, int_fast16_t iz0, int_fast16_t ix1, int_fast16_t iy1, int_fast16_t iz1, double l0, const BondType& type ){
    int_fast32_t id1 = insertNode( ix0, iy0, iz0 );
    int_fast32_t id2 = insertNode( ix1, iy1, iz1 );
    return insertBond( id1, id2, l0, type );
}

Bond& TrussBuilder::insertBond( int_fast16_t ix0, int_fast16_t iy0, int_fast16_t iz0, int_fast16_t ix1, int_fast16_t iy1, int_fast16_t iz1, const BondType& type ){
    double l0  = dist2( {ix0, iy0, iz0}, { ix1,  iy1,  iz1}  );
    return insertBond( ix0, iy0, iz0, ix1, iy1, iz1, l0, type );
}

bool TrussBuilder::removeBond( int_fast32_t i, int_fast32_t j ){
    //if( i > j ) { int_fast32_t sw = i; i=j; j=sw; };  // To be sure we canput it also here
    int_fast64_t key = xy2id( i, j );
    return bonds.erase ( key ) > 0;
}

int_fast64_t TrussBuilder::getBondKey( int_fast16_t ix0, int_fast16_t iy0, int_fast16_t iz0, int_fast16_t ix1, int_fast16_t iy1, int_fast16_t iz1 ){
    int_fast32_t i = getNodeIndex( ix0, iy0, iz0 );
    int_fast32_t j = getNodeIndex( ix1, iy1, iz1 );
    if( (i<0) || (j<0) ) return -1;
    if( i > j ) { int_fast32_t sw = i; i=j; j=sw; };
    return xy2id( i, j );
}

bool TrussBuilder::removeBond( int_fast16_t ix0, int_fast16_t iy0, int_fast16_t iz0, int_fast16_t ix1, int_fast16_t iy1, int_fast16_t iz1 ){
    int_fast64_t key = getBondKey( ix0, iy0, iz0, ix1, iy1, iz1 );
    if( key > 0 ){
        return bonds.erase ( key ) > 0;
    }
}


bool TrussBuilder::removeNodesWithoutBond( ){
    int  * permut = new int[nodes.size()];
    for(int i=0; i<nodes.size(); i++ ){
        permut[i] = -1;
    }
    for( auto it : bonds ){
        Bond& bond = it.second;
        permut[bond.i] = 1;
        permut[bond.j] = 1;
    }
    int n=0;
    for(int i=0; i<nodes.size(); i++ ){
        if( permut[i]>=0 ){
            permut[i]=n;
            nodes[n] = nodes[i];
            n++;
        }
    }
    for( auto it : bonds ){
        Bond& bond = it.second;
        bond.i = permut[bond.i];
        bond.j = permut[bond.j];
    }
    delete permut;
}

// ============= IO functions

void TrussBuilder::init( int nNodesGuess, int nBondsGuess, int nBondsTypesGuess ){
    bonds .clear(); nodeIs.clear(); // nodes .clear(); bondTypes.clear();
    bonds .reserve ( nBondsGuess );
    nodeIs.reserve ( nNodesGuess );
    nodes .reserve ( nNodesGuess );
    bondTypes.reserve(nBondsTypesGuess);
    setScaling( {1.0d,1.0d,1.0d} );
    pos0.set  ( -ioff, -ioff, -ioff );
    for( int i=0; i<bondTypes.size(); i++ ){
        bondTypes[i].id = i;
    }
};

void TrussBuilder::toFile( char * fname ){
    std::ofstream fout;
    fout.open( fname );
    fout << nodes.size() << "\n";
    for( int i=0; i<nodes.size(); i++ ){
        GridNode& node = nodes[i];
        fout << (node.ix-ioff) << " " << (node.iy-ioff) << " " << (node.iz-ioff) << "\n";
    }
    fout << bonds.size() << "\n";
    for( auto it : bonds ){
        Bond& bond = it.second;
        //fout << bond.i << " " << bond.j << " " << bond.type.id << "\n";
        fout << bond.i << " " << bond.j << " " << 0 << "\n";
    }
    fout.close();
}

void TrussBuilder::fromFile( char * fname ){
    std::ifstream fin;
    fin.open( fname );
    if ( fin.is_open() ){
        int nnodes,nbonds;
        fin >> nnodes;     printf( "nnodes %i \n", nnodes );
        nodes.clear(); nodeIs.clear(); bonds.clear();
        nodes .reserve( nnodes );
        nodeIs.reserve( nnodes );
        for(int i=0; i<nnodes; i++){
            int ix,iy,iz;
            fin >> ix >> iy >> iz;
            //fin >> iz >> iy >> ix;
            printf( " %i %i %i \n", ix, iy, iz );
            ix+=ioff; iy+=ioff; iz+=ioff;
            insertNode( ix, iy, iz );
        }
        fin >> nbonds;     printf( "nbonds %i \n", nbonds );
        for(int i=0; i<nbonds; i++){
            int inod, jnod, itype;
            fin >> inod >> jnod >> itype;
            //fin >> itype >> jnod >> inod;
            printf( " %i %i %i \n", inod, jnod, itype );
            insertBond( nodes[inod].ix, nodes[inod].iy, nodes[inod].iz, nodes[jnod].ix, nodes[jnod].iy, nodes[jnod].iz, bondTypes[itype] );
        }
    }
    fin.close();
}

void TrussBuilder::toSoftBody( SoftBody& truss ){
    truss.allocate( nodes.size(), bonds.size(), 3, NULL, NULL, NULL, NULL );
    /*
    for( auto it : nodeIs ){
        int_fast32_t node = it.second;
        //truss.points[ node.id ] = node.pos; // do we need this at all ?
        ID64 index;
        index.id = it.first;
        index2pos( {index.x,index.y,index.z}, truss.points[ node.id ] );
    }
    */
    for( int i=0; i<nodes.size(); i++ ){
        //truss.points[ node.id ] = node.pos; // do we need this at all ?
        GridNode& node = nodes[i];
        index2pos( {node.ix,node.iy,node.iz}, truss.points[ i ] );
    }
    for( auto it : bonds ){
        Bond& bond = it.second;
        truss.bonds[ bond.id ] = bond;
    }
}

