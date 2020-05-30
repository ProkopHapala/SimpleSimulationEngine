
#ifndef  AABBTree3D_h
#define  AABBTree3D_h

#include "geom3D.h"

// References
//  - http://www.randygaul.net/2013/08/06/dynamic-aabb-tree/
//  - http://allenchou.net/2014/02/game-physics-broadphase-dynamic-aabb-tree/
//  - https://github.com/lohedges/aabbcc
//  - https://www.azurefromthetrenches.com/introductory-guide-to-aabb-tree-collision-detection/


template<typename O, int NSUBMAX, int NOBJMAX>
class AABBNode3D{
    static constexpr int nsubmax = NSUBMAX;
    static constexpr int nobjmax = NOBJMAX;
    int nobj=0,nsub=0;
    Box span;
    O           objs[NOBJMAX];
    AABBNode3D* subs[NSUBMAX];

    bool insert( O o, Box box){
        if( span.pointIn(o) ){
            // ToDo: Branching !!!
            if(nobj<nobjmax-1){
                objs[nobj]=o;
                nobj++;
            }
        }
        return false;
    }

};

struct Box3f{
    Vec3f a,b;
};


/*
IDEA:
- Node can occupy up to 6 leafs
- When number of leafs reach 6 we split it to 2 sub-nodes, each containing 3 sub nodes
    - this is not strict criterium, we may also prefer to minimize total volume of the sub nodes
- Since we now have free nodes

*/

struct AABBNode6{
    Box3f    span;     // 6x float32
    uint32_t subs[6];  // 6x int32 ... shared for branches and leafs, if 6th occupied it is leaf

    inline bool isLeaf(){ return subs[5]; };
};

class AABBTree6{
    std::vector<AABBNode6> nodes;

    inline void sort2(AABBNode6& nd ){
        int i=0;
        uint32_t s;
        do{
            s = nd.subs[i];
            nodes[s];
            i++;
        }while( s );
    }

};



// See Reference Here:
// - https://github.com/erincatto/Box2D/blob/master/Box2D/Collision/b2DynamicTree.cpp
// - https://github.com/bulletphysics/bullet3/blob/master/src/Bullet3Collision/BroadPhaseCollision/b3DynamicBvh.cpp
// Balancing by Tree Rotation AVL
// - https://www.geeksforgeeks.org/avl-tree-set-1-insertion/


#define NULL_NODE -1

struct AABBCellBin{   // NOTE: this is data cell in array, not Node class
    Box   span;    // TODO Box3f to save some cache
    int   subs;    // pointer to branches or leafs data block, if top bit set it is pointer to body
    int   parrent; // index of parrent
    inline float insertCostBranch(const Box3f& box,Box3f& tmpBox1) const { return 0.; };
    inline float insertCostLeaf  (const Box3f& box,Box3f& tmpBox1) const { return 0.; };
};


// Data Structure according Box2
// home/prokop/Dropbox/KnowDev/CODEs/Box2D_Bullet_AABBTree
// https://github.com/erincatto/Box2D/blob/master/Box2D/Collision/b2DynamicTree.cpp
struct AABBNodeBin{
	Box   aabb;   // TODO Box3f to save some cache
	void* userData;
	union{
		int parent;
		int next;
	};
	int child1;
	int child2;
	int height;  // leaf = 0, free node = -1

    inline bool isLeaf() const{ return child1 == NULL_NODE; }
};

class AABBTreeBin{
    //std::vector<AABBCellBin> cells;
    //Box3f tmpBox1,tmpBox2; // insertCostBranch and insertCostBranch calculate new bbox anyway, we need to reuse it when inserting

    std::vector<AABBNodeBin> nodes;
    int root;

    inline float costFunc( const Box& box ){
        return box.surfArea();
        //return box.volume();
    };


    /*
    int insertSubs(int iold, int i1, int o, const Box3f& box, const Box3f& tmpBox){
        AABBCellBin& a = cells[i1];
        if(i1<0){   // inserting to leaf
            cells.push_back( cells[iold] );   // https://arne-mertz.de/2016/02/modern-c-features-in-place-construction/
            cells[iold] = a;
            return i1;
        }
        // inserting to branch
        if(1){ // push down cost
            /// .....
        }
        int i2=i1+1; // we know that second sub-node is allocated right next to the first one - this saves some cache misses
        AABBCellBin& b = cells[i2];
        float cost1,cost2;
        //Box newBox1,newBox2;  insertCostBranch and insertCostBranch calculate new bbox anyway, we need to reuse it when inserting
        if( a.subs>0 ){ cost1 = a.insertCostBranch(box,tmpBox1); }else{ cost1 = a.insertCostLeaf(box,tmpBox1); };
        if( b.subs>0 ){ cost2 = b.insertCostBranch(box,tmpBox2); }else{ cost2 = b.insertCostLeaf(box,tmpBox2); };
        if( cost1 < cost2 ){ insertSubs( i1, a.subs, o, box, tmpBox1); }else{ insertSubs( i2, b.subs, o, box, tmpBox2); }
    }
    */

    // ==================    InsertLeaf !!!!!!!!!!!!!!!!!!!

inline void deleteNode( int i){

}

inline float evalInsertCost( int index, const Box& leafAABB ){
    Box aabb;
    aabb.combine( leafAABB, nodes[index].aabb );
    if ( nodes[index].isLeaf() ){
        return costFunc( aabb );
    }else{
        float oldArea = costFunc( nodes[index].aabb );
        float newArea = costFunc( aabb );
        return newArea - oldArea;
    }
}

int selectInsertLeaf( int index, const Box& leafAABB ){
	while ( nodes[index].isLeaf() == false) {

        index = selectInsertLeaf( index, leafAABB );

        int child1 = nodes[index].child1;
        int child2 = nodes[index].child2;

        // make combined AABB and measure area
        Box combinedAABB;
        combinedAABB.combine( nodes[index].aabb, leafAABB );
        float cost            =        2.0f * costFunc( combinedAABB      );
        float inheritanceCost = cost + 2.0f * costFunc( nodes[index].aabb );  // Minimum cost of pushing the leaf further down the tree

        float cost1 = evalInsertCost( child1, leafAABB ) + inheritanceCost;
        float cost2 = evalInsertCost( child2, leafAABB ) + inheritanceCost;
        // Descend according to the minimum cost.
        if (cost < cost1 && cost < cost2){ break; }  // no Descend
        if (cost1 < cost2){ index = child1; }        // Descend child1
        else              { index = child2; }        // Descend child2

	}
    return index;
}

void fixBoundUpwards(int index){
	// Walk back up the tree fixing heights and AABBs
	while (index != NULL_NODE ){
		index = balance(index);  // Balance  - this could be very costly
		int child1 = nodes[index].child1;
		int child2 = nodes[index].child2;
		nodes[index].height = 1 + _max( nodes[child1].height, nodes[child2].height );
		nodes[index].aabb.combine( nodes[child1].aabb,  nodes[child2].aabb );
		index = nodes[index].parent;
	}
}

// ==================    InsertLeaf !!!!!!!!!!!!!!!!!!!

void insertLeaf(int leaf){
	//++m_insertionCount;
	if (root == NULL_NODE ) {
		root = leaf;
		nodes[root].parent = NULL_NODE;
		return;
	}
	// Find the best sibling for this node

    Box leafAABB = nodes[leaf].aabb;
	int sibling = selectInsertLeaf( root, leafAABB );

	// Create a new parent.
	int oldParent = nodes[sibling].parent;

	//int newParent = AllocateNode();
	nodes.push_back( AABBNodeBin() );
	int newParent = nodes.size()-1;

	nodes[newParent].parent   = oldParent;
	nodes[newParent].userData = nullptr;
	nodes[newParent].height   = nodes[sibling].height + 1;
	nodes[newParent].aabb.combine(leafAABB, nodes[sibling].aabb);

	if (oldParent != NULL_NODE ){ // The sibling was not the root.
		if ( nodes[oldParent].child1 == sibling){ nodes[oldParent].child1 = newParent; }
		else                                    { nodes[oldParent].child2 = newParent; }
		nodes[newParent].child1 = sibling;
		nodes[newParent].child2 = leaf;
		nodes[sibling  ].parent = newParent;
		nodes[leaf     ].parent = newParent;
	}else{ // The sibling was the root.
		nodes[newParent].child1 = sibling;
		nodes[newParent].child2 = leaf;
		nodes[sibling  ].parent = newParent;
		nodes[leaf     ].parent = newParent;
		root = newParent;
	}
    fixBoundUpwards( nodes[leaf].parent );
	//Validate();
}

// ==================    RemoveLeaf !!!!!!!!!!!!!!!!!!!

void removeLeaf(int leaf){
	if (leaf == root){
		root = NULL_NODE;
		return;
	}
	int parent      = nodes[leaf].parent;
	int grandParent = nodes[parent].parent;
	int sibling;
	if (nodes[parent].child1 == leaf){ sibling = nodes[parent].child2; }
	else                               { sibling = nodes[parent].child1; }
	if (grandParent != NULL_NODE ){
		// Destroy parent and connect sibling to grandParent.
		if (nodes[grandParent].child1 == parent){ nodes[grandParent].child1 = sibling; }
		else                                    { nodes[grandParent].child2 = sibling; }
		nodes[sibling].parent = grandParent;
		deleteNode(parent);
		fixBoundUpwards( grandParent );
	}else{
		root = sibling;
		nodes[sibling].parent = NULL_NODE;
		deleteNode(parent);
	}
	//Validate();
}

// ==================    Balance  !!!!!!!!!!!!!!!!!!!

// Perform a left or right rotation if node A is imbalanced.
// Returns the new root index.
// https://www.geeksforgeeks.org/avl-tree-set-1-insertion/

inline void rotate( int iA, int iB, int iC, int iF, int iG ){
    nodes[iG].parent = iA;

    nodes[iA].child2 = iG;
    nodes[iA].aabb.combine( nodes[iB].aabb, nodes[iG].aabb );
    nodes[iA].height = 1 + _max( nodes[iB].height, nodes[iG].height);

    nodes[iC].child2 = iF;
    nodes[iC].aabb.combine( nodes[iA].aabb, nodes[iF].aabb );
    nodes[iC].height = 1 + _max( nodes[iA].height, nodes[iF].height);
}

inline int rotateUp( int iA, int iC, int iB ){
    int iF = nodes[iC].child1;
    int iG = nodes[iC].child2;
    // Swap A and C
    nodes[iC].child1 = iA;
    nodes[iC].parent = nodes[iA].parent;
    nodes[iA].parent = iC;
    // A's old parent should point to C
    int iCpar = nodes[iC].parent;
    if ( iCpar !=NULL_NODE ){
        if ( nodes[iCpar].child1 == iA){  nodes[iCpar].child1 = iC; }
        else                           {  nodes[iCpar].child2 = iC; }
    }else{
        root = iC;
    }
    // Rotate
    if ( nodes[iF].height > nodes[iG].height ){ rotate(iA, iB, iC, iF, iG ); }
    else                                      { rotate(iA, iB, iC, iG, iF ); }
    return iC;
}

int balance(int iA){
	if ( nodes[iA].isLeaf() || nodes[iA].height < 2){ return iA; }
	int iB = nodes[iA].child1;
	int iC = nodes[iA].child2;
	int balance = nodes[iC].height - nodes[iB].height;
	// Rotate C up
	if         (balance > 1){ // C.height >  B.height   =>   Rotate C up
        rotateUp( iA, iC, iB );
		return iC;
	} else 	if (balance < -1){  // C.height <  B.height   =>  Rotate B up
        rotateUp( iA, iB, iC );
		return iB;
	}
	return iA;
}

};




#endif

