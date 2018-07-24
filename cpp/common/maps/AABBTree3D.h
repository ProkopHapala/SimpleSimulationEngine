
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
        if( span ){
        }
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

struct AABBCellBin{   // NOTE: this is data cell in array, not Node class
    Box3f     span;    //
    int  subs;    // pointer to branches or leafs data block, if top bit set it is pointer to body
    int  parrent; // index of parrent

    inline float insertCostBranch(const Box3f& box,Box3f& tmpBox1) const { };
    inline float insertCostLeaf  (const Box3f& box,Box3f& tmpBox1) const { };

};

class AABBTreeBin{
    std::vector<AABBCellBin> cells;
    Box3f tmpBox1,tmpBox2; // insertCostBranch and insertCostBranch calculate new bbox anyway, we need to reuse it when inserting
    
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

};




#endif

