
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

#endif

