
#ifndef  HashMap3D_h
#define  HashMap3D_h

#include <vector>

#include "fastmath.h"
#include "Vec3.h"
#include "geom3D.h"

#include "grids3D.h"
#include "HashMap64.h"
#include "HashMapT64.h"  // this is for very big scenes

typedef unsigned short  UHALF;
//typedef short  UHALF;

const UHALF MAP_OFFSET = 0x7FFF;

class HashMap3D : public HashMap64{ public:
    // this works only for maps size   NX < 2^21 = 2,097,152

    CubeGridRulerUnbound  ruler;
    HashMap64             hmap;

    int        nIndTmpMax = 256;
    uint64_t * tmpCells = 0;
    int      * tmpInds  = 0;
    bool     * isOld    = 0;  // This array keeps track if particular object was already touched

    inline void insertInds(int o, int n, uint64_t* inds ){
        for(int i=0; i<n; i++){
            hmap.insert( o, inds[i] );
        }
    }

    inline void insert( int o, const Box& b ){ insertInds( o, ruler.overlap_BBox(b.a,b.b,tmpCells), tmpCells ); };

    inline int collide( int o, Box* boxes, Vec2i* colPairs, bool bSelf ){
        const Box& b = boxes[o];
        int ncell = ruler.overlap_BBox(b.a,b.b,tmpCells);
        int nfound=0;
        for(int i=0; i<ncell; i++){
            nfound += hmap.getAllInBoxOnce( tmpCells[i], tmpInds+nfound, isOld );
        }
        int ncol = 0;
        for(int ii=0; ii<nfound; ii++){
            int i    =           tmpInds[ii];
            int oi    = (uint32_t)hmap.store[i];
            isOld[oi] = false;
            if( (o==oi) & bSelf ) continue;
            if( b.overlap(boxes[oi]) ){
                colPairs[ncol] = {o,oi};
                ncol++;
            }
        }
        return ncol;
    }

};

#endif

