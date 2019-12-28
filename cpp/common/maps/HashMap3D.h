
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

    bool bDEBUG = true;
    //int DEBUG_ninserts = 0;

    int   nMaxCell   = 256;
    int   nMaxFound  = 4096;
    int   nbodies    = 0;

    uint64_t * tmpCells = 0;
    int      * tmpFound = 0;
    bool     * isOld    = 0;  // This array keeps track if particular object was already touched

    int realloc( int n ){
        _realloc(tmpCells,nMaxCell);   //printf( " DEBUG 1 \n"  );
        _realloc(tmpFound,nMaxFound);  //printf( " DEBUG 2 \n"  );
        return hmap.reserve( n );      //printf( " DEBUG 3 \n" );
    }

    void setBodyBuff(int nbodies_){
        nbodies = nbodies_;
        _realloc(isOld,nbodies);
        for(int i=0; i<nbodies; i++){ isOld[i] = false; }
    }

    inline void insertInds(int o, int n, uint64_t* inds ){
        //printf( "insertInds o %i n %i \n", o, n );
        for(int i=0; i<n; i++){
            //DEBUG_ninserts++;
            hmap.insert( o, inds[i] );
        }
    }

    inline void insert( int o, const Box& b ){ insertInds( o, ruler.overlap_BBox(b.a,b.b,tmpCells), tmpCells ); };

    inline int collide( int o, const Box& b, Box* boxes, Vec2i* colPairs ){
        int ncell = ruler.overlap_BBox(b.a,b.b,tmpCells);

        if(bDEBUG){ if(ncell>nMaxCell){ printf("ncell(%i) > nMaxCell(%i)\n", ncell, nMaxCell); exit(0); }   };

        //printf( "collide ncell %i \n", ncell );
        int nfound=0;
        for(int i=0; i<ncell; i++){
            nfound += hmap.getAllInBoxOnce( tmpCells[i], tmpFound+nfound, isOld );
            //printf( "collide cell[%i] : %li nfound %i \n", i, tmpCells[i], nfound );
        }
        if(bDEBUG){ if(nfound>nMaxFound){ printf("ncell(%i) > nMaxCell(%i)\n", nfound, nMaxFound); exit(0); }   };

        //printf( "collide nfound %i \n", nfound );
        int ncol = 0;
        for(int ii=0; ii<nfound; ii++){
            int i    =           tmpFound[ii];
            int oi    = (uint32_t)hmap.store[i];
            //printf( "%i store[%i] %li  oi %i  o %i\n", ii, i, oi, hmap.store[i], o );
            isOld[oi] = false;
            if( o==oi ) continue;
            if( b.overlap(boxes[oi]) ){
                //printf( "col[%i] (%g,%g,%g) (%g,%g,%g) j %i \n", ncol, boxes[oi].a.x, boxes[oi].a.y, boxes[oi].a.z,     boxes[oi].b.x, boxes[oi].b.y, boxes[oi].b.z, oi );
                colPairs[ncol] = {o,oi};
                ncol++;
            }
        }
        return ncol;
    }

    //inline int nbytes(){ return sizeof(*this) + ( capacity*(  sizeof(uint64_t)*2 +  sizeof(uint32_t) + sizeof(uint16_t)  ) ); }

};

#endif

