
#ifndef broadPhaseCollision_h
#define broadPhaseCollision_h

// http://buildnewgames.com/broad-phase-collision-detection/

#include <unordered_map>   // http://www.cplusplus.com/reference/unordered_map/unordered_multimap/
//#include <unordered_set> 

//#include "Projectile3D.h"
//#include "Gun3D.h"
//#include "Object3D.h"

#include "geom3D.h"
#include "grids3D.h"

class BroadSpaceMapHash{
    int nIndTmpMax = 256;
    CubeGridRuler ruler;
    std::unordered_multimap<int,int> map;
    bool* isNew = nullptr;



    inline void insertInds(int o, int n, int* inds ){
        for(int i=0; i<n; i++){
            map.insert( {i,o} );
        }
    }

    inline void insert( int o, const Vec3d&      p ){ map.insert( {ruler.icell(p),o} ); };
    inline void insert( int o, const Sphere&     s ){ int inds[nIndTmpMax]; insertInds( o, ruler.overlap_Sphere  (s.p,s.r    ,inds), inds ); };
    inline void insert( int o, const Box&        b ){ int inds[nIndTmpMax]; insertInds( o, ruler.overlap_BBox    (b.a,b.b    ,inds), inds ); };
    inline void insert( int o, const Line&       l ){ int inds[nIndTmpMax]; insertInds( o, ruler.overlap_Line    (l.a,l.b    ,inds), inds ); };
    inline void insert( int o, const Triangle3D& t ){ int inds[nIndTmpMax]; insertInds( o, ruler.overlap_Triangle(t.a,t.b,t.c,inds), inds ); };

    inline int cell2array( int icell, int* inds ){ // TODO: if we use own hashmap map, we do not have to do this
        int io = 0;
        auto range = map.equal_range( icell );
        for(auto it = range.first; it != range.second; it++){
            inds[io] = it->second;
            io++;
        }
        return io;
    }

    inline int find( const Vec3d& p, int* inds){
        return cell2array( ruler.icell(p), inds );
    }

    /*
    inline int findInCells_naive( int n, int* cells, int* objects ){
        // problem - what if we have same object in two different buckets? we will insert it multiple times
        // We can perhaps "paint" objects which are already found
        int io=0;
        for(int i=0; i<n; i++){                                    // go over cells
            auto range = map.equal_range( cells[i] );             // list all objects in std::unordered_multimap<int,int> map;
            for(auto it = range.first; it != range.second; it++){  // 
                outInds[io]=it->second;
                io++;
            }
        }
        return io;
    }

    inline int findInCells_set( int n, int* cells, int* objects ){
        std::unordered_set<int> found;
        for(int i=0; i<n; i++){                                   // go over cells
            auto range = map.equal_range( cells[i] );             // list all objects in std::unordered_multimap<int,int> map;
            for(auto it = range.first; it != range.second; it++){ // 
                found.insert(it->second);
            }
        }
        int io=0;
        for(int o : found){ objects[io]=o; }
        return io;
    }
    */

    inline int findInCells_paint( int n, int* cells, int* objects ){
        int io=0;
        for(int i=0; i<n; i++){                                    // go over cells
            auto range = map.equal_range( cells[i] );             // list all objects in std::unordered_multimap<int,int> map;
            for(auto it = range.first; it != range.second; it++){  // 
                int o = it->second;
                if( isNew[o] ){
                    objects[io]=o;
                    isNew[o]=false;
                    io++;
                }
            }
        }
        for(int i=0; i<io; i++){ isNew[objects[io]]=true; };
        return io;
    }

    inline int find( const Sphere&     s, int* inds){ int cells[nIndTmpMax]; return findInCells_paint( ruler.overlap_Sphere  (s.p,s.r    ,cells), cells, inds ); }
    inline int find( const Box&        b, int* inds){ int cells[nIndTmpMax]; return findInCells_paint( ruler.overlap_BBox    (b.a,b.b    ,cells), cells, inds ); }
    inline int find( const Line&       l, int* inds){ int cells[nIndTmpMax]; return findInCells_paint( ruler.overlap_Line    (l.a,l.b    ,cells), cells, inds ); }
    inline int find( const Triangle3D& t, int* inds){ int cells[nIndTmpMax]; return findInCells_paint( ruler.overlap_Triangle(t.a,t.b,t.c,cells), cells, inds ); }

};


#endif

