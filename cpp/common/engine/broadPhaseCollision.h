
#ifndef broadPhaseCollision_h
#define broadPhaseCollision_h

// http://buildnewgames.com/broad-phase-collision-detection/
// http://www.metanetsoftware.com/technique/tutorialB.html
// http://www.cs.yorku.ca/~amana/research/grid.pdf
// http://www.flipcode.com/archives/Raytracing_Topics_Techniques-Part_4_Spatial_Subdivisions.shtml
// http://www.bulletphysics.org/mediawiki-1.5.8/index.php/Broadphase
// http://www.cs.jhu.edu/~cohen/Publications/icollide.pdf



// Sweep-and-prune inside larger boxes ?

#include <unordered_map>   // http://www.cplusplus.com/reference/unordered_map/unordered_multimap/
//#include <unordered_set>
//#include <set>
#include <map>

//#include "Projectile3D.h"
//#include "Gun3D.h"
//#include "Object3D.h"

#include "geom3D.h"
#include "grids3D.h"

#include "arrayAlgs.h"

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
    inline void insert( int o, const Line3D&       l ){ int inds[nIndTmpMax]; insertInds( o, ruler.overlap_Line    (l.a,l.b    ,inds), inds ); };
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

    inline int overlap( const Sphere&     s, int* inds){ int cells[nIndTmpMax]; return findInCells_paint( ruler.overlap_Sphere  (s.p,s.r    ,cells), cells, inds ); }
    inline int overlap( const Box&        b, int* inds){ int cells[nIndTmpMax]; return findInCells_paint( ruler.overlap_BBox    (b.a,b.b    ,cells), cells, inds ); }
    inline int overlap( const Line3D&     l, int* inds){ int cells[nIndTmpMax]; return findInCells_paint( ruler.overlap_Line    (l.a,l.b    ,cells), cells, inds ); }
    inline int overlap( const Triangle3D& t, int* inds){ int cells[nIndTmpMax]; return findInCells_paint( ruler.overlap_Triangle(t.a,t.b,t.c,cells), cells, inds ); }

};





// =========================
// ====== Sweep And Prune
// =========================
//  - http://www.cs.jhu.edu/~cohen/Publications/icollide.pdf
//  - https://ieeexplore.ieee.org/document/7935404/   Adaptive Collision Culling for Massive Simulations by a Parallel and Context-Aware Sweep and Prune Algorithm
//  - http://www.codercorner.com/SAP.pdf

//https://www.gamedev.net/forums/topic/473529-sweep-and-prune-sort-and-sweep/

// Neraly Sorted Array https://www.toptal.com/developers/sorting-algorithms/nearly-sorted-initial-order
// https://rosettacode.org/wiki/Sorting_algorithms/Insertion_sort#C
// https://rosettacode.org/wiki/Sorting_algorithms/Shell_sort#C

int sweepAndPrune_simple( int n, Box* boxes, Vec2i* collisions ){
    // sort by box.a.array[ax];
    int ncol = 0;
    for( int i=0; i<n; i++){
        Box& bi = boxes[i];
        double xbi = bi.b.x;
        for( int j=i+1; j<n; j++){
            Box& bj = boxes[j];
            if ( bj.a.x > xbi ) break; // cannot colide, and the later j also not
            if ( ( bj.a.x > bi.b.y ) || ( bj.b.y < bi.a.y ) ) continue;
            if ( ( bj.a.z > bi.b.z ) || ( bj.b.z < bi.a.z ) ) continue;
            collisions[ncol] = {i,j};
            ncol++;
        }
    }
    return ncol;
}

int sweepAndPrune_simple( int n, Box* boxes, Vec3i axes, Vec2i* collisions ){
    // sort by box.a.array[ax];
    int ax=axes.x;
    int ncol = 0;
    for( int i=0; i<n; i++){
        Box& bi = boxes[i];
        double xbi = bi.b.array[ax];
        for( int j=i+1; j<n; j++){
            Box& bj = boxes[j];
            if ( bj.a.array[ax] > xbi ) break; // cannot colide, and the later j also not
            if ( ( bj.a.array[axes.b] > bi.b.array[axes.b] ) || ( bj.b.array[axes.b] < bi.a.array[axes.b] ) ) continue;
            if ( ( bj.a.array[axes.c] > bi.b.array[axes.c] ) || ( bj.b.array[axes.c] < bi.a.array[axes.c] ) ) continue;
            collisions[ncol] = {i,j};
            ncol++;
        }
    }
    return ncol;
}






/*
int sweepAndPrune_2set( int n, Box* boxes, Box* boxes, Vec3i axes, Vec2i* collisions ){
    // sort by box.a.array[ax];
    int ax=axes.x;
    int ncol = 0;
    for( int i=0; i<n; i++){
        Box& bi = boxes[i];
        double xbi = bi.b.array[ax];
        for( int j=i+1; j<n; j++){
            Box& bj = boxes[j];
            if ( bj.a.array[ax] > xi ) break; // cannot colide, and the later j also not
            if ( ( bj.a.array[axes.b] > bi.b.array[axes.b] ) || ( bj.b.array[axes.b] < bi.a.array[axes.b] ) ) continue;
            if ( ( bj.a.array[axes.c] > bi.b.array[axes.c] ) || ( bj.b.array[axes.c] < bi.a.array[axes.c] ) ) continue;
            collisions[ncol] = {i,j};
            ncol++;
        }
    }
    return ncol;
}
*/





class AOO{ public:
    int o;          // object
    //double center;
    //double size;
    Box bbox;
};

// see https://stackoverflow.com/questions/2620862/using-custom-stdset-comparator
//auto comp_bbox_x = [](const AOO& a, const AOO& b){ return a.bbox.a.x > b.bbox.a.x; };
//std::map<int,, decltype(comp_bbox_x)> sweep(comp_bbox_x);

class SweepAndPrune{  public:
    Vec3d dir;
    double maxSize=10.0;
    //std::set<int, decltype(comp_bbox_x)>; // sweep(comp_bbox_x);
    std::map<double,AOO> sweep;

    int insert( int o, const Box& bbox ){
        //spanAlongDir();
        sweep.insert( { bbox.a.x, (AOO){o,bbox}} );
    }

    int overlap( const Box& bbox, int *inds ){
        auto xitlow = sweep.lower_bound( bbox.a.x -maxSize );
        auto xittup = sweep.upper_bound( bbox.b.x +maxSize );
        //auto xitlow = sweep.lower_bound( dir.dot( bbox.a ) -maxSize );
        //auto xittup = sweep.upper_bound( dir.dot( bbox.b ) +maxSize );
        int io=0;
        for( auto it = xitlow; it!=xittup; ++it ){
            AOO& aoo = it->second;
            if( bbox.overlap(aoo.bbox) ){
                inds[io] = aoo.o;
                io++;
            };
        }
        return io;
    }

};



/*
class AOO{
    int o;          // object
    double center;
    double size;
}

class PruneAndSweep{
    Vec3i axesOrder = axesOrder{0,1,2};
    double maxSize=10.0;
    //int n;
    //std::vector<AOO> *xos=0,*yos=0,*zos=0;
    std::vector<AOO> axis[3];

    // insertion search should be better https://en.wikipedia.org/wiki/Insertion_sort
    void sortAxis( std::vector<AOO>& axos ){ std::sort( axos.begin(), axos.end(), []bool(const AOO& a, const AOO& a) { return a.center > b.center; }); }

    sort(){
        sortAxis( axis[0] );
        sortAxis( axis[1] );
        sortAxis( axis[2] );
    }

    int overlap( const Box& bbox, int *inds ){
        double xmin = bbox.a.array[axesOrder.x];
        double ymin = bbox.b.array[axesOrder.x];
        for(AOO& aoo: axis[axesOrder.x] ){

        }
    }

}
*/


#endif

