
#ifndef broadPhaseCollision_h
#define broadPhaseCollision_h

// http://buildnewgames.com/broad-phase-collision-detection/
// http://www.metanetsoftware.com/technique/tutorialB.html
// http://www.cs.yorku.ca/~amana/research/grid.pdf
// http://www.flipcode.com/archives/Raytracing_Topics_Techniques-Part_4_Spatial_Subdivisions.shtml
// http://www.bulletphysics.org/mediawiki-1.5.8/index.php/Broadphase
// http://www.cs.jhu.edu/~cohen/Publications/icollide.pdf


// IDEA - Fast Update Of Polar Basis For Point Clusters
//  - instead of diagonalization of 3x3 matrix, we may just do few(1-5) steps of PowerIteration
//  - bounds are calculated for previous orientation => does not have to screen points twice


// Sweep-and-prune inside larger boxes ?

#include <unordered_map>   // http://www.cplusplus.com/reference/unordered_map/unordered_multimap/
//#include <unordered_set>
//#include <set>
#include <map>
#include <algorithm>

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
    inline void insert( int o, const Sphere3d&   s ){ int inds[nIndTmpMax]; insertInds( o, ruler.overlap_Sphere  (s.p,s.r    ,inds), inds ); };
    inline void insert( int o, const Box&        b ){ int inds[nIndTmpMax]; insertInds( o, ruler.overlap_BBox    (b.a,b.b    ,inds), inds ); };
    inline void insert( int o, const Line3d&     l ){ int inds[nIndTmpMax]; insertInds( o, ruler.overlap_Line    (l.a,l.b    ,inds), inds ); };
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

    inline int overlap( const Sphere3d&   s, int* inds){ int cells[nIndTmpMax]; return findInCells_paint( ruler.overlap_Sphere  (s.p,s.r    ,cells), cells, inds ); }
    inline int overlap( const Box&        b, int* inds){ int cells[nIndTmpMax]; return findInCells_paint( ruler.overlap_BBox    (b.a,b.b    ,cells), cells, inds ); }
    inline int overlap( const Line3d&     l, int* inds){ int cells[nIndTmpMax]; return findInCells_paint( ruler.overlap_Line    (l.a,l.b    ,cells), cells, inds ); }
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

struct SAPitem{
    // TODO: perhaps it can be more efficient use independnet arrays per each axis
    void* o;
    Vec2f span;
    inline bool operator < (const SAPitem& o) const{ return (span.x < o.span.x); }
};


static int collideCrossSAP_Objects( std::vector<SAPitem>& sweepi, std::vector<SAPitem>& sweepj, bool inv=false ){
    int ni = sweepi.size();
    int nj = sweepj.size();
    int ncol=0;
    for(int i=0; i<ni; i++ ){
        const SAPitem& itemi = sweepi[i];
        Object3d* oi = (Object3d*)itemi.o;
        const Vec2f& spani = itemi.span;
        int j0 = 0;
        while( sweepj[j0].span.x<spani.x ){ j0++; if(j0>=nj) return ncol;  }
        int j = j0;
        while( sweepj[j].span.x<spani.y ){
            Object3d* oj = (Object3d*)sweepj[j].o;
            if(inv){ oj->collide(oi); }else{ oi->collide(oj); };
            //println( "add "+i+" "+j+" "+collisions.size()+" ("+pi.x+","+pi.y+") "+" ("+set2[j].x+","+set2[j].y+") " );
            j++;
            if(j>=nj) break;
        }
    }
    return ncol;
}

class SAPbuff{ public:
    // TODO:
    //  - This is general version, we may define specialized versions (e.g. along cartesian axis) later
    //  - We may cast generatized version to this version later
    //  - We sould store permutation, so that we can update particles
    int    ndirs;
    Vec3d* dirs = 0;
    //std::vector<int>[ndirs] sweeps;
    std::vector<std::vector<SAPitem>> sweeps; // TODO: should be later optimized to 2D array?

    void init( int ndirs_, Vec3d* dirs_, int nItemGeuss ){
        ndirs=ndirs_;
        dirs = dirs_;
        sweeps.resize(ndirs);
        int i=0;
        for( int i=0; i<ndirs; i++ ){
            sweeps[i].reserve(nItemGeuss);
        }
    }

    inline void addSphere(const void* o, const Vec3d& p, double R){
        // TODO: more effcient can be first add all objects to one axis-sweep, than loop over each axis per object
        for(int i=0; i<ndirs; i++){
            double x = dirs[i].dot(p);
            sweeps[i].push_back( { (void*)o, (float)(x-R), (float)(x+R) } ); // we can better use insertion sort or something like that
        }
    }

    inline void addObject( Object3d* o ){
        for(int idir=0; idir<ndirs; idir++){
            sweeps[idir].push_back( { (void*)o, o->spanAlongDir(dirs[idir]) } ); // we can better use insertion sort or something like that
        }
    }

    void addObjects(int n, Object3d** objs){ for(int i=0; i<n; i++){ addObject( objs[i] ); } }

    void updateObjects(){
        for(int idir=0; idir<ndirs; idir++){
            Vec3d dir = dirs[idir];
            std::vector<SAPitem>& sweep = sweeps[idir];
            for(SAPitem& item : sweep){
                //const Object3d& o = *(Object3d*)item.o;
                //double x = dir.dot( o.pos );
                //item.span.set( x-o.R,x+o.R );
                item.span = ((Object3d*)item.o)->spanAlongDir(dir);
            }
        }
    }

    int collideSelfObjects(int idir){
        //printf("collideSelfObjects \n");
        std::vector<SAPitem>& sweep = sweeps[idir];
        int n = sweep.size();
        int ncol=0,ncomp=0,nhit=0;

        //for(SAPitem& item : sweep){
        //    printf( "-- span.x=%g \n", item.span.x  );
        //}
        for( int i=0; i<n; i++){
            //printf( "i=%i \n", i  );
            //printf( "i=%i span.x=%g \n", i, sweep[i].span.x  );
            const SAPitem& itemi = sweep[i];
            Object3d* oi = (Object3d*)itemi.o;
            const Vec2f& spani = itemi.span;
            double xmax = spani.y;
            //printf( "%i(%g,%g)\n", i, spani.x, spani.y  );
            for( int j=i+1; j<n; j++){
                const SAPitem& itemj = sweep[j];
                //printf( "(%i,%i) i.%g >? j.%g  \n", i,j, xmax, itemj.span.x ); ncomp++;
                if ( itemj.span.x > xmax ) break; // cannot colide, and the later j also not
                Object3d* oj = (Object3d*)itemj.o;
                nhit+=
                oi->collide(oj);
                ncol++;
            }
        }
        //printf( "collideSelfObjects : ncomps %i ncadidates : %i nhits %i | nBrute %i \n", ncomp, ncol, nhit, n*(n-1)/2 );
        return ncol;
    }

    int collideCrossObjects(int idir, std::vector<SAPitem>& sweepj ){
        std::vector<SAPitem>& sweepi = sweeps[idir];
        int ncol = 0;
        ncol+=collideCrossSAP_Objects( sweepi, sweepj, false );
        ncol+=collideCrossSAP_Objects( sweepj, sweepi, true  );
        return ncol;
    }

    void sort(){
        // TODO: std::sort (quicksort) may not be optimal, std::stable_sort (with insertion sort) may be faster for almost-sorted array
        //http://www.cplusplus.com/reference/algorithm/stable_sort/
        //https://stackoverflow.com/questions/23985891/what-is-the-difference-between-stdsort-and-stdstable-sort
        int ncomp=0;
        //for( auto& sweep : sweeps ){
        for( std::vector<SAPitem>& sweep: sweeps ){
        std::sort( sweep.begin(), sweep.end() );
        //    std::sort( sweep.begin(), sweep.end(), [](const SAPitem& a, const SAPitem& b){ return a.span.x<b.span.x; } );
        /*
        std::sort( sweep.begin(), sweep.end(), [&](const SAPitem& a, const SAPitem& b){
            ncomp++;
            bool out  = a.span.x<b.span.x;
            printf( "compare %i %i a(%g,%g) b(%g,%g) \n", ncomp, out, a.span.x, a.span.y, b.span.x, b.span.y );
            return out; }
        );
        */
        }
    }
};

template<class A,class B,int NAXIS, int NAMAX, int NBMAX>
class CrossSAP{
    static constexpr int naxis  = NAXIS;
    static constexpr int namax  = NAMAX;
    static constexpr int nbmax  = NBMAX;
    Vec3d dirs[NAXIS];

    int na=0,nb=0;
    A*    aos [NAXIS][NAMAX];
    float amin[NAXIS][NAMAX];
    float amax[NAXIS][NAMAX];

    B*    bos [NAXIS][NBMAX];
    float bmin[NAXIS][NBMAX];
    float bmax[NAXIS][NBMAX];

};


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
    std::map<double,AOO> sweep; // unordered (?)

    void insert( int o, const Box& bbox ){
        //spanAlongDir();
        sweep.insert( { bbox.a.x, (AOO){o,bbox}} ).first;
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

