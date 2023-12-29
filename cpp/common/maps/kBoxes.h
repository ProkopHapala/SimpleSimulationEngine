
#ifndef  kBoxes_h
#define  kBoxes_h

//#include "geom3D.h";

#include <vector>
#include "geom3D.h"

#include "sweep.h"
#include "arrayAlgs.h"

// Inspired by k-Means algorithm and AABB-Tree
//  Pseudocode:
//  a. construction
//      1. pick K pivots out of N box-bound bodies and make pivot bounding box B_k out of it;        NOTE: k~sqrt(N), at random or according better heuristic)
//      2. for each body o_i i<N find bounding box B_kmin with minimal cost( combinedBox(B_kmin,o_i) ) ) and insert it into the box
//  b. collisions
//      1. find pairs B_i, B_j which overlaps ( K^2 comparisions )
//      2. for each such pair B_i B_j find which overlaps
//  c. updates
//      - check size and overlap of all boxes, disolve M worst boxes and create M new out of them
//  d. balancing
//      - find N boxes with most overlap and disolve them, add points to remaining boxes
//      - find N boxes with most particles/larges and split them to two halves
//
//  - It is actually very flat tree
//  - Wa can make it recursively for splitting further the larger boxes (=> N-tree)

// !!!!!!!!!!!!!!!!!!!!!
// !!!!!   IDEA :   !!!!
// !!!!!!!!!!!!!!!!!!!!!
//  We can do SweepAndPrune for those boxes
//   - this means that optimal K would be slightly higher than sqrt(N)


// ======================================
// ========= Free Functions =============
// ======================================

template<typename T>
void applyPermut( int n, int* permut, T* dataIn, T* dataOut ){
    for(int i=0; i<n; i++ ){ dataOut[i] = dataIn[permut[i]]; }
}

template<typename T>
void applyPermutInv( int n, int* permut, T* dataIn, T* dataOut ){
    for(int i=0; i<n; i++ ){ dataOut[permut[i]] = dataIn[i]; }
}

template<typename T>
void applyPermutTmp( int n, int* permut, T* data, T* tmp ){
    for(int i=0; i<n; i++ ){ tmp[i] = data[i]; }
    applyPermut( n, permut, tmp, data );
}

template<typename T>
void applyPermut_relloc( int np, T*& ps, int* c2o ){
    T* ps_ = ps;
    ps     = new T[np];
    applyPermut( np,c2o, ps_, ps );
    delete [] ps_;
}

void pickKofN( int K, int n, int* picked ){
    //_realloc( bounds, K );
    //_realloc( ranges, K );
    //int occupied[K];
    int i=0;
     // https://stackoverflow.com/questions/196017/unique-non-repeating-random-numbers-in-o1
    // https://en.wikipedia.org/wiki/Linear-feedback_shift_register
    // https://en.wikipedia.org/wiki/Maximum_length_sequence
    // random numbers should not repeate, but after modulo N they can
    while(i<K){
        int ii = rand()%n;
        bool isNew = true;
        for(int j=0;j<i;j++){ // check if new - can we do it in less than O(K^2) steps ?
            if( picked[j] == ii ){ isNew = false; break; };
            // TODO: alternatively we can check some heuristics to minimize overlap with previous
        }
        if(isNew){
            picked[i] = ii;
            i++;
        }
    }
}

// Universal function for putting objects into cells
void cellPermut( int nc, int no, Vec2i* cni, const int* o2c, int* c2o ){
    // cni[nc](work) ... (Vec2i){a=ni,b=i0}; a=ni=number of objects in cell;  b=i0=starting index
    // o2c[no](in  ) ... object->cell
    // c2o[no](out ) ... array of object indexes assigned to each cell (i.e. continuous segments [i0:i0+ni] of the array belong to each cell)
    // count number of objets in each cell
    for(int ic=0; ic<nc; ic++ ){ cni[ic].a=0;        };
    for(int io=0; io<no; io++ ){ cni[ o2c[io] ].a++; };
    //DEBUG
    int ntot=0;
    //int nmax=0;
    // find starting index for each cell by cumulative-sum of points in each cell
    //DEBUG
    for(int ic=0; ic<nc; ic++ ){
        Vec2i& c = cni[ic];
        c.b=ntot;
        ntot+=c.a;
        //printf( "[%i] ntot %i ni %i \n", ic, ntot, c.a );
        //if( c.a>nmax ) nmax=c.a;
        c.a=0;
    }
    //exit(0);
    //DEBUG
    for(int io=0; io<no; io++){
        int ic   = o2c[io];
        Vec2i& c = cni[ic];
        int jo   = c.b + c.a;
        //printf( "c2o[%i|%i] | ni %i i0 %i \n", jo, no, c.a, c.b );
        c2o[jo]   = io;
        c.a++;
    }
    //DEBUG
    //exit(0);
}

int findClosest(const Vec3d& p, int np, const Vec3d* ps, double r2max ){
    double r2min=r2max;
    int    ibest=-1;
    for(int i=0; i<np; i++){
        Vec3d  d  = ps[i]; d.sub(p);
        double r2 = d.norm2();
        if(r2<r2min){ ibest=i; r2min=r2; };
    }
    return ibest;
}

// ========================================
// ============= Class KPivots ============
// ========================================


class KPivots{ public:
    int    np;     // number of points
    Vec3d* ps;     // points
    int    K;
    //int*   pivots;

    Vec3d* pivots;

    Vec2i* cellNIs;
    int*   o2c;
    int*   c2o;

    //KPivots(int K_, int np_, Vec3d* ps_):K(K_),np(np_),ps(ps_){cells = new Vec2i[K];}
    KPivots(int K_):K(K_){ cellNIs = new Vec2i[K]; np=0;o2c=0;c2o=0; pivots=0; }
    ~KPivots(){ delete [] cellNIs; if(pivots) delete [] pivots; }

    void findPivots( bool bStorePos ){
        for(int i=0; i<K; i++){
            int ip = (rand()%(np-i))+i;
            if(bStorePos){
                //printf( "piv[%i] ip %i\n", i, ip );
                pivots[i] = ps[ip];
            }else{
                _swap(ps[i]    ,ps[ip]    );
                //_swap(permut[i],permut[ip]);
                o2c[i]=i;
            }
        }
    }

    void repelPivots( double rmin, double rmax ){
        if(rmax<0){ rmax=1e+300; }
        for(int i=0; i<K; i++){
            Vec3d& pi=pivots[i];
            for(int j=0; j<i; j++){
                Vec3d& pj=pivots[j];
                Vec3d d; d.set_sub(pj,pi);
                double r = d.norm();
                if(r<rmin){
                    double rnew = 0.5*(rmin-r);
                    if(rnew>rmax){ rnew=rmax; }
                    d.mul(  rnew/r );
                    pj.add(d);
                    pi.sub(d);
                }
            }
        }
    }

    void insert_distance(){
        // WARRNING : we assume first K points in ps are pivots
        for(int i=K;i<np; i++){
            int k = findClosest(ps[i], K, ps, 1e+300 );
            //printf( "o2c[%i] %i \n", i, k );
            o2c[i]=k;
        }
    }

    void insert_distance_piv( bool bCOG ){
        if(bCOG){ for(int k=0;k<K; k++){ cellNIs[k].a=0; } };
        for(int i=0;i<np; i++){
            int k = findClosest(ps[i], K, pivots, 1e+300 );
            o2c[i]=k;
            if(bCOG){
                int& ni   = cellNIs[k].a;
                double f  = ni/(ni+1.);
                Vec3d& pk = pivots[k];
                pk.mul(f);
                pk.add_mul( ps[i], 1-f );
                ni++;
            }
        }
    }

    int* build( int np_, Vec3d* ps_, int* o2c_=0, int* c2o_=0, bool bStorePos=false, bool bCOG=true, double rmin=-1, double rmax=1e+300 ){
        ps=ps_; //bool new_o2c = (o2c==0);
        //DEBUG
        if(np_!=np){
            np=np_;
            if(o2c_ ){ o2c=o2c_; }else{ _realloc(o2c,np); };
            if(c2o_ ){ c2o=c2o_; }else{ _realloc(c2o,np); };
        }
        if(bStorePos && !pivots ){ pivots=new Vec3d[K]; };
        //DEBUG
        //o2c   =new int[np];
        //cells = new int[K ];
        findPivots     (bStorePos);  //DEBUG
        //for(int i=0; i<K; i++){ printf( "pivot[%i] (%g,%g,%g)\n", i, pivots[i].x,pivots[i].y,pivots[i].z ); }
        if(bStorePos){
            if(rmin>0)repelPivots( rmin, rmax );
            insert_distance_piv( bCOG );
        }else{
            insert_distance();           //DEBUG
        }
        cellPermut( K, np, cellNIs, o2c, c2o );  //DEBUG
        return c2o;
    }

};

// ========================================
// ============= Class HierarchicalKpivot ============
// ========================================

class HierarchicalKpivot{ public:
    int np;
    Vec3d *ps,*ps_;
    int   *o2c,*c2o;
    std::vector<Vec3d> pivots;
    std::vector<Vec2i> cells;
    //std::vector<Vec2i> cells;

    int perCellMax;
    int kmax;
    int levMax;

    HierarchicalKpivot( int np_, Vec3d* ps__ ){
        np=np_; ps=ps__;
        ps_ = new Vec3d[np];
        o2c = new int[np];
        c2o = new int[np];
    }

    void run( int perCellMax_, int kmax_, int levMax_ ){
        perCellMax=perCellMax_; kmax=kmax_; levMax=levMax_;
        recur( np, ps, ps_, o2c, c2o, 0, 0 );
    }

    void recur( int np, Vec3d* ps, Vec3d* ps_, int* o2c, int* c2o, int lev, int i0 ){
        //if(lev>levMax)return;
        int K = _min( 1+np/perCellMax, kmax );
        if(K<=1)return;
        printf( "%0*cK %i np %i \n", lev+1, '+', K, np );
        //if(np<100)for(int i=0; i<np; i++){ printf( "ps[%i] (%g,%g,%g)\n", i, ps[i].x,ps[i].y,ps[i].z ); }
        KPivots kpiv(K);
        kpiv.build( np, ps, o2c, c2o, true, false, 1.0, 2.0 );
        applyPermutTmp( np, c2o, ps, ps_ );
        //if(np<100)for(int i=0; i<np; i++){ printf( "ps_[%i] c2o[] %i (%g,%g,%g) \n", i, c2o[i], ps[i].x,ps[i].y,ps[i].z ); }
        //return;
        int nsum=0;
        for(int k=0; k<kpiv.K; k++){
            int ni = kpiv.cellNIs[k].a;
            nsum+=ni;
            printf( "%0*c[%i] ni %i \n", lev+2, '.', k, ni );
            if( (ni>perCellMax)&&(lev<levMax) ){
                int j0=kpiv.cellNIs[k].b;
                recur( ni, ps+j0, ps_+j0, o2c+j0, c2o+j0, lev+1, i0+j0 );
            }else{
                int npiv = pivots.size();
                pivots.push_back( kpiv.pivots[k] );
                kpiv.cellNIs[k].b+=i0;
                cells .push_back( kpiv.cellNIs[k] );
                int* o2c_ = o2c+i0;
                //for(int i=0; i<ni; i++){ o2c[i]+=npiv; }
            }
        }
        //printf( "%0*cnsum %i np %i \n", lev, '-', nsum, np );
    }

};

// ========================================
// ============= Class KBoxes  ============
// ========================================


/**
    * @brief Represents axis-aligned bouding box in 3D space and range of indexes pointing to bodies contained in the box which are stored in a permutation array.
    * 
    * This struct contains information about a box in a k-d tree, including the number of bodies in the block,
    * the index pointing to the beginning of the block in the permutation array, and the bounding box of the block.
    */
struct KBox{
    int n;    // number of bodies in the block
    int i0;   // pointing to beginning of block in permutation array
    Box span; // bounding box 
    //float cost; // this may be also outside
};

/**
 * @file kBoxes.h
 * @brief Defines the KBoxes class for spatial partitioning using K-d trees.
 * 
 * The KBoxes class is used for spatial partitioning using K-d trees. It is used for collision detection and for finding the nearest neighbour of a point and ray-tracing.
 */
class KBoxes{ public:
    std::vector<KBox>    branches; 
    // ToDo : if we split KBox to (Box)span and (Vec2i){n,i0} we can use general function   cellPermut()
    //std::vector<Box>    branches;  // we should probably split this
    //std::vector<Vec2i>  branche_ni;
    std::vector<Vec2i>  collisionPairs; // to output collision pairs for narrow-phase 
    //Box   * bounds;
    //Vec2i * ranges;
    int   * permut      = 0;    // indexes of bodies sorted by the branches  ( should call it branch2body )
    int   * body2branch = 0;    // inverse of permut

    bool   bakePermut = false;  // if true we will store permutated bodies in bodiesPremuted array
    float  growFactor = 1.6;    // how much we grow the arrays when we need more space
    int    nbodies =0;          // number of bodies
    int    nBodyMax=0;          // maximum number of bodies
    int    nBinKmax = 0;        // maximum number of bodies in a single branch
    Box   * bodies         = 0; // bounding boxes of individual bodies     (if we permute also this array we can optimize some cache access)
    Box   * bodiesPremuted = 0; // bounding boxes of individual bodies sorted by the branches (if bakePermut==true), this is optimized for cache access

    double costDisolveTrashold = 1e+6;  // if cost of the box (given by volume and number of bodies) is larger than this we disolve it   (We should think how to set this number properly)

    // ====== Sweep and Prune temporary variables
    bool bSweep = false;
    int*         Kpermut     = 0;  // permutation of branches
    sweep::Span* Kintervals  = 0;  // intervals of branches (KBoxes)
    sweep::Span* Bintervals  = 0;  // intervals of bodies   (Boxes)
    Int2*        Kcols       = 0;  // collision pairs of branches
    Int2*        colPairsTmp = 0;  // temporary storage for collision pairs
    Box   globBBox;                // global bounding box of all bodies
    float sweepXStep,invSweeXStep; // step of line-sweep and its inverse


    // ========== Functions

    /**
     * @brief Reserves memory for the bodies in the kBoxes class.
     * 
     * This function reserves memory for the bodies in the kBoxes class based on the given number of bodies (n) and a flag indicating whether permutation is required (bPermut).
     * If the given number of bodies is greater than the current maximum number of bodies (nBodyMax), the memory is reallocated to accommodate the new number of bodies.
     * 
     * @param n The number of bodies to reserve memory for.
     * @param bPermut A flag indicating whether permutation is required.
     */
    void reserveBodies( int n, bool bPermut ){
        //branches.reserve(K);
        //branches.resize(K);
        if(n>nBodyMax){
            nBodyMax = n*growFactor;
            _realloc( permut,      nBodyMax );
            _realloc( body2branch, nBodyMax );
            if(bPermut){
                _realloc( bodiesPremuted, nBodyMax);
            }
        }
    }

    void realocSweep( int nk, int n ){
        _realloc(Kpermut    ,nk );
        _realloc(Kintervals ,nk );
        _realloc(Bintervals ,n  );
        for(int i=0; i<nk; i++){ Kpermut[i]=i; };
    }

    void findGlobalBounds(){
        globBBox.a.set(+1e+100,+1e+100,+1e+100);
        globBBox.b.set(-1e+100,-1e+100,-1e+100);
        for(int i=0; i<nbodies; i++){
            globBBox.enclose( bodies[i] );
        }
    };

    /**
     * Sketches branches sweep.
     * 
     * This function divides the global bounding box into nk equal branches along the x-axis
     * and assigns each branch a span based on the sweep step size.
     * 
     * @param nk The number of branches to create.
     */
    void sketchBrachesSweep( int nk ){
        sweepXStep   = (globBBox.b.x-globBBox.a.x)/nk;
        invSweeXStep = 1.0/sweepXStep;
        for(int i=0; i<nk; i++){
            branches[i].n = 0;
            branches[i].span     = globBBox;
            branches[i].span.a.x = globBBox.a.x + sweepXStep*(i  );
            branches[i].span.b.x = globBBox.a.x + sweepXStep*(i+1);
        }
    };

    inline double insertCost( int i, const Box& box ){
        Box b;
        b.combine( box, branches[i].span );
        //return b.volume() + branches[i].n;  // TODO: think more about this function
        //return b.volume();  // TODO: think more about this function
        //return b.surfArea();  // TODO: think more about this function
        return box.center().dist2( branches[i].span.center() );
        //return sq( box.center().x - branches[i].span.center().x );    // Good in combination with line-sweep
        //return (b.dimensions() - branches[i].span.dimensions()).norm2();
    }

    inline double insertCostSweep( int i, const Box& box ){
        Box b;
        float a,cost = 0.0;
        a = branches[i].span.a.x - box.a.x;  if( a > 0 ) cost+=a;
        a = box.b.x - branches[i].span.b.x;  if( a > 0 ) cost+=a;
        //cost += branches[i].n*sweepXStep;
        //printf( "cost %i %g\n", i, cost );
        //   TODO : we should consider cost of N-elements in the set
        return cost;
    }

    /**
     * Inserts a box into the kBoxes data structure. It finds the branch with the lowest cost of inserting the box and inserts the box into that branch.
     * 
     * @param ib  The index of the box to be inserted.
     * @param box The box to be inserted.
     */
    void insert( int ib, const Box& box ){
        int    imin;
        double minCost = +1e+300;
        for(int i=0; i<branches.size(); i++){
            double cost = insertCost( i, box );
            if(cost<minCost){ minCost = cost; imin=i; }
        }
        branches   [imin].span.enclose( box );
        body2branch[ib] = imin;
    }

    /**
     * Inserts a box into the sweep structure at the specified index. It finds the branch with the lowest cost of inserting the box and inserts the box into that branch.
     * 
     * @param ib The index of the box to be inserted.
     * @param box The box to be inserted.
     */
    void insertSweep( int ib, const Box& box ){
        int    imin;
        double minCost = +1e+300;
        int ibeg = (int)( (box.a.x - globBBox.a.x)*invSweeXStep );
        int iend = (int)( (box.b.x - globBBox.a.x)*invSweeXStep )+1;
        //int nk = branches.size(); if(imax>nk) imax=nk;
        for(int i=ibeg; i<iend; i++){
            double cost = insertCostSweep( i, box );
            if(cost<minCost){ minCost = cost; imin=i; }
        }
        branches   [imin].span.enclose( box );
        body2branch[ib] = imin;
    }

    /**
     * Updates the permutation array based on new mapping of bodies to branches (body2branch array).
     * This function should be called after each rearrangement of bodies.
     */
    void updatePermut(){
        //printf( "updatePermut nbodies %i bSweep %i\n", nbodies, bSweep );
        // --- count number of bodies in each branch
        for(int i=0; i<nbodies; i++ ){ branches[ body2branch[i] ].n++; };
        // --- find starting index for each branch, also find maximum number of bodies in a single branch (nBinKmax)
        int ntot=0;
        for(int k=0; k<branches.size(); k++ ){
            branches[k].i0=ntot;
            int& ni =branches[k].n;
            ntot+=ni;
            if( ni>nBinKmax ) nBinKmax=ni;
            ni=0;
        }
        // --- fill the permutation array (permut) (inverse of body2branch)
        for(int i=0; i<nbodies; i++){
            int k =  body2branch[i];
            int j =  branches[k].i0 + branches[k].n;
            permut[j] = i;
            branches[k].n++;
        }
        // ToDo : use cellPermut here - problem is {ni,i0} is inside KBox
        //cellPermut( nc, no, cni, permut );
        if(bSweep){ realocTemp(); }
    }

    /**
     * Updates the sweep intervals for stable sorting of branches and bodies.
     * 
     * @param kSorted Flag indicating that the Kintervals array is almost sorted, therefore insertion sort is used rather than quicksort.
     * @param bSorted Flag indicating that the Bintervals array is almost sorted, therefore insertion sort is used rather than quicksort.
     */
    void updateSweepStable( bool kSorted=true, bool bSorted=true ){
        int nk = branches.size();
        // --- update Kintervals ( temporary array with x-intervals of branches used for sorting )
        for(int i=0; i<nk; i++){
            Box& span =  branches[i].span;
            Kintervals[i] = (sweep::Span){(float)span.a.x,(float)span.b.x};
        }
        // --- update Bintervals ( temporary array with x-intervals of bodies used for sorting )
        for(int i=0; i<nbodies; i++){ // TODO: we can merge this with findGlobalBounds() ?????
            Box& span =  bodies[i];
            Bintervals[i] = (sweep::Span){(float)span.a.y,(float)span.b.y};
        }
        // --- sort Kintervals  ( sort branches by x-coordinate )
        sort_permut( nk, Kpermut, Kintervals, false, kSorted ); //printf( "insertSort N: %i niters: %i \n",  nk, niter );
        // --- sort Bintervals  ( sort bodies inside each branch by x-coordinate )
        for(int i=0; i<nk; i++){
            KBox& B =  branches[i];
            sort_permut( B.n, permut+B.i0, Bintervals, false, bSorted );
        }
        //printf("updateSweep  1.3 \n");
    }

    /**
     * @brief Rebuilds the sweep structure.
     * This function updates the sweep structure by performing the following steps:
     * 1. Calculates the number of branches (nk) in the structure.
     * 2. Finds the global bounds of the structure.
     * 3. Sketches initial layout of intervals (branches) equaly spaced along the x-axis.
     * 4. Inserts the bodies into the sweep structure.
     */
    void rebuildSweep(){
        int nk = branches.size();
        findGlobalBounds();
        sketchBrachesSweep(nk);
        for(int i=0; i<nbodies; i++ ){ insertSweep( i, bodies[i] ); }
    }

    /**
     * @brief Updates the sweep by rebuilding it, updating the permutation, and updating the sweep by sorting branches and bodies.
     * This function updates the sweep structure by performing the following steps:
     * 1. Rebuilds the sweep structure.
     * 2. Updates the permutation.
     * 3. Updates the sweep by sorting branches and bodies.
     */
    void updateSweep(){
        rebuildSweep();   // find bounds and sketch branches equally spaced along x-axis
        updatePermut();   // update permutation array mapping bodies to branches
        //updateSweepStable(false,true);
        updateSweepStable(true,false);  // update sweep by sorting branches and bodies
        //updateSweepStable(true,true);
    }

    void build( int K, int nbodies_, bool bakePermut_, Box* bodies_, bool bSweep_=false ){
        bakePermut=bakePermut_;
        nbodies=nbodies_;
        bodies = bodies_;
        bSweep = bSweep_;
        reserveBodies( nbodies, bakePermut );
        //printf( "KBoxes::build \n" );
        //pickPivots( K, n );
        // insert pick pivots
        //   - TODO: maybe would be more efficient to use K-Means to pick pivots?
        branches.clear();
        branches.resize(K); // it will not free memory //https://stackoverflow.com/questions/1155693/stdvector-resize-downward
        //printf( "n %i K %i size %i \n", n, K, branches.size()  );
        int picked[K];
        //printf( "KBoxes::build 1.1  K %i  nb %i\n", K, nbodies );
        if(bSweep){
            //findGlobalBounds();
            //sketchBrachesSweep(K);
            //for(int i=0; i<nbodies; i++ ){ insertSweep( i, bodies[i] ); }
            rebuildSweep();
        }else{
            pickKofN( K, nbodies, picked );   // pick K bodies at random as pivots
            for(int k=0; k<K; k++ ){
                branches[k].span = bodies[ picked[k] ];
                branches[k].n = 0;
            }
            for(int i=0; i<nbodies; i++ ){ insert( i, bodies[i] ); }
        }
        //printf( "KBoxes::build 1.1.1\n" );
        //printf( "KBoxes::build 1.2\n" );
        // insert bodies
        updatePermut();
        if(bakePermut){
            applyPermut( nbodies, permut, bodies, bodiesPremuted );
            // TODO: we probably do not need to keep both arrays in memory
        }
        //printf( "KBoxes::build 1.3\n" );
        if(bSweep){
            realocSweep( K, nbodies );
            updateSweepStable( false );
            //updateSweep( true );
        }
        //printf(" KBoxes::build DONE! \n");
    }

    void updateStable(){
        // update without changing connectivity
        for(int i=0; i<branches.size(); i++){
            KBox& Bi = branches[i];
            for( int ii=Bi.i0; ii<Bi.i0+Bi.n; ii++ ){
                int ib = permut[ii];
                Bi.span.enclose( bodies[ib] );
            }
        }
        if(bSweep) updateSweepStable( true );
        //if(bakePermut){ applyPermut( n, permut, bodies, bodiesPremuted ); }
    }

    void updateStablePermuted(){
        applyPermut( nbodies, permut, bodies, bodiesPremuted );
        for(int i=0; i<branches.size(); i++){
            KBox& Bi = branches[i];
            for( int ii=Bi.i0; ii<Bi.i0+Bi.n; ii++ ){
                Bi.span.enclose( bodiesPremuted[ii] );
            }
        }
    }

    bool checkDisolve( const KBox& B ){
        double cost = B.span.volume() + B.n;   // TODO: think more about this function
        if( cost > costDisolveTrashold ){
            return true;
        }
        // we can do more complex tests - like checking overlap with others
        return false;
    }

    /**
     * Updates the kBoxes data structure with dissolving and generating new branches.
     * 1. Updates the bounding box of each branch based on the bodies it contains.
     * 2. Dissolves branches that meet certain criteria and stores the bodies and branch indexes.
     * 3. Generates new branches based on the dissolved branches and randomly picks pivots from the free bodies.
     * 4. Reinserts the free bodies into the new branches.
     * 5. Updates the permutation array.
     * 6. Updates the sweep structure if enabled.
     */
    void update(){ 
        int *freeBodies = new int[nbodies];   // to store indexes of bodies which are not in any branch
        int freeBranches[branches.size()];    // to store indexes of branches which are disolved
        int nfreeBranches = 0;                // number of disolved branches
        int nfreeBodies   = 0;                // number of bodies which are not in any branch
        for(int i=0; i<branches.size(); i++){
            KBox& Bi = branches[i];
            // --- update branch bounding box
            for( int ii=Bi.i0; ii<Bi.i0+Bi.n; ii++ ){
                int ib = permut[ii];
                Bi.span.enclose( bodies[ib] );
            }
            // --- try to disolve branch
            if( checkDisolve(Bi) ){
                for( int ii=Bi.i0; ii<Bi.i0+Bi.n; ii++ ){
                    freeBodies[ii] = permut[ii];
                }
                freeBranches[nfreeBranches]=i;
                nfreeBodies+=Bi.n;
            }
        }
        // --- Now we have to generate new branches ( What strategy? )
        int picked[nfreeBranches];
        pickKofN( nfreeBranches, nfreeBodies, picked );  // pick k pivots from free-bodies at random 
        for(int i=0; i<nfreeBranches;i++ ){              // setup new branches by picked pivots
            int k  = picked[i];
            int ib = freeBodies[k];
            branches[i].span = bodies[ib];
        }
        // --- reinsert free bodies to new branches
        for(int i=0; i<nfreeBodies; i++){
            int ib = freeBodies[i];
            insert( ib, bodies[ib] );
        }
        delete [] freeBodies;
        updatePermut();
        if(bSweep) updateSweepStable( true );
    }

    void collideBranch( const KBox& Bi ){
        int i0=Bi.i0;
        int n =Bi.n;
        for( int jj=i0; jj<i0+n; jj++ ){
            int jb  = permut[jj];
            const Box& bj = bodies[jb];
            for( int ii=jj+1; ii<i0+n; ii++ ){
                int ib  = permut[ii];
                if( bj.overlap( bodies[ib] ) ){
                    collisionPairs.push_back({ib,jb});
                };
            }
        }
    }

    void collideBranches( const KBox& Bi, const KBox& Bj ){
        // TODO : optimization where we fist check shared BBox
        for( int jj=Bj.i0; jj<Bj.i0+Bj.n; jj++ ){
            int jb  = permut[jj];
            const Box& bj = bodies[jb];
            if( Bi.span.overlap(bj) ){
                for( int ii=Bi.i0; ii<Bi.i0+Bi.n; ii++ ){
                    int ib  = permut[ii];
                    if( bj.overlap( bodies[ib] ) ){
                        collisionPairs.push_back({ib,jb});
                        // callback
                        // generate NarrowPhase pair ?
                    };
                }
            }
        }
    }

    inline int selectOverlap( const KBox& Bi, const KBox& Bj, int* sel ){
        int nsel = 0;
        for( int ii=Bi.i0; ii<Bi.i0+Bi.n; ii++ ){
            int ib = permut[ii];
            if( Bj.span.overlap( bodies[ib] ) ){ sel[nsel]=ib; nsel++; };
        }
        return nsel;
    }

    void collideBranchesNew( const KBox& Bi, const KBox& Bj ){
        // TODO : optimization where we fist check shared BBox
        //printf(" collideBranchesNew \n");
        int isel[Bi.n];  // selection of bodies in Bi which overlap with Bj
        int nsel = selectOverlap( Bi, Bj, isel );
        // --- select bodies in Bi which overlap with Bj
        //int nsel = 0;    // number of selected bodies
        //for( int ii=Bi.i0; ii<Bi.i0+Bi.n; ii++ ){
        //    int ib = permut[ii];
        //    if( Bj.span.overlap( bodies[ib] ) ){ isel[nsel]=ib; nsel++; };
        // }
        // --- check overlap of selected bodies with bodies in Bj
        for( int jj=Bj.i0; jj<Bj.i0+Bj.n; jj++ ){
            int jb        = permut[jj];
            const Box& bj = bodies[jb];
            if( !Bi.span.overlap(bj) ) continue;
            for( int ii=0; ii<nsel; ii++ ){
                int ib  = isel[ii];
                if( bj.overlap( bodies[ib] ) ){
                    collisionPairs.push_back({ib,jb});
                    // ToDo: callback
                    // ToDo: generate NarrowPhase pair ?
                };
                //printf( "-- ib %i jb %i \n", ib, jb );
            }
        }
    }

    void collideBranchesNew2( const KBox& Bi, const KBox& Bj ){
        // TODO : optimization where we fist check shared BBox
        //printf(" collideBranchesNew \n");
        int isel[Bi.n];  // selection of bodies in Bi which overlap with Bj
        int jsel[Bj.n];  // selection of bodies in Bj which overlap with Bi
        int nj = selectOverlap( Bj, Bi, jsel );
        int ni = selectOverlap( Bi, Bj, isel );
        // --- check overlap of selected bodies with bodies in Bj
        for( int jj=0; jj<nj; jj++ ){
            int jb        = jsel[jj];
            const Box& bj = bodies[jb];
            if( !Bi.span.overlap(bj) ) continue;
            for( int ii=0; ii<ni; ii++ ){
                int ib  = isel[ii];
                if( bj.overlap( bodies[ib] ) ){
                    collisionPairs.push_back({ib,jb});
                    // ToDo: callback
                    // ToDo: generate NarrowPhase pair ?
                };
                //printf( "-- ib %i jb %i \n", ib, jb );
            }
        }
    }

    void collideSelf(){
        collisionPairs.clear();
        int K = branches.size();
        for( int i=0; i<K; i++ ){
            const KBox& Bi = branches[i];
            collideBranch( Bi );
            for( int j=i+1; j<K; j++ ){
                if( Bi.span.overlap( branches[j].span ) ){
                    collideBranches( Bi, branches[j] );
                }
            }
        }
    }

    void collideSelfNew(){
        collisionPairs.clear();
        int K = branches.size();
        for( int i=0; i<K; i++ ){
            const KBox& Bi = branches[i];
            collideBranch( Bi );
            for( int j=i+1; j<K; j++ ){
                if( Bi.span.overlap( branches[j].span ) ){
                    collideBranchesNew(  Bi, branches[j] );
                }
            }
        }
    }

    void collideBranches_permuted( const KBox& Bi, const KBox& Bj ){
        for( int jj=Bj.i0; jj<Bj.i0+Bj.n; jj++ ){
            const Box& bj = bodiesPremuted[jj];
            for( int ii=Bi.i0; ii<Bi.i0+Bi.n; ii++ ){
                if( bj.overlap( bodiesPremuted[ii] ) ){
                    //collisionPairs[].push_back({ib,jb});
                    // callback
                    // generate NarrowPhase pair ?
                };
            }
        }
    }

    void collideSelf_permuted(){
        int K = branches.size();
        for( int i=0; i<K; i++ ){
            const KBox& Bi = branches[i];
            for( int j=0; j<K; j++ ){
                if( Bi.span.overlap( branches[j].span ) ){
                    collideBranches_permuted( branches[j], Bi );
                }
            }
        }
    }

    inline void narrowColPairs( int nc, Int2* colPairsTmp ){
        //printf( "narrowColPairs \n" );
        for(int j=0; j<nc; j++){
            const Int2& p = colPairsTmp[j];
            if( bodies[p.i].overlap( bodies[p.j] ) ){
                collisionPairs.push_back({p.i,p.j});
            }
        }
    }

    void realocTemp(){
        //Kcols       = Int2[nk*(nk-1)/2];
        //colPairsTmp = Int2[ nBinKmax*(nBinKmax-1)/2 ];
        int nk = branches.size();
        //printf( "realocTemp : nk %i(%i) nBinKmax %i (%i) \n", nk, nk*(nk-1)/2, nBinKmax, nBinKmax*nBinKmax );
        _realloc( Kcols,       nk*(nk-1)/2 );
        _realloc( colPairsTmp, nBinKmax*nBinKmax );
    }

    int collideSelfSweep(){
        int nk = branches.size();
        //bDEBUG = true;
        //printf("DEBUG 1 \n" );
        int nkcol  = sweep::collideSelf( nk, Kpermut, Kintervals, Kcols );
        //bDEBUG = false;
        //Int2* colPairs = (Int2*)collisionPairs.data();
        int ncol=0;
        //printf("DEBUG 2 \n" );
        for(int i=0; i<nk; i++){
            int nc = sweep::collideSelf(  branches[i].n, permut+branches[i].i0, Bintervals, colPairsTmp );
            narrowColPairs(nc,colPairsTmp);
        }
        //printf("DEBUG 3 \n" );
        for(int i=0; i<nkcol; i++){
            KBox& Bi = branches[Kcols[i].i];
            KBox& Bj = branches[Kcols[i].j];
            //if( anyOf(Kcols[i].i,0,3) && anyOf(Kcols[i].j,0,3) ){ printf( "DEBUG branch-pair: (%i,%i) \n", Kcols[i].i, Kcols[i].j  ); }
            //bool bOverlap = !( (Bi.span.a.x>Bj.span.b.x) || (Bj.span.a.x>Bi.span.b.x) ); // NOT OVERLAP
            //printf( "collideSelfSweep branch %i (%i,%i)   %i   \n", i, Kcols[i].i, Kcols[i].j, bOverlap );
            // TODO: this can be perhaps still optimized if we test againts overlap(Bi,Bj)   ?????
            //printf( " ppermuti %i ppermutj %i \n", permut+Bi.i0+Bi.n, permut+Bj.i0+Bj.n );
            //printf( "icol %i i,j(%i,%i) i0,j0(%i,%i) ni,nj(%i,%i) *permuti %i *permutj %i \n", i, Kcols[i].i, Kcols[i].j,   Bi.i0, Bj.i0, Bi.n, Bj.n, permut[Bi.i0], permut[Bj.i0] );
            int nc = sweep::collideCross(  Bi.n, permut+Bi.i0, Bintervals, Bj.n, permut+Bj.i0, Bintervals, colPairsTmp );
            narrowColPairs(nc,colPairsTmp);
        };
        //collisionPairs.resize(ncol);
        return ncol;
    };

    int findBody( int i, int& ibranch ){
        ibranch = body2branch[i];
        int ii0=branches[ibranch].i0;
        for(int ii=ii0; ii<ii0+branches[ibranch].n; ii++){
            if( permut[ii] == i ) return ii;
        }
        return -1;
    }

    void printBranchInfo(){
        int nzero = 0;
        for(int i=0; i<branches.size(); i++){
            const KBox& Bi = branches[i];
            if( Bi.n == 0 ) nzero++;
            printf( "branchs[%i] n %i i0 %i pbeg %i pend %i \n", i, Bi.n, Bi.i0, permut[Bi.i0], permut[Bi.i0+Bi.n] );
        }
        printf( "%i of %i branches (%g) has zero elements  \n", nzero, branches.size(), nzero/((float)branches.size()) );
    }

};

template<typename Lambda>
void collideKBoxes( KBoxes& kA, KBoxes& kB, Lambda callback ){
    int Ka = kA.branches.size();
    int Kb = kB.branches.size();
    for( int i=0; i<Ka; i++ ){
        KBox& Bi = kA.branches[i];
        for( int j=0; j<Kb; j++ ){
            if( Bi.span.overlap( kB.branches[j].span ) ){
                KBox& Bj = kB.branches[j];
                for( int jj=Bj.i0; jj<Bj.i0+Bj.n; jj++ ){
                    int jb = kB.permut[jj];
                    for( int ii=Bi.i0; ii<Bi.i0+Bi.n; ii++ ){
                        // We actually do not need body-boxes here
                        int ia = kB.permut[ii];
                        callback(ia,jb); // Narrow space collision
                    }
                }
            }
        }
    }
}

template<typename Lambda>
void collideKBoxes_permuted( KBoxes& kA, KBoxes& kB, Lambda callback ){
    int Ka = kA.branches.size();
    int Kb = kB.branches.size();
    for( int i=0; i<Ka; i++ ){
        KBox& Bi = kA.branches[i];
        for( int j=0; j<Kb; j++ ){
            if( Bi.span.overlap( kB.branches[j].span ) ){
                KBox& Bj = kB.branches[j];
                for( int jj=Bj.i0; jj<Bj.i0+Bj.n; jj++ ){
                    Box& bj = kB.bodiesPremuted[jj];
                    for( int ii=Bi.i0; ii<Bi.i0+Bi.n; ii++ ){
                        if( bj.overlap( kA.bodiesPremuted[ii] ) ){
                            callback(ii,jj);
                        }
                    }
                }
            }
        }
    }
}

#endif

