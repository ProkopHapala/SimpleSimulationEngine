
#ifndef  kBoxes_h
#define  kBoxes_h

//#include "geom3D.h";

#include <vector>
#include "geom3D.h"

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
//  - It is actually very flat tree
//  - Wa can make it recursively for splitting further the larger boxes (=> N-tree)

// !!!!!!!!!!!!!!!!!!!!!
// !!!!!   IDEA :   !!!!
// !!!!!!!!!!!!!!!!!!!!!
//  We can do SweepAndPrune for those boxes
//   - this means that optimal K would be slightly higher than sqrt(N)

struct KBox{
    // we should probably split this
    int n;  // 
    int i0; // pointg to beginning of block in permutation array
    Box span;
    //float cost; // this may be also outside
};

template<typename T> 
void applyPermut( int n, int* permut, T* dataIn, T* dataOut ){
    for(int i=0; i<n; i++ ){ dataOut[i] = dataIn[permut[i]]; }
}

template<typename T> 
void applyPermutInv( int n, int* permut, T* dataIn, T* dataOut ){
    for(int i=0; i<n; i++ ){ dataOut[permut[i]] = dataIn[i]; }
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
        }
    }
}

class KBoxes{ public:
    std::vector<KBox> branches;  // we should probably split this 
    //Box   * bounds;
    //Vec2i * ranges;
    int   * permut      = 0;
    int   * body2branch = 0; // inv permut
    
    bool   bakePermut = false;
    float  growFactor = 1.6;
    int    nbodies =0;
    int    nBodyMax=0;
    Box   * bodies         = 0;    // if we permute also this array we can optimize some cache
    Box   * bodiesPremuted = 0;
    
    double costDisolveTrashold = 1e+6; // We should think how to set this number properly
    
    // ========== Functions
    
    void reserveBodies( int n, bool bPermut ){
        //branches.reserve(K);
        //branches.resize(K);
        if(n>nBodyMax){
            nBodyMax = n*growFactor;
            _realloc( permut, nBodyMax);
            if(bPermut){
                _realloc( body2branch,    nBodyMax);
                _realloc( bodiesPremuted, nBodyMax);
            }
        }
    }
    
    /*
    void pickPivots( int K, int n, int* picked ){
        branches.reserve(K);
        branches.erase();
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
                if( occupied[j] == ii ){ isNew = false; break; };
                // TODO: alternatively we can check some heuristics to minimize overlap with previous
            }
            if(isNew){
                occupied[i] = ii;
                branches[i].span = bodies[ii];
            }
        }
    }
    */
    
    /*
    void pickPivotsKMeans( int K ){
        Vec3d cpos[K];
        int   cns [K];
        for(int i=0; i<K; i++){
            cpos.set(randf(),randf(),randf());
        }
        for(int i=0; i<nbodies; i++ ){
            int rmin=1e+300;
            int imin;
            for(int k=0; k<K; k++ ){
                cpos[k].;
                
            };
        };
    }
    */
    
    /*
    void applyKMeans(int K){
        for(int i=0; i<nbodies; i++ ){
            int rmin=1e+300;
            int imin;
            for(int k=0; k<K; k++ ){
                branches[k].span.a.add();
                
            };
        };
        span.
    }
    */
    
    inline double insertCost( int i, const Box& box ){
        Box b; 
        b.combine( box, branches[i].span );
        return b.volume() + branches[i].n;  // TODO: think more about this function
    }
    
    void insert( int ib, const Box& box ){
        // find best bounding-box to insert
        int    imin;
        double minCost = 1e+300;
        for(int i=0; i<branches.size(); i++){
            double cost = insertCost( i, box );
            if(cost<minCost){ minCost = cost; imin=i; }
        }
        // insert
        branches   [imin].span.enclose( box );
        body2branch[ib] = imin;
        //branches[imin].n++;
    }
    
    void updatePermut(){
        // we need to call this after each rearrangement of bodies
        for(int i=0; i<nbodies; i++ ){ branches[ body2branch[i] ].n++; };
        int ntot=0;
        for(int k=0; k<branches.size(); k++ ){ 
            branches[k].i0=ntot; 
            ntot+=branches[k].n; 
            branches[k].n=0;  
        }
        for(int i=0; i<nbodies; i++){
            int k =  body2branch[i];
            int j =  branches[k].i0 + branches[k].n; 
            permut[j] = i;
            branches[k].n++;
        }
    }
    
    void build( int K, int n, bool bakePermut_, Box* bodies_ ){
        bakePermut=bakePermut_;
        bodies = bodies_;
        reserveBodies( n, bakePermut );
        //pickPivots( K, n );
        // insert pick pivots
        //   - TODO: maybe would be more efficient to use K-Means to pick pivots?
        branches.reserve(K);
        branches.clear();
        int picked[K];
        pickKofN( K, n, picked );
        for(int k=0; k<K; k++ ){ branches[k].span = bodies[ picked[k] ]; }
        // insert bodies
        for(int i=0; i<n; i++ ){ insert( i, bodies[i] ); }
        updatePermut();
        if(bakePermut){
            applyPermut( n, permut, bodies, bodiesPremuted );
            // TODO: we probably do not need to keep both arrays in memory
        }
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
        double cost = B.span.volume() + B.n;
        if( cost > costDisolveTrashold ){ 
            return true;
        }
        // we can do more complex tests - like checking overlap with others
        return false;
    }
    
    void update(){
        // update with changing connectivity
        //int freeBodies  [nBodies]; // is this stack-allocation too much? 
        int *freeBodies = new int[nbodies];
        int freeBranches[branches.size()];
        int nfreeBranches = 0;
        int nfreeBodies   = 0;
        for(int i=0; i<branches.size(); i++){ 
            KBox& Bi = branches[i];
            for( int ii=Bi.i0; ii<Bi.i0+Bi.n; ii++ ){
                int ib = permut[ii];
                Bi.span.enclose( bodies[ib] );
            }
            if( checkDisolve(Bi) ){
                for( int ii=Bi.i0; ii<Bi.i0+Bi.n; ii++ ){
                    freeBodies[ii] = permut[ii];
                }
                freeBranches[nfreeBranches]=i;
                nfreeBodies+=Bi.n;
            }
        }
        // Now we have to generate new branches ... What strategy?
        int picked[nfreeBranches];
        pickKofN( nfreeBranches, nfreeBodies, picked );
        for(int i=0; i<nfreeBranches;i++ ){ 
            int k  = picked[nfreeBranches];
            int ib = freeBodies[i]; 
            branches[k].span = bodies[ib]; 
        }
        // reinsert free free bodies
        for(int i=0; i<nfreeBodies; i++){
            int ib = freeBodies[i]; 
            insert( ib, bodies[ib] );
        }
        delete [] freeBodies;
        updatePermut();
    }
    
    int collideBranches( const KBox& Bi, const KBox& Bj ){
        for( int jj=Bj.i0; jj<Bj.i0+Bj.n; jj++ ){
            int jb  = permut[jj];  
            const Box& bj = bodies[jb];
            for( int ii=Bi.i0; ii<Bi.i0+Bi.n; ii++ ){
                int ib  = permut[ii];
                if( bj.overlap( bodies[ib] ) ){
                    // callback
                    // generate NarrowPhase pair ?
                };
            }
        }
    }
    
    int collideSelf(){
        int K = branches.size();
        for( int i=0; i<K; i++ ){
            const KBox& Bi = branches[i];
            for( int j=0; j<K; j++ ){
                if( Bi.span.overlap( branches[j].span ) ){
                    collideBranches( branches[j], Bi );
                }
            }
        }
    }
        
    int collideBranches_permuted( const KBox& Bi, const KBox& Bj ){
        for( int jj=Bj.i0; jj<Bj.i0+Bj.n; jj++ ){
            const Box& bj = bodiesPremuted[jj];
            for( int ii=Bi.i0; ii<Bi.i0+Bi.n; ii++ ){
                if( bj.overlap( bodiesPremuted[ii] ) ){
                    // callback
                    // generate NarrowPhase pair ?
                };
            }
        }
    }
    
    int collideSelf_permuted(){
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
    
};

template<typename Lambda>
int collideKBoxes( KBoxes& kA, KBoxes& kB, Lambda callback ){
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
int collideKBoxes_permuted( KBoxes& kA, KBoxes& kB, Lambda callback ){
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

