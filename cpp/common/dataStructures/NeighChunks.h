
#ifndef NeighChunks_h
#define NeighChunks_h

/// @file NeighChunks.h
/// @brief Cache-line-aligned neighbor lists with overflow chaining — designed for dynamic mesh editing.
///
/// The problem: in a dynamic mesh, vertex degree changes as edges are added/removed.
/// A fixed-size array per vertex wastes space for low-degree vertices and overflows for high-degree ones.
/// A heap-allocated vector per vertex causes cache misses and fragmentation.
///
/// The solution: each vertex gets a primary chunk of exactly 1 cache line (64 bytes = 16 int32).
/// Chunk = [count, n0..n13, link]. 14 neighbors fit inline. If a vertex has more (rare but possible
/// at high-valence vertices in refined meshes), an extension chunk is allocated at the END of the
/// flat data array and linked. This keeps the common case (≤14 neighbors) at zero indirection,
/// one cache line, no heap allocation — while still handling the degenerate case gracefully.
///
/// The count+link design (vs. sentinel-only) was chosen because we need O(1) "how many neighbors"
/// for iteration bounds, and O(1) "is this chunk full" for overflow decisions. The sentinel (-1)
/// is still used within slots to mark unused entries after removal (compact() cleans them).

#include "macroUtils.h"

// Variable-length neighbor lists stored in fixed-size chunks with overflow linking.
// Each owner (e.g. vertex) gets a primary chunk of `slotSize` int32s.
// Chunk layout: [count, n0, n1, ..., n_{nMaxNeigh-1}, link]
//   count = number of valid neighbors in this chunk (not including extensions)
//   n_i  = neighbor index, or -1 if unused
//   link = index of extension chunk in data[], or -1 if none
// slotSize = nMaxNeigh + 2  (count + neighbors + link)
// Default nMaxNeigh=14 => slotSize=16 => 64 bytes = 1 cache line

class NeighChunks{ public:
    int nMaxNeigh = 14;
    int slotSize  = 16;  // nMaxNeigh + 2
    int nOwner    = 0;   // number of primary owners (e.g. vertices)
    int nTot      = 0;   // total number of chunks (primary + extension) currently allocated
    int nCap      = 0;   // capacity in number of chunks (data has nCap*slotSize ints)
    int* data     = 0;   // flat array of chunks, each slotSize long

    NeighChunks() = default;
    ~NeighChunks(){ _dealloc(data); }

    void setMaxNeigh(int n){
        nMaxNeigh = n;
        slotSize  = n + 2;
    }

    void realloc(int nOwner_){
        nOwner = nOwner_;
        nTot   = nOwner;  // initially just primary chunks, no extensions
        nCap   = nOwner + nOwner/4 + 16;  // 25% headroom for extensions
        _realloc(data, nCap * slotSize);
        for(int i=0; i<nCap*slotSize; i++) data[i] = -1;  // init all to -1
        for(int i=0; i<nOwner; i++) data[i*slotSize] = 0; // count=0 for primary chunks
    }

    void ensureCap(int nNeeded){
        if(nNeeded <= nCap) return;
        int nCapNew = nNeeded + nNeeded/4 + 16;
        int* dataNew = new int[nCapNew * slotSize];
        for(int i=0; i<nTot*slotSize; i++) dataNew[i] = data[i];
        for(int i=nTot*slotSize; i<nCapNew*slotSize; i++) dataNew[i] = -1;
        delete[] data;
        data = dataNew;
        nCap = nCapNew;
    }

    // get pointer to chunk at slot index
    inline int* chunk(int iSlot){ return data + iSlot*slotSize; }
    inline const int* chunk(int iSlot) const { return data + iSlot*slotSize; }

    // count = chunk[0]
    inline int count(int iSlot) const { return data[iSlot*slotSize]; }

    // total neighbor count across all extension chunks
    int countTotal(int iSlot) const {
        int n = 0;
        while(iSlot >= 0){
            n += data[iSlot*slotSize];
            iSlot = data[iSlot*slotSize + slotSize - 1]; // link
        }
        return n;
    }

    // add neighbor `neigh` to owner's list starting at primary slot `iSlot`
    // returns true on success
    bool add(int iSlot, int neigh){
        while(true){
            int* ch = chunk(iSlot);
            int& cnt = ch[0];
            if(cnt < nMaxNeigh){
                ch[1 + cnt] = neigh;
                cnt++;
                return true;
            }
            // this chunk is full, follow or create extension
            int& link = ch[slotSize - 1];
            if(link < 0){
                // allocate extension chunk
                ensureCap(nTot + 1);
                int iExt = nTot++;
                int* ext = chunk(iExt);
                ext[0] = 0;  // count=0
                ext[slotSize-1] = -1;  // no further extension
                link = iExt;
            }
            iSlot = link;
        }
    }

    // remove neighbor `neigh` from owner's list starting at primary slot `iSlot`
    // returns true if found and removed
    bool remove(int iSlot, int neigh){
        while(iSlot >= 0){
            int* ch = chunk(iSlot);
            int& cnt = ch[0];
            for(int i=0; i<cnt; i++){
                if(ch[1+i] == neigh){
                    // shift last element into this slot
                    cnt--;
                    ch[1+i] = ch[1+cnt];
                    ch[1+cnt] = -1;
                    return true;
                }
            }
            iSlot = ch[slotSize - 1]; // link
        }
        return false;
    }

    // find neighbor `neigh` in owner's list, returns true if found
    bool find(int iSlot, int neigh) const {
        while(iSlot >= 0){
            const int* ch = chunk(iSlot);
            int cnt = ch[0];
            for(int i=0; i<cnt; i++){
                if(ch[1+i] == neigh) return true;
            }
            iSlot = ch[slotSize - 1];
        }
        return false;
    }

    // iterate: copy all neighbors of owner into out[], return count
    int getAll(int iSlot, int* out) const {
        int n = 0;
        while(iSlot >= 0){
            const int* ch = chunk(iSlot);
            int cnt = ch[0];
            for(int i=0; i<cnt; i++) out[n++] = ch[1+i];
            iSlot = ch[slotSize - 1];
        }
        return n;
    }

    // remove dead neighbors (value == -1 or matching predicate) from all chunks of owner
    // deadVal: entries with this value are removed. Use -1 for default sentinel.
    void compact(int iSlot, int deadVal=-1){
        while(iSlot >= 0){
            int* ch = chunk(iSlot);
            int& cnt = ch[0];
            int w = 0;
            for(int r=0; r<cnt; r++){
                int v = ch[1+r];
                if(v != deadVal && v >= 0){
                    ch[1+w] = v;
                    w++;
                }
            }
            // clear remaining slots
            for(int i=w; i<cnt; i++) ch[1+i] = -1;
            cnt = w;
            iSlot = ch[slotSize - 1];
        }
    }

    // print one owner's neighbor list
    void print(int iSlot) const {
        printf("[");
        while(iSlot >= 0){
            const int* ch = chunk(iSlot);
            int cnt = ch[0];
            for(int i=0; i<cnt; i++) printf(" %i", ch[1+i]);
            iSlot = ch[slotSize - 1];
            if(iSlot >= 0) printf(" |");
        }
        printf(" ]\n");
    }

    void printStats() const {
        printf("NeighChunks: nOwner=%i nTot=%i nCap=%i slotSize=%i nMaxNeigh=%i\n", nOwner, nTot, nCap, slotSize, nMaxNeigh);
        int nExt = nTot - nOwner;
        int nOver = 0;
        for(int i=0; i<nOwner; i++){
            if(data[i*slotSize + slotSize-1] >= 0) nOver++;
        }
        printf("  owners with overflow: %i (%.1f%%), extension chunks: %i\n", nOver, 100.0*nOver/_max(1,nOwner), nExt);
    }
};

#endif
