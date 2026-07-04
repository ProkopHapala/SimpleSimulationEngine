
#ifndef Slots_h
#define Slots_h

/// @file Slots.h
/// @brief Fixed-capacity stack array for the specific case where each element has a small, bounded number of associations.
///
/// Born from the edge-loop sorting problem: when sorting N edges into a loop, each vertex
/// has at most 2 incident edges from that loop. Using a heap-allocated vector<int> per vertex
/// would be wildly wasteful (heap allocation, cache misses, fragmentation). Slots<int,2>
/// stores exactly 2 ints on the stack, with sentinel-based emptiness (-1).
///
/// The swap-with-last remove strategy preserves contiguity — important when iterating
/// the slots in order. This is a classic "small size optimization" applied to associative storage.

template <typename T, int N, T empty_=-1>
class Slots{ public:
    static const T   empty = empty_;
    static const int size  = N;
    T data[N]; // Raw array to store the slots

    Slots(){ fill(empty); }
    Slots( T v ){ fill(v); }

    void fill( T v){ for(int i=0; i<N; i++){data[i]=v;} }

    int count( ) const {
        for(int i=0; i<N;++i){ if(data[i]==empty)return i; }
        return N;
    }

    int find( T v ) const {
        for(int i=0; i<N;++i){ T vi=data[i]; if(vi==v){return i;}else if(vi==empty){return -1;} }
        return -1;
    }

    bool add( T v) {
        int i = count();
        if(i>=N){ return false; }
        data[i] = v;
        return true;
    }

    bool remove( T v){
        int i = find(v);
        if(i<0){ return false; }
        if(i<N-1){
            int j;
            for(j=i+1; j<N;++j){ if(data[j]==empty){break;} } // find the last non-empty slot
            data[i] = data[j];
            data[j] = empty;
        }else{
            data[i] = empty;
        }
        return true;
    }

};


// void addPointToEdge(LoopDict& map, int point, int edge) {
//     auto& slots = map[point];
//     bool success = slots.add(edge);
//     if (!success) {
//         std::cerr << "Warning: All slots are filled for point " << point << ". Edge " << edge << " was not added." << std::endl;
//     }
// }

// Define the unordered_map with std::string as the key and SlotArray<NUM_SLOTS> as the value
//using MyDataStructure = std::unordered_map<std::string, SlotArray<NUM_SLOTS>>;

#endif
