
#ifndef  allocation_manager_h
#define  allocation_manager_h

//#include <string>
//#include <vector>
#include <unordered_map>

struct AllocationRecord{
    void* ptr=0;
    int   n;
    int   bytes;
    
    void dealoc(){ if(n==1){ delete ptr; }else if(n>1){ delete[] ptr; } }
    void print(){ printf( "AllocationRecord[%p] size=%i func=%s\n", ptr, n, bytes ); }
};

std::unordered_map<void*,AllocationRecord> allocationDict;

inline void registerPointer_( void* o, int perItem, int n=1 ){
    if(o==0){ printf( "ERROR: registerPointer_() pointer is null\n" ); exit(0); }    
    AllocationRecord rec;
    rec.ptr   = o;
    rec.n     = n;
    rec.bytes = n*perItem;
    allocationDict.insert( {o,rec} );
}

template<typename T> void registerPointer( T*& o,  int n=1 ){ registerPointer_( (void*)o, sizeof(T), n ); }

/**
 * @brief Allocates memory for an object or an array of objects and keeps track of the allocation.
 * 
 * @tparam T The type of the object(s) to allocate.
 * @param o A reference to a pointer that will store the allocated memory.
 * @param n The number of objects to allocate (default is 1 meaning a single object not an array ).
 */
template<typename T> void allocMan( T*& o, int n=1 ){      
    if( n>1 ){ o = new T[n]; }else{ o = new T; }
    registerPointer_( o, sizeof(T), n );
    //printf( "allocMan[%p] size=%i line=%i file=%s func=%s\n", o, rec.size ); 
};

/**
 * Tries to deallocate the memory pointed to by the given pointer.
 * 
 * @param o The pointer to the memory to deallocate.
 * @return True if the deallocation was successful, false otherwise.
 */
inline bool tryDeallocMan( void*& o ){      
    bool ret;
    auto it = allocationDict.find(o);
    if( it == allocationDict.end() ){  
        ret = false; }
    else{
        ret = true; 
        AllocationRecord& rec = it->second;
        rec.dealoc();
        allocationDict.erase(o);        
    }
    o = 0;
    return ret;
    //printf( "deallocMan[%p] size=%i line=%i file=%s func=%s\n", o, rec.size ); 
};

#endif
