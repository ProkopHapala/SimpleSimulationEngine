
#ifndef  synchronize_h
#define  synchronize_h

#include <vector>
#include "macroUtils.h"

template<typename T>
struct Ptr2<T>{
    union{
        struct{T *x,*y};
        struct{T *a,*b};
        struct{T *i,*j};
        T* arr[2];
    }
}

template<typename T>
class Synchronizer(){
    std::vector<Ptr2> links;
    read (){ for(Ptr2& p:links){ p.x=p.y; } };
    write(){ for(Ptr2& p:links){ p.y=p.x; } };
}

#endif
