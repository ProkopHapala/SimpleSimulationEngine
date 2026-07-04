
#ifndef  synchronize_h
#define  synchronize_h

/// @file synchronize.h
/// @brief Double-buffer synchronizer for ping-pong time-stepped simulations.
///
/// Ptr2<T> holds paired read/write pointers; Synchronizer maintains a list of such pairs
/// and copies y→x on read(), x→y on write(). This decouples the "which buffer is current"
/// logic from the simulation code — the solver always reads from x and writes to y,
/// then calls read() to swap. (Note: currently a prototype — syntax errors present.)

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
