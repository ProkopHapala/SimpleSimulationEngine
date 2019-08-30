
#ifndef  binarySwitch_h
#define  binarySwitch_h

#include <cstdint>

/*
template <class T, uint8_t npow, size_t n >
size_t binarySwitch(T c, T* levels){
    size_t i = 0;
    for(uint8_t ipow=npow;ipow>=0;ipow--){
        const size_t bit = (1<<ipow);
        const size_t i_  = i|bit;
        printf( "binarySwitch[%i]  %i   ?[%i]: %i<%i \n", ipow,  bit,   i_, c, levels[i_] );
        if( (i_<n)&&(c>levels[i_]) ){ i=i_; }
    }
    return i;
};
*/

size_t binarySearch(int c, int8_t npow, size_t n, int* levels){
    size_t i = 0;
    for(int8_t ipow=npow;ipow>=0;ipow--){
        const size_t bit = (1<<ipow);
        const size_t i_  = i|bit;
        //printf( "binarySwitch[%i]  %i   ?[%i]: %i<%i \n", ipow,  bit,   i_, c, levels[i_] );
        if( (i_<n)&&(c>levels[i_]) ){ i=i_; }
    }
    return i;
};


int binarySwitch(const int c, const int8_t npow, const int n, const int* levels){
    int i = 0;
    for(int8_t ipow=npow;ipow>=0;ipow--){
        const int bit = (1<<ipow);
        const int i_  = i|bit;
        //printf( "binarySwitch[%i]  %i   ?[%i]: %i<%i \n", ipow,  bit,   i_, c, levels[i_] );
        const int l =levels[i_];
        if( c>=l ){
            if( c==l )return i_;
            i=i_;
        }
    }
    return -1;
}


template <class T, int8_t npow, int n >
int binarySwitchT(const T c,const T* levels){
    //printf( "binarySwitchT\n" );
    int i = 0;
    for(int8_t ipow=npow;ipow>=0;ipow--){
        //const int bit = (1<<ipow);
        //const int i_  = i|bit;
        const int i_  = i|(1<<ipow);
        //printf( "binarySwitch[%i]  %i   ?[%i]: %i<%i \n", ipow,  bit,   i_, c, levels[i_] );
        //printf( "binarySwitch[%i]    ?[%i]: %c<%c \n", ipow,   i_, c, levels[i_] );
        const int l   = levels[i_];
        if( c>=l ){
            if( c==l )return i_;
            i=i_;
        }
    }
    return -1;
}




template <class Tx,class Ty>
void bakeSwitchTable( int nmax, int ncase, const Tx* xs, const Ty* ys, Ty* table, Tx xoff, Ty yoff, Ty zero ){
    for(int i=0;i<nmax;i++){ table[i]=zero; }
    for(int i=0;i<nmax;i++){
        Tx x = i+xoff;
        for(int j=0;j<ncase;j++){
            if( xs[j] == x ){ table[i]=(ys[j]-yoff); break; };
        }
    }
}


#endif
