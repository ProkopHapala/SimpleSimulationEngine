
#ifndef  sweep_h
#define  sweep_h

struct Int2{
    int i,j;
};

namespace sweep{

struct Span{
    float xmin,xmax;
    inline bool operator<( const Span& rhs){ return xmin < rhs.xmin; }
    inline bool operator>( const Span& rhs){ return xmax < rhs.xmax; }
};

int collideSelf(int n, int* permut, Span* bounds, Int2* cols){
    int ncol=0;
    for( int ii=0; ii<n; ii++){
        int i= permut[ii];
        const Span& spani = bounds[i];
        for( int jj=ii+1; jj<n; jj++){
            int j = permut[jj];
            if ( bounds[j].xmin > spani.xmax ) break; // cannot colide, and the later j also not
            cols[ncol] = {i,j};
            ncol++;
        }
    }
    //printf( "collideSelfObjects : ncomps %i ncadidates : %i nhits %i | nBrute %i \n", ncomp, ncol, nhit, n*(n-1)/2 );
    return ncol;
}

int collideCross( int ni, int* permuti, Span* boundsi,  int nj, int* permutj, Span* boundsj, Int2* cols,  bool inv=false ){
    int ncol=0;
    for(int ii=0; ii<ni; ii++ ){
        int i= permuti[ii];
        const Span& spani = boundsi[i];
        int jj0 = 0;
        while( boundsj[ permutj[jj0] ].xmin<spani.xmin ){ jj0++; if(jj0>=nj) return ncol;  }
        int jj = jj0;
        while( boundsj[ permutj[jj] ].xmin<spani.xmax ){
            if(inv){ cols[ncol] = {i,permutj[jj]}; }else{ cols[ncol] = {permutj[jj],i}; };
            ncol++;
            jj++;
            if(jj>=nj) break;
        }
    }
    return ncol;
}

};

#endif
