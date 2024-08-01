
#ifndef  sweep_h
#define  sweep_h

struct Int2{
    int i,j;
};

inline bool anyOf( int i, int ia, int ib ){ return (i==ia)||(i==ib); };

/*
int DEBUG_ib1=;
int DEBUG_ip1=;
int DEBUG_ib2=;
int DEBUG_ip2=;
*/

//bool bDEBUG = false;

namespace sweep{

struct Span{
    float xmin,xmax;
    inline bool operator<( const Span& rhs) const { return xmin < rhs.xmin; }
    inline bool operator>( const Span& rhs) const { return xmin > rhs.xmin; }
    //inline bool operator>( const Span& rhs){ return xmax > rhs.xmax; }
};

int collideSelf(int n, int* permut, Span* bounds, Int2* cols){
    int ncol=0;
    for( int ii=0; ii<n; ii++){
        //printf( "sweep::collideSelf ii %i \n", ii );
        int i = permut[ii];
        //printf( "---sweep::collideSelf ii %i \n", ii );
        const Span& spani = bounds[i];
        //printf( "--sweep::collideSelf ii %i \n", ii );
        for( int jj=ii+1; jj<n; jj++){
            //printf( "sweep::collideSelf ii %i jj %i \n", ii, jj );
            int j = permut[jj];

            //if(bDEBUG){
            //    printf( "DEBUG branch-pair: i,j(%i,%i)  |  ii,jj (%i,%i) \n", i, j,  ii, jj );
            //    //if( anyOf(i,0,3) && anyOf(j,0,3) ){ printf( "DEBUG found branch-pair: (%i,%i) \n", i, j  ); }
            //    if( anyOf(i,12,58) && anyOf(j,12,58) ){ printf( "DEBUG found branch-pair: (%i,%i) \n", i, j  ); }
            //}

            if ( bounds[j].xmin > spani.xmax ) break; // cannot colide, and the later j also not
            //printf( "sweep::collideSelf i,j %i,%i \n", i, j );
            cols[ncol] = {i,j};
            ncol++;
        }
    }
    //printf( "collideSelfObjects : ncomps %i ncadidates : %i nhits %i | nBrute %i \n", ncomp, ncol, nhit, n*(n-1)/2 );
    return ncol;
}

int collideCrossOneWay( int ni, int* permuti, Span* boundsi,  int nj, int* permutj, Span* boundsj, Int2* cols,  bool inv=false ){
    //printf( "DEBUG collideCrossOneWay %i \n", inv );
    if(nj<=0)return 0;
    int ncol = 0;
    int jj0  = 0;
    for(int ii=0; ii<ni; ii++ ){
        //printf( "DEBUG ii %i  permuti+ii %i \n", permuti+ii );
        int i= permuti[ii];
        const Span& spani = boundsi[i];
        //printf( "DEBUG ii %i jj0 %i permutj[jj0] %i \n", ii, jj0, permutj[jj0] );
        while( boundsj[ permutj[jj0] ].xmin<spani.xmin ){
            jj0++;
            //printf( "DEBUG jj0 %i  permutj+jj0 %i \n", permutj+jj0 );
            if(jj0>=nj) return ncol;
        }
        int jj = jj0;
        while( boundsj[ permutj[jj] ].xmin<spani.xmax ){
            //if(bDEBUG){
            //    printf( "DEBUG i,j(%i,%i)  |  ii,jj (%i,%i) \n", i, permutj[jj],  ii, jj );
            //}
            //printf( "DEBUG  ncol %i  :  i,j(%i,%i) \n", ncol, i,permutj[jj] );
            if(inv){ cols[ncol] = {i,permutj[jj]}; }else{ cols[ncol] = {permutj[jj],i}; };
            ncol++;
            jj++;
            if(jj>=nj) break;
        }
    }
    return ncol;
}

int collideCross( int ni, int* permuti, Span* boundsi,  int nj, int* permutj, Span* boundsj, Int2* cols ){
    int ncol;
    ncol  =collideCrossOneWay( ni, permuti, boundsi,  nj, permutj, boundsj, cols     , false );
    ncol +=collideCrossOneWay( nj, permutj, boundsj,  ni, permuti, boundsi, cols+ncol, true  );
    return ncol;
}

};

#endif
