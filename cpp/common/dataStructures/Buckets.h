
#ifndef  Buckets_h
#define  Buckets_h

#include "macroUtils.h"

class Buckets{ public:
    int ncell;
    int* cellNs=0;   //[ncell]
    int* cellI0s=0;  //[ncell]
    int* cell2obj=0; //[nobj]

    int  maxInBucket=0;
    
    int  nobjSize=-1;
    int  nobj_bind=-1;
    int* obj2cell_bind=0; 

    inline void clean(){ for(int k=0; k<ncell; k++ ){ cellNs[k]=0; } }
    inline void count( int nobj, int* obj2cell ){ 
        for(int i=0; i<nobj; i++ ){ cellNs[ obj2cell[i] ]++;  } 
    }
    inline void updateOffsets(){   
        int ntot=0;
        maxInBucket=0;
        for(int k=0; k<ncell; k++ ){
            cellI0s[k] = ntot;
            int ni    = cellNs[k];
            ntot      += ni;
            if( ni>maxInBucket ) maxInBucket=ni;     // ? check most occupied cell
            cellNs[k]=0;
        }
    }

    inline void objectsToCells( int nobj, int* obj2cell ){ 
        for(int i=0; i<nobj; i++){
            //printf( "[%i] ", i );
            int k =  obj2cell[i];
            int j =  cellI0s[k] + cellNs[k];
            //printf( " k %i j %i | nobj %i \n", k, j, nobj );
            cell2obj[j] = i;
            cellNs[k]++;
        }
    }

    inline bool resizeCells( int ncell_ ){ bool b=(ncell_>ncell); if(b){ ncell=ncell_; _realloc(cellNs,ncell_); _realloc(cellI0s,ncell_); }; return b; };
    inline bool resizeObjs ( int nobj_, bool bO2C=false ){ bool b=(nobj_>nobjSize); nobj_bind=nobj_; if(b){ nobjSize =nobj_; _realloc(cell2obj,nobjSize); if(bO2C)_realloc(obj2cell_bind,nobj); }; return b; };
    inline void bindObjs( int nobj_, int* obj2cell_ ){ resizeObjs ( nobj_, false ); obj2cell_bind=obj2cell_; }

    inline int getInCell( int icell, int* out ){    
        const int i0 = cellI0s[icell];
        const int n  = cellNs [icell];
        const int i1 = i0+n;
        for(int i=i0; i<i1; i++){
            out[i] = cell2obj[i];
        }
        return n;
    }

    inline void updateCells( int nobj=-1, int* obj2cell=0 ){
        if( obj2cell==0 ) obj2cell=obj2cell_bind;
        if( nobj<0      ) nobj=nobj_bind;
        clean();
        count( nobj, obj2cell );
        updateOffsets();
        objectsToCells( nobj, obj2cell );
    }


};

#endif
