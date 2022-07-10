
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
        for(int i=0; i<nobj; i++ ){ 
            int ic = obj2cell[i]; //printf( "obj[%i ]-> cell %i \n", i, ic );
            cellNs[ ic ]++;  
        } 
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
    inline bool resizeObjs ( int nobj_, bool bO2C=false ){ bool b=(nobj_>nobjSize); nobj_bind=nobj_; if(b){ nobjSize =nobj_; _realloc(cell2obj,nobjSize); if(bO2C)_realloc(obj2cell_bind,nobj_); }; return b; };
    inline void bindObjs( int nobj_, int* obj2cell_ ){ resizeObjs ( nobj_, false ); obj2cell_bind=obj2cell_; }

    inline int getInCell( int icell, int* out ){    
        const int i0 = cellI0s[icell];
        const int n  = cellNs [icell];
        const int i1 = i0+n;
        for(int i=0; i<n; i++){
            //printf("getInCell(%i) i %i(%i|%i) \n", icell, i, i0,i1);
            //printf("getInCell[%i] out[%i]=%i \n", icell, i, cell2obj[i] );
            out[i] = cell2obj[i+i0];
        }
        return n;
    }

    inline void updateCells( int nobj=-1, int* obj2cell=0 ){
        if( obj2cell==0 ) obj2cell=obj2cell_bind;
        if( nobj<0      ) nobj=nobj_bind;
        clean();                   //printf( "DEBUG updateCells.clean() \n"  );
        count( nobj, obj2cell );   //printf( "DEBUG updateCells.count() \n"  );
        updateOffsets();           //printf( "DEBUG updateCells.updateOffsets() \n"  );
        objectsToCells( nobj, obj2cell ); //printf( "DEBUG updateCells.objectsToCells() \n"  );
    }

    inline void printCells(int verb=0){
        for(int i=0; i<ncell; i++){
            printf( "cell[%i] n %i i0 %i \n", i, cellNs[i], cellI0s[i] );
            if(verb>0){
                int i0=cellI0s[i];
                int ni=cellNs [i];
                int i1=i0+ni;
                for(int j=0; j<ni; j++){ printf( "cell2obj[%i|%i,%i] = %i \n", i0+j, i,j, cell2obj[i0+j] ); }
            }
        }
    }


};

#endif
