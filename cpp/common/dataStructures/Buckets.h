
#ifndef  Buckets_h
#define  Buckets_h

#include "macroUtils.h"

class Buckets{ public:
    int ncell,nobj;
    int* cellNs=0;   //[ncell] number of objects contained in given cell, it may be not necessary as it can be computed as cellI0s[i+1]-cellI0s[i], but it is more convenient to have it here
    int* cellI0s=0;  //[ncell] index of first object contained given cell in the array cell2obj
    int* cell2obj=0; //[nobj]  indices of objects contained in given cell

    int  maxInBucket=0;
    
    int  nobjSize=-1;
    int* obj2cell=0; 

    inline void clean(){ for(int k=0; k<ncell; k++ ){ cellNs[k]=0; } }
    inline void cleanO2C( int icell=-1 ){ for(int i=0; i<nobj; i++ ){ obj2cell[i]=icell; } }

    /**
     * Counts the number of objects in each cell based on the given mapping. It stores the result in the array cellNs.
     *
     * @param nobj The number of objects.
     * @param obj2cell An array mapping each object to its corresponding cell.
     */
    inline void count( int nobj, int* obj2cell ){ 
        for(int i=0; i<nobj; i++ ){ 
            int ic = obj2cell[i]; printf( "obj[%i ]-> cell %i \n", i, ic );
            cellNs[ ic ]++;  
        } 
    }

    /**
     * @brief Updates the offsets of the buckets (cellI0s) based on the number of elements in each bucket (cellNs). It also determines the maximum number of elements in a bucket (maxInBucket).
     * 
     * This function updates the offsets of the buckets based on the number of elements in each bucket.
     * It also determines the maximum number of elements in a bucket.
     */
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

    /**
     * Assigns objects to cells based on the given mapping. The result is stored in the array cell2obj.
     * 
     * @param nobj The number of objects.
     * @param obj2cell An array representing the mapping of objects to cells.
     */
    inline void objectsToCells( int nobj, int* obj2cell ){ 
        for(int i=0; i<nobj; i++){
            //printf( "[%i] ", i );
            int k =  obj2cell[i];
            int& ni = cellNs[k];
            int j =  cellI0s[k] + ni;
            //printf( " k %i j %i | nobj %i \n", k, j, nobj );
            cell2obj[j] = i;
            ni++;
        }
    }


    /**
     * Updates the cells based on the given object-to-cell mapping.
     * This function performs the following steps:
     * 1. Cleans the cells.
     * 2. Counts the number of objects in each cell.
     * 3. Updates the offsets of the cells.
     * 4. Assigns objects to their respective cells.
     *
     * @param nobj_     number of objects.           If not specified, it uses the default number of objects.
     * @param obj2cell_ object-to-cell mapping.  If not specified, it uses the default mapping. (which have to be allocated or binded before).
     */
    inline void updateCells( int nobj_=-1, int* obj2cell_=0 ){
        if( obj2cell_==0 ) obj2cell_=obj2cell;
        if( nobj_<0      ) nobj_    =nobj;
        clean         ();                      printf("updateCells.clean() \n"  );
        count         ( nobj_, obj2cell_ );    printf("updateCells.count() \n"  );
        updateOffsets ();                      printf("updateCells.updateOffsets() \n"  );
        objectsToCells( nobj_, obj2cell_ );    printf("updateCells.objectsToCells() \n"  );
    }

    inline bool resizeCells( int ncell_                 ){ bool b=(ncell_>ncell); if(b){ ncell=ncell_; _realloc(cellNs,ncell_); _realloc(cellI0s,ncell_); }; return b; };
    inline bool resizeObjs ( int nobj_, bool bO2C=false ){ bool b=(nobj_>nobjSize); nobj=nobj_; if(b){ nobjSize =nobj_; _realloc(cell2obj,nobjSize); if(bO2C)_realloc(obj2cell,nobj_); }; return b; };
    inline void bindObjs   ( int nobj_, int* obj2cell_  ){ resizeObjs ( nobj_, false ); obj2cell=obj2cell_; }

    inline void realloc( int ncell_, int nobj_, bool bO2C=false ){ 
        ncell=-1; resizeCells(ncell_); 
        nobj=-1;  resizeObjs(nobj_,bO2C); 
    }

    inline int getInCell( int icell, int* out ){    
        // why we copy? we can just return pointer to cell2obj[ cellI0s[icell] ]
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
    inline int* inCell( int icell ){ return cell2obj + cellI0s[icell]; }

    inline void printObjCellMaping(int i0=-1, int i1=-1, bool bBoundary=true){
        if(i0<0)i0=0; 
        if(i1<0)i1=nobj;
        printf( "Buckets::printObjCellMaping() ncell %i nobj %i \n", ncell, nobj );
        int ib=0;
        for(int i=i0; i<i1; i++){
            if( bBoundary & (i==cellI0s[ib])){ printf( "---- start cell[%i] \n", i, ib ); ib++; }
            printf( "[%i] o2c %i c2o %i \n", i, obj2cell[i], cell2obj[i] );
        }
    }

    inline void printCells(int verb=0){
        printf( "Buckets::printCells() ncell %i nobj %i \n", ncell, nobj );
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
