package Common;


///////////////////////////////
//  class      CellSort
//  ==========================
//  accelerate local interaction by mapping multiple ants into grid cells
///////////////////////////////

import java.lang.System;

public class CellSort {

  int nid;
  int ncell;
  int [] id2cell;
  int [] cellNs;
  int [] cellIs;
  int [] cell2id;
  
public CellSort(int ncell_, int nid_ ){
    System.out.println("ncell "+ncell_+" nid "+nid_);
    id2cell=new int[nid_  ];
    cellNs =new int[ncell_];
    cellIs =new int[ncell_];
    cell2id=new int[nid_];
    ncell=ncell_;
    //nid  =nid_; 
    //nid=0;
    clean();
  }

public void clean(){
    nid=0;
}

public final void insert( int ic ){
    id2cell[nid]=ic;
    nid++;
}
  
public final void sort(){
    // prepare cell regions
    for(int i=0; i<nid; i++){
      int ic = id2cell[i];
      cellNs[ic]++;
    }
    int n=0;
    for(int i=0; i<ncell; i++){
      cellIs[i]=n;
      int ni = cellNs[i];
      n+=ni;
      System.out.println( i+": i0 "+cellIs[i]+" ni "+ni );
      cellNs[i] = 0;
    }
    // put to cells
    for(int i=0; i<nid; i++){
      int ic = id2cell[i];
      int ni = cellNs[ic];
      if((cellIs[ic] + ni)>ncell)System.out.println( " i0 "+cellIs[ic]+" ni "+ni );
      cell2id[ cellIs[ic] + ni ] = i;
      cellNs[ic]=ni+1;
    }
  }
  
public final void printCells(){
    for(int i=0; i<ncell; i++){
      if(cellNs[i]>0){
      System.out.print("["+i+"]");
      int i0=cellIs[i];
      for(int j=0;j<cellNs[i];j++){
        System.out.print(" "+cell2id[i0+j]);
      }
      System.out.println();
      }
    }
  }
  
}; // class CellSort{    
