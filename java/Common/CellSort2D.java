package Common;


public class CellSort2D<T extends GetXY> extends CellSort {
  double dcell;
  double step;
  int nx,ny;
  T[] army;
  //Ant[] tmp;   // temporary array for loading
  int nl0;
  
  CellSort2D( T[] army_, int w,int h, int dcell_ ){ super( (w/dcell_)*(h/dcell_), army_.length ); dcell=dcell_; step=1./dcell; nx=(w/dcell_); ny=(h/dcell_); army=army_;  }
  
  /*
  void drawCells(){
    noFill();
    float d=1./step;
    for(int iy=0; iy<ny; iy++){
      for(int ix=0; ix<nx; ix++){
        int ic = iy*nx+ix;
        if(cellNs[ic]>0){
          rect(ix*d,iy*d,(ix+1)*d,(iy+1)*d);
          int i0=cellIs[ic];
          for(int j=0;j<cellNs[ic];j++){
            //print(" "+cell2id[i0+j]);
            int id = cell2id[i0+j];
            line( (ix+.5)*d, (iy+.5)*d, army[id].getX(), army[id].getY() );
          }
          //println();
        }
      }
    }
  }
*/
  
  void army2cells( ){
    nid=0;
    for(int i=0; i<army.length; i++){
      int ix = (int)(army[i].getX()*step);
      int iy = (int)(army[i].getY()*step);
      id2cell[i] = nx*iy + ix;
      nid++;
    }
  }

  int toArray( int il0, T[] loaded, int ix, int iy ){
    int ic = iy*nx+ix;
    if( (ix>0)&&(ix<nx) && (iy>0)&&(iy<ny) ){
      int i0 = cellIs[ic];
      int nl = cellNs[ic];
      for(int i=0; i<nl; i++ ){
        int id = cell2id[i0+i];
        loaded[il0+i]=army[id];
      }
      return nl;
    }else{
      return 0;
    }
  }
  
  int loadNeighs( T[] tmp, int ix, int iy, boolean bAlways ){
    nl0=toArray( 0, tmp, ix-1,iy-1 );
    int nl=nl0;
    if( bAlways && (nl0>0) ){
      nl+=toArray( nl, tmp, ix-1,iy-1 );
      nl+=toArray( nl, tmp, ix  ,iy-1 );
      nl+=toArray( nl, tmp, ix+1,iy-1 );
      
      nl+=toArray( nl, tmp, ix-1,iy-1 );
      nl+=toArray( nl, tmp, ix+1,iy-1 );
      
      nl+=toArray( nl, tmp, ix-1,iy+1 );
      nl+=toArray( nl, tmp, ix  ,iy+1 );
      nl+=toArray( nl, tmp, ix+1,iy+1 );
    }
    return nl;
  }
  
}; // class CellSort2D{