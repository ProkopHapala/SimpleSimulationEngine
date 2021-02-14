
package Common;

public class CellSort2D_ind extends CellSort {
double step =1.0;
double invStep=1/step;
  int nx,ny;
  //int[] IDs;
  int nl0;
  
  CellSort2D_ind( int w,int h, int step_, int nid_  ){ 
      super( (w/step_)*(h/step_), nid_ ); 
      setStep( step_ ); nx=(w/step_); ny=(h/step_); 
  }
  
  public final void setStep( double step_ ){ step=step_; invStep=1/step; }
  
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
  
void insert( int i, double x, double y ){
    int ix = (int)(x*invStep);
    int iy = (int)(y*invStep);
    id2cell[i] = nx*iy + ix;
}
  
  void army2cells( GetXY [] ps ){
    nid=0;
    for(int i=0; i<ps.length; i++){
      int ix = (int)(ps[i].getX()*step);
      int iy = (int)(ps[i].getY()*step);
      id2cell[i] = nx*iy + ix;
      nid++;
    }
  }

  int toArray( int il0, int[] loaded, int ix, int iy ){
    int ic = iy*nx+ix;
    if( (ix>0)&&(ix<nx) && (iy>0)&&(iy<ny) ){
      int i0 = cellIs[ic];
      int nl = cellNs[ic];
      for(int i=0; i<nl; i++ ){
        int id = cell2id[i0+i];
        loaded[il0+i]=id;
      }
      return nl;
    }else{
      return 0;
    }
  }
  
  int loadNeighs( int[] tmp, int ix, int iy, boolean bAlways ){
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
