
package Common;

public class CellSort2D_ind extends CellSort {
double step =1.0;
double invStep=1/step;
  int nx,ny;
  //int[] IDs;
  int nl0;
  
public CellSort2D_ind( int w,int h, int step_, int nid_  ){ 
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

public final int getIx(double x){ return (int)(x*invStep); }
public final int getIy(double y){ return (int)(y*invStep); }

public final void insert( double x, double y ){
    int ix = getIx(x);
    int iy = getIy(y);
    int ic = nx*iy + ix;
    if( ic>=ncell ) System.out.println("ix,iy "+ix+","+iy+" ic "+ic);
    //id2cell[i] =
    insert(ic);
}
  
public final void army2cells( GetXY [] ps ){
    nid=0;
    for(int i=0; i<ps.length; i++){
      int ix = getIx( ps[i].getX() );
      int iy = getIy( ps[i].getY() ); 
      id2cell[i] = nx*iy + ix;
      nid++;
    }
  }

public final int loadCell( int il0, int[] loaded, int ix, int iy ){
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
  
public final int loadNeighs( int ix, int iy, int[] tmp ){
    nl0=loadCell( 0, tmp, ix,iy );
    int nl=nl0;

    nl+=loadCell( nl, tmp, ix-1,iy-1 );
    nl+=loadCell( nl, tmp, ix  ,iy-1 );
    nl+=loadCell( nl, tmp, ix+1,iy-1 );

    nl+=loadCell( nl, tmp, ix-1,iy   );
    nl+=loadCell( nl, tmp, ix+1,iy   );

    nl+=loadCell( nl, tmp, ix-1,iy+1 );
    nl+=loadCell( nl, tmp, ix  ,iy+1 );
    nl+=loadCell( nl, tmp, ix+1,iy+1 );

    return nl;
}

public final int loadNeighs( double x, double y, int[] tmp ){
    int ix = getIx(x);
    int iy = getIx(y);
    return loadNeighs( ix, iy, tmp );
}
  
}; // class CellSort2D{
