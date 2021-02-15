
package Common;

public class CellSort2D extends CellSort implements Map2D{
double step =1.0;
double invStep=1/step;
  int nx,ny;
  //int[] IDs;
  int nl0;
  
public CellSort2D( int w,int h, int step_, int nid_  ){ 
      super( (w/step_)*(h/step_), nid_ ); 
      setStep( step_ ); nx=(w/step_); ny=(h/step_); 
  }
  
public final void setStep( double step_ ){ step=step_; invStep=1/step; }
  
public final int getIx(double x){ return (int)(x*invStep); }
public final int getIy(double y){ return (int)(y*invStep); }

@Override
public final void insert( double x, double y ){
    int ix = getIx(x);
    int iy = getIy(y);
    int ic = nx*iy + ix;
    if( ic>=ncell ) System.out.println("ix,iy "+ix+","+iy+" ic "+ic);
    //id2cell[i] =
    insert(ic);
}
  
@Override
public final void insert( GetXY [] ps ){
    nid=0;
    for(int i=0; i<ps.length; i++){
      int ix = getIx( ps[i].getX() );
      int iy = getIy( ps[i].getY() ); 
      id2cell[i] = nx*iy + ix;
      nid++;
    }
  }

@Override
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
  
@Override
public final int loadNeighs( int ix, int iy, int[] tmp, boolean bAlways ){
           nl0      =loadCell( 0,  tmp, ix,   iy  );
    if( (nl0==0) && (!bAlways) ) return 0;
    int nl=nl0;
    int mx=nx-1;
    if(iy>0){
        if(ix>0) nl+=loadCell( nl, tmp, ix-1,iy-1 );
                 nl+=loadCell( nl, tmp, ix  ,iy-1 );
        if(ix<mx)nl+=loadCell( nl, tmp, ix+1,iy-1 );
    }
        if(ix>0 )nl+=loadCell( nl, tmp, ix-1,iy   );
    //           nl0=loadCell( 0,  tmp, ix,   iy  ); 
        if(ix<mx)nl+=loadCell( nl, tmp, ix+1,iy   );
    if(iy<(nx-1)){
        if(ix>0) nl+=loadCell( nl, tmp, ix-1,iy+1 );
                 nl+=loadCell( nl, tmp, ix  ,iy+1 );
        if(ix<mx)nl+=loadCell( nl, tmp, ix+1,iy+1 );
    }
    return nl;
}

@Override
public final int loadNeighs( double x, double y, int[] tmp, boolean bAlways ){
    int ix = getIx(x);
    int iy = getIx(y);
    return loadNeighs( ix, iy, tmp, bAlways );
}

public int interactNeighs( Interactor inter, int[] tmp ){
    int nint=0;
    for(int iy=0; iy<ny; iy++){
        for(int ix=0; ix<ny; ix++){
            int n = loadNeighs( ix, iy, tmp, false );
            nint += inter.interactNeighs_inds( nl0, n, tmp );
        }
    }
    return nint;
}
  
}; // class CellSort2D{
