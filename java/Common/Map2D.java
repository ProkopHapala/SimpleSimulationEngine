
package Common;

public interface Map2D{
    
//public void setStep( double step_ );
//public int getIx(double x;
//public int getIy(double y);

public void insert    ( double x, double y );
public void insert    ( GetXY [] ps );
public int  loadCell  ( int il0, int[] loaded, int ix, int iy );
public int  loadNeighs( int ix, int iy, int[] tmp, boolean bAlways );
public int  loadNeighs( double x, double y, int[] tmp, boolean bAlways );

};