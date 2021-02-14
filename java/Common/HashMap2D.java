package Common;


//
// 2D map using Open-indexing Hashmap
//

public class HashMap2D extends HashMapMulti{
    double step =1.0;
    double invStep=1/step;
    //int i0=-100000.0;
    //int i0=-100000.0;
    
    public HashMap2D( int power_ ){
        super(power_);
    }
    
    public final void setStep( double step_ ){ step=step_; invStep=1/step; }
    
    public final long ibox2d( double x, double y ){
        //int ix =((int)(x*invStep))  +  2147483647;
        //int iy =((int)(y*invStep))  +  2147483647;
        //return (((long)iy)<<32) + ix;
        int ix =((int)(x*invStep))  + 32767;
        int iy =((int)(y*invStep))  + 32767; 
        //int ix =((int)(x*invStep));
        //int iy =((int)(y*invStep)); 
        return (((long)iy)<<16) + ix;
    }
    
    public final int hash( double x, double y ){
        return hash( ibox2d( x, y ) );
    }
    
    public final int insert( Object p, double x, double y ){
        long ibox = ibox2d( x, y );
        return insert( p, ibox );
    };
    
    public final int getAllInBox( double x, double y, Object[] outi ){
        long ibox = ibox2d( x, y );
        return getAllInBox( ibox, outi );
    }
    
}
