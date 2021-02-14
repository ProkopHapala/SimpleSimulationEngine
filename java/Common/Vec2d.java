package Common;


final public class Vec2d implements GetXY{
    public double x,y;
    
    final public double getX(){ return x; };
    final public double getY(){ return y; };
    
    Vec2d(){};
    Vec2d(double f           ){x=f ; y=f ;};
    Vec2d(double fx,double fy){x=fx; y=fy; };
    Vec2d(Vec2d v            ){x=v.x;y=v.y;};
    
    double norm2(){ return            x*x + y*y;   };
    double norm (){ return Math.sqrt( x*x + y*y ); };
    
    final public void   set( double f ){ x =f; y =f; };
    final public void   add( double f ){ x+=f; y+=f; };
    final public void   mul( double f ){ x*=f; y*=f; };
    
    final public double dot( double fx, double fy ){ return x*fx+y*fy; };
    final public void   set( double fx, double fy ){ x =fx; y =fy;  };
    final public void   add( double fx, double fy ){ x+=fx; y+=fy;  };
    final public void   mul( double fx, double fy ){ x*=fx; y*=fy;  };
    
    final public double dot( Vec2d v ){ return x*v.x+y*v.y; };
    final public void   add( Vec2d v ){ x+=v.x; y+=v.y; };
    final public void   mul( Vec2d v ){ x*=v.x; y*=v.y; };
    final public void   sub( Vec2d v ){ x-=v.x; y-=v.y; };
    final public void   div( Vec2d v ){ x/=v.x; y/=v.y; };
    
    final public void   fma( Vec2d v, double f ){ x+=v.x*f; y+=v.y*f; };
    
    final public void set_sub( Vec2d a, Vec2d b ){ x=a.x-b.x; y=a.y-b.y; };
    
    
    //static Vec2d ofX( float x ){ return new Vec2d };
    //static Vec2d ofY( float y ){ };
    //static Vec2d  of ( float x, float y, float z ){ return  };
    final public static Vec2d add( Vec2d a, Vec2d b ){ return new Vec2d( a.x+b.x, a.y+b.y ); };
    final public static Vec2d sub( Vec2d a, Vec2d b ){ return new Vec2d( a.x-b.x, a.y-b.y ); };
    final public static Vec2d mul( Vec2d a, Vec2d b ){ return new Vec2d( a.x*b.x, a.y*b.y ); };
    final public static Vec2d div( Vec2d a, Vec2d b ){ return new Vec2d( a.x/b.x, a.y/b.y ); };
    
    final public static Vec2d add( Vec2d a, double f ){ return new Vec2d( a.x+f, a.y+f ); };
    final public static Vec2d mul( Vec2d a, double f ){ return new Vec2d( a.x*f, a.y*f ); };
    
}
