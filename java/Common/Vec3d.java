package Common;

final public class Vec3d{
    public double x,y,z;
    
    double norm2(){ return            x*x + y*y + z*z  ; };
    double norm (){ return Math.sqrt( x*x + y*y + z*z ); };
    
    final public void   set( double f ){ x =f; y =f; z =f; };
    final public void   add( double f ){ x+=f; y+=f; z+=f; };
    final public void   mul( double f ){ x*=f; y*=f; z*=f; };
    
    final public double dot( double fx, double fy, double fz ){ return x*fx+y*fy+z*fz; };
    final public void   set( double fx, double fy, double fz ){ x =fx; y =fy; z =fz; };
    final public void   add( double fx, double fy, double fz ){ x+=fx; y+=fy; z+=fz; };
    final public void   mul( double fx, double fy, double fz ){ x*=fx; y*=fy; z*=fz; };
    
    final public double dot( Vec3d v ){ return x*v.x+y*v.y+z*v.z; };
    final public void   add( Vec3d v ){ x+=v.x; y+=v.y; z+=v.z; };
    final public void   mul( Vec3d v ){ x*=v.x; y*=v.y; z*=v.z; };
    final public void   sub( Vec3d v ){ x-=v.x; y-=v.y; z-=v.z; };
    final public void   div( Vec3d v ){ x/=v.x; y/=v.y; z/=v.z; };
    
    final public void   fma( Vec3d v, double f ){ x+=v.x*f; y+=v.y*f; z+=v.z*f; };
    
}
