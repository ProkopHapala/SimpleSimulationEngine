package Common;

public class NBody_array {
    public double [] ps;
    public double [] vs;
    public double [] fs;
    public double m = -1.0;
    
    public NBody_array(int n){
        ps = new double[n*2];
        vs = new double[n*2];
        fs = new double[n*2];
    }
    
    public void evalForce(){
        int n = ps.length/2;
        for(int i=0; i<n; i++){
            final int ix=i*2;
            final int iy=ix+1;
            double pix = ps[ix];
            double piy = ps[iy];
            for(int j=i+1; j<n; j++){
                final int jx=j*2;
                final int jy=jx+1;
                double dx=ps[jx]-pix;
                double dy=ps[jy]-piy;
                double r  = Math.sqrt( dx*dx + dy*dy );
                //System.out.println( "ps["+i+","+j+"] "+r );
                double fr = m/( 1e-8 + r*r*r);
                    dx*=fr;     dy*=fr;
                fs[jx]+=dx; fs[jy]+=dy;
                fs[ix]-=dx; fs[iy]-=dy;
            }
        }
    }
    
    public void move_leapFrog(double dt){
        int n = ps.length/2;
        for(int i=0; i<n; i++){
            final int ix=i*2;
            final int iy=ix+1;
            double vx = vs[ix];
            double vy = vs[iy];
            vx += fs[ix]*dt;
            vy += fs[iy]*dt;
            ps[ix] += vx* dt;
            ps[iy] += vy* dt;
            vs[ix]=vx;
            vs[iy]=vy;
        }
    }
    
    public void update( double dt){
        for(int i=0; i<ps.length; i++){ fs[i]=0; }
        evalForce();
        move_leapFrog(dt);
    }
    
}
