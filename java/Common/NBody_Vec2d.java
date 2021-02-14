package Common;

final public class NBody_Vec2d {
    public Vec2d [] ps;
    public Vec2d [] vs;
    public Vec2d [] fs;
    public double m = -1.0;
    
    public NBody_Vec2d(int n){
        ps = new Vec2d[n];
        vs = new Vec2d[n];
        fs = new Vec2d[n];
        for(int i=0; i<n; i++){ ps[i]=new Vec2d(); };
        for(int i=0; i<n; i++){ vs[i]=new Vec2d(); };
        for(int i=0; i<n; i++){ fs[i]=new Vec2d(); };
    }
    
    public void evalForce(){
        Vec2d d = new Vec2d();
        for(int i=0; i<ps.length; i++){
            Vec2d pi = ps[i];
            Vec2d fi = fs[i];
            for(int j=i+1; j<ps.length; j++){
                d.set_sub( ps[j], pi );
                double r  = d.norm();
                //System.out.println( "ps["+i+","+j+"] "+r );
                double fr = m/( 1e-8 + r*r*r);
                d.mul(fr);
                fs[j].add(d);
                fi   .sub(d);
            }
        }
    }
    
    public void evalForce_Naive(){
        for(int i=0; i<ps.length; i++){
            Vec2d pi = ps[i];
            Vec2d fi = fs[i];
            for(int j=i+1; j<ps.length; j++){
                Vec2d  d  = Vec2d.sub( ps[j], pi );
                double r  = d.norm();
                double fr = m/( 1e-8 + r*r*r);
                //d.mul(fr);  fs[j].add(d); fi   .sub(d);
                Vec2d f = Vec2d.mul( d, fr ); fs[j].add(f); fi   .sub(f);
            }
        }
    }
    
    public void move_leapFrog(double dt){
        for(int i=0; i<ps.length; i++){
            Vec2d v = vs[i];
            //v.mul(1-drag*dt);
            v    .fma( fs[i], dt );
            ps[i].fma(     v, dt );
        }
    }
    
    public void update( double dt, boolean bNaive ){
        for(int i=0; i<ps.length; i++){ fs[i].set(0); }
        if(bNaive){
            evalForce_Naive();
        }else{
            evalForce();
        }
        move_leapFrog(dt);
    }
    
}
