package Common;

import java.util.Set;



final public class NBody implements Interactor {
    public Vec2d [] ps;
    public Vec2d [] vs;
    public Vec2d [] fs;
    public double m = -1.0;
    
    TempPool pool;
    public double collisionRadius=20.0;
    public double collisionStiffness=0.1;
    public boolean bDEBUG=false;
    
    public NBody(int n, TempPool pool_){
        pool=pool_;
        ps = new Vec2d[n];
        vs = new Vec2d[n];
        fs = new Vec2d[n];
        for(int i=0; i<n; i++){ ps[i]=new Vec2d(); };
        for(int i=0; i<n; i++){ vs[i]=new Vec2d(); };
        for(int i=0; i<n; i++){ fs[i]=new Vec2d(); };
    }
    
    public void cleanForce(){
        for(Vec2d f: fs){ f.set(0,0); }
    }
    
    @Override
    public int interactNeighs_inds( int n0, int nng, int[] inds ){    
        int nint=0;
        //System.out.println( "interactNeighs_inds n0 "+n0+" nng "+nng );
        Vec2d d = (Vec2d)pool.borrow();
        double R2 = collisionRadius*collisionRadius;
        for(int i=0; i<nng; i++){
            int ip=inds[i];
            Vec2d pi=ps[ip];
            Vec2d fi=fs[ip];
            //Vec2d vi=vs[i]p;
            for(int j=0; j<n0; j++){
                int jp=inds[j];
                if(j<=i) continue;
                Vec2d pj=ps[jp];
                Vec2d fj=fs[jp];
                //Vec2d vj=vs[jp];
                d.set_sub( pj, pi );
                double r2  = d.norm2();
                if(r2<R2){
                    double r  = Math.sqrt(r2);
                    //System.out.println( "ps["+i+","+j+"] "+r );
                    double fr = -collisionStiffness * (r-collisionRadius)/(r+1e-3);
                    //if(bDEBUG)System.out.println( "fr."+i+"."+j+" "+fr );
                    d.mul(fr);
                    fj.add(d);
                    fi.sub(d);
                    nint++;
                }
            }
        }
        pool.payback(d);
        return nint;
    }
    
    public int interactNeighs_naive( ){    
        int nint=0;
        //System.out.println( "interactNeighs_inds n0 "+n0+" nng "+nng );
        Vec2d d = (Vec2d)pool.borrow();
        int n=ps.length;
        double R2 = collisionRadius*collisionRadius;
        for(int i=0; i<n; i++){
            Vec2d pi=ps[i];
            Vec2d fi=fs[i];
            //Vec2d vi=vs[i]p;
            for(int j=0; j<i; j++){
                if(j==i) continue;
                Vec2d pj=ps[j];
                Vec2d fj=fs[j];
                //Vec2d vj=vs[jp];
                d.set_sub( pj, pi );
                double r2  = d.norm2();
                if(r2<R2){
                    double r  = Math.sqrt(r2);
                    //System.out.println( "ps["+i+","+j+"] "+r );
                    double fr = -collisionStiffness * (r-collisionRadius)/(r+1e-3);
                    //if(bDEBUG)System.out.println( "fr."+i+"."+j+" "+fr );
                    d.mul(fr);
                    fj.add(d);
                    fi.sub(d);
                    nint++;
                }
            }
        }
        pool.payback(d);
        return nint;
    }
    
    public void project2map( Map2D map ){
        for (Vec2d p : ps) { map.insert( p.x, p.y); }
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
    
    public void move_GD(double dt){
        for(int i=0; i<ps.length; i++){
            ps[i].fma( fs[i], dt );
        }
    }
        
    public void move_leapFrog(double dt, double damping){
        double damp = 1-damping*dt;
        if(damp<0)damp=0;
        for(int i=0; i<ps.length; i++){
            Vec2d v = vs[i];
            v.mul(damp);
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
        move_leapFrog(dt, 0.0 );
    }
    
}
