
package Tests;


import Common.HashMap2D;
import Common.Vec2d;
import Common.NBody;
import Common.CellSort2D;
import Common.TempPool;
import Utils.PDraw.PDraw;


import processing.core.PApplet;
import processing.core.PConstants;

public class Test_CollisionGrids extends PApplet implements PConstants  {

    float zoomRate = 1.2f;
    float zoom     = 100.0f;
    int perFrame=1;
    
    NBody       nbody1; //= new NBody_Vec2d(nbody);
    HashMap2D   hmap;
    CellSort2D  cmap;

    TempPool pool_Vec2d;
    PDraw pg;
    
    Object [] tmpOuts;
    int    [] outIDs;
        
    @Override
    public void settings() {
        size(800, 800 );
        //size(800, 800, P2D);
        //size(800, 800, P3D);
    }
    
    @Override
    public void setup() {
        colorMode(HSB, 1.0f );
        
        pg = new PDraw(this);
        
        pool_Vec2d = new TempPool( 10, Vec2d.class );
        Vec2d d = (Vec2d)pool_Vec2d.borrow();
        System.out.println( "pool[0]: "+d.x+" "+ d.y );
        pool_Vec2d.payback(d);
        
        int nbody    = 1000;
        nbody1 = new NBody(nbody, pool_Vec2d );
        hmap   = new HashMap2D( 5 ); 
        
        int step = 50;
        hmap.setStep( 50.0 );
        
        cmap = new CellSort2D( width, height, step, nbody );
                
        //double px=400;
        //double py=400;
        float speed=30.0f;
        for(int i=0; i<nbody; i++){
            double px=random(0,width);
            double py=random(0,height);
            //px+=random(-speed,+speed);
            //py+=random(-speed,+speed);
            //if(px<0){ px=width -1; }else if(px>width  ){ px=1; }
            //if(py<0){ py=height-1; }else if(py>height ){ py=1; }
            Vec2d p = nbody1.ps[i];
            p.set(px,py);
            int im = hmap.insert( p, p.x, p.y );
            
            cmap.insert( p.x, p.y );
            
            //int ip = hmap.hash( p.x, p.y );
            //System.out.println(  "p."+i+" -> "+ip  );
        }
        //hmap.insert( nbody1.ps[0], nbody1.ps[0].x, nbody1.ps[0].y );
        tmpOuts = new Object[100];       
        outIDs  = new int   [100];
        
        
        cmap.sort();
        
        int ntot=0;
        for(int i=0; i<hmap.capacity; i++){
            if(hmap.hits[i]>0){
                ntot+=hmap.hits[i];
                //System.out.println( "hmap.hits."+i+" : "+hmap.hits[i]  );
            }
        }
        if(ntot!=nbody1.ps.length){
            System.out.println( "ntot!=nbody1.ps.length | "+ntot+" != "+nbody1.ps.length  );
            System.exit(0); 
        }
        
        int n = hmap.getAllInBox( nbody1.ps[0].x, nbody1.ps[0].y, tmpOuts );
        
        nbody1.bDEBUG = true;
        nbody1.cleanForce();
        n = cmap.interactNeighs(nbody1, outIDs);
        //System.out.println( "n-interactions "+n );
        pg.vecsInPos( nbody1.ps.length, nbody1.ps, nbody1.fs, 10.0f );
        nbody1.bDEBUG = false;        
        
        System.out.println( " *** SETUP DONE *** " );
    }
    
    @Override
    public void draw() {
        background(0.5f);
        
        // ======== Simulation
        perFrame=1;
        for(int i=0; i<perFrame; i++){
            nbody1.cleanForce();

            cmap.clean();
            cmap.insert( nbody1.ps );
            int n = cmap.interactNeighs(nbody1, outIDs);
            //System.out.println( "n-interactions "+n );
            pg.vecsInPos( nbody1.ps.length, nbody1.ps, nbody1.fs, 10.0f );
            //vecsInPos(nbody1.ps.length, nbody1.ps, nbody1.fs, (float)10.0);

            nbody1.move_GD(0.1);
        }
    
        // ======== Drawing
        
        strokeWeight(2);
        
        for(int i=0; i<nbody1.ps.length; i++){
            Vec2d p = nbody1.ps[i];
            int ip = hmap.hash( p.x, p.y );
            stroke( (float)(ip*16.5454)%1, 1.f, 1.f  );
            point( (float)p.x,(float)p.y );
            //point( (float)p.x*zoom+width*0.5f,(float)p.y*zoom+height*0.5f );
            //System.out.println(  i+" "+p.x+" "+p.y  );
        }
        
        stroke(1);
        strokeWeight(1);
        
        //int ip = hmap.hash( mouseX, mouseY );
        //long ip = hmap.ibox2d( mouseX, mouseY );
        //int n = hmap.getAllInBox( mouseX, mouseY, tmpOuts );
        int n = cmap.loadNeighs( (float)mouseX, (float)mouseY, outIDs, true );
        //if(n>0)System.out.println(  "nfpund " + n  );
        for(int i=0; i<n; i++){
            //Vec2d p = (Vec2d)tmpOuts[i];
            Vec2d p = nbody1.ps[ outIDs[i] ];
            line( mouseX, mouseY, (float)p.x, (float)p.y );
            
        }
        
        //System.out.println(  "mouse: "+((int)mouseX)+" "+((int)mouseY)+" -> "+ip  );

    }
    
    @Override
    public void keyPressed( ){
        switch(key){
            case '+': zoom*=zoomRate; break;
            case '-': zoom/=zoomRate; break;
        }
    }
    
    public static void main (String... args) {
        Test_CollisionGrids pt = new Test_CollisionGrids();
        PApplet.runSketch(new String[]{"test_CollisionGrid"}, pt);
    }
}

