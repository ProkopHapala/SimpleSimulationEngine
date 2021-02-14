
package Tests;


import Common.HashMap2D;
import Common.Vec2d;
import Common.NBody_Vec2d;


import processing.core.PApplet;
import processing.core.PConstants;

public class Test_CollisionGrids extends PApplet implements PConstants  {

    float zoomRate = 1.2f;
    float zoom     = 100.0f;
    NBody_Vec2d nbody1; //= new NBody_Vec2d(nbody);
    HashMap2D   hmap;
    
    Object [] tmpOuts;
        
    @Override
    public void settings() {
        size(800, 800 );
        //size(800, 800, P2D);
        //size(800, 800, P3D);
    }
    
    @Override
    public void setup() {
        colorMode(HSB, 1.0f );
        int nbody    = 1000;
        nbody1 = new NBody_Vec2d(nbody);
        hmap   = new HashMap2D( 5 ); 
        hmap.setStep( 50.0 );
        
        double px=400;
        double py=400;
        float speed=30.0f;
        for(int i=0; i<nbody; i++){
            //double px=random(0,width);
            //double py=random(0,height);
            px+=random(-speed,+speed);
            py+=random(-speed,+speed);
            if(px<0){ px=width -1; }else if(px>width  ){ px=1; }
            if(py<0){ py=height-1; }else if(py>height ){ py=1; }
            Vec2d p = nbody1.ps[i];
            p.set(px,py);
            int im = hmap.insert( p, p.x, p.y );
            
            //int ip = hmap.hash( p.x, p.y );
            //System.out.println(  "p."+i+" -> "+ip  );
        }
        //hmap.insert( nbody1.ps[0], nbody1.ps[0].x, nbody1.ps[0].y );
        tmpOuts = new Object[100];       
        
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
        
    }

    @Override
    public void draw() {
        background(0.5f);
        
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
        int n = hmap.getAllInBox( mouseX, mouseY, tmpOuts );
        //if(n>0)System.out.println(  "nfpund " + n  );
        for(int i=0; i<n; i++){
            Vec2d p = (Vec2d)tmpOuts[i];
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

