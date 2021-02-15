package Tests;



import Common.NBody;
import Common.Vec2d;
import Common.NBody_arr;
import Common.TempPool;
import processing.core.PApplet;
import processing.core.PConstants;

public class Sketch_1 extends PApplet implements PConstants  {

    float zoomRate = 1.2f;
    float zoom     = 100.0f;
    NBody nbody1; //= new NBody_Vec2d(nbody);
    NBody_arr nbody2; //= new NBody_array(nbody);
    TempPool pool_Vec2d;
        
    @Override
    public void settings() {
        size(800, 800 );
        //size(800, 800, P2D);
        //size(800, 800, P3D);
        
    }
    
    @Override
    public void setup() {
        colorMode(RGB, 1.0f );
        
        pool_Vec2d = new TempPool( 10, Vec2d.class );
        
        int nbody    = 5000;
        nbody1 = new NBody(nbody, pool_Vec2d );
        nbody2 = new NBody_arr(nbody);
        for(int i=0; i<nbody; i++){
            double px=Math.random();
            double py=Math.random();
            nbody1.ps[i].set(px,py);
            nbody2.ps[i*2  ]=px;
            nbody2.ps[i*2+1]=py;
        }
    }

    @Override
    public void draw() {
        
        background(1.0f);
        
        int perFrame = 1;
        double dt = 0.0001;
        
        //long nop  = (nbody*nbody/2 + nbody)*perFrame; 
        
        stroke( 0, 0, 0 );
        strokeWeight(2.0f);
        
        for(int i=0; i<perFrame; i++)nbody1.update(dt,false);
        for(int i=0; i<nbody1.ps.length; i++){
            Vec2d p = nbody1.ps[i];
            point( (float)p.x*zoom+width*0.5f,(float)p.y*zoom+height*0.5f );
            //System.out.println(  i+" "+p.x+" "+p.y  );
        }
        
        /*
        for(int i=0; i<perFrame; i++)nbody2.update(dt);
        for(int i=0; i<nbody2.ps.length/2; i++){
            float px = (float)nbody2.ps[i*2  ];
            float py = (float)nbody2.ps[i*2+1];
            point( px*zoom+width*0.5f,py*zoom+height*0.5f );
            //System.out.println(  i+" "+p.x+" "+p.y  );
        }
        */
        //line(mouseX, mouseY, width / 2, height / 2);
        
    }
    
    @Override
    public void keyPressed( ){
        switch(key){
            case '+': zoom*=zoomRate; break;
            case '-': zoom/=zoomRate; break;
        }
    }
    
    public static void main (String... args) {
        Sketch_1 pt = new Sketch_1();
        PApplet.runSketch(new String[]{"Sketch_1"}, pt);
    }

}

