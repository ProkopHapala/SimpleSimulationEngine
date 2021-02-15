/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Utils.PDraw;

import Common.Vec2d;

import processing.core.PApplet;
import processing.core.PConstants;

public class PDraw{
    
    PApplet app;
    
    public PDraw(PApplet app_){
        app=app_;
    }
    
    public void vecsInPos(int n, Vec2d[] ps, Vec2d[] vs, float sc){
        for(int i=0; i<n; i++){
            Vec2d p = ps[i];
            Vec2d v = vs[i];
            app.line( (float)p.x, (float)p.y, (float)(p.x+v.x*sc),  (float)(p.y+v.y*sc) );
        }
    }
    
}
