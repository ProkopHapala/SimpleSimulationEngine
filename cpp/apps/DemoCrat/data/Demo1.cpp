
#include <stdlib.h>
#include <stdio.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include <GL/glu.h>

#include "../Demo.h"

class Demo1:  public Demo{  public:

    int n;
    Vec2f * ps;
    
    void setup(){ 
        printf("Initializing Chain \n"); 
        n  = 30;
        ps = new Vec2f[n];
        for(int i=0; i<n; i++ ){ ps[i]=(Vec2f){0.0,i*1.0f}; };
    };
    
    void draw(){
        glBegin(GL_LINE_STRIP);
        for(int i=0; i<n; i++ ){
            glVertex3f( ps[i].x, ps[i].y, 0.0 );
        }
        glEnd();
    }

    void move( float x, float y ){
        ps[0] =  (Vec2f){x,y};
        //printf( " (%f,%f) (%f,%f) \n", ps[0].x, ps[0].y, ps[1].x, ps[1].y );
        for(int i=1; i<n; i++ ){
            Vec2f d = ps[i] - ps[i-1];
            float renorm = 0.1/d.norm();
            ps[i] = ps[i-1] + (d*renorm);
        }
    }
    
    void draw( float x, float y ){
        //printf( " plDrawXY \n" );
        move( x, y );
        draw();
    }

};

Demo* MakeDemo(){
    return (Demo*) ( new Demo1() );
};
