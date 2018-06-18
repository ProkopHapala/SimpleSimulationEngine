
#include <stdlib.h>
#include <stdio.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include <GL/glu.h>

#include "../Demo.h"
#include "Solids.h"
#include "CMesh.h"
#include "Draw3D.h"

class Demo1:  public Demo{  public:
    int ogl;

    void setup(){
        ogl = glGenLists(1);
        glNewList(ogl, GL_COMPILE);
        glColor3f(0.0,0.0,1.0);
        Draw3D::drawLines( Solids::Tetrahedron.nedge, (int*)Solids::Tetrahedron.edges, Solids::Tetrahedron.verts );
        glColor3f(1.0,0.0,0.0);
        for( int i=0; i<Solids::Tetrahedron.nvert; i++ ){
            //Draw3D::drawPointCross( Solids::Tetrahedron.verts[i], 0.1 );
            Vec3d pav = Vec3dZero;
            for( int j=0; j<Solids::Tetrahedron.nvert; j++ ){
                if(i!=j) pav.add( Solids::Tetrahedron.verts[j] );
            }
            pav.mul( (2.0/(Solids::Tetrahedron.nvert-1)) );
            pav.sub( Solids::Tetrahedron.verts[i] );
            for( int j=0; j<Solids::Tetrahedron.nvert; j++ ){
                if(i!=j) Draw3D::drawLine(pav,Solids::Tetrahedron.verts[j]);
            }
            //Draw3D::drawPointCross( pav, 0.1 );
        };

        glEndList();
    };

    void draw(){
        glCallList(ogl);
    }

    void move( float x, float y ){
    }

    void onMouse( float x, float y, uint8_t buttons ){
    }

};

Demo* CreateDemo(){
    return (Demo*) ( new Demo1() );
};
