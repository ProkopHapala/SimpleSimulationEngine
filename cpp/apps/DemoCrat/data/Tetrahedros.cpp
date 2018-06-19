
#include <stdlib.h>
#include <stdio.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include <GL/glu.h>

#include "../Demo.h"
#include "Solids.h"
#include "CMesh.h"
#include "Draw.h"
#include "Draw3D.h"

void reflect( int pivot, Vec3d from[], Vec3d to[] ){
    int j=0;
    Vec3d cog; cog.set(0.0);
    for( int i=0; i<4; i++ ){
        if(i!=pivot){
            to[j]=from[i];
            cog.add(from[i]);
            //print(from[i]); printf("\n");
            j++;
        }
    }
    cog.mul(1/3.0);
    //print(cog); printf("\n");
    //Draw3D::drawPointCross(cog,0.1);
    to[3]=(cog*2)-from[pivot];
}

void recur( Vec3d from[], int level ){
    Draw3D::drawLine( from[0], from[3] );
    Draw3D::drawLine( from[1], from[3] );
    Draw3D::drawLine( from[2], from[3] );
    if(level<=0) return;
    level--;
    //Draw::color_of_hash((level+44877)*454545+1144);
    Draw::color_of_hash( level+15448878 );
    Vec3d t1[4]; reflect( 0, from, t1 ); recur( t1, level );
    Vec3d t2[4]; reflect( 1, from, t2 ); recur( t2, level );
    Vec3d t3[4]; reflect( 2, from, t3 ); recur( t3, level );
    //Vec3d t4[4]; reflect( 1, from, t4 ); //recur( t1, level );
    //Vec3d t5[4]; reflect( 1, from, t5 ); //recur( t2, level );
    //Vec3d t6[4]; reflect( 1, from, t6 ); //recur( t3, level );
    //Vec3d t4[4]; reflect( 1, from, t4 ); recur( t4, level );
}

class Demo1:  public Demo{  public:
    int ogl;

    void setup(){
        ogl = glGenLists(1);
        glNewList(ogl, GL_COMPILE);
        glColor3f(0.0,0.0,1.0);

        const CMesh& msh = Solids::Tetrahedron;

        Draw3D::drawLines( msh.nedge, (int*)msh.edges, msh.verts );
        glColor3f(1.0,0.0,0.0);

        Vec3d cog = Vec3d::average( msh.nvert, msh.verts );
        Mat3d rot; rot.fromDirUp( (msh.verts[0]-msh.verts[1]).normalized(), msh.verts[2]-msh.verts[3] );
        Draw3D::drawMatInPos(rot,cog);

        //recur( msh.verts, 3 );

        int nlevels = 1;
        Vec3d t1[4]; reflect( 0, msh.verts, t1 ); recur( t1, nlevels );
        Vec3d t2[4]; reflect( 1, msh.verts, t2 ); recur( t2, nlevels );
        Vec3d t3[4]; reflect( 2, msh.verts, t3 ); recur( t3, nlevels );
        Vec3d t4[4]; reflect( 3, msh.verts, t4 ); recur( t4, nlevels );

        /*
        for( int i=0; i<msh.nvert; i++ ){
            //Draw3D::drawPointCross( Solids::Tetrahedron.verts[i], 0.1 );
            Vec3d pav = Vec3dZero;
            for( int j=0; j<msh.nvert; j++ ){
                if(i!=j) pav.add( msh.verts[j] );
            }
            pav.mul( (2.0/(msh.nvert-1)) );
            pav.sub( msh.verts[i] );
            for( int j=0; j<msh.nvert; j++ ){
                if(i!=j) Draw3D::drawLine(pav,msh.verts[j]);
            }
            //Draw3D::drawPointCross( pav, 0.1 );
        };
        */

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
