
#ifndef  MeshBuilder2Draw_h
#define  MeshBuilder2Draw_h

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw2D.h"
#include "Draw3D.h"
//#include "Draw3D_Surf.h"
//#include "SDL_utils.h"

#include "MeshBuilder2.h"

namespace Mesh {

void drawPointLabels( const Mesh::Builder2& mesh, float sz=0.02 ){
    for(int i=0; i<mesh.verts.size(); i++){
        Draw3D::drawInt( mesh.verts[i].pos, i, fontTex, sz );
    }
}

void drawEdgeLabels( const Mesh::Builder2& mesh, float sz=0.02 ){
    for(int i=0; i<mesh.edges.size(); i++){
        Vec2i e =  mesh.edges[i].lo;
        Vec3d p = (mesh.verts[e.x].pos + mesh.verts[e.y].pos)*0.5;
        Draw3D::drawInt( p, i, fontTex, sz );
    }
}

void drawSelectedEdgeLabels( const Mesh::Builder2& mesh, float sz=0.02 ){
    for(int ie: mesh.selset){
        Vec2i e =  mesh.edges[ie].lo;
        Vec3d p = (mesh.verts[e.x].pos + mesh.verts[e.y].pos)*0.5;
        Draw3D::drawInt( p, ie, fontTex, sz );
    }
}

void drawPolygonBorder( const Mesh::Builder2& mesh, int ich ){
    Quat4i ch = mesh.chunks[ich];
    glBegin(GL_LINE_LOOP);
    for(int i=0;i<ch.z;i++){
        int iv = mesh.strips[ch.x+i];
        //printf( " %i ", i, ch.x+i, iv );
        Draw3D::vertex( mesh.verts[iv].pos );
    }
    glEnd();
}


void drawFaces( const Mesh::Builder2& mesh, bool bFlat=true ){
    const int face_typ = (int)Mesh::Builder2::ChunkType::face;
    for( Quat4i ch: mesh.chunks ){
        //printf( "drawFaces() chunk(%i,%i,%i,%i)\n", ch.x, ch.y, ch.z, ch.w );
        if(ch.w==face_typ ){
            glBegin(GL_TRIANGLE_FAN);
            Vec3d nr;
            if(bFlat){
                const Vec3d& a = mesh.verts[mesh.strips[ch.x  ]].pos;
                const Vec3d& b = mesh.verts[mesh.strips[ch.x+1]].pos;
                const Vec3d& c = mesh.verts[mesh.strips[ch.x+2]].pos;
                nr = cross(b-a, c-a); nr.normalize(); 
            }
            for(int i=0;i<ch.z;i++){
                int iv = mesh.strips[ch.x+i];
                //printf( "i: %i iv: %i\n", i, iv, mesh.verts.size() );
                if(!bFlat){ nr = mesh.verts[iv].nor; }
                Draw3D::normal( nr );
                Draw3D::vertex( mesh.verts[iv].pos );
            }
            glEnd();
        }
    }
}

void drawSelectedEdges( const Mesh::Builder2& mesh ){
    glBegin(GL_LINES);
    for(int ie: mesh.selset){
        Vec2i e = mesh.edges[ie].lo;
        Draw3D::drawLine( mesh.verts[e.i].pos, mesh.verts[e.j].pos );
    }
    glEnd();
}

void drawEdges( const Mesh::Builder2& mesh ){
    glBegin(GL_LINES);
    for(int i=0;i<mesh.edges.size();i++){
        Vec2i e = mesh.edges[i].lo;
        Draw3D::drawLine( mesh.verts[e.i].pos, mesh.verts[e.j].pos );
    }
    glEnd();
}

}; // namespace Mesh

#endif
