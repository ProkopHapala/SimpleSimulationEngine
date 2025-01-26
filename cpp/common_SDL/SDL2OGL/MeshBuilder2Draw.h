
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

void drawFaceLabels( const Mesh::Builder2& mesh, float sz=0.02 ){
    for(int i=0; i<mesh.chunks.size(); i++){
        Quat4i ch = mesh.chunks[i];
        if( ch.w != (int)Mesh::Builder2::ChunkType::face ) continue;
        Vec3d p=Vec3dZero;
        for( int j=0; j<ch.z; j++ ){
            int iv = mesh.strips[ch.x+j];
            p.add( mesh.verts[iv].pos );
        }
        p.mul( 1./ch.z );
        Draw3D::drawInt( p, i, fontTex, sz );
    }
}

void drawTriLabels( const Mesh::Builder2& mesh, float sz=0.02 ){
    for(int i=0; i<mesh.tris.size(); i++){
        Quat4i t = mesh.tris[i];
        Vec3d p=Vec3dZero;
        for( int j=0; j<3; j++ ){
            p.add( mesh.verts[t.array[j]].pos );
        }
        p.mul( 1./3 );
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

void drawSelectedVertLabels( const Mesh::Builder2& mesh, float sz=0.02 ){
    for(int iv: mesh.selset){
        Draw3D::drawInt( mesh.verts[iv].pos, iv, fontTex, sz );
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

void drawFace( const Mesh::Builder2& mesh, Quat4i ch, bool bFlat=true, bool bNormal=false ){
    Vec3d nr;
    //Vec3d pcog=Vec3dZero;
    if(bFlat){
        const Vec3d& a = mesh.verts[mesh.strips[ch.x  ]].pos;
        const Vec3d& b = mesh.verts[mesh.strips[ch.x+1]].pos;
        const Vec3d& c = mesh.verts[mesh.strips[ch.x+2]].pos;
        nr = cross(b-a, c-a); nr.normalize(); 
    }
    glBegin(GL_TRIANGLE_FAN);
    for(int i=0;i<ch.z;i++){
        int iv = mesh.strips[ch.x+i];
        //printf( "i: %i iv: %i\n", i, iv, mesh.verts.size() );
        if(!bFlat){ nr = mesh.verts[iv].nor; }
        Draw3D::normal( nr );
        Draw3D::vertex( mesh.verts[iv].pos );
        //pcog.add( mesh.verts[iv].pos );
    }
    glEnd();
    // if(bNormal){ 
    //     glColor3f(0.f,1.f,0.f);
    //     pcog.mul( 1./ch.z );
    //     Draw3D::drawLine( pcog, pcog+nr );
    // }
}

void drawFaces( const Mesh::Builder2& mesh, bool bFlat=true, bool bNormals=false ){
    const int face_typ = (int)Mesh::Builder2::ChunkType::face;
    for( Quat4i ch: mesh.chunks ){
        //printf( "drawFaces() chunk(%i,%i,%i,%i)\n", ch.x, ch.y, ch.z, ch.w );
        if(ch.w==face_typ ){
            drawFace( mesh, ch, bFlat );
        }
    }
}

void drawFaceNormals( const Mesh::Builder2& mesh){
    const int face_typ = (int)Mesh::Builder2::ChunkType::face;
    glBegin(GL_LINES);
    for( Quat4i ch: mesh.chunks ){
        //printf( "drawFaces() chunk(%i,%i,%i,%i)\n", ch.x, ch.y, ch.z, ch.w );
        Vec3d nr;
        Vec3d pcog=Vec3dZero;
        if(ch.w!=face_typ ) continue;
        const Vec3d& a = mesh.verts[mesh.strips[ch.x  ]].pos;
        const Vec3d& b = mesh.verts[mesh.strips[ch.x+1]].pos;
        const Vec3d& c = mesh.verts[mesh.strips[ch.x+2]].pos;
        nr = cross(b-a, c-a); nr.normalize(); 
        for(int i=0;i<ch.z;i++){
            int iv = mesh.strips[ch.x+i];
            pcog.add( mesh.verts[iv].pos );
        }
        //glColor3f(0.f,1.f,0.f);
        pcog.mul( 1./ch.z );
        //Draw3D::drawLine( pcog, pcog+nr );
        Draw3D::vertex( pcog );
        Draw3D::vertex( pcog+nr );
    }
    glEnd();
}

void drawVertLoop( const Mesh::Builder2& mesh, int n, int* ivs ){
    for(int i=0;i<n;i++){
        Draw3D::vertex( mesh.verts[ivs[i]].pos );
    }
}

void drawVertLoop( const Mesh::Builder2& mesh, int n, int* ivs, int drawAs=GL_LINE_LOOP ){
    if(drawAs!=0){ glBegin(drawAs); }
    for(int i=0;i<n;i++){
        Draw3D::vertex( mesh.verts[ivs[i]].pos );
    }
    if(drawAs!=0){ glEnd(); }
}

void drawSelectedEdges( const Mesh::Builder2& mesh ){
    glBegin(GL_LINES);
    for(int ie: mesh.selset){
        Vec2i e = mesh.edges[ie].lo;
        Draw3D::drawLine( mesh.verts[e.i].pos, mesh.verts[e.j].pos );
    }
    glEnd();
}

void drawSelectedVerts( const Mesh::Builder2& mesh ){
    glBegin(GL_POINTS);
    for(int iv: mesh.selection){
        Draw3D::vertex( mesh.verts[iv].pos );
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
