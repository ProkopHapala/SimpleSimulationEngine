
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
    //printf( "drawEdgeLabels()\n" );
    for(int i=0; i<mesh.edges.size(); i++){
        Vec2i e =  mesh.edges[i].lo;
        Vec3d p = (mesh.verts[e.x].pos + mesh.verts[e.y].pos)*0.5;
        Draw3D::drawInt( p, i, fontTex, sz );
    }
}

void drawTriagles( const Mesh::Builder2& mesh ){
    glBegin(GL_TRIANGLES);
    for(int i=0; i<mesh.tris.size(); i++){
        //Vec3d p=Vec3dZero;
        const Quat4i& t = mesh.tris[i];
        Draw3D::vertex( mesh.verts[t.x].pos );
        Draw3D::vertex( mesh.verts[t.y].pos );
        Draw3D::vertex( mesh.verts[t.z].pos );
        // p.mul( 1./3 );
        // Draw3D::drawInt( p, i, fontTex, sz );
    }
    glEnd();
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
    for(int ie: mesh.curSelection->vec){
        Vec2i e =  mesh.edges[ie].lo;
        Vec3d p = (mesh.verts[e.x].pos + mesh.verts[e.y].pos)*0.5;
        Draw3D::drawInt( p, ie, fontTex, sz );
    }
}

void drawSelectedVertLabels( const Mesh::Builder2& mesh, float sz=0.02, bool bOrder=false ){
    std::vector<int>& selection = mesh.curSelection->vec;
    if( bOrder ){ for(int i=0; i<selection.size(); i++){ Draw3D::drawInt( mesh.verts[selection[i] ].pos, i,  fontTex, sz ); } }
    else        { for(int iv: selection               ){ Draw3D::drawInt( mesh.verts[iv           ].pos, iv, fontTex, sz ); } }
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


static void drawVerts(const Mesh::Builder2& mesh) {
    glBegin(GL_POINTS);
    for (const auto& v : mesh.verts) {
        glVertex3dv(v.pos.array);
    }
    glEnd();
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
    for(int ie: mesh.curSelection->vec){
        Vec2i e = mesh.edges[ie].lo;
        Draw3D::drawLine( mesh.verts[e.x].pos, mesh.verts[e.y].pos );
    }
    glEnd();
}

void drawSelectedVerts( const Mesh::Builder2& mesh ){
    glBegin(GL_POINTS);
    for(int iv: mesh.curSelection->vec){
        //printf( "drawSelectedVerts() iv=%i \n", iv );
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

class Renderer{
    public:
    Mesh::Builder2* mesh=0;
    bool bViewVerts        = false;
    bool bViewEdges        = true;
    bool bViewTris         = true;
    bool bViewFaces        = true;
    bool bViewFaceNormals  = false;
    bool bViewPointLabels  = false;
    bool bViewEdgeLabels   = false;
    bool bViewFaceLabels   = false;
    bool bViewTriLabels    = false;

    Renderer(Mesh::Builder2& mesh): mesh(&mesh){}

    void draw(){
        glColor3f( 1.0,1.0,1.0 );
        // Enable polygon offset to push the solid faces back slightly in the depth buffer.
        // This prevents "z-fighting" where the wireframe edges would be hidden by the faces.
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset( 1.0, 1.0 );
        if(bViewFaces ){ drawFaces( *mesh );    }
        if(bViewTris  ){ drawTriagles( *mesh ); }
        glDisable(GL_POLYGON_OFFSET_FILL);
        if(bViewFaceNormals) {
            glColor3f( 0.0,0.5,1.0 );
            drawFaceNormals( *mesh );
        }
        glDisable(GL_LIGHTING);
        glDisable(GL_CULL_FACE);
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
        if(bViewEdges) {
            // Draw All Edges
            glColor3f(0.0,0.0,0.0);
            glLineWidth(1.0);
            drawEdges( *mesh );
        }
    }
    
    void drawLabels(){
        if(bViewPointLabels) {
            glColor3f(0.f,0.f,0.f);
            drawPointLabels( *mesh, 0.02 );
        }
        if(bViewEdgeLabels) {
            glColor3f(0.f,.5f,0.f);
            drawEdgeLabels( *mesh, 0.02 );
        }
        if(bViewFaceLabels) {
            glColor3f(1.f,0.f,0.f);
            drawFaceLabels( *mesh, 0.02 );
        }
        if(bViewTriLabels) {
            glColor3f(0.7f,0.5f,0.f);
            drawTriLabels( *mesh, 0.02 );
        }
    }

    void drawSelection( int ipick=-1 ){
        if( mesh->selection_mode == (int)Mesh::Builder2::SelectionMode::face ){
        if(ipick>=0){
                glLineWidth(5.0);
                glColor3f(0.0,0.7,0.0);
                drawPolygonBorder( *mesh, ipick );
            }
        }else if( mesh->selection_mode==(int)Mesh::Builder2::SelectionMode::edge ){
            glColor3f(0.0,0.7,0.0);
            glLineWidth(5.0);
            drawSelectedEdges( *mesh );
            drawSelectedEdgeLabels( *mesh, 0.02 );
            glLineWidth(1.0);
        }else if( mesh->selection_mode==(int)Mesh::Builder2::SelectionMode::vert ){
            glColor3f(0.0,1.0,0.0);
            glPointSize(8.0);
            drawSelectedVerts( *mesh );
            glColor3f(0.0,0.5,0.0);
            drawSelectedVertLabels( *mesh, 0.02, true );
            glPointSize(1.0);
        }
    }

}; // class MeshRenderer

}; // namespace Mesh

#endif
