module glObjects;

import std.stdio;
import std.conv;

//import derelict.sdl2.sdl;
import derelict.opengl3.gl3;
import glFunctions;

import Vector;

//public import derelict.util.exception;

final class GLObject{
    union{
        //struct{ GLuint vpos=0,vnor=0,vcol=0,vUVs=0; }; // union does not work with initializer
        struct{ GLuint vpos,vnor,vcol,vUVs; };
        GLuint [4] buffs;
    };
    GLuint inds=0;
    int nVerts =0;
    int nInds  =0;
    GLenum draw_mode=GL_TRIANGLES;

    void init( int nVerts_, int nInds_, const uint[] c_inds, const float[] c_vpos, const float[] c_vnor, const float[] c_vcol, const float[] c_vUVs, GLenum usage=GL_STATIC_DRAW ){
        nVerts = nVerts_;
        nInds  = nInds_;
        if(nInds  ){ inds = newElementBuffer( c_inds, usage ); }else{ inds=0; }
        if(c_vpos ){ vpos = newArrayBuffer( c_vpos, usage ); }else{ vpos=0; }
        if(c_vnor ){ vnor = newArrayBuffer( c_vnor, usage ); }else{ vnor=0; }
        if(c_vcol ){ vcol = newArrayBuffer( c_vcol, usage ); }else{ vcol=0; }
        if(c_vUVs ){ vUVs = newArrayBuffer( c_vUVs, usage ); }else{ vUVs=0; }
    };

    void init_d( int nVerts_, int nInds_, const uint[] c_inds, const double[] c_vpos, const double[] c_vnor, const double[] c_vcol, const double[] c_vUVs, GLenum usage=GL_STATIC_DRAW ){
        nVerts = nVerts_;
        nInds  = nInds_;
        //float * fs = new float[4*nVerts];
        //float [] fs = new float[4*nVerts];
        if(nInds  ){ inds = newElementBuffer( to!(GLuint[])(c_inds), usage ); }else{ inds=0; }
        if(c_vpos ){ vpos = newArrayBuffer( to!(GLfloat[])(c_vpos), usage ); }else{ vpos=0; }
        if(c_vnor ){ vnor = newArrayBuffer( to!(GLfloat[])(c_vnor), usage ); }else{ vnor=0; }
        if(c_vcol ){ vcol = newArrayBuffer( to!(GLfloat[])(c_vcol), usage ); }else{ vcol=0; }
        if(c_vUVs ){ vUVs = newArrayBuffer( to!(GLfloat[])(c_vUVs), usage ); }else{ vUVs=0; }
        //delete[] fs;
    };


/*
    void init_wireframe( const CMesh& msh ){
        draw_mode = GL_LINES;
        init_d( msh.nvert, msh.nedge*2, (int*)msh.edges, (double*)msh.verts, NULL, NULL, NULL );
    }

    void init_hardface( const CMesh& msh ){
        draw_mode = GL_TRIANGLES;
        int nVerts = countVerts( msh.nfaces, msh.ngons );
        Vec3f * model_vpos = new Vec3f[nVerts];
        Vec3f * model_vnor = new Vec3f[nVerts];
        hardFace( msh.nfaces, msh.ngons, msh.faces, msh.verts, (GLfloat*)model_vpos, (GLfloat*)model_vnor );
        init( nVerts, 0, nullptr, (float*)model_vpos, (float*)model_vnor, nullptr, nullptr );
        delete [] model_vpos;
        delete [] model_vnor;
    }
*/

    void deleteBuffs(){ for(int i=0; i<4; i++){ if( buffs[i] ){ glDeleteBuffers(1, &buffs[i] ); } } }

    int preDraw ()const{
        int iarg = 0;
        if(vpos){ bindVertexAttribPointer( iarg, vpos, 3, GL_FLOAT, GL_FALSE );  iarg++; }
        if(vnor){ bindVertexAttribPointer( iarg, vnor, 3, GL_FLOAT, GL_FALSE );  iarg++; }
        if(vcol){ bindVertexAttribPointer( iarg, vcol, 3, GL_FLOAT, GL_FALSE );  iarg++; }
        if(vUVs){ bindVertexAttribPointer( iarg, vUVs, 2, GL_FLOAT, GL_FALSE );  iarg++; }
        return iarg;
    };
    void drawRaw( GLenum draw_mode )const{
        if(nInds){ drawElements( draw_mode, inds, nInds ); }else{ glDrawArrays( draw_mode, 0, nVerts); }
    }
    void postDraw(int narg)const{ for( int i=0; i<narg; i++ ){ glDisableVertexAttribArray(i); } }
    void draw( GLenum draw_mode )const{ int narg = preDraw(); drawRaw( draw_mode ); postDraw(narg); }
    void drawRaw()const{ drawRaw( draw_mode ); }
    void draw   ()const{ draw   ( draw_mode ); }
    void drawPointsRaw( float sz )const{ glPointSize(sz); glDrawArrays( GL_POINTS, 0, nVerts); } // mostly for debugging - ignores indexes
    void drawPoints   ( float sz )const{ int narg = preDraw(); drawPointsRaw(sz); postDraw(narg);}
};

