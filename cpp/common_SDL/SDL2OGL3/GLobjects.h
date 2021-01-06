
#ifndef  GLobjects_h
#define  GLobjects_h

#include <GL/glew.h>
#include "GLfunctions.h"
//#include <SDL2/SDL.h>

#include "CMesh.h"
#include "GLUtils.h"

//float * double2float( int n, double * ds ){ float*fs=float[n]; for(int i=0; i<n; i++){ fs[i]=(float)ds[i]; }; return fs; }
void double2float( int n, const double * ds, float * fs ){ for(int i=0; i<n; i++){ fs[i]=(float)ds[i]; }; }

// ==============================
// ========== GLMesh
// ==============================

static const float DEFAULT_Bilboard_verts[] = {
    0.0f,0.0f,0.0f,   1.0f,0.0f,0.0f,   0.0f,1.0f,0.0f,
    1.0f,1.0f,0.0f,   1.0f,0.0f,0.0f,   0.0f,1.0f,0.0f
};

static const float DEFAULT_Bilboard_UVs[] = {
    0.0f,0.0f,   1.0f,0.0f,   0.0f,1.0f,
    1.0f,1.0f,   1.0f,0.0f,   0.0f,1.0f
};

static const float DEFAULT_Bilboard_UVs_2x2[] = {
    0.0f,0.0f,   2.0f,0.0f,   0.0f,2.0f,
    2.0f,2.0f,   2.0f,0.0f,   0.0f,2.0f
};

class GLMesh{ public:
    union{
        //struct{ GLuint vpos=0,vnor=0,vcol=0,vUVs=0; }; // union does not work with initializer
        struct{ GLuint vpos,vnor,vcol,vUVs; };
        GLuint buffs[4];
    };
    GLuint inds=0;
    int nVerts =0;
    int nInds  =0;
    GLenum draw_mode=GL_TRIANGLES;

    void init( int nVerts_, int nInds_, const void * c_inds, const void * c_vpos,  const void * c_vnor, const void * c_vcol, const void * c_vUVs, GLenum usage=GL_STATIC_DRAW ){
        nVerts = nVerts_;
        nInds  = nInds_;
        if(nInds  ){ newElementBuffer( inds, nInds   *sizeof(GLuint),  c_inds, usage ); }else{ inds=0; };
        if(c_vpos ){ newArrayBuffer  ( vpos, nVerts*3*sizeof(GLfloat), c_vpos, usage ); }else{ vpos=0; };
        if(c_vnor ){ newArrayBuffer  ( vnor, nVerts*3*sizeof(GLfloat), c_vnor, usage ); }else{ vnor=0; };
        if(c_vcol ){ newArrayBuffer  ( vcol, nVerts*3*sizeof(GLfloat), c_vcol, usage ); }else{ vcol=0; };
        if(c_vUVs ){ newArrayBuffer  ( vUVs, nVerts*2*sizeof(GLfloat), c_vUVs, usage ); }else{ vUVs=0; };
    };

    void init_d( int nVerts_, int nInds_, const int * c_inds, const double * c_vpos, const double * c_vnor, const double * c_vcol, const double * c_vUVs, GLenum usage=GL_STATIC_DRAW ){
        nVerts = nVerts_;
        nInds  = nInds_;
        float * fs = new float[4*nVerts];
        if(nInds  ){ newElementBuffer( inds, nInds   *sizeof(GLuint),  c_inds, usage ); }else{ inds=0; };
        if(c_vpos ){ double2float( nVerts*3, c_vpos,fs); newArrayBuffer( vpos, nVerts*3*sizeof(GLfloat), fs, usage ); }else{ vpos=0; };
        if(c_vnor ){ double2float( nVerts*3, c_vnor,fs); newArrayBuffer( vnor, nVerts*3*sizeof(GLfloat), fs, usage ); }else{ vnor=0; };
        if(c_vcol ){ double2float( nVerts*3, c_vcol,fs); newArrayBuffer( vcol, nVerts*3*sizeof(GLfloat), fs, usage ); }else{ vcol=0; };
        if(c_vUVs ){ double2float( nVerts*2, c_vUVs,fs); newArrayBuffer( vUVs, nVerts*2*sizeof(GLfloat), fs, usage ); }else{ vUVs=0; };
        delete[] fs;
    };

    GLMesh() = default;
    GLMesh(int nVerts_, int nInds_, const void * c_inds, const void * c_vpos,  const void * c_vnor, const void * c_vcol, const void * c_vUVs, GLenum usage=GL_STATIC_DRAW){
        init(nVerts_,nInds_,c_inds,c_vpos,c_vnor,c_vcol,c_vUVs,usage );
    };
    GLMesh( int nVerts_, int nInds_, const int * c_inds, const double * c_vpos, const double * c_vnor, const double * c_vcol, const double * c_vUVs, GLenum usage=GL_STATIC_DRAW ){
        init_d(nVerts_,nInds_,c_inds,c_vpos,c_vnor,c_vcol,c_vUVs, usage );
    };

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

    void deleteBuffs(){ for(int i=0; i<4; i++){ if( buffs[i] ){ glDeleteBuffers(1, &buffs[i] ); } }; }

    int preDraw ()const{
        int iarg = 0;
        if(inds){ glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, inds ); };
        if(vpos){ bindVertexAttribPointer( iarg, vpos, 3, GL_FLOAT, GL_FALSE );  iarg++; }
        if(vnor){ bindVertexAttribPointer( iarg, vnor, 3, GL_FLOAT, GL_FALSE );  iarg++; }
        if(vcol){ bindVertexAttribPointer( iarg, vcol, 3, GL_FLOAT, GL_FALSE );  iarg++; }
        if(vUVs){ bindVertexAttribPointer( iarg, vUVs, 2, GL_FLOAT, GL_FALSE );  iarg++; }
        return iarg;
    };
    void drawRaw( GLenum draw_mode )const{
        //if(nInds){ drawElements( draw_mode, inds, nInds ); }else{ glDrawArrays( draw_mode, 0, nVerts); };
        if(nInds){ glDrawElements( draw_mode,nInds,GL_UNSIGNED_INT,(void*)0 ); }else{ glDrawArrays( draw_mode, 0, nVerts); };
    }
    void drawRaw( GLenum draw_mode, int i0, int n )const{
        // https://stackoverflow.com/questions/9431923/using-an-offset-with-vbos-in-opengl
        //if(nInds){ drawElements( draw_mode, inds, nInds ); }else{ glDrawArrays( draw_mode, 0, nVerts); };
        if(nInds){ glDrawElements( draw_mode,n,GL_UNSIGNED_INT, (void*)(i0 * sizeof(GLuint)) ); }else{ glDrawArrays( draw_mode, i0, n); };
    }
    void postDraw(int narg)const{ for( int i=0; i<narg; i++ ){ glDisableVertexAttribArray(i); } };
    void draw( GLenum draw_mode )const{ int narg = preDraw(); drawRaw( draw_mode ); postDraw(narg); };
    void drawRaw()const{ drawRaw( draw_mode ); };
    void draw   ()const{ draw   ( draw_mode ); };
    void drawPointsRaw( float sz )const{ glPointSize(sz); glDrawArrays( GL_POINTS, 0, nVerts); } // mostly for debugging - ignores indexes
    void drawPoints   ( float sz )const{ int narg = preDraw(); drawPointsRaw(sz); postDraw(narg);}

    //inline GLMesh(){ vpos=0,vnor=0,vcol=0,vUVs=0; };
};


void makeBilboard( GLMesh*& mesh ){
    //mesh = new GLMesh();
    //mesh->init( 6, 0,  NULL, DEFAULT_Bilboard_verts, NULL, NULL, DEFAULT_Bilboard_UVs);
    mesh = new GLMesh( 6, 0,  NULL, DEFAULT_Bilboard_verts, NULL, NULL, DEFAULT_Bilboard_UVs);
}


// ==============================
// ========== GLMesh
// ==============================

class FrameBuffer{ public:
    GLuint buff =0;
    GLuint buffZ=0, buffRGB=0;
    GLuint texZ =0, texRGB=0;
    int W=0,H=0;

    inline void init( GLuint texRGB_, GLuint texZ_, int W_, int H_, bool bFloat=false ){

        texZ   =  texZ_;
        texRGB =  texRGB_;
        W=W_;H=H_;
        glGenFramebuffers(1, &buff);              GL_DEBUG;
        glBindFramebuffer(GL_FRAMEBUFFER, buff);  GL_DEBUG;

        if( bFloat ){
        // See: https://community.khronos.org/t/full-precision-renderbuffer/63716/5
        // The color buffer
        glGenRenderbuffers       (1, &buffRGB);                           GL_DEBUG;
        glBindRenderbuffer       (GL_RENDERBUFFER, buffRGB );             GL_DEBUG;
        glRenderbufferStorage    (GL_RENDERBUFFER, GL_RGBA16F, W, H );    GL_DEBUG;
        glFramebufferRenderbuffer(GL_FRAMEBUFFER,  GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, buffRGB ); GL_DEBUG;
        }

        // The depth buffer
        glGenRenderbuffers       (1, &buffZ);                                                     GL_DEBUG;
        glBindRenderbuffer       (GL_RENDERBUFFER, buffZ );                                       GL_DEBUG;
        glRenderbufferStorage    (GL_RENDERBUFFER, GL_DEPTH_COMPONENT, W, H );                    GL_DEBUG;
        glFramebufferRenderbuffer(GL_FRAMEBUFFER,  GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, buffZ ); GL_DEBUG;

        //printf( "texZ %i texRGB %i \n", texZ, texRGB );
        glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,  texZ,   0 ); GL_DEBUG;
        glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, texRGB, 0 ); GL_DEBUG;
        // Depth buffer should not be used with glDrawBuffers !!!!
        //https://stackoverflow.com/questions/25587359/how-is-gldrawbuffers-associated-to-drawing-to-a-depth-texture
        //GLenum DrawBuffers[2] = {GL_DEPTH_ATTACHMENT, GL_COLOR_ATTACHMENT0};
        //glDrawBuffers(2, DrawBuffers); GL_DEBUG;
        GLenum DrawBuffers[1] = {GL_COLOR_ATTACHMENT0};
        glDrawBuffers(1, DrawBuffers); GL_DEBUG;
        if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE){
            printf(" problem in FBO ! \n ");
            checkFramebufferStatus();
        }
        glBindFramebuffer(GL_FRAMEBUFFER, 0); GL_DEBUG;
    }

    void init( int W, int H, bool bFloat=false ){
        //newTexture2D( texRGB, W, H, NULL, GL_RGB,             GL_UNSIGNED_BYTE );
        //newTexture2D( texRGB, W, H, NULL, GL_RGBA,             GL_UNSIGNED_BYTE );
        if(bFloat){ newTexture2D( texRGB, W, H, NULL, GL_RGBA, GL_FLOAT         ); }
        else      { newTexture2D( texRGB, W, H, NULL, GL_RGBA, GL_UNSIGNED_BYTE ); }
        newTexture2D( texZ  , W, H, NULL, GL_DEPTH_COMPONENT,  GL_FLOAT  );
        init( texRGB, texZ, W, H );
    }

    void bind(){
        glBindFramebuffer(GL_FRAMEBUFFER,buff);
        glViewport(0,0,W,H);
        glScissor (0,0,W,H);
    }

    FrameBuffer() = default;
    FrameBuffer(int W, int H, bool bFloat=false){  init( W, H, bFloat ); };

};

#endif
