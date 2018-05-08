
#ifndef  GLInstances_h
#define  GLInstances_h

#include <GL/glew.h>
//#include <SDL2/SDL.h>
#include "GLfunctions.h"

// =====================
// ===== GLInstances
// =====================

class GLInstances{ public:
    GLuint model_vpos=0;
    GLuint model_vnor=0;
    GLuint pose_pos  =0;
    GLuint pose_dir  =0;
    GLuint pose_Up   =0;
    GLuint pose_sc   =0;
    int nInstances=0, nVerts=0;

    void init( int nInstances_, int nVerts_, const void * c_model_vpos, const void * c_model_vnor, const void * c_pose_pos, const void * c_pose_dir, const void * c_pose_Up, const void * c_pose_sc ){
        nInstances = nInstances_;
        nVerts     = nVerts_;
        newArrayBuffer( model_vpos, nVerts    *3*sizeof(GLfloat),  c_model_vpos, GL_STATIC_DRAW );
        newArrayBuffer( model_vnor, nVerts    *3*sizeof(GLfloat),  c_model_vpos, GL_STREAM_DRAW );
        newArrayBuffer( pose_pos,   nInstances*3*sizeof(GLfloat),  c_pose_pos,   GL_STREAM_DRAW );
        newArrayBuffer( pose_dir,   nInstances*3*sizeof(GLfloat),  c_pose_dir,   GL_STREAM_DRAW );
        newArrayBuffer( pose_Up,    nInstances*3*sizeof(GLfloat),  c_pose_Up,    GL_STREAM_DRAW );
        newArrayBuffer( pose_sc,    nInstances*3*sizeof(GLfloat),  c_pose_sc,    GL_STREAM_DRAW );
    };

    void draw( GLenum mode ){
        bindVertexAttribPointer( 0, model_vpos, 3, GL_FLOAT, GL_FALSE );
        bindVertexAttribPointer( 1, model_vnor, 3, GL_FLOAT, GL_FALSE );
        bindVertexAttribPointer( 2, pose_pos  , 3, GL_FLOAT, GL_FALSE );
        bindVertexAttribPointer( 3, pose_dir  , 3, GL_FLOAT, GL_FALSE );
        bindVertexAttribPointer( 4, pose_Up   , 3, GL_FLOAT, GL_FALSE );
        bindVertexAttribPointer( 5, pose_sc   , 3, GL_FLOAT, GL_FALSE );
        glVertexAttribDivisor(0, 0);
        glVertexAttribDivisor(1, 0);
        glVertexAttribDivisor(2, 1);
        glVertexAttribDivisor(3, 1);
        glVertexAttribDivisor(4, 1);
        glVertexAttribDivisor(5, 1);
        glDrawArraysInstanced( mode, 0, nVerts, nInstances);
        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(1);
        glDisableVertexAttribArray(2);
        glDisableVertexAttribArray(3);
        glDisableVertexAttribArray(4);
        glDisableVertexAttribArray(5);
    };
};

// =====================
// ===== GLInstances
// =====================

//static const float DEFAULT_Bilboard_UVs[] = {
//    0.0f,0.0f,   1.0f,0.0f,   0.0f,1.0f,
//    1.0f,1.0f,   1.0f,0.0f,   0.0f,1.0f
//};

class GLBillboards{ public:
    GLuint model_UVs =0;
    GLuint pose_pos  =0;
    GLuint pose_sc   =0;
    int nInstances=0, nVerts=0;

    void init( int nInstances_, int nVerts_, const void * c_model_UVs, const void * c_pose_pos, const void * c_pose_sc ){
        nInstances = nInstances_;
        nVerts     = nVerts_;
        newArrayBuffer( model_UVs,  nVerts    *2*sizeof(GLfloat),  c_model_UVs,  GL_STATIC_DRAW );
        newArrayBuffer( pose_pos,   nInstances*3*sizeof(GLfloat),  c_pose_pos,   GL_STREAM_DRAW );
        newArrayBuffer( pose_sc,    nInstances*2*sizeof(GLfloat),  c_pose_sc,    GL_STREAM_DRAW );
    };

    void draw( GLenum mode ){
        bindVertexAttribPointer( 0, model_UVs , 2, GL_FLOAT, GL_FALSE );
        bindVertexAttribPointer( 1, pose_pos  , 3, GL_FLOAT, GL_FALSE );
        bindVertexAttribPointer( 2, pose_sc   , 2, GL_FLOAT, GL_FALSE );
        glVertexAttribDivisor(0, 0);
        glVertexAttribDivisor(1, 1);
        glVertexAttribDivisor(2, 1);
        glDrawArraysInstanced( mode, 0, nVerts, nInstances);
        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(1);
        glDisableVertexAttribArray(2);
    };

};

// =====================
// ===== GLParticles
// =====================

class GLParticles{ public:
    GLuint buff_verts =0;
    GLuint buff_pos   =0;
    GLuint buff_color =0;
    int nInstances=0, nVerts=0;

    void init( int nInstances_, int nVerts_, void * cbuff_verts,  void * cbuff_pos, void * cbuff_color ){
        nInstances = nInstances_;
        nVerts     = nVerts_;

        glGenBuffers(1, &buff_verts);
        glBindBuffer(GL_ARRAY_BUFFER, buff_verts );
        glBufferData(GL_ARRAY_BUFFER, nVerts*3*sizeof(GLfloat),     cbuff_verts, GL_STATIC_DRAW);

        glGenBuffers(1, &buff_pos);
        glBindBuffer(GL_ARRAY_BUFFER, buff_pos   );
        glBufferData(GL_ARRAY_BUFFER, nInstances*4*sizeof(GLfloat), cbuff_pos, GL_STREAM_DRAW);

        glGenBuffers(1, &buff_color);
        glBindBuffer(GL_ARRAY_BUFFER, buff_color );
        glBufferData(GL_ARRAY_BUFFER, nInstances*4*sizeof(GLubyte), cbuff_color, GL_STREAM_DRAW);
    };

    void upload_pos( int n, void * data){
        glBindBuffer(GL_ARRAY_BUFFER, buff_pos);
        glBufferData(GL_ARRAY_BUFFER,       n * 4 * sizeof(GLfloat), NULL, GL_STREAM_DRAW); // Buffer orphaning, a common way to improve streaming perf. See above link for details.
        glBufferSubData(GL_ARRAY_BUFFER, 0, n * 4 * sizeof(GLfloat), data);
    };

    void upload_colors( int n, void * data ){
        glBindBuffer(GL_ARRAY_BUFFER, buff_color);
        glBufferData(GL_ARRAY_BUFFER,       n * 4 * sizeof(GLubyte), NULL, GL_STREAM_DRAW); // Buffer orphaning, a common way to improve streaming perf. See above link for details.
        glBufferSubData(GL_ARRAY_BUFFER, 0, n * 4 * sizeof(GLubyte), data);
    };

    void draw(){
        glEnableVertexAttribArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, buff_verts);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0,(void*)0 );

        glEnableVertexAttribArray(1);
        glBindBuffer(GL_ARRAY_BUFFER, buff_pos );
        glVertexAttribPointer(1,4,GL_FLOAT,GL_FALSE,0,(void*)0);

        glEnableVertexAttribArray(2);
        glBindBuffer(GL_ARRAY_BUFFER, buff_color);
        glVertexAttribPointer(2,4,GL_UNSIGNED_BYTE,GL_TRUE,0,(void*)0);

        glVertexAttribDivisor(0, 0);
        glVertexAttribDivisor(1, 1);
        glVertexAttribDivisor(2, 1);
        //glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, nVerts, ParticlesCount);
        //glDrawArraysInstanced(GL_TRIANGLE_FAN, 0, nVerts, ParticlesCount);
        glLineWidth(5); glDrawArraysInstanced(GL_LINE_STRIP, 0, nVerts, nInstances);
        //glDrawArraysInstanced(GL_TRIANGLES, 0, nVerts, ParticlesCount);
        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(1);
        glDisableVertexAttribArray(2);
    };

};

#endif
