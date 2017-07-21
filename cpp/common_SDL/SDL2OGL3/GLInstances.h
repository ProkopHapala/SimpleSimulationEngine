
#ifndef  GLInstances_h
#define  GLInstances_h

#include <GL/glew.h>
//#include <SDL2/SDL.h>

class GLInstances{ public:
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
