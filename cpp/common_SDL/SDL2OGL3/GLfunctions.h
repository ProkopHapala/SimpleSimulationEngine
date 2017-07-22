
#ifndef  GLfunctions_h
#define  GLfunctions_h

#include <GL/glew.h>
//#include <SDL2/SDL.h>

// =========== GL_ARRAY_BUFFER

inline void newArrayBuffer( GLuint& buff, GLuint sz,  void * c_buff, GLenum usage ){
    glGenBuffers(1, &buff);
    glBindBuffer(GL_ARRAY_BUFFER, buff );
    glBufferData(GL_ARRAY_BUFFER, sz, c_buff, usage);
};

inline void uploadArrayBuffer( GLuint buff, GLuint sz,  void * c_buff ){
    glBindBuffer   (GL_ARRAY_BUFFER, buff);
    glBufferData   (GL_ARRAY_BUFFER,    sz, NULL, GL_STREAM_DRAW); // Buffer orphaning, a common way to improve streaming perf. See above link for details.
    glBufferSubData(GL_ARRAY_BUFFER, 0, sz, c_buff );
};

inline void bindVertexAttribPointer( GLuint id, GLuint buff, GLint sz, GLenum type, GLboolean normalized  ){
    glEnableVertexAttribArray( id );
    glBindBuffer(GL_ARRAY_BUFFER, buff );
    glVertexAttribPointer(id, sz, type, normalized, 0, (void*)0 );
};

// =========== GL_ELEMENT_ARRAY_BUFFER

inline void newElementBuffer( GLuint& buff, GLuint sz,  void * c_buff, GLenum usage ){
    glGenBuffers(1, &buff);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buff);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sz, c_buff, usage );
}

inline void drawElements( GLenum draw_mode, GLuint buff, GLuint n ){
    glBindBuffer  ( GL_ELEMENT_ARRAY_BUFFER, buff );
    glDrawElements( draw_mode,n,GL_UNSIGNED_INT,(void*)0 );
}

#endif
