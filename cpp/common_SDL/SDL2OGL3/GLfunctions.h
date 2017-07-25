
#ifndef  GLfunctions_h
#define  GLfunctions_h

#include <GL/glew.h>
//#include <SDL2/SDL.h>

// =========== GL_ARRAY_BUFFER

inline void newArrayBuffer( GLuint& buff, GLuint sz, const void * c_buff, GLenum usage ){
    glGenBuffers(1, &buff);
    glBindBuffer(GL_ARRAY_BUFFER, buff );
    glBufferData(GL_ARRAY_BUFFER, sz, c_buff, usage);
};

inline void uploadArrayBuffer( GLuint buff, GLuint sz, const void * c_buff ){
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

inline void newElementBuffer( GLuint& buff, GLuint sz, const void * c_buff, GLenum usage ){
    glGenBuffers(1, &buff);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buff);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sz, c_buff, usage );
}

inline void drawElements( GLenum draw_mode, GLuint buff, GLuint n ){
    glBindBuffer  ( GL_ELEMENT_ARRAY_BUFFER, buff );
    glDrawElements( draw_mode,n,GL_UNSIGNED_INT,(void*)0 );
}

// =========== GL_ELEMENT_ARRAY_BUFFER

inline void newTexture2D( GLuint& textureID, int W, int H, const void * cbuff, GLint format, GLenum type ){
    // example:
    // newTexture2D( textureID, 800, 600, imgData, GL_RGBA, GL_UNSIGNED_BYTE );
    // newTexture2D( textureID, 800, 600, imgData, GL_R, GL_FLOAT );
    glGenTextures(1, &textureID);    // Create one OpenGL texture
    glBindTexture  (GL_TEXTURE_2D, textureID);
    glTexImage2D   (GL_TEXTURE_2D, 0, format, W, H, 0, format, type, cbuff );
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
}

inline void bindTexture( GLuint slot, GLuint textureID, GLuint uloc ){
    // https://www.opengl.org/discussion_boards/showthread.php/163092-Passing-Multiple-Textures-from-OpenGL-to-GLSL-shader
    // https://www.khronos.org/registry/OpenGL-Refpages/es1.1/xhtml/glActiveTexture.xml
    glActiveTexture(GL_TEXTURE0 + slot );
    glBindTexture(GL_TEXTURE_2D, textureID);
    glUniform1i( uloc, slot );
}

// ========== Frame Buffer

bool checkFramebufferStatus(){
    // check FBO status
    GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
    switch(status){
        case GL_FRAMEBUFFER_COMPLETE:                      printf( "Framebuffer complete.\n" );                                              return true;
        case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT:         printf( "[ERROR] Framebuffer incomplete: Attachment is NOT complete.\n" );        return false;
        case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT: printf( "[ERROR] Framebuffer incomplete: No image is attached to FBO.\n" );       return false;
        case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER:        printf( "[ERROR] Framebuffer incomplete: Draw buffer.\n" );                       return false;
        case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER:        printf( "[ERROR] Framebuffer incomplete: Read buffer.\n" );                       return false;
        case GL_FRAMEBUFFER_UNSUPPORTED:                   printf( "[ERROR] Framebuffer incomplete: Unsupported by FBO implementation.\n" ); return false;
        default:                                           printf( "[ERROR] Framebuffer incomplete: Unknown error.\n" );                     return false;
    }
}

#endif
