
#ifndef  GLfunctions_h
#define  GLfunctions_h

#include <GL/glew.h>
//#include <SDL2/SDL.h>

bool checkFramebufferStatus(){
    // check FBO status
    GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
    switch(status){
        case GL_FRAMEBUFFER_COMPLETE:                      printf( "GL_FRAMEBUFFER_COMPLETE.\n" );                              return true;
        case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT:         printf( "[ERROR] GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT\n" );          return false;
        case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT: printf( "[ERROR] GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT\n" );  return false;
        case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER:        printf( "[ERROR] GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER\n" );         return false;
        case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER:        printf( "[ERROR] GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER\n" );         return false;
        case GL_FRAMEBUFFER_UNSUPPORTED:                   printf( "[ERROR] GL_FRAMEBUFFER_UNSUPPORTED\n" );                    return false;
        default:                                           printf( "[ERROR] Framebuffer incomplete: Unknown error.\n" );        return false;
    }
}

bool checkOpenGLError(){
    GLenum status = glGetError();
    switch(status){
        case GL_NO_ERROR:        return true;
        case GL_INVALID_ENUM:    printf( "[ERROR] GL_INVALID_ENUM\n"    );           return false;
        case GL_INVALID_VALUE:   printf( "[ERROR] GL_INVALID_VALUE\n"   );           return false;
        case GL_INVALID_FRAMEBUFFER_OPERATION: printf( "[ERROR] GL_INVALID_FRAMEBUFFER_OPERATION\n" ); return false;
        case GL_OUT_OF_MEMORY:   printf( "[ERROR] GL_OUT_OF_MEMORY\n"   );           return false;
        case GL_STACK_UNDERFLOW: printf( "[ERROR] GL_STACK_UNDERFLOW\n" );           return false;
        case GL_STACK_OVERFLOW:  printf( "[ERROR] GL_STACK_OVERFLOW\n"  );           return false;
        default:                 printf( "[ERROR] glGetError(): Unknown error.\n" ); return false;
    }
}


#define GL_DEBUG { if(!checkOpenGLError()){ printf( "@ %i %s %s \n", __LINE__, __FUNCTION__, __FILE__ ); exit(0); } }


// =========== GL

inline void newUniformBuffer( GLuint& buff, GLuint sz ){
    glGenBuffers( 1, &buff );
    glBindBuffer( GL_UNIFORM_BUFFER, buff );
    glBufferData( GL_UNIFORM_BUFFER, sz, NULL, GL_STREAM_DRAW );
    glBindBuffer( GL_UNIFORM_BUFFER, 0 );
}


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
    glGenTextures(1, &textureID);                GL_DEBUG;              // Create one OpenGL texture
    glBindTexture  (GL_TEXTURE_2D, textureID);   GL_DEBUG;
    //
    // see : https://learnopengl.com/Advanced-Lighting/HDR
    if(format==GL_RGBA && type==GL_FLOAT){
        //glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA32F, W, H, 0, GL_RGBA, GL_FLOAT, cbuff );
        glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA16F, W, H, 0, GL_RGBA, GL_FLOAT, cbuff );
    }else{
        glTexImage2D( GL_TEXTURE_2D, 0, format,     W, H, 0, format,  type,     cbuff );
    } GL_DEBUG;

    //glTexImage2D   ( GL_TEXTURE_2D, 0, GL_RGB16F, W, H, 0, GL_RGB, GL_FLOAT, cbuff );
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);  GL_DEBUG;
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    // how to properly minimap:http://www.opengl-tutorial.org/beginners-tutorials/tutorial-5-a-textured-cube/
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR); GL_DEBUG;
    //glActiveTexture(GL_TEXTURE0 );
    glGenerateMipmap(GL_TEXTURE_2D); GL_DEBUG;
}

inline void bindTexture( GLuint slot, GLuint textureID, GLuint uloc ){
    // https://www.opengl.org/discussion_boards/showthread.php/163092-Passing-Multiple-Textures-from-OpenGL-to-GLSL-shader
    // https://www.khronos.org/registry/OpenGL-Refpages/es1.1/xhtml/glActiveTexture.xml
    glActiveTexture(GL_TEXTURE0 + slot );
    glBindTexture(GL_TEXTURE_2D, textureID);
    glUniform1i( uloc, slot );
}


inline void makeRandomTexture( GLuint& textureID, int W, int H, GLenum type=GL_UNSIGNED_BYTE ){
    float* cbuff  = new float [4*W*H];
    int i=0;
    for(int iy=0; iy<H; iy++){
        for(int ix=0; ix<W; ix++){
            cbuff[i+0] = randf();
            cbuff[i+1] = randf();
            cbuff[i+2] = randf();
            cbuff[i+3] = randf();
            i+=4;
        }
    }
    newTexture2D( textureID, W, H, cbuff, GL_RGBA, type );
    delete [] cbuff;
}

// ========== Frame Buffer

/*
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
*/


#endif
