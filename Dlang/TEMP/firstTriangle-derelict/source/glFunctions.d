module glFunctions;

import std.stdio;
import derelict.opengl3.gl3;

// =========== GL_ARRAY_BUFFER

GLuint newArrayBuffer( const GLfloat[] c_buff, GLenum usage ){
    GLuint buff;
    glGenBuffers(1, &buff);
    glBindBuffer(GL_ARRAY_BUFFER, buff );
    glBufferData(GL_ARRAY_BUFFER, c_buff.sizeof, c_buff.ptr, usage);
    //glBindBuffer(GL_ARRAY_BUFFER, vertexBufferID);
    //glBufferData(GL_ARRAY_BUFFER, vertices.length * Vec3f.sizeof, vertices.ptr, GL_STATIC_DRAW );
    return buff;
}

void uploadArrayBuffer( GLuint buff, GLuint sz, const void * c_buff ){
    glBindBuffer   (GL_ARRAY_BUFFER, buff);
    glBufferData   (GL_ARRAY_BUFFER,    sz, null, GL_STREAM_DRAW); // Buffer orphaning, a common way to improve streaming perf. See above link for details.
    glBufferSubData(GL_ARRAY_BUFFER, 0, sz, c_buff );
}

void bindVertexAttribPointer( GLuint id, GLuint buff, GLint sz, GLenum type, GLboolean normalized  ){
    glEnableVertexAttribArray( id );
    glBindBuffer(GL_ARRAY_BUFFER, buff );
    glVertexAttribPointer(id, sz, type, normalized, 0, cast(void*)0 );
}

void drawVertexArray( uint n, GLuint vbID, GLenum mode=GL_TRIANGLES, uint vertSize=3, GLboolean normalized=GL_FALSE ){
    glEnableVertexAttribArray(0);
    glBindBuffer( GL_ARRAY_BUFFER, vbID );
    glVertexAttribPointer( 0, vertSize, GL_FLOAT, GL_FALSE, 0, cast(void*)0 );
    //glUseProgram(programID);
    glDrawArrays( GL_TRIANGLES, 0, n );
}

// =========== GL_ELEMENT_ARRAY_BUFFER

GLuint newElementBuffer( const GLuint[] c_buff, GLenum usage ){
    GLuint buff;
    glGenBuffers(1, &buff);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buff );
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, c_buff.sizeof, c_buff.ptr, usage );
    return buff;
}

void drawElements( GLenum draw_mode, GLuint buff, GLuint n ){
    glBindBuffer  ( GL_ELEMENT_ARRAY_BUFFER, buff );
    glDrawElements( draw_mode,n,GL_UNSIGNED_INT, cast(void*)0 );
}

// =========== GL_ELEMENT_ARRAY_BUFFER

GLuint newTexture2D( int W, int H, const void * cbuff, GLint format, GLenum type ){
    GLuint textureID;
    // example:
    // newTexture2D( textureID, 800, 600, imgData, GL_RGBA, GL_UNSIGNED_BYTE );
    // newTexture2D( textureID, 800, 600, imgData, GL_R, GL_FLOAT );
    glGenTextures(1, &textureID);    // Create one OpenGL texture
    glBindTexture  (GL_TEXTURE_2D, textureID);
    glTexImage2D   (GL_TEXTURE_2D, 0, format, W, H, 0, format, type, cbuff );
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    // how to properly minimap:http://www.opengl-tutorial.org/beginners-tutorials/tutorial-5-a-textured-cube/
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    //glActiveTexture(GL_TEXTURE0 );
    glGenerateMipmap(GL_TEXTURE_2D);
    return textureID;
}

 void bindTexture( GLuint slot, GLuint textureID, GLuint uloc ){
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
        case GL_FRAMEBUFFER_COMPLETE:                      writeln( "Framebuffer complete." );                                              return true;
        case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT:         writeln( "[ERROR] Framebuffer incomplete: Attachment is NOT complete." );        return false;
        case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT: writeln( "[ERROR] Framebuffer incomplete: No image is attached to FBO." );       return false;
        case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER:        writeln( "[ERROR] Framebuffer incomplete: Draw buffer." );                       return false;
        case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER:        writeln( "[ERROR] Framebuffer incomplete: Read buffer." );                       return false;
        case GL_FRAMEBUFFER_UNSUPPORTED:                   writeln( "[ERROR] Framebuffer incomplete: Unsupported by FBO implementation." ); return false;
        default:                                           writeln( "[ERROR] Framebuffer incomplete: Unknown error." );                     return false;
    }
}

