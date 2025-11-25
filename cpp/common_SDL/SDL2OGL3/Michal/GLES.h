
//#define DEBUG_GLES

#ifndef _GLES_H_
#define _GLES_H_

#include <GLES3/gl3.h>
#include "Vec2.h"

template<typename T> class CameraT;

namespace GLES{
    extern CameraT<float>* active_camera;
    extern GLuint currentGL_ARRAY_BUFFER;
    extern Vec2i screen_size;
    extern GLuint activeProgram;
    extern bool context;

    void checkError(const char* file, int line);
    void pushFramebuffer(GLuint handle);
    void popFramebuffer(GLuint handle);

    inline void useProgram(GLuint program){
        if (program != activeProgram) {glUseProgram(program); activeProgram = program;}
    }
}

#define GL_CHECK_ERROR() GLES::checkError(__FILE__, __LINE__)

#define glBindBuffer(target, buffer) do { \
    glBindBuffer(target, buffer); \
    if (target == GL_ARRAY_BUFFER) GLES::currentGL_ARRAY_BUFFER = buffer; \
} while (0)


#define glEnable(cap) do { \
    glEnable(cap); \
    GL_CHECK_ERROR(); \
} while(0)



#ifdef DEBUG_GLES
#include <cstdlib>
#include <execinfo.h>
#include <cstdio>

#define glActiveTexture(texture) do { \
    glActiveTexture(texture); \
    GL_CHECK_ERROR(); \
} while(0)

#define glAttachShader(program, shader) do { \
    glAttachShader(program, shader); \
    GL_CHECK_ERROR(); \
} while(0)

#define glBindAttribLocation(program, index, name) do { \
    glBindAttribLocation(program, index, name); \
    GL_CHECK_ERROR(); \
} while(0)

#define glBindFramebuffer(target, framebuffer) do { \
    glBindFramebuffer(target, framebuffer); \
    GL_CHECK_ERROR(); \
} while(0)

#define glBindRenderbuffer(target, renderbuffer) do { \
    glBindRenderbuffer(target, renderbuffer); \
    GL_CHECK_ERROR(); \
} while(0)

#define glBindTexture(target, texture) do { \
    glBindTexture(target, texture); \
    GL_CHECK_ERROR(); \
} while(0)

#define glBlendColor(red, green, blue, alpha) do { \
    glBlendColor(red, green, blue, alpha); \
    GL_CHECK_ERROR(); \
} while(0)

#define glBlendEquation(mode) do { \
    glBlendEquation(mode); \
    GL_CHECK_ERROR(); \
} while(0)

#define glBlendEquationSeparate(modeRGB, modeAlpha) do { \
    glBlendEquationSeparate(modeRGB, modeAlpha); \
    GL_CHECK_ERROR(); \
} while(0)

#define glBlendFunc(sfactor, dfactor) do { \
    glBlendFunc(sfactor, dfactor); \
    GL_CHECK_ERROR(); \
} while(0)

#define glBlendFuncSeparate(srcRGB, dstRGB, srcAlpha, dstAlpha) do { \
    glBlendFuncSeparate(srcRGB, dstRGB, srcAlpha, dstAlpha); \
    GL_CHECK_ERROR(); \
} while(0)

#define glBufferData(target, size, data, usage) do { \
    glBufferData(target, size, data, usage); \
    GL_CHECK_ERROR(); \
} while(0)

#define glBufferSubData(target, offset, size, data) do { \
    glBufferSubData(target, offset, size, data); \
    GL_CHECK_ERROR(); \
} while(0)

#define glCheckFramebufferStatus(target) ({ \
    GLenum status = glCheckFramebufferStatus(target); \
    GL_CHECK_ERROR(); \
    status; \
})

#define glClear(mask) do { \
    glClear(mask); \
    GL_CHECK_ERROR(); \
} while(0)

#define glClearColor(red, green, blue, alpha) do { \
    glClearColor(red, green, blue, alpha); \
    GL_CHECK_ERROR(); \
} while(0)

#define glClearDepthf(depth) do { \
    glClearDepthf(depth); \
    GL_CHECK_ERROR(); \
} while(0)

#define glClearStencil(s) do { \
    glClearStencil(s); \
    GL_CHECK_ERROR(); \
} while(0)

#define glColorMask(red, green, blue, alpha) do { \
    glColorMask(red, green, blue, alpha); \
    GL_CHECK_ERROR(); \
} while(0)

#define glCompileShader(shader) do { \
    glCompileShader(shader); \
    GL_CHECK_ERROR(); \
} while(0)

#define glCompressedTexImage2D(target, level, internalformat, width, height, border, imageSize, data) do { \
    glCompressedTexImage2D(target, level, internalformat, width, height, border, imageSize, data); \
    GL_CHECK_ERROR(); \
} while(0)

#define glCompressedTexSubImage2D(target, level, xoffset, yoffset, width, height, format, imageSize, data) do { \
    glCompressedTexSubImage2D(target, level, xoffset, yoffset, width, height, format, imageSize, data); \
    GL_CHECK_ERROR(); \
} while(0)

#define glCopyTexImage2D(target, level, internalformat, x, y, width, height, border) do { \
    glCopyTexImage2D(target, level, internalformat, x, y, width, height, border); \
    GL_CHECK_ERROR(); \
} while(0)

#define glCopyTexSubImage2D(target, level, xoffset, yoffset, x, y, width, height) do { \
    glCopyTexSubImage2D(target, level, xoffset, yoffset, x, y, width, height); \
    GL_CHECK_ERROR(); \
} while(0)

#define glCreateProgram() ({ \
    GLuint program = glCreateProgram(); \
    GL_CHECK_ERROR(); \
    program; \
})

#define glCreateShader(type) ({ \
    GLuint shader = glCreateShader(type); \
    GL_CHECK_ERROR(); \
    shader; \
})

#define glCullFace(mode) do { \
    glCullFace(mode); \
    GL_CHECK_ERROR(); \
} while(0)

#define glDeleteBuffers(n, buffers) do { \
    glDeleteBuffers(n, buffers); \
    GL_CHECK_ERROR(); \
} while(0)

#define glDeleteFramebuffers(n, framebuffers) do { \
    glDeleteFramebuffers(n, framebuffers); \
    GL_CHECK_ERROR(); \
} while(0)

#define glDeleteProgram(program) do { \
    glDeleteProgram(program); \
    GL_CHECK_ERROR(); \
} while(0)

#define glDeleteRenderbuffers(n, renderbuffers) do { \
    glDeleteRenderbuffers(n, renderbuffers); \
    GL_CHECK_ERROR(); \
} while(0)

#define glDeleteShader(shader) do { \
    glDeleteShader(shader); \
    GL_CHECK_ERROR(); \
} while(0)

#define glDeleteTextures(n, textures) do { \
    glDeleteTextures(n, textures); \
    GL_CHECK_ERROR(); \
} while(0)

#define glDepthFunc(func) do { \
    glDepthFunc(func); \
    GL_CHECK_ERROR(); \
} while(0)

#define glDepthMask(flag) do { \
    glDepthMask(flag); \
    GL_CHECK_ERROR(); \
} while(0)

#define glDepthRangef(zNear, zFar) do { \
    glDepthRangef(zNear, zFar); \
    GL_CHECK_ERROR(); \
} while(0)

#define glDetachShader(program, shader) do { \
    glDetachShader(program, shader); \
    GL_CHECK_ERROR(); \
} while(0)

//#define glDisable(cap) do { \
//    glDisable(cap); \
//    GL_CHECK_ERROR(); \
//} while(0)

#define glDisableVertexAttribArray(index) do { \
    glDisableVertexAttribArray(index); \
    GL_CHECK_ERROR(); \
} while(0)

#define glDrawArrays(mode, first, count) do { \
    glDrawArrays(mode, first, count); \
    GL_CHECK_ERROR(); \
} while(0)

#define glDrawElements(mode, count, type, indices) do { \
    glDrawElements(mode, count, type, indices); \
    GL_CHECK_ERROR(); \
} while(0)

//#define glEnable(cap) do { \
//    glEnable(cap); \
//    GL_CHECK_ERROR(); \
//} while(0)

#define glEnableVertexAttribArray(index) do { \
    glEnableVertexAttribArray(index); \
    GL_CHECK_ERROR(); \
} while(0)

#define glFinish() do { \
    glFinish(); \
    GL_CHECK_ERROR(); \
} while(0)

#define glFlush() do { \
    glFlush(); \
    GL_CHECK_ERROR(); \
} while(0)

#define glFramebufferRenderbuffer(target, attachment, renderbuffertarget, renderbuffer) do { \
    glFramebufferRenderbuffer(target, attachment, renderbuffertarget, renderbuffer); \
    GL_CHECK_ERROR(); \
} while(0)

#define glFramebufferTexture2D(target, attachment, textarget, texture, level) do { \
    glFramebufferTexture2D(target, attachment, textarget, texture, level); \
    GL_CHECK_ERROR(); \
} while(0)

#define glFrontFace(mode) do { \
    glFrontFace(mode); \
    GL_CHECK_ERROR(); \
} while(0)

#define glGenBuffers(n, buffers) do { \
    glGenBuffers(n, buffers); \
    GL_CHECK_ERROR(); \
} while(0)

#define glGenerateMipmap(target) do { \
    glGenerateMipmap(target); \
    GL_CHECK_ERROR(); \
} while(0)

#define glGenFramebuffers(n, framebuffers) do { \
    glGenFramebuffers(n, framebuffers); \
    GL_CHECK_ERROR(); \
} while(0)

#define glGenRenderbuffers(n, renderbuffers) do { \
    glGenRenderbuffers(n, renderbuffers); \
    GL_CHECK_ERROR(); \
} while(0)

#define glGenTextures(n, textures) do { \
    glGenTextures(n, textures); \
    GL_CHECK_ERROR(); \
} while(0)

#define glGetActiveAttrib(program, index, bufSize, length, size, type, name) do { \
    glGetActiveAttrib(program, index, bufSize, length, size, type, name); \
    GL_CHECK_ERROR(); \
} while(0)

#define glGetActiveUniform(program, index, bufSize, length, size, type, name) do { \
    glGetActiveUniform(program, index, bufSize, length, size, type, name); \
    GL_CHECK_ERROR(); \
} while(0)

#define glGetAttachedShaders(program, maxCount, count, shaders) do { \
    glGetAttachedShaders(program, maxCount, count, shaders); \
    GL_CHECK_ERROR(); \
} while(0)

#define glGetAttribLocation(program, name) ({ \
    GLint location = glGetAttribLocation(program, name); \
    GL_CHECK_ERROR(); \
    location; \
})

#define glGetBooleanv(pname, params) do { \
    glGetBooleanv(pname, params); \
    GL_CHECK_ERROR(); \
} while(0)

#define glGetBufferParameteriv(target, pname, params) do { \
    glGetBufferParameteriv(target, pname, params); \
    GL_CHECK_ERROR(); \
} while(0)

#define glGetFloatv(pname, params) do { \
    glGetFloatv(pname, params); \
    GL_CHECK_ERROR(); \
} while(0)

#define glGetFramebufferAttachmentParameteriv(target, attachment, pname, params) do { \
    glGetFramebufferAttachmentParameteriv(target, attachment, pname, params); \
    GL_CHECK_ERROR(); \
} while(0)

#define glGetIntegerv(pname, params) do { \
    glGetIntegerv(pname, params); \
    GL_CHECK_ERROR(); \
} while(0)

#define glGetProgramiv(program, pname, params) do { \
    glGetProgramiv(program, pname, params); \
    GL_CHECK_ERROR(); \
} while(0)

#define glGetProgramInfoLog(program, bufSize, length, infoLog) do { \
    glGetProgramInfoLog(program, bufSize, length, infoLog); \
    GL_CHECK_ERROR(); \
} while(0)

#define glGetRenderbufferParameteriv(target, pname, params) do { \
    glGetRenderbufferParameteriv(target, pname, params); \
    GL_CHECK_ERROR(); \
} while(0)

#define glGetShaderiv(shader, pname, params) do { \
    glGetShaderiv(shader, pname, params); \
    GL_CHECK_ERROR(); \
} while(0)

#define glGetShaderInfoLog(shader, bufSize, length, infoLog) do { \
    glGetShaderInfoLog(shader, bufSize, length, infoLog); \
    GL_CHECK_ERROR(); \
} while(0)

#define glGetShaderPrecisionFormat(shadertype, precisiontype, range, precision) do { \
    glGetShaderPrecisionFormat(shadertype, precisiontype, range, precision); \
    GL_CHECK_ERROR(); \
} while(0)

#define glGetShaderSource(shader, bufSize, length, source) do { \
    glGetShaderSource(shader, bufSize, length, source); \
    GL_CHECK_ERROR(); \
} while(0)

#define glGetString(name) ({ \
    const GLubyte* result = glGetString(name); \
    GL_CHECK_ERROR(); \
    result; \
})

#define glVertexAttribPointer(index, size, type, normalized, stride, pointer) do { \
    glVertexAttribPointer(index, size, type, normalized, stride, pointer); \
    GL_CHECK_ERROR(); \
} while(0)

#define glUniformMatrix4fv(location, count, transpose, value) do { \
    glUniformMatrix4fv(location, count, transpose, value); \
    GL_CHECK_ERROR(); \
} while(0)

#define glGetUniformLocation(program, name) ({ \
    GLint location = glGetUniformLocation(program, name); \
    GL_CHECK_ERROR(); \
    location; \
})

#define glUniform3f(location, v0, v1, v2) do { \
    glUniform3f(location, v0, v1, v2); \
    GL_CHECK_ERROR(); \
} while(0)

#endif
#endif
