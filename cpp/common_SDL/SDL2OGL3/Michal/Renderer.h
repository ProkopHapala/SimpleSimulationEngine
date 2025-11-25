#ifndef _Renderer_H_
#define _Renderer_H_

#include "GLES.h"
#include <math.h>
#include <cstdlib>
#include <stdint.h>
#include <vector>

#include "Vec3.h"
#include "Mat4.h"

#include <functional>
#include "GLMesh.h"

#ifndef __SKIP_DEFS__
#include "SDL_opengl_defs.h" // type and constants for old OpenGL
#endif


class OpenGL1Renderer {
    public:
        Vec3f color = {1, 1, 1};
        float color_alpha = 1;

        Vec3f normal = {0, 0, 0};

        bool begun = false;
        GLenum mode = -1;

        GLMesh<MPOS,MNORMAL,MCOLOR> current_mesh = GLMesh<MPOS,MNORMAL,MCOLOR>(GL_TRIANGLES, GL_STREAM_DRAW, defaultShader<MPOS,MNORMAL,MCOLOR>);

        GLenum MatrixMode = GL_MODELVIEW;


        std::vector<Mat4f> mvMatStack;
        std::vector<Mat4f> projMatStack;
        std::vector<Mat4f> textureMatStack;
        std::vector<Mat4f> colorMatStack;

        std::vector<std::vector<std::function<void()>>> callLists;

        std::vector<std::function<void(void)>> current_call_list_builder;
        int current_callList = 0;
        GLenum current_callListMode = GL_COMPILE;

    public:
        OpenGL1Renderer();

        // ===== OpenGL 1.x functions =====
        void begin(GLenum mode);
        void end();
        void enable(GLenum cap);
        void disable(GLenum cap);
        void flush();
        void finish();

        void color3f(float r, float g, float b);
        void color3d(double r, double g, double b);
        void color4f(float r, float g, float b, float a);

        void vertex2d(double x, double y);
        void vertex3f(float x, float y, float z);
        void vertex3d(double x, double y, double z);
        
        void normal3f(float x, float y, float z);
        void normal3d(double x, double y, double z);

        void newList(GLuint list, GLenum mode);
        GLuint genLists(GLsizei range);
        void deleteLists(GLuint list, GLsizei range);
        void endList();
        void callList(GLuint list);
        
        void clear(GLbitfield mask);
        void clearColor(float r, float g, float b, float a);

        void pushMatrix();
        void popMatrix();
        void loadMatrixf(const GLfloat *m);
        void matrixMode(GLenum mode);
        void multMatrixf(const GLfloat *m);

        void depthFunc(GLenum func);
        void blendFunc(GLenum sfactor, GLenum dfactor);

        void translatef(float x, float y, float z);
        void readPixels(GLint x, GLint y, GLsizei width, GLsizei height, GLenum format, GLenum type, GLvoid *pixels);
        void genTextures(GLsizei n, GLuint *textures);
        void shadeModel(GLenum mode);
        void viewport(GLint x, GLint y, GLsizei width, GLsizei height);
        void bindTexture(GLenum target, GLuint texture);
        void texParameteri(GLenum target, GLenum pname, GLint param);
        void lineWidth(GLfloat width);
        void pixelStorei(GLenum pname, GLint param);
        void loadIdentity();
        void frustum(GLdouble left, GLdouble right, GLdouble bottom, GLdouble top, GLdouble zNear, GLdouble zFar);
        void ortho(GLdouble left, GLdouble right, GLdouble bottom, GLdouble top, GLdouble zNear, GLdouble zFar);
        void getFloatv(GLuint pname, GLfloat *params);
        void rotatef(GLfloat angle, GLfloat x, GLfloat y, GLfloat z);
        void scalef(GLfloat x, GLfloat y, GLfloat z);
        void pointSize(GLfloat size);
        void scissor(GLint x, GLint y, GLsizei width, GLsizei height);
        void frontFace(GLenum mode);
        void materialfv(GLenum face, GLenum pname, const GLfloat *params);
        void lightModeli(GLenum pname, GLint param);
        void polygonMode(GLenum face, GLenum mode);
        void lightfv(GLenum light, GLenum pname, const GLfloat *params);
        void texImage2D(GLenum target, GLint level, GLint internalformat, GLsizei width, GLsizei height, GLint border, GLenum format, GLenum type, const GLvoid *pixels);
        void texCoord2f(GLfloat s, GLfloat t);
        void hint(GLenum target, GLenum mode);
        void lightModelfv(GLenum pname, const GLfloat *params);

        void TEST();

        const GLubyte* getString(GLenum name);

        void copyTexImage2D(GLenum target, GLint level,
            GLenum internalformat,
            GLint x, GLint y,
            GLsizei width, GLsizei height,
            GLint border);

};

extern OpenGL1Renderer opengl1renderer;




#endif // _Renderer_H_
