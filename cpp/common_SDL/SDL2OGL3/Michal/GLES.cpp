#include "GLES.h"
#include <cstdlib>
#include <cstdio>

CameraT<float>* GLES::active_camera = nullptr;
GLuint GLES::currentGL_ARRAY_BUFFER = 0;
Vec2i GLES::screen_size = {1820, 980};
GLuint GLES::activeProgram = -1;
bool GLES::context = false;


static const char* GLErrorString(GLenum error) {
    switch (error) {
        case GL_NO_ERROR: return "GL_NO_ERROR";
        case GL_INVALID_ENUM: return "GL_INVALID_ENUM";
        case GL_INVALID_VALUE: return "GL_INVALID_VALUE";
        case GL_INVALID_OPERATION: return "GL_INVALID_OPERATION";
        case GL_INVALID_FRAMEBUFFER_OPERATION: return "GL_INVALID_FRAMEBUFFER_OPERATION";
        case GL_OUT_OF_MEMORY: return "GL_OUT_OF_MEMORY";
        default: return "UNKNOWN_GL_ERROR";
    }
}

void GLES::checkError(const char* file, int line){
    GLenum err = glGetError();
    if (err == GL_NO_ERROR) return;
    
    printf("\033[1m\033[31m GL error: %s at %s:%i\033[0m\n", GLErrorString(err), file, line);
}

static std::vector<GLuint> framebufferStack;
void GLES::pushFramebuffer(GLuint handle){
    framebufferStack.push_back(handle);
    glBindFramebuffer(GL_FRAMEBUFFER, handle);
}
void GLES::popFramebuffer(GLuint handle){
    if (framebufferStack.empty() || framebufferStack.back() != handle) {
        printf("\033[1m\033[31m GL error: popFramebuffer called with wrong handle\033[0m\n");
        exit(1);
    }
    framebufferStack.pop_back();
    if (framebufferStack.empty()) {
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
    } else {
        glBindFramebuffer(GL_FRAMEBUFFER, framebufferStack.back());
    }
}
