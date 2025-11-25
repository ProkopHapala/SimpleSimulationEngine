
#include "GLES.h"
#include "GLTexture.h"
#include <cstdio>
#include <cstdlib>

class GLFramebuffer{
public:
    GLTexture colorBuffer = GLTexture(GLES::screen_size, GL_RGBA);
    GLTexture depthBuffer = GLTexture(GLES::screen_size, GL_DEPTH_COMPONENT, GL_DEPTH_COMPONENT24, GL_UNSIGNED_INT);

private:
    GLuint handle = 0;
    GLuint depthBufferHandle = 0;
    bool active = false;
    bool paused = false;

    void ensure_handle(){
        if (handle) return;

        // adapted from https://www.opengl-tutorial.org/intermediate-tutorials/tutorial-14-render-to-texture
        glGenFramebuffers(1, &handle);
        glBindFramebuffer(GL_FRAMEBUFFER, handle);

        colorBuffer.setMagFilter(GL_NEAREST);
        colorBuffer.setMinFilter(GL_NEAREST);

        depthBuffer.setMagFilter(GL_NEAREST);
        depthBuffer.setMinFilter(GL_NEAREST);

        //glGenRenderbuffers(1, &depthBufferHandle);
        //glBindRenderbuffer(GL_RENDERBUFFER, depthBufferHandle);
        //glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, GLES::screen_size.x, GLES::screen_size.y);
        //glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthBufferHandle);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT , GL_TEXTURE_2D, depthBuffer.getHandle(), 0);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, colorBuffer.getHandle(), 0);

        //Check framebuffer status
        GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
        if (status != GL_FRAMEBUFFER_COMPLETE) {
            char* status_str;
            switch(status) {
                case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT:         status_str = "GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT"; break;
                case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT: status_str = "GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT"; break;
                case GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS: status_str = "GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS"; break;
                case GL_FRAMEBUFFER_UNSUPPORTED: status_str = "GL_FRAMEBUFFER_UNSUPPORTED"; break;
                default: status_str = "Unknown error"; break;
            }

            printf("Framebuffer is not complete: %s\n", status_str);
            exit(0);
        }
    }

public:
    GLFramebuffer(){};
    ~GLFramebuffer(){ if (handle) glDeleteFramebuffers(1, &handle); }

    GLuint getHandle(){
        ensure_handle();
        return handle;
    }

    void begin(){
        if (active) {
            printf("ERROR: framebuffer is already active\n");
            exit(1);
        }
        colorBuffer.resize(GLES::screen_size);
        depthBuffer.resize(GLES::screen_size);

        ensure_handle();
        GLES::pushFramebuffer(handle);
        glViewport(0, 0, GLES::screen_size.x, GLES::screen_size.y);
        glClearColor(0, 0, 0, 0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        active = true;
    }

    void end(){
        if (paused){
            printf("ERROR: framebuffer is paused\n");
            exit(1);
        }
        if (!active){
            printf("ERROR: framebuffer is not active\n");
            exit(1);
        }
        GLES::popFramebuffer(handle);
        active = false;
    }

    void pause(){
        if (!active){
            printf("ERROR: framebuffer is not active\n");
            exit(1);
        }
        if (paused){
            printf("ERROR: framebuffer is already paused\n");
            exit(1);
        }
        GLES::popFramebuffer(handle);
        paused = true;
    }

    void unpause(){
        if (!active){
            printf("ERROR: framebuffer is not active\n");
            exit(1);
        }
        if (!paused){
            printf("ERROR: framebuffer is not paused\n");
            exit(1);
        }
        GLES::pushFramebuffer(handle);
        paused = false;
    }
};

void mergeFramebuffers(GLFramebuffer& fb1, GLFramebuffer& fb2);
