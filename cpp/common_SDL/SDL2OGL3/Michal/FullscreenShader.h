#include "GLES.h"
#include "GLMesh.h"
#include "GLTexture.h"
#include "Framebuffer.h"

const static char* vertexShaderSource = R"(#version 300 es

in vec4 vPosition;
out vec2 fUV;

void main(){
    gl_Position = vPosition;
    fUV = (vPosition.xy + vec2(1, 1)) / 2.0;
}

)";

class FullscreenShader{
private:
    Shader shader;
    GLMeshBase<MPOS> mesh;
    GLFramebuffer framebuffer_back; // everything rendered between "begin()" and "end()" will be rendered to this framebuffer
    
public:
    GLFramebuffer out_framebuffer; // after "end()", the result of the shader is stored in this framebuffer
    
    FullscreenShader(const char* fragShaderSource) : shader(Shader(vertexShaderSource, fragShaderSource))
    {
        printf("%s\n", fragShaderSource);
        mesh = GLMeshBase<MPOS>(GL_TRIANGLES, GL_STATIC_DRAW, &shader);
        mesh.addVertex({-1, -1, -1});
        mesh.addVertex({ 3, -1, -1});
        mesh.addVertex({-1,  3, -1});
        mesh.setUniformTex("uTexture", &framebuffer_back.colorBuffer);
        mesh.setUniformTex("uDepth"  , &framebuffer_back.depthBuffer);
    }

    void begin(){
        framebuffer_back.begin();
    }

    void end(){
        framebuffer_back.end();
        out_framebuffer.begin();
        
        glEnable(GL_DEPTH_TEST);
        mesh.draw();

        out_framebuffer.end();
    }

    inline void pause(){ framebuffer_back.pause(); }
    inline void unpause(){ framebuffer_back.unpause(); }
};
