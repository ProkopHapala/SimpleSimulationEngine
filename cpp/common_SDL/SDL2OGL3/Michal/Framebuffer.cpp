#include "Framebuffer.h"
#include "GLES.h"
#include "GLTexture.h"
#include "GLMesh.h"

const static char* vertexShaderSource = R"(#version 300 es

in mediump vec4 vPosition;
out mediump vec2 fUV;

void main(){
    gl_Position = vPosition;
    fUV = (vPosition.xy + vec2(1, 1)) / 2.0;
}

)";

static const char* fragmentShaderSource = R"(#version 300 es

uniform sampler2D uTexture1;
uniform sampler2D uDepth1;

uniform sampler2D uTexture2;
uniform sampler2D uDepth2;

in mediump vec2 fUV;
layout(location=0) out mediump vec4 FragColor;

void main(){
    highp float depth1 = texture(uDepth1, fUV).x;
    highp float depth2 = texture(uDepth2, fUV).x;

    mediump vec4 color1 = texture(uTexture1, fUV);
    mediump vec4 color2 = texture(uTexture2, fUV);

    if (depth1 < depth2){
        FragColor = vec4((color2.rgb*color2.a).rgb, color2.a);
        FragColor = vec4((color1.rgb*color1.a + (1.0-color1.a)*FragColor.rgb).rgb, 1.0 - (1.0-color1.a)*(1.0-FragColor.a));
        gl_FragDepth = depth1;
    }else{
        FragColor = vec4((color1.rgb*color1.a).rgb, color1.a);
        FragColor = vec4((color2.rgb*color2.a + (1.0-color2.a)*FragColor.rgb).rgb, 1.0 - (1.0-color2.a)*(1.0-FragColor.a));
        gl_FragDepth = depth2;
    }
}

)";

static Shader shader = Shader(vertexShaderSource, fragmentShaderSource);

static inline GLMeshBase<MPOS>makeMergeMesh(){
    GLMeshBase<MPOS> m = GLMeshBase<MPOS>(GL_TRIANGLES, GL_STATIC_DRAW, &shader);
    m.addVertex({-1, -1, -1});
    m.addVertex({ 3, -1, -1});
    m.addVertex({-1,  3, -1});
    return m;
}
static GLMeshBase<MPOS> mergeMesh = makeMergeMesh();

void mergeFramebuffers(GLFramebuffer& fb1, GLFramebuffer& fb2){
    mergeMesh.setUniformTex("uTexture1", &fb1.colorBuffer);
    mergeMesh.setUniformTex("uTexture2", &fb2.colorBuffer);
    mergeMesh.setUniformTex("uDepth1"  , &fb1.depthBuffer);
    mergeMesh.setUniformTex("uDepth2"  , &fb2.depthBuffer);
    
    glEnable(GL_DEPTH_TEST);
    mergeMesh.draw();
}

