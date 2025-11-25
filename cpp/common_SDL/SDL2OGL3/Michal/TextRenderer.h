#ifndef _TEXT_RENDERER_H_
#define _TEXT_RENDERER_H_

#include <string>
#include "GLMesh.h"
#include "MeshBuilder.h"

extern GLTexture text_font;

class TextRenderer{
    GLMesh<MPOS,MUV> mesh = GLMesh<MPOS,MUV>(defaultTexColorShader<MPOS,MUV>);
    std::string text = "";

public:
    
    TextRenderer(const std::string txt = "MISSING_TEXT") {
        mesh.uniforms.setTex("uTexture", &text_font);
        set(txt);
    }

    void set(const std::string& txt, int maxWidth=0xFFFF, int maxHeight = 0xFFFF){
        if (txt == text) return;
        text = txt;

        mesh.verts->clear();

        Vec2f pos = {0,0};
        for (int i=0; i<text.size(); i++){
            if (text[i] == '\n' || pos.x >= maxWidth){
                pos.x = 0;
                pos.y -= 2;
                if (text[i] == '\n') continue;
            }
            if (pos.y <= -2*maxHeight) break;

            MeshBuilder::addChar( *mesh.verts, text[i], pos, {1, 2});
            pos.x += 1;
        }
    }
    void setf(const char* fmt, ...);

    void draw2D(Vec3f pos, float size, Vec3f color=COLOR_BLACK){
        text_font.getHandle();
        
        mesh.setUniform3f("uColor", color);
        mesh.draw2D(pos, size);
    }
    void draw3D(Vec3f pos, Vec3f color=COLOR_BLACK);
};



#endif // _TEXT_RENDERER_H_