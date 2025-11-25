#include "TextRenderer.h"

static GLTexture makeTextFont(){
    GLTexture tex = GLTexture("common_resources/dejvu_sans_mono_RGBA_pix-UpDown.bmp");
    tex.setMagFilter(GL_NEAREST);
    tex.setMinFilter(GL_NEAREST);
    return tex;
}

GLTexture text_font = makeTextFont();
