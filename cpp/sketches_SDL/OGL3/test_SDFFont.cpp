/**
 * test_SDFFont.cpp — Demo of SDF (Signed Distance Field) font rendering.
 *
 * Shows crisp text at any zoom level using a low-res SDF texture atlas.
 * Based on Chris Green (Valve), SIGGRAPH 2007.
 *
 * Controls:
 *   Mouse right-drag: rotate camera
 *   WASD: move camera
 *   +/-: zoom in/out
 *   F: toggle between SDF and binary (non-SDF) rendering for comparison
 */

#include <stdlib.h>
#include <stdio.h>
#include <cmath>

#include <vector>
#include <GL/glew.h>

#include <fastmath.h>
#include <Vec2.h>
#include <Vec3.h>
#include <Mat3.h>
#include <quaternion.h>
#include <raytrace.h>

#include "Shader.h"
#include "GLObject.h"
#include "SceneOGL3.h"
#include "ScreenSDL2OGL3.h"
#include "AppSDL2OGL3.h"
#include "GLfunctions.h"
#include "GLobjects.h"

#include "SDFFont.h"

// ========== Test App

class TestAppSDFFont: public AppSDL2OGL3, public SceneOGL3 { public:

    SDFFontRenderer sdfFont;
    GLuint dummyTex = 0;

    bool bUseSDF = true;

    TestAppSDFFont():AppSDL2OGL3(800,600),SceneOGL3(){
        for( ScreenSDL2OGL3* screen: screens ) screen->scenes.push_back( this );

        // Generate SDF font atlas from real TTF font
        int nGlyphs    = 95;   // printable ASCII 32..126
        int perGlyphW  = 64;   // SDF texture width per glyph
        int perGlyphH  = 64;   // SDF texture height
        const char* fontPath = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf";
        dummyTex = makeSDFFontAtlas(fontPath, nGlyphs, 32, perGlyphW, perGlyphH, 8.0f);

        // Init SDF font renderer
        sdfFont.init(dummyTex, nGlyphs, 32, perGlyphW, perGlyphH, 0.2f, 0.2f);

        Camera& cam = screens[0]->cam;
        cam.zmin = 1.0; cam.zmax = 1000.0; cam.zoom = 5.0f;
        cam.aspect = screens[0]->WIDTH/(float)screens[0]->HEIGHT;  // width/height
        cam.pos.set(0.0f, 4.0f, -30.0f);  // looking at text from distance 30
    }

    virtual void eventHandling(const SDL_Event &event) {
        AppSDL2OGL3::eventHandling(event);
        if(event.type == SDL_KEYDOWN) {
            if(event.key.keysym.sym == SDLK_f) {
                bUseSDF = !bUseSDF;
                printf("SDF rendering: %s\n", bUseSDF ? "ON" : "OFF (binary)");
            }
        }
    }

    int frameCount2 = 0;

    virtual void draw(Camera& cam) {
        glClearColor(0.15f, 0.15f, 0.2f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glDisable(GL_DEPTH_TEST); // text is overlay

        // Add labels
        sdfFont.clear();
        sdfFont.addLabel({0.0f, 0.0f, 0.0f},   "HELLO SDF WORLD");
        sdfFont.addLabel({0.0f, -2.0f, 0.0f},  "CRISP AT ANY ZOOM");
        sdfFont.addLabel({0.0f, -4.0f, 0.0f},  "ABCDEFGHIJKLMNOPQRSTUVWXYZ");
        sdfFont.addLabel({0.0f, -6.0f, 0.0f},  "abcdefghijklmnopqrstuvwxyz");
        sdfFont.addLabel({0.0f, -8.0f, 0.0f},  "0123456789 !@#$%^&*()");

        frameCount2++;

        sdfFont.draw(cam);

        glEnable(GL_DEPTH_TEST);
    }
};

// ================== main

TestAppSDFFont * app;

int main(int argc, char *argv[]) {
    app = new TestAppSDFFont();
    app->loop(1000000);
    app->quit();
    return 0;
}
