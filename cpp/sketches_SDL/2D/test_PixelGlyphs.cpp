/// @file @brief This program demonstrates pixel-based glyph rendering for displaying text. It loads a font texture and renders characters on screen, showcasing basic text rendering capabilities and potentially exploring effects like anti-aliasing or color variations.

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "globals.h"

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"


#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw2D.h"
#include "AppSDL2OGL.h"
#include "SDL_utils.h"

#include "Plot2D.h"
#include "PlotScreen2D.h"

// ======================  TestApp


// Ref See:
// Chapter 9. Drawing Pixels, Bitmaps, Fonts, and Images -  http://www.dei.isep.ipp.pt/~matos/cg/docs/OpenGL_PG/ch09.html#id35226
// glBitmap     -  https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/glBitmap.xml
// glDrawPixels -  https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/glDrawPixels.xml

/*
How to generate pixel fonts form Xterm:
- Ctrl + RMB opens VT font menu  https://wiki.archlinux.org/index.php/Xterm
- >> xterm -geometry 100x24   // to setup window size
*/


class TestAppPixelGlyphs : public AppSDL2OGL{
	public:

    Plot2D plot1;
    //int fontTex;

    SDL_Surface * surf;

	virtual void draw   ();
	virtual void drawHUD();
    //virtual void eventHandling( const SDL_Event& event );

	TestAppPixelGlyphs( int& id, int WIDTH_, int HEIGHT_ );

	void demonstrateTextOffsetEffect( const char* caption, const char* testStr, int sz, float dx, float dy ){
        Draw::drawText( caption, fontTex, sz, 0 );
        glTranslatef(dx,dy-sz*2,0.0);
        Draw::drawText( testStr, fontTex, sz, 0 );
        glTranslatef(-dx,-dy+sz*2,0.0);
    }
};

TestAppPixelGlyphs::TestAppPixelGlyphs( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    //fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    //fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );
    fontTex = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );
    //fontTex = makeTexture( "common_resources/dejvu_sans_mono_A_pix.bmp" );

    //SDL_Surface * surf

    //surf = SDL_LoadBMP( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );
    //surf = SDL_LoadBMP( "common_resources/dejvu_sans_mono_RGBA_pix-UpDown.bmp" );
    surf = SDL_LoadBMP( "common_resources/dejvu_sans_mono_A_pix-UpDown.bmp" );
    surf->pixels;

    for(int i=0; i<95; i++){ printf("%c", i+32 ); }; printf("\n");


    /*
    for(int ix=0;ix<surf->w;ix++){
        for(int iy=0;iy<surf->h;iy++){
            int i = iy*surf->w + ix;
            printf( "%03i ", ((uint8_t*)surf->pixels)[i] );
        }
        printf("\n");
    }
    */

    //SDL_FreeSurface( surf );

     /*
    Func1d myFunc = &sin;

    plot1.init();
    plot1.fontTex = fontTex; !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~

    DataLine2D * line1 = new DataLine2D(100);
    line1->linspan(-3*M_PI,2*M_PI);
    line1->yfunc = myFunc;

    DataLine2D * line2 = new DataLine2D(400);
    line2->linspan(-M_PI,3*M_PI);
    for(int i=0; i<line2->n; i++){
        double x     = line2->xs[i];
        line2->ys[i] = sin(x*x);
    }
    line2->clr = 0xFF00FF00;

    plot1.lines.push_back( line1 );
    plot1.lines.push_back( line2 );
    plot1.render();
    */

}

void TestAppPixelGlyphs::draw(){
    glClearColor( 1.0f, 1.0f, 1.0f, 0.0f );
    //glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );


	/*
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glColor3f(1.0,0.0,0.0);
    //glDrawPixels(664,  14, GL_RGBA,  GL_UNSIGNED_BYTE,  surf->pixels );
    //glDrawPixels(664,  14, GL_ALPHA,  GL_UNSIGNED_BYTE,  surf->pixels );
    //glDrawPixels(664,  14, GL_R,  GL_UNSIGNED_BYTE,  surf->pixels );
    //glBitmap (664, 14, 0.0, 0.0, 0, 0.0, (const GLubyte*)surf->pixels );
    */


    //Draw2D::drawText("sdfhsdfhjegjfgj544615464*/-/*;\;'",(Vec2d){0.0,0.0},(Vec2d){50.0,50.0},fontTex,0.5);

	//plot1.drawAxes();
    //plot1.view();
};

void TestAppPixelGlyphs::drawHUD(){
    //printf( "drawGUI \n");
    //Draw2D::drawText( );
    const char testStr[] = " !\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";
    int sz = 7;

    const char testTextBlock[] = R"(
    class TestAppPixelGlyphs : public AppSDL2OGL{ public:
    Plot2D plot1;
    int fontTex;
    SDL_Surface * surf;
	virtual void draw   ();
	virtual void drawHUD();
    //virtual void eventHandling( const SDL_Event& event );
	TestAppPixelGlyphs( int& id, int WIDTH_, int HEIGHT_ );
        void demonstrateTextOffsetEffect( const char* caption, const char* testStr, int sz, float dx, float dy ){
            Draw::drawText( caption, fontTex, sz, 0 );
            glTranslatef(dx,dy-sz*2,0.0);
            Draw::drawText( testStr, fontTex, sz, 0 );
            glTranslatef(-dx,-dy+sz*2,0.0);
        }
    };)";


    glColor3f(0.0,0.0,0.0);
    glTranslatef(10.0,500.0,0.0);
    demonstrateTextOffsetEffect( "noOffset", testStr, sz, 0.0,0.0 ); glTranslatef(0.0,sz*-4,0.0);
    demonstrateTextOffsetEffect( "x0.5",     testStr, sz, 0.5,0.0 ); glTranslatef(0.0,sz*-4,0.0);
    demonstrateTextOffsetEffect( "y0.5",     testStr, sz, 0.0,0.5 ); glTranslatef(0.0,sz*-4,0.0);
    demonstrateTextOffsetEffect( "x0.5y0.5", testStr, sz, 0.5,0.5 ); glTranslatef(0.0,sz*-4,0.0);

    Draw::drawText( testTextBlock, fontTex,7,{30,60}  );

};

// ===================== MAIN

TestAppPixelGlyphs * testApp;

int main(int argc, char *argv[]){

	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestAppPixelGlyphs( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}
