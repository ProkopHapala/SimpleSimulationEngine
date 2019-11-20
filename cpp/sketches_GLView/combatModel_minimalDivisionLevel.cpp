
#include <cmath>
#include <cstdio>
#include "GLView.h"

#include "Draw.h"
#include "Draw2D.h"
#include "Plot2D.h"
#include "PlotScreen2D.h"
#include "testUtils.h"

#include "SDL_utils.h"
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "MinimalDivisionLevel.h"

// ===== Globals

Plot2D plot1;
int fontTex;

void my_draw(){

    //printf(" my_draw ! %i %i %g \n", plot1.lines.size(), plot1.lines.back()->n, plot1.lines.back()->ys[2]  );
    glClearColor( 1.0f, 1.0f, 1.0f, 0.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glDisable( GL_DEPTH_TEST );
    //plot1.drawAxes();

    Draw2D::drawLine({-1,0},{1,0});

    plot1.drawAxes();
    plot1.view();

}

void setup(){

fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );

    
}

int main(){

    setup();

/*
    init( 640, 480 );
    set_draw_function( my_draw );
    setup();
    run_Nframes(5000);
*/

}
