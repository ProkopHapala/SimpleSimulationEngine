
#include <cmath>
#include <cstdio>
#include "GLView.h"

#include <vector>
#include <unordered_map>

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

//std::unordered_map<std::string,WeaponType*> weaponTypes;
//std::unordered_map<std::string,UnitType*> unitTypes;

WeaponType kar98k;
WeaponType mg42;
WeaponType gun37;

UnitType inf1;
UnitType tank1;



void my_draw(){

    //printf(" my_draw ! %i %i %g \n", plot1.lines.size(), plot1.lines.back()->n, plot1.lines.back()->ys[2]  );
    glClearColor( 1.0f, 1.0f, 1.0f, 0.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glDisable( GL_DEPTH_TEST );
    //plot1.drawAxes();

    //Draw2D::drawLine({-1,0},{1,0});
    //plot1.drawAxes();
    //plot1.view();

}

void setup(){

    //fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );

    char stmp[8000];

    //               reach(dir, indir),  supression, movePenalty, supply, firepower[]
    kar98k.fromString("400.    700. .1 .25 1.     2. 0. 0. 0. ");
    mg42.fromString  ("500.   1000. .1 .25 100.  60. 0. 0. 0. ");
    gun37.fromString ("1500.  5000. .2 .25 200.  10. 2. .1 0. ");

    kar98k.info(stmp); printf( "kar98k : \n %s\n", stmp );
    mg42.  info(stmp); printf( "mg42   : \n %s\n", stmp );
    gun37. info(stmp); printf( "gun37  : \n %s\n", stmp );

    //              armorClass,    size,  weight, supply, baseCost,    speed_transit, speed_combat   
    inf1.fromString ("0   1.  0.1   1.   0.  1000. 15.  2.");
    inf1.primary   = &kar98k;
    inf1.secondary = 0;
    tank1.fromString("2  10. 11.0   4. 100. 20000. 10. 6.");
    tank1.primary   = &gun37;
    tank1.secondary = &mg42;
    inf1. info(stmp,true); printf( "inf1  : \n %s\n", stmp );
    tank1.info(stmp,true); printf( "tank1 : \n %s\n", stmp );

    exit(0);
    
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
