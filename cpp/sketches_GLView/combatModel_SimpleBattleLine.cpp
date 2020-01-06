
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
#include "SimpleBattleLine.h"

// ===== Globals

Plot2D plot1;
int fontTex;

//std::unordered_map<std::string,WeaponType*> weaponTypes;
//std::unordered_map<std::string,UnitType*> unitTypes;

WeaponType kar98k;
WeaponType m67;
WeaponType rpg17;
WeaponType mg42;
WeaponType gun37;

UnitType inf1;
UnitType infATG;
UnitType infATR;
UnitType tank1;

UnitState unit_inf;
UnitState unit_tank;
UnitState unit_ATG;
UnitState unit_ATR;

Combat combat1;

void my_draw(){

    //printf(" my_draw ! %i %i %g \n", plot1.lines.size(), plot1.lines.back()->n, plot1.lines.back()->ys[2]  );
    glClearColor( 1.0f, 1.0f, 1.0f, 0.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glDisable( GL_DEPTH_TEST );
    //plot1.drawAxes();

    //Draw2D::drawLine({-1,0},{1,0});
    plot1.drawAxes();
    plot1.view();

    glColor3f(0.,0.,0.);
    Vec2d mouse_pos;
    getMousePos(&mouse_pos.x,&mouse_pos.y);
    char str[256];
    //sprintf(str,"mpos(%g %g)", mouse_pos.x, mouse_pos.y );
    //sprintf(str,"%3.3f %3.3f", mouse_pos.x, mouse_pos.y );
    sprintf(str,"%3.3f %3.3f", pow(10,mouse_pos.x), mouse_pos.y );
    Draw2D::drawText(str, mouse_pos, {10,10}, fontTex, 0.1);

}

void setup(){

    //printf("\033c");
    printf( "======================\n");
    printf( "======================\n");
    printf( "======================\n");
    printf( "======================\n");
    printf( "setup() \n");

    fontTex = makeTexture( "../common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    plot1.init();
    plot1.fontTex = fontTex;

    char stmp[8000];

    //               reach(dir, indir),  supression, movePenalty, supply, firepower[]
    kar98k.fromString(" 400.   700. .1 .25 1.     2. 0. 0.  0. "); // rifle
    m67.fromString   ("  10.    30. .1 .25 5.    10. 1. 0.1 0."); // granade
    rpg17.fromString ("  50.   250. .1 .25 10.   10. 3. .5 0. ");
    mg42.fromString  (" 500.  1000. .1 .25 100.  60. 0. 0. 0. "); // machine gun
    gun37.fromString ("1500.  5000. .2 .25 200.  10. 2. .1 0. "); // AT cannon light
    

    //kar98k.info(stmp); printf( "kar98k : \n %s\n", stmp );
    //mg42.  info(stmp); printf( "mg42   : \n %s\n", stmp );
    //gun37. info(stmp); printf( "gun37  : \n %s\n", stmp );

    //              armorClass, size, weight,    supply, baseCost,    speed_transit, speed_combat   
    inf1.fromString ("0   1.  0.1   1. 1000. 15.  2.");
    inf1.primary   = &kar98k;
    inf1.secondary = &m67;
    //inf1. info(stmp,true); printf( "inf1  : \n %s\n", stmp );

    infATG.fromString ("0   1.  0.1   1. 1000. 15.  2.");
    infATG.primary   = &gun37;
    infATG.secondary = &kar98k;

    infATR.fromString ("0   1.  0.1   1. 1000. 15.  2.");
    infATR.primary   = &rpg17;
    infATR.secondary = &kar98k;

    tank1.fromString("2  10. 11.0   4. 20000. 10. 6.");
    tank1.primary   = &gun37;
    tank1.secondary = &mg42;
    //tank1.info(stmp,true); printf( "tank1 : \n %s\n", stmp );

    unit_inf .init(&inf1  ,10000);
    unit_ATG .init(&infATG,1000 );
    unit_ATR .init(&infATR,1000 );
    unit_tank.init(&tank1 ,100  );

    combat1.dist=40.0;
    combat1.conds.maxCamo=1.0;

    unit_tank.type->info(stmp,true); printf( "attacker  unit_tank  : \n %s\n", stmp );
    combat1.attacker.composition.units.push_back(&unit_tank);

    unit_inf.type->info(stmp,true); printf( "defender unit_inf  : \n %s\n", stmp );
    combat1.defender.composition.units.push_back(&unit_inf);
    //combat1.defender.composition.units.push_back(&unit_ATG);
    unit_ATR.type->info(stmp,true); printf( "defender unit_ATR  : \n %s\n", stmp );
    combat1.defender.composition.units.push_back(&unit_ATR);

    combat1.attacker.composition.reset();
    combat1.defender.composition.reset();
    combat1.start( 20.0 );
    iDEBUG = 1;
    combat1.round(0,0);


    int nsamp = 20;
    DataLine2D * l_vis1  = plot1.add( new DataLine2D(nsamp,0,1.0,      0xFF0000FF, "visibility tank" )  );
    DataLine2D * l_vis2  = plot1.add( new DataLine2D(nsamp,l_vis1->xs, 0xFFFF0000, "visibility inf"  )  );

    iDEBUG = 0;
    double zoom = 100;
    combat1.conds.viewRange = 500;
    for(int i=0; i<nsamp; i++){
        double dist = 1<<i;
        l_vis1 ->xs[i]=log10(dist);
        l_vis1 ->ys[i]=unit_tank.visibilityFunction(dist,zoom,combat1.conds);
        l_vis2 ->ys[i]=unit_inf .visibilityFunction(dist,zoom,combat1.conds);
        printf( "dist %g vis %g | %g %g \n", dist, l_vis1 ->ys[i], l_vis1 ->xs[i], l_vis2 ->ys[i] );
    }
    plot1.logX = true;
    plot1.logY = false;

    plot1.clrGrid   = 0xFFE0E0E0;
    plot1.clrTicksX = 0xFFA0A0A0;
    plot1.clrTicksY = 0xFFA0A0A0;
    plot1.xlabel="log dist[m]";
    plot1.ylabel="log got_shot";
    plot1.render(true);

    //exit(0);
    
}

int main(){

    //init( 640, 480 );
    init( 1024, 768 );
    set_draw_function( my_draw );
    setup();
    run_Nframes(500000);

}
