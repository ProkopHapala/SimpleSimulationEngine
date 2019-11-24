
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

}

void setup(){

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

    iDEBUG = 1;
    combat1.dist=10.0;
    combat1.round();

    combat1.dist=1000.0;
    combat1.round();

    int nsamp = 15;
    DataLine2D * l_tank = plot1.add( new DataLine2D(nsamp,0,1.0,      0xFF0000FF ) );
    DataLine2D * l_inf  = plot1.add( new DataLine2D(nsamp,l_tank->xs, 0xFFFF0000 ) );
    DataLine2D * l_ATR  = plot1.add( new DataLine2D(nsamp,l_tank->xs, 0xFFFF8000 ) );

    iDEBUG = 0;
    for(int i=0; i<nsamp; i++){
        //double dist = 0.1;
        combat1.dist= 10*i*i;
        //l_tank->xs[i] = log10(combat1.dist);
        l_tank->xs[i] = i;
        //unit_inf .debug_got_shot=0;
        //unit_tank.debug_got_shot=0;
        //unit_ATR.debug_got_shot=0;
        //printf( "dist %g ", combat1.dist );
        combat1.round();
        l_tank->ys[i]=log10(unit_tank.debug_got_shot);
        l_inf ->ys[i]=log10(unit_inf.debug_got_shot);
        l_ATR ->ys[i]=log10(unit_ATR.debug_got_shot);
        printf( "dist %g fp12 %g %g | fp21 %g \n", combat1.dist, unit_inf.debug_got_shot, unit_ATR.debug_got_shot, unit_tank.debug_got_shot );
        //printf( "dist %g fp12 %g %g | fp21 %g \n", dist, l_inf ->ys[i]), log10(l_ATR ->ys[i]), log10(l_tank->ys[i]) );
    }

    plot1.render();

    //exit(0);
    
}

int main(){

    setup();

    init( 640, 480 );
    set_draw_function( my_draw );
    setup();
    run_Nframes(5000);

}
