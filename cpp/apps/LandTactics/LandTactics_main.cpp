
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include <SDL2/SDL_image.h>
//#include <SDL2/SDL_ttf.h>
//#include "Texture.h"
#include "Draw.h"
#include "Draw2D.h"
#include "Draw3D.h"


#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"
#include "AppSDL2OGL.h"

#include "testUtils.h"
#include "SDL_utils.h"

#include "TerrainCubic.h"
#include "TiledView.h"

#include "LTcommon.h"

#include "LTUnitType.h"
#include "LTUnit.h"
#include "LTShelter.h"
#include "LTFaction.h"
#include "LTWorld.h"

#include "LTrender.h"

#include "GUI.h"
#include "Plot2D.h"

#include "AirCombatModel.h"

// font rendering:
//  http://www.willusher.io/sdl2%20tutorials/2013/12/18/lesson-6-true-type-fonts-with-sdl_ttf
//  http://stackoverflow.com/questions/28880562/rendering-text-with-sdl2-and-opengl

/*
//   Units should have 2 scales of resolution
// - Company scale (individual units are abstracted out- not rendered, not evaluated)
// - Detailed scale ( individual units are rendered and evaluated )
*/

int   default_font_texture;
char strBuf[0x10000];

#include "LTdraw.h"

class FormationTacticsApp : public AppSDL2OGL {
	public:
    LTWorld world;

    int formation_view_mode = 0;
    LTSquad   * currentSquad      = NULL;
    LTFaction * currentFaction   = NULL;
    int       ifaction = 0;

    bool bDrawing  = false;
    bool bDrawGoal = true;

    GLuint oglTerrain=0;

    //double xsc,ysc;

    //GLuint       itex;

    // ==== function declaration
    void printASCItable( int imin, int imax  );

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling( const SDL_Event& event  );
	void debug_buffinsert( );

	FormationTacticsApp( int& id, int WIDTH_, int HEIGHT_ );

	void renderStaticTerrain();

};

void FormationTacticsApp::printASCItable( int imin, int imax  ){
    int len = imax-imin;
    char str[len];
    for ( int i=0; i<len; i++ ){
        str[i] = (char)(i+imin);
    }
    printf("%s\n", str );
};

void FormationTacticsApp::renderStaticTerrain(){

    glDisable    ( GL_LIGHTING   );
    glDisable    ( GL_DEPTH_TEST );
    glShadeModel ( GL_SMOOTH     );

    const uint32_t colorsBW[2]{ 0xFF000000, 0xFFFFFFF};
    //drawMap( SimplexRuler& ruler, world.ground,   int ncol=ncolors, const uint32_t * colors=&colors_rainbow[0] );
    //drawMap( world.ruler, world.ground, 0, 1000 );
    drawMap( world.ruler, world.ground, 0, 1000, 2, colorsBW );
    //return;

    // TODO !!!! This function Draw2D::drawTriaglePatch<>() is somehow lost !!!!
	//glPushMatrix();
	//glScalef(world.ruler.step,world.ruler.step,1.0);
	//Draw2D::drawTriaglePatch<cmapHeight>( {0,0}, {128,128}, world.ruler.na, world.ground, 0.0, world.maxHeight );
	//Draw2D::drawTriaglePatch<cmapHeight>( {0,0}, {128,128}, world.ruler.na, world.pathFinder.moveCosts, 0.0, world.maxHeight );
	//Draw2D::drawTriaglePatchBas( {0,0}, {128,128}, world.ruler.na, world.pathFinder.toBasin, 0.0, world.maxHeight );
	//glPopMatrix();

	glColor3f(1.0,1.0,1.0);

	for( Vec2i& ip : world.pathFinder.centers ){
        Vec2d p;
        world.ruler.nodePoint( ip, p );
        Draw2D::drawPointCross_d(p, 50.0);
        Draw2D::drawCircle_d(p, 50.0, 8, false );
	}

	glDisable(GL_DEPTH_TEST);
    for( auto& item : world.pathFinder.pass ){
        //PathFinder
        Vec2d p1,p2;
        Vec2i ip1 = world.pathFinder.i2ip( item.second.x);
        Vec2i ip2 = world.pathFinder.i2ip( item.second.y);
        world.ruler.nodePoint( ip1, p1   );
        world.ruler.nodePoint( ip2, p2   );
        //printf( " point %i %i %g %g\n", ip1.x, ip2.x, p.x, p.y );
        //Draw2D::z_layer = 1.0;
        //Draw2D::drawPointCross_d( p, 10.0 );
        Draw2D::drawLine_d( p1, p2 );
    }

    for(Way* w: world.pathFinder.paths ){
        drawPath( world.ruler, *w );
    }

    for( River* river: world.hydraulics.rivers ){
        Draw::color_of_hash( river->path[0]+54877 );
        drawRiver( world.ruler, *river );
    }

	//printf( "world.objects.size %i \n", world.objects.size() );
	glColor3f(0.7,0.7,0.7);
	for( LTStaticObject& o : world.objects ){
        //Draw2D::drawShape( o.pos, o.dir, o.type->glo );
        o.view();
        //o.type->render( o.pos, o.dir );
    }

    glColor3f(0.7,0.7,0.7);
    for( LTLinearObject& o : world.linObjects ){
        //printf( " (%f,%f) (%f,%f) \n", o.p1.x, o.p1.y, o.p2.x, o.p2.y );
        Draw2D::drawLine_d( o.p1, o.p2 );
    }
    //exit(0)
}


void evalAeroCraft( const CombatAirCraft& aero1 ){
    double CL, speed, accel, Rturn, sa, vy, ta;
    speed = aero1.maxSpeed_simple ();                                   printf( " speed : %g[km/h]  \n", speed*3.6 );
    double omega = aero1.trunRate_simple ( CL, speed, accel, Rturn  );  printf( " turn  : Time %g[s] R %g[m] v %g[km/h] accel %g[g] \n", 2*M_PI/omega, Rturn, speed*3.6, accel/const_GravAccel );
    vy    = aero1.climbRate_simple( sa );
    ta    = sa/sqrt(1-sa*sa);
    printf( " climb    :  %g[m/s] @ %g [km/h] tan %g alpha %g[deg] \n", vy, (vy/ta)*3.6, ta, asin(sa)*180/M_PI );
    vy    = aero1.climbRate_CLmax ( sa );
    ta    = sa/sqrt(1-sa*sa);
    printf( " climb_CL :  %g[m/s] @ %g [km/h] tan %g alpha %g[deg] \n", vy, (vy/ta)*3.6, ta, asin(sa)*180/M_PI );
}

FormationTacticsApp::FormationTacticsApp( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    CombatAirCraft aero1;

    printf( "### I-153 \n" );
    // https://en.wikipedia.org/wiki/Polikarpov_I-153#Specifications_(I-153_(M-62))
    aero1.mass      = 1960;      // [kg]
    aero1.power     = 0.597e+6;  // [W]   engine power
    aero1.area_wing = 22.14;     // [m^2] wing area - makes
    aero1.area_hull = 10.0*0;    // [m^2] hull area - makes only drag
    aero1.aspect    = 4.5;       // [1]   aspec ratio of wing - affetcs predominantly lift-induced drag (Oswald efficiency factor included in)
    aero1.CD0       = 0.05;      // [1    min Drag Coef at zero lift
    aero1.CLmax     = 1.0;       // [1]   max Lift coef at stall
    aero1.accelMax  = 5*const_GravAccel;    // [m/s^2]  maximum acceleration
    evalAeroCraft( aero1 );


    printf( "### Bf 109 G6 \n" );
    // https://en.wikipedia.org/wiki/Messerschmitt_Bf_109#Specifications_(Bf_109G-6)
    aero1.mass      = 3400;     // [kg]
    aero1.power     = 1.085e+6; // [W]   engine power
    aero1.area_wing = 16.05;    // [m^2] wing area - makes
    aero1.area_hull = 10.0*0;     // [m^2] hull area - makes only drag
    aero1.aspect    = 6.14;      // [1]   aspec ratio of wing - affetcs predominantly lift-induced drag (Oswald efficiency factor included in)
    aero1.CD0       = 0.03;      // [1    min Drag Coef at zero lift
    aero1.CLmax     = 1.0;       // [1]   max Lift coef at stall
    aero1.accelMax  = 5*const_GravAccel;    // [m/s^2]  maximum acceleration
    evalAeroCraft( aero1 );

    printf( "### F4U-4 \n" );
    // https://en.wikipedia.org/wiki/Republic_P-47_Thunderbolt#Specifications_(P-47D-40_Thunderbolt)
    aero1.mass      = 6592;     // [kg]
    aero1.power     = 1.770e+6; // [W]   engine power
    aero1.area_wing = 29.17;    // [m^2] wing area - makes
    aero1.area_hull = 10.0*0;   // [m^2] hull area - makes only drag
    aero1.aspect    = 5.35;     // [1]   aspec ratio of wing - affetcs predominantly lift-induced drag (Oswald efficiency factor included in)
    aero1.CD0       = 0.03;     // [1    min Drag Coef at zero lift
    aero1.CLmax     = 1.0;      // [1]   max Lift coef at stall
    aero1.accelMax  = 5*const_GravAccel;    // [m/s^2]  maximum acceleration
    evalAeroCraft( aero1 );

    //exit(0);

    printASCItable( 33, 127  );

    world.init();

    camX0 = world.map_center.x;
    camY0 = world.map_center.y;

    currentFaction = world.factions[0];  printf( "currentFaction: %s\n", currentFaction->name );
    currentSquad   = currentFaction->squads[0];

    //TiledView::init( 6, 6 );
    //tiles    = new int[ nxy ];
    //TiledView::renderAll( -10, -10, 10, 10 );

    default_font_texture     = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );
    GUI_fontTex = default_font_texture;


    //printf( " pass.size:  %i \n" , world.pathFinder.pass.size() );
    //for( auto& item : world.pathFinder.pass ){
    //    printf( "pass %i %i %i %i \n", item.first&0xFFFF, item.first>>32, item.second.x, item.second.y );
    //}

    zoom = 1000.00;

    oglTerrain=Draw::list(oglTerrain);
    renderStaticTerrain();
    glEndList();

    printf( "default_font_texture :  %i \n", default_font_texture );

}



void FormationTacticsApp::draw(){
    //long tTot = getCPUticks();
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

    glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

	glCallList( oglTerrain );

    const LTFaction* testFaction = world.factions[0];
    const LTUnit*    testUnit    = &testFaction->squads[0]->units[0];
    Vec2d pmouse = { mouse_begin_x, mouse_begin_y };
    world.prepareSurroundings( testFaction, pmouse, 100.0 );
    drawMap( {20,20}, pmouse-(Vec2d){10.0,10.0}, {1.0,1.0},
        [&](Vec2d p){ return world.tmpSur.unitPosFittness( testUnit, p ); },
        -1.0, 1.0, Draw::ncolors, Draw::colors_RWB
    );

    /*
	if(bDrawing){
        Vec2i ind;Vec2d dind;
        //int i = world.ruler.simplexIndex({mouse_begin_x,mouse_begin_y},ind,dind)>>1;
        world.ruler.simplexIndexBare({mouse_begin_x,mouse_begin_y},ind);
        ind = world.ruler.wrap_index(ind);
        int i = world.ruler.ip2i(ind);
        //printf("i %i (%i,%i) \n",i, ind.x, ind.y );
        world.ground[i] = randf(0.0,1.0);
	}
	*/


    float camMargin = ( camXmax - camXmin )*0.1;
    //float camMargin = 0;
    //TiledView::draw(  camXmin-camMargin, camYmin-camMargin, camXmax+camMargin, camYmax+camMargin  );
    //printf( " camRect  %f %f %f %f \n", camXmin-camMargin, camYmin-camMargin, camXmax+camMargin, camYmax+camMargin );
    //long tComp = getCPUticks();
    world.update( );
    //tComp = getCPUticks() - tComp;

    long tDraw = getCPUticks();
    int i=0;
    for( LTSquad* u : world.squads ){
        if( (u!= NULL) ){
            // TODO : check if on screen
            //printf( "squad %i \n", i ); i++;
            if  ( u == currentSquad ){ render( *u, u->faction->color&0x1fFFFFFF, 1, bDrawGoal ); }
            else                     { render( *u, u->faction->color&0x1fFFFFFF, 1, false     ); }
        }
    }
    tDraw = getCPUticks() - tDraw;

    if( currentSquad != 0 ){
        //glColor3f(1.0,0.0,1.0);
        //glColor3f(0.0,1.0,0.0);
        glColor4f( 0.0,1.0,0.0, 0.2 );
        //Draw2D::drawCircle_d   ( currentSquad->pos, 0.5, 16, false );
        //currentSquad->renderJob( currentSquad->faction->color );
        drawVisibilityIsolines ( world, currentSquad->pos, 5, 50, 0, 2*M_PI, -0.1, +0.1, 500.0, true );
    }

    //world.tmpSur.bConstr=false;
    world.tmpSur.bConstr=true;
    world.tmpSur.ConstrPos = currentSquad->goal;
    world.tmpSur.ConstrRad = currentSquad->goalRadius;
    world.tmpSur.ConstrE   = -1.0;

    world.tmpSur.clear();
    world.getSurroundings( world.tmpSur, currentSquad->faction, currentSquad->pos, 50.0 );
    //plotSurrounding( world.tmpSur, {mouse_begin_x,mouse_begin_y} );
    plotSiteFittness( world.tmpSur, currentSquad->units[0], {mouse_begin_x,mouse_begin_y} );

    world.optimizeDeployment( currentSquad, 50.0, 5, 5, true );
    //for( LTUnit& u: currentSquad->units ){ u.pos = u.goal_pos; }


};

void FormationTacticsApp::drawHUD(){

    if(currentSquad){
        glColor3f(1.0,1.0,1.0);
        glPushMatrix();
        //printf( "currentSquad.type %i %s \n", currentSquad->type, currentSquad->type->name.c_str() );
        //Draw::drawText( "abcdefghijklmnopqrstuvwxyz \n0123456789 \nABCDEFGHIJKLMNOPQRTSTUVWXYZ \nxvfgfgdfgdfgdfgdfgdfg", fontTex, 8, {10,5} );
        glTranslatef( 10.0,HEIGHT-20,0.0  ); currentSquad->type->toStrCaptioned(strBuf,true); Draw::drawText( strBuf, default_font_texture, fontSizeDef, {80,50} );
        //glTranslatef( 300.0,      0.0,0.0 ); glColor3f(0.0,0.5,0.0); currentFormation->reportStatus(strBuf); Draw::drawText( strBuf, fontTex, 8, {80,15} );
        //glTranslatef( 300.0,      0.0,0.0 ); glColor3f(0.0,0.5,0.8); currentFormation->soldiers[0].type->toStrCaptioned(strBuf); Draw::drawText( strBuf, fontTex, 8, {80,15} );
        glPopMatrix();
    }

};

void FormationTacticsApp::eventHandling ( const SDL_Event& event  ){
    //printf( "NBodyWorldApp::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_0:  formation_view_mode = 0;            printf( "view : default\n" ); break;
                //case SDLK_1:  formation_view_mode = VIEW_INJURY;  printf( "view : injury\n"  ); break;
                //case SDLK_2:  formation_view_mode = VIEW_STAMINA; printf( "view : stamina\n" ); break;
                //case SDLK_3:  formation_view_mode = VIEW_CHARGE;  printf( "view : charge\n"  ); break;
                //case SDLK_4:  formation_view_mode = VIEW_MORAL;   printf( "view : moral\n"   ); break;
                case   SDLK_n: ifaction++; if(ifaction>=world.factions.size()) ifaction=0; currentFaction = world.factions[ifaction]; printf("ifaction %i\n",ifaction); break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    if( currentFaction != NULL ) currentSquad = currentFaction->getUnitAt( { mouse_begin_x, mouse_begin_y } );
                    //bDrawing=true;
                break;
                case SDL_BUTTON_RIGHT:
                    //printf( "left button pressed !!!! " );
                    if( currentSquad != NULL ){
                        int imin = world.getUnitAt( { mouse_begin_x, mouse_begin_y }, currentFaction );
                        if( imin > -1 ) {
                            printf( "target selected %i %i\n", imin, world.squads[imin] );
                            currentSquad->setOpponent( world.squads[imin] );
                        }else{
                            printf( "goal selected (%3.3f,%3.3f)\n", mouse_begin_x, mouse_begin_y );
                            currentSquad->setGoal  ( { mouse_begin_x, mouse_begin_y } );
                        }
                    }
                break;
            }
            break;
            /*
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    world.picked = NULL;
                    break;
            }
            break;
            */
    };
    AppSDL2OGL::eventHandling( event );
    camStep = zoom*0.05;
}

// ===================== MAIN

FormationTacticsApp * thisApp;

int main(int argc, char *argv[]){
	// SDL_Init(SDL_INIT_VIDEO);
	// SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
    // SDL_DisplayMode dm;
    // SDL_GetDesktopDisplayMode(0, &dm);
    // //printf( "Display size W: %i H:%i \n", dm.w, dm.h );

    setbuf(stdout, NULL);                 // disable stdout buffering
    //setvbuf(stdout, NULL, _IONBF, 0);  //Or use the more flexible setvbuf:
    // example: use like : ./spaceCraftEditor -s data/ship_ICF_interceptor_1.lua
    printf( "argc %i \n", argc );
    SDL_DisplayMode dm = initSDLOGL( 8 );
    int junk;
    thisApp = new FormationTacticsApp( junk, dm.w-150, dm.h-100 ); 
    
	thisApp->zoom = 30;
	thisApp->loop( 1000000 );
	return 0;
}
















