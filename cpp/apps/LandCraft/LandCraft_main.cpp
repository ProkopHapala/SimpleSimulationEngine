
/*

#####################
##### LandCraft #####
#####################

* Econmy-Building game with focus on use of Geography, and distribution of natural resources
* Build civilization in given natrual conditions ( Terrain, Rivers, distribution of minral )
* solve problem of transport of resources
    * it is more efficient transport it by river (down-the-flow), by train (problematic in hills), by road (laborious,low capacity, high energy consumption)
* Modify terrain
    * better trains
    * build Damns
    * watter channels ( Suez, Volga-Baltic )


### ToDo :

* Path-finding Roads - minimum energy path between two points with respect to given vehicle ( slope vs distance )
    * pathFinder
        * /cpp/common/maps/PathFinder.h
        * /cpp/apps/LandTactics/LandTactics_main.cpp
    * Roads are asymmetric
        * down-hill vs. up-hill and both-direction compromise (different energy cost)
        * rivers down-the-flow vs up-the-flow
            * raft work only down-the-flow
* Use new Commodity-Network solver (from different files)
    * /cpp/common/dynamics/CommodityNetwork.h
    * /cpp/sketches_SDL/2D/test_CommodityNetwork.cpp
*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <unordered_map>
#include <math.h>
#include <assert.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include <SDL2/SDL_image.h>
//#include <SDL2/SDL_ttf.h>
//#include "Texture.h"
#include "Draw.h"
#include "Draw2D.h"

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"
#include "AppSDL2OGL.h"

#include "testUtils.h"
#include "SDL_utils.h"

#include "TerrainCubic.h"
#include "TiledView.h"

#include "SimplexRuler.h"
#include "Ruler2DFast.h"
#include "TerrainHydraulics.h"
#include "hydraulics1D.h"
#include "PathFinder.h"    // ToDo : Use this for searching roads

#include "macroUtils.h"
#include "IO_utils.h"

#include "CommandParser.h"

#include "Economy.h"
//#include "CommodityNetwork.h"  // ToDo : Use This instead of Economy.h
#include "Roads.h"

#include "GUI.h"
#include "Plot2D.h"

#include <algorithm>

#include "LandCraftWorld.h"

// C API for LandCraft (used by Lua bindings and some app helpers)
//#include "../../libs/LandCraft/LandCraftLib.h"

LandCraftWorld W; // shared world instance

#ifdef LUA
#include "LuaHelpers.h"
#include "LandCraftLua.h"
static lua_State* gLua = nullptr;
#endif

#include "argparse.h"


// font rendering:
//  http://www.willusher.io/sdl2%20tutorials/2013/12/18/lesson-6-true-type-fonts-with-sdl_ttf
//  http://stackoverflow.com/questions/28880562/rendering-text-with-sdl2-and-opengl

int   default_font_texture;

void cmapHeight(double g){
    /// Terrain Color-map
    //glColor3f(0.2+0.8*g*g,0.2+0.3*(1-g)*g,0.2);
    //double snow = g*g*g*g;
    //glColor3f(g*0.8+0.2,0.5+0.5*snow,0.2+0.8*snow);
    glColor3f(g,g,g);
    //return g;
}

FILE* file_relax_debug = 0;

class LandCraftApp : public AppSDL2OGL { public:
    int      fontTex = 0;

    enum class Actions{ generateTerrain, relaxWater, relaxWaterHex, save, load, outflow, inflow, findRivers, terrainViewMode, traceDroplet };

    // ---- World propertis which are being moved to LandCraftWorld W;
    // ---- Map Geography & 2D Hydraulics
    // Vec2d  map_center;
    // double maxHeight = 500.0;
    SimplexRuler     ruler;
    Ruler2DFast      square_ruler;
    // HydraulicGrid2D  hydraulics;
    // double * ground   = NULL;
    // double * water    = NULL;
    double drawHeight = 0;
    // ---- Hydraulics 1D
    Hydraulics1D hydro1d;
    const int   nTraceMax = 256;
    int         nTrace    = 0;
    int       * trace     = NULL;
    // std::vector<int> river;
    // std::vector<int> feeders;
    // // ------- Economy
    // std::unordered_map<std::string,Commodity*>  commodities;
    // std::unordered_map<std::string,Technology*> technologies;
    // // ------- Roads & Vehicles
    // RoadBuilder roadBuilder;
    // std::vector<RoadVehicleType*>  vehicleTypes;
    // std::vector<Road*>             roads;
    // std::vector<RoadVehicle*>      vehicles;

    // ------- GUI
    GUI gui;
    DropDownList* riverList = 0;
    CheckBoxList  chklDrawLayers;
    CommandList   cmdList1;

    Commander     keybinds;
    CommandParser cmdPars;

    Plot2D riverProfile;
    Plot2D roadProfile;
    int iBuildStartHex = -1;

    // ------- Terrain Rendering
    bool bDrawing = false;
    int terrainViewMode  = 1;

    // ------- Update & Draw Swithces
    int  doDrain=0;
    bool bRunHydroRelax=0;
    bool bRunHydro1D=0;
    bool bDrawRivers=0;
    bool bDrawTraceDroplet=0;
    bool bDrawRoads=0;
    bool bRunVehicles=0;
    bool bDrawVehicles=0;
    bool bRoadProfile=0;
    bool bRiverProfile=0;

    // =================================
    // ========= Functions  ============
    // =================================

    void printASCItable( int imin, int imax  );
    //GLuint makeTexture( char * fname );
    //GLuint renderImage( GLuint itex, const Rect2d& rec );
    //void drawString( char * str, int imin, int imax, float x, float y, float sz, int itex );
    //void drawString( char * str, float x, float y, float sz, int itex );

	virtual void draw   ();
	virtual void drawHUD();
	virtual void mouseHandling( );
	virtual void eventHandling( const SDL_Event& event  );
	void debug_buffinsert( );
	//void pickParticle( Particle2D*& picked );
	//virtual int tileToList( float x0, float y0, float x1, float y1 );

	void registerCommands();
    void registerDrawLayers();

    // void generateTerrain();
    // void loadTechnologies( const char* fname);
    // void makeMap( int sz, double step, bool newMap );
    // void makeRivers();
    // void makeRoads();
    // void makeVehicles();

    void hydro1D_update();
	void hydroRelaxUpdate();
    void updateVehicles();

	void terrainColor( int i );
	void drawTerrain( Vec2i i0, Vec2i n, int NX );
	void drawRoad( Road* road );
	void drawRiver( River* road );
    void drawRivers();
    void drawDropletTrace();
    void drawSinks();

    void drawVehicles();

    int addRoadStright( Vec2i p1, Vec2i p2 ){
        int id = W.makeRoadStraight(p1,p2);
        if(id>=0 && !W.roads.empty()) addRoadPlot(W.roads.back());
        return id;
    }
    void addRoadPlot(Road* road);
    void addRiverPlot(int iRiver, float scFlow);

    void commandDispatch( int icommand );

    LandCraftApp( int& id, int WIDTH_, int HEIGHT_ );

};

LandCraftApp::LandCraftApp( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {
    //default_font_texture   = makeTexture    ( "common_resources/dejvu_sans_mono.bmp" );
    default_font_texture = makeTexture    ( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    fontTex                = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );

    riverProfile.fontTex = fontTex;
    roadProfile.fontTex  = fontTex;
    GUI_fontTex = fontTex;
    printASCItable( 33, 127  );

    registerCommands();
    registerDrawLayers();

    // TODO(modular-init): Config-driven initialization; see data/world.ini
    // The old explicit sequence is kept commented below for reference.
    W.loadConfig("data/world.ini");
    // Ensure the hex ruler matches the world grid, otherwise drawTerrain indices go out of bounds
    if(W.ready){
        ruler.setSize(W.hydro.n.x, W.hydro.n.y);
        ruler.setStep(50.0);
    }
    // W.loadTechnologies("data/Technologies.txt");
    // W.makeMap(128,50,false);
    // W.generateTerrain();
    // W.makeRivers();
    // W.makeRoadStraight({10,15}, {55,38});
    // W.makeVehicles();
    cmdPars.execFile( "data/comands.ini" );


    //hydro1d.realloc(256);
    //hydro1d.clear();
    //bisectNoise1D(8,hydro1d.ground,-1.0,0.0);
    hydro1d.realloc(512);
    hydro1d.clear();
    bisectNoise1D(9,hydro1d.ground,-1.0,0.0);
    //hydro1d.realloc(16);
    //hydro1d.clear();
    //bisectNoise1D(4,hydro1d.ground,-1.0,0.0);
    VecN::set(hydro1d.n,5.0,hydro1d.water);
    // allocate droplet trace buffer now that map exists
    trace = new int[nTraceMax];


    // Fill random water on top of ground to visualize
    if(W.ready){
        for(int i=0; i<W.hydro.ntot; i++) W.hydro.water[i] = W.hydro.ground[i] + randf(0,20.0);
        W.relaxHex(36,39);
    }

    gui.layoutRow(10,10);

#ifdef LUA
    if(!gLua){
        gLua = luaL_newstate();
        luaL_openlibs(gLua);
        LandCraftLua::register_api(gLua);
        // Optional demo script auto-run
        if( fileExist("data/landcraft_demo.lua") ){
            Lua::dofile(gLua, "data/landcraft_demo.lua");
        }
    }
#endif
}

// ------------------------------------------------------------------

void LandCraftApp::draw(){
    //long tTot = getCPUticks();
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );
    glShadeModel ( GL_SMOOTH );

    // ----- Terrain
    glPushMatrix();
	glScalef(ruler.step,ruler.step,1.0);
	drawTerrain( {0,0}, {128,128}, ruler.na );
	glPopMatrix();

	// ----- Hydraulics 2D
    if(bRunHydroRelax) W.hydroRelaxStep();
    //if(bRunOutflow)    hydraulics.outflow_step();
	//if( doDrain==1 ){ hydraulics.outflow_step(); }else if (doDrain==-1){ hydraulics.inflow_step();  }

	// ----- Hydraulics 1D
    if(bRunHydro1D)       hydro1D_update  ( );
    if(bDrawRivers)       drawRivers      ( );
    if(bDrawTraceDroplet) drawDropletTrace( );

	// ----- Transport (Roads, Vehicles)
    if(bDrawRoads   ){ glColor3f( 1.0, 1.0, 1.0 ); if(!W.roads.empty()) drawRoad( W.roads[0] ); }
    if(bRunVehicles ) W.vehiclesStep(2.0);
    if(bDrawVehicles)drawVehicles();

    // ----- Mouse  & Misc (debug?)
    Draw2D::drawPointCross({mouse_begin_x,mouse_begin_y}, 100.0);

};

// --------------------------------------------------------------

void LandCraftApp::drawHUD(){
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);

    glPushMatrix();
    glTranslatef(WIDTH-300,0,200.0);
    if(bRoadProfile){
        roadProfile.view();
        glTranslatef(100,0.0,0.0);
    }
    if(bRiverProfile)riverProfile.view();
    glPopMatrix();

    gui.draw();
}

// --------------------------------------------------------------

void LandCraftApp::commandDispatch( int icommand ){
    switch( (Actions)icommand ){
        case Actions::generateTerrain : W.generateTerrain(); break;
        case Actions::relaxWater      : W.hydro.relaxWater(); break;
        //case myKeyBind: // does not work
        case Actions::relaxWaterHex   :{
            int ihex = ruler.hexIndex({mouse_begin_x,mouse_begin_y});
            W.hydro.relaxWater( ruler.i2ip(ihex) );
            }break;
        case Actions::save :
            saveBin( "data/ground.bin", sizeof(double)*W.hydro.ntot, (char*)W.hydro.ground );
            saveBin( "data/water.bin",  sizeof(double)*W.hydro.ntot,  (char*)W.hydro.water  );
            break;
        case Actions::load :{
            bool ok = true;
            if( fileExist("data/ground.bin") ){
                if( loadBin( "data/ground.bin", sizeof(double)*W.hydro.ntot, (char*)W.hydro.ground, /*bExitOnErr*/false ) != 0 ) ok = false;
            }else{ ok = false; }
            if( fileExist("data/water.bin") ){
                if( loadBin( "data/water.bin", sizeof(double)*W.hydro.ntot,  (char*)W.hydro.water,  /*bExitOnErr*/false ) != 0 ) ok = false;
            }else{ ok = false; }
            if(!ok){
                printf("[LandCraft] Load requested but terrain files missing/corrupt -> generating new terrain.\n");
                W.generateTerrain();
            }
            terrainViewMode = 1;
            }break;
        case Actions::outflow :{
            for(int i=0; i<W.hydro.ntot; i++){ W.hydro.known[i]=false; }
            int ihex = ruler.hexIndex({mouse_begin_x,mouse_begin_y});
            W.hydro.contour2[0] = ihex;
            W.hydro.nContour++;
            printf( "idrain  %i \n", ihex  );
            W.hydro.water[ihex]  = W.hydro.ground[ihex];
            W.hydro.isOutflow    = true;
            terrainViewMode = 1;
            //doDrain = 1;
            }break;
        case Actions::inflow :{
            for(int i=0; i<W.hydro.ntot; i++){ W.hydro.known[i]=false; }
            int ihex = ruler.hexIndex({mouse_begin_x,mouse_begin_y});
            W.hydro.contour2[0] = ihex;
            W.hydro.nContour++;
            printf( "idrain  %i \n", ihex   );
            W.hydro.water[ihex]  = W.hydro.ground[ihex] + 10.0;
            W.hydro.isOutflow    = false;
            terrainViewMode = 1;
            //doDrain = -1;
            }break;
        case Actions::findRivers :
            W.hydro.gatherRain( 100.0 );
            terrainViewMode = 2;
            //val=0.0; ihex=0;
            //for( int isink : W.hydro.sinks ){ double w=W.hydro.water[isink]; if(w>val){val=w; ihex=isink; } } // find largest river
            //W.hydro.trackRiver( ihex, 50.0, river, feeders );
            W.hydro.findAllRivers( 50.0 );
            //for(int i=0; i<W.hydro.ntot; i++){ W.hydro.known[i]=false; }
            //W.hydro.trackRiverRecursive( ihex, 50.0, NULL );
            break;
        case Actions::terrainViewMode :
            terrainViewMode=(terrainViewMode%2)+1;
            break;
        case Actions::traceDroplet :{
            int ihex   = ruler.hexIndex({mouse_begin_x,mouse_begin_y});
            nTrace = W.hydro.traceDroplet( {ihex%W.hydro.n.x, ihex/W.hydro.n.x}, nTraceMax, trace );
            }break;
    }
    roadProfile.clear();
    addRoadPlot(W.roads[0]);
    roadProfile.update();
    roadProfile.autoAxes(0.5,0.2);
    roadProfile.render();
}

// --------------------------------------------------------------

void LandCraftApp::eventHandling ( const SDL_Event& event  ){
    //printf( "NBodyWorldApp::eventHandling() \n" );
    int ihex;
    double val;

    //SDL_EventType keyT = SDLK_e;
    int myKeyBind = SDLK_e;

    GUIAbstractPanel* activeGUIPanel = gui.onEvent(mouseX,mouseY,event);
    switch( event.type ){
        case SDL_KEYDOWN :{
            if( !cmdList1.getKeyb( event.key.keysym.sym ) ){
                auto got = keybinds.keymap.find( event.key.keysym.sym );
                if( got!=keybinds.keymap.end() ){
                    commandDispatch( got->second );
                }
            }
            } break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    if( activeGUIPanel == 0 ){
                        ihex = ruler.hexIndex( {mouse_begin_x, mouse_begin_y} );
                        if( (ihex>=0) && (ihex<W.hydro.ntot) ){
                            drawHeight = W.hydro.ground[ihex];
                            printf( "[LandCraft] pick drawHeight=%f at ihex=%d\n", drawHeight, ihex );
                        }else{
                            printf( "[LandCraft] left-click outside map: ihex=%d ntot=%d (ignoring)\n", ihex, W.hydro.ntot );
                        }
                    }
                    else if( activeGUIPanel == riverList ){
                        riverProfile.clear();
                        //iCurrentRiver = riverList.iSelected;
                        addRiverPlot(riverList->iSelected,  -1 );
                        addRiverPlot(riverList->iSelected,  0.1);
                    }
                    break;
                case SDL_BUTTON_RIGHT:
                    //printf( "left button pressed !!!! " );
                    iBuildStartHex = ruler.hexIndex( {mouse_begin_x, mouse_begin_y} );
                    if( (iBuildStartHex<0) || (iBuildStartHex>=W.hydro.ntot) ){
                        printf( "[LandCraft] right-click start outside map: ihex=%d ntot=%d\n", iBuildStartHex, W.hydro.ntot );
                        iBuildStartHex = -1;
                    }
                    break;
            }
            break;

        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_RIGHT:
                    if( iBuildStartHex>=0 ){
                        int ihex = ruler.hexIndex( {mouse_begin_x, mouse_begin_y} );
                        if( (ihex>=0) && (ihex<W.hydro.ntot) ){
                            for(Road* road : W.roads){ delete road; }
                            W.roads.clear();
                            addRoadStright( ruler.i2ip(iBuildStartHex), ruler.i2ip(ihex) );

                            if(!W.roads.empty()){
                                roadProfile.clear();
                                addRoadPlot(W.roads[0]);
                                roadProfile.update();
                                roadProfile.autoAxes(0.5,0.2);
                                roadProfile.render();
                            }
                        }else{
                            printf( "[LandCraft] right-click end outside map: ihex=%d ntot=%d (ignoring)\n", ihex, W.hydro.ntot );
                        }
                    }
                    iBuildStartHex=-1;
                    break;
            }
            break;

        case SDL_MOUSEWHEEL:
            //if( event.type == SDL_MOUSEWHEEL ){
            //printf( " SDL_MOUSEWHEEL \n" );
            break;
    };
    AppSDL2OGL::eventHandling( event );
    camStep = zoom*0.05;
}


// ------------------------------------------------------------------

void LandCraftApp::mouseHandling( ){
    uint32_t buttons = SDL_GetMouseState( &mouseX, &mouseY );
    mouseY=HEIGHT-mouseY;
    defaultMouseHandling( mouseX, mouseY );
    mouse_begin_y = mouseUp_(mouseY) + camY0;
    //printf( "mouse %g %g \n", mouse_begin_x, mouse_begin_y );
    //mouse_begin_y = -mouse_begin_y;
    //mouse_begin_y = HEIGHT - mouse_begin_y;
    if( buttons & SDL_BUTTON(SDL_BUTTON_LEFT) ){
        int ihex = ruler.hexIndex( {mouse_begin_x, mouse_begin_y} );
        if( (ihex>=0) && (ihex<W.hydro.ntot) ){
            W.hydro.ground[ihex] = drawHeight;
        }else{
            printf("[LandCraft] mouse paint outside map: ihex=%d ntot=%d\n", ihex, W.hydro.ntot);
            // Click outside map; ignore safely
            //printf("[LandCraft] mouse paint outside map: ihex=%d ntot=%d\n", ihex, W.hydro.ntot);
        }
    }
};

// ------------------------------------------------------------------
// ------------------------------------------------------------------
//                    IMPLEMENTATION DETIALS
// ------------------------------------------------------------------
// ------------------------------------------------------------------

void LandCraftApp::registerCommands(){
    #define _bindKey(key,name) _.add( (int)Actions::name, (int)key, #name );
    {WITH(keybinds)
        _bindKey( SDLK_n, generateTerrain );
        _bindKey( SDLK_r, relaxWater      );
        _bindKey( SDLK_e, relaxWaterHex   );
        _bindKey( SDLK_s, save            );
        _bindKey( SDLK_l, load            );
        _bindKey( SDLK_o, outflow         );
        _bindKey( SDLK_i, inflow          );
        _bindKey( SDLK_g, findRivers      );
        _bindKey( SDLK_m, terrainViewMode );
        _bindKey( SDLK_t, traceDroplet    );
    }
    #undef _bindKey
    for( Command& cmd : keybinds.commands ){
        printf( "key %i %c %s\n", cmd.key, (char)cmd.key, cmd.name.c_str() );
    };
    cmdList1.initCommandList(300,100,450);
    cmdList1.commands = &keybinds;
    gui.addPanel(&cmdList1);
    //cmdList1.commandDispatch = &LandCraftApp::commandDispatch;
    //cmdList1.commandDispatch = &(this->commandDispatch);
    cmdList1.commandDispatch = [&](int i){ commandDispatch(i); };
    //cmdList1.commandDispatch = std::bind( &LandCraftApp::commandDispatch, this  );
    cmdList1.bDispatch=true;
}

void LandCraftApp::registerDrawLayers(){
    chklDrawLayers.initCheckBoxList( 100,100, 200);
    {WITH(chklDrawLayers)
        _addBox(bRunHydroRelax);
        _addBox(bRunHydro1D);
        _addBox(bDrawRivers);
        _addBox(bDrawTraceDroplet);
        _addBox(bDrawRoads);
        _addBox(bRunVehicles);
        _addBox(bDrawVehicles);
        _addBox(bRoadProfile);
        _addBox(bRiverProfile);
    }
    gui.addPanel(&chklDrawLayers);
}

void LandCraftApp::printASCItable( int imin, int imax  ){
    /// print ASCITable for debugging
    assert(imax > imin); // DEBUG: ensure positive length
    int len = imax - imin;
    char str[len + 1]; // +1 for null terminator
    for ( int i=0; i<len; i++ ){
        str[i] = (char)(i+imin);
    }
    str[len] = '\0';
    printf("%s\n", str );
};

// ------------------------------------------------------------------

void LandCraftApp::addRiverPlot(int iRiver, float scFlow){
    River* thisRiver = W.hydro.rivers[iRiver];
    //DataLine2D * riverH1 = new DataLine2D(thisRiver->path.size());
    //riverProfile.lines.push_back( riverH1 );
    //riverH1->linspan(0,thisRiver->path.size());
    DataLine2D * riverH1 = riverProfile.add( new DataLine2D( thisRiver->path.size(),0,1, hash_Wang( (iRiver+15464)*5646 ) | (0xFF*(scFlow<0)), "ground") );
    int ii=0;
    if(scFlow>0){
        for( int i : thisRiver->path ){
            //W.hydro.i2ip(i);
            //riverH1->ys[ii] = W.hydro.water[i] * scFlow;
            riverH1->ys[ii] = thisRiver->flow[ii] * scFlow;
            //riverH1->ys[ii] = thisRiver->flow[ii];
            //printf( "road[%i] h flow %g \n", ii, riverH1->ys[ii] );
            ii++;
        };
    }else{
        for( int i : thisRiver->path ){
            riverH1->ys[ii] = W.hydro.ground[i];
            //printf( "road[%i] h height %g \n", ii, riverH1->ys[ii] );
            ii++;
        };
    }
    //exit(0);
    //riverH1->clr = hash_Wang( (iRiver+15464)*5646 );
    //DataLine2D * lDrag = new DataLine2D(nsamp); mainWingL.lines.push_back( lDrag ); lDrag->linspan(-phiRange,phiRange); lDrag->clr = 0xFF0000ff;
    riverProfile.update();
    riverProfile.autoAxes(0.5,0.2);
    riverProfile.render();
}

// ------------------------------------------------------------------

void LandCraftApp::addRoadPlot(Road* road){
    int n = road->n;

    //DataLine2D * roadH1 = new DataLine2D(n);
    //roadProfile.lines.push_back( roadH1 );
    //roadH1->linspan(0,n);

    DataLine2D * roadH1 = roadProfile.add( new DataLine2D(n,0,1, 0xFFff0000, "ground") );
    DataLine2D * roadH2 = roadProfile.add( new DataLine2D(n,0,1, 0xFF0000ff, "water" ) );

    //DataLine2D * roadH2 = new DataLine2D(n);roadProfile.lines.push_back( roadH2 );roadH2->linspan(0,n);

    for( int ii=0; ii<n; ii++ ){
        RoadTile& p = road->path[ii];
        int i = W.hydro.ip2i({p.ia,p.ib});
        double h = W.hydro.ground[i];
        double w = W.hydro.water [i];
        roadH1->ys[ii] = w;
        roadH2->ys[ii] = h;
        //roadH2->ys[ii] = h + w;
        //printf( "road[%i] h %g w %g w+h %g p.h %g \n", ii, h,w,h+w, p.height );
    };
    //exit(0);
    //roadH2->clr = 0xFF0000ff;
    //roadH1->clr = 0xFFff0000;
    //DataLine2D * lDrag = new DataLine2D(nsamp); mainWingL.lines.push_back( lDrag ); lDrag->linspan(-phiRange,phiRange); lDrag->clr = 0xFF0000ff;
}


// ------------------------------------------------------------------

void LandCraftApp::terrainColor( int i ){
    double maxDepth = 50.0;
    double depth;
    double w = W.hydro.water [i];
    double g = W.hydro.ground[i];
    switch( terrainViewMode ){
        case 1:
            depth = clamp( (w-g)/maxDepth, 0.0, 1.0 );
            //w-=g; if(w)depth=clamp( sqrt(w-g)*0.05, 0.0, 1.0 );
            g /= W.maxHeight;
            glColor3f( (1-depth)*g, (1-depth)*(1-g), depth );
            break;
        case 2:
            //w *= 0.1;
            //w *= log(w*5.0)*0.1;
            w=sqrt(w)*0.05;
            if( W.hydro.known[i] ){ glColor3f( 1, 0, 0 ); }else{ glColor3f( w, w, w ); }
            break;
    }
    //glColor3f( 1.0,1.0,1.0 );
}

// ------------------------------------------------------------------

void LandCraftApp::drawTerrain( Vec2i i0, Vec2i n, int NX ){
    Vec2f a,b,p;
    a.set( 1.0, 0.0           ); //a.mul(scale);
    b.set( 0.5, 0.86602540378 ); //b.mul(scale);
    //glDisable(GL_SMOOTH);
    //int ii = 0;
    //double renorm=1.0d/(vmax-vmin);
    for (int iy=0; iy<n.y-1; iy++){
        glBegin( GL_TRIANGLE_STRIP );
        int ii = (i0.y+iy)*NX + i0.x;
        for (int ix=0; ix<n.x; ix++){
            p.set( ix*a.x+iy*b.x, ix*a.y+iy*b.y );
            terrainColor( ii    ); glVertex3f( p.x    , p.y    , 0 );
            terrainColor( ii+NX ); glVertex3f( p.x+b.x, p.y+b.y, 0 );
            ii++;
        }
        glEnd();
    }
}

// ------------------------------------------------------------------

void LandCraftApp::drawRiver( River* river ){
    //printf( " path size: %i \n", river->path.size() );
    //Draw::color_of_hash( river->path[0]+54877 );
    //glLineWidth(10.0);
    glBegin(GL_LINES);
    int ii=0;
    Vec2d p;
    ruler.nodePoint( river->path[0],p);
    glBegin(GL_LINES);
    for(int i : river->path ){
        glLineWidth( river->flow[ii]*0.005 );
        glBegin(GL_LINES);
        glVertex3f( p.x, p.y, 100.0 );
        ruler.nodePoint(i,p);
        glVertex3f( p.x, p.y, 100.0 );
        glEnd();
        ii++;
    }
    //glEnd();
    /*
    //int ii=0;
    glBegin(GL_LINE_STRIP);
    Vec2d p;
    for(int i : river->path ){
        //glLineWidth( river->flow[ii]*0.01 );
        ruler.nodePoint(i,p);
        glVertex3f( p.x, p.y, 100.0 );
        //ii++;
    }
    glEnd();
    */
}

void LandCraftApp::drawSinks(){
    // draw sinks
    glColor3f( 0.0, 1.0, 1.0 );
    for(int isink : W.hydro.sinks ){
        Vec2d p;
        ruler.nodePoint( isink,p);
        Draw2D::drawPointCross_d( p, 100.0 );
    }
}

// --------------------------------------------

void LandCraftApp::drawDropletTrace(){
    // draw droplet trace
    glBegin(GL_LINE_STRIP);
    glColor3f( 1.0, 0.0, 1.0 );
    for(int ii=0; ii<nTrace; ii++){
        Vec2d p;
        ruler.nodePoint( trace[ii],p);
        glVertex3f( p.x, p.y, 100.0 );
    }
    glEnd();
}

// ----------------------------------------

void LandCraftApp::drawRivers( ){
    // draw all rivers from world
    for( River* rv: W.hydro.rivers ){
        Draw::color_of_hash( rv->path[0]+54877 );
        drawRiver( rv );
    }
    // highlight selected
    if(riverList && riverList->iSelected>=0 && riverList->iSelected<(int)W.hydro.rivers.size()){
        glLineWidth( 3 );
        glColor3f( 0.0,1.0,0.0 );
        drawRiver( W.hydro.rivers[riverList->iSelected] );
        glLineWidth( 1 );
    }
}

// ------------------------------------------------------------------

void LandCraftApp::drawRoad( Road* road ){
    glBegin( GL_LINE_STRIP );
    RoadTile* path = road->path;
    //printf( "road->n %i\n", road->n );
    for(int i=0; i<road->n; i++){
        Vec2d p;
        ruler.nodePoint( {path[i].ia,path[i].ib}, p );
        glVertex3f( (float)p.x, (float)p.y, 100.0f );
        //printf( "(%i,%i)  %f %f \n", path[i].ia,path[i].ib, p.x, p.y );
        //glVertex3f();
    }
    glEnd();
    //exit(0);
}

// ------------------------------------------------------------------

void LandCraftApp::     updateVehicles   (){
	for( RoadVehicle* veh : W.vehicles ){
        //veh->moveStep( 1.0 );
        veh->move( 2.0 );
        //if( !veh->onWay ){ veh->idir*=-1; veh->onWay=true; }
        if( !veh->onWay ) veh->depart();
	}
}
// ------------------------------------------------------------------

void LandCraftApp::     drawVehicles    (){
	//if(frameCount>100) exit(0);
	for( RoadVehicle* veh : W.vehicles ){
        Vec2d p;
        //veh->ipath
        RoadTile& tile = veh->road->path[veh->ipath];
        ruler.nodePoint( {tile.ia,tile.ib}, p );
        Draw2D::drawCircle_d( p, 25, 16, false );
	}
}

// ------------------------------------------------------------------

void LandCraftApp::hydro1D_update(){
    printf( "frame %i \n", frameCount );
    for(int i=0; i<10; i++)
    //hydro1d.step( 9.0, 10.0 );
    hydro1d.step( 0.0, 0.0 );
    //hydro1d.step( 0.008, 0.01 );
    //hydro1d.step( 0.08, 0.1 );

    double wtot1 = VecN::sum(hydro1d.n,hydro1d.water);
    //hydro1d.deepAccel(0.1);
    double wtot2 = VecN::sum(hydro1d.n,hydro1d.water);
    printf( "water %g -> %g \n", wtot1, wtot2 );
    //Draw2D::plot( , );
    double dx = 3.0;
    glColor3f(0,0,1);
    glBegin(GL_LINE_STRIP);
    for(int i=0;i<hydro1d.n;i++){
        double h = hydro1d.ground[i] + hydro1d.water[i];
        //double g = hydro1d.ground[i];
        glVertex3f( i*dx, h, 0.0 );
    }
    glEnd();
    glColor3f(1,0,0);
    Draw2D::plot(hydro1d.n,dx, hydro1d.ground);
}

// ------------------------------------------------------------------

void LandCraftApp::hydroRelaxUpdate(){
    // ============ // Hydraulic relaxation
    if(!file_relax_debug) file_relax_debug = fopen( "data/relax_debug.log", "w" );
    terrainViewMode = 1;
    double DEBUG_wsum_before=0.0;
    double DEBUG_E_before   =0.0;
    W.hydro.sumWater( DEBUG_wsum_before, DEBUG_E_before );
    // TODO(deprecated): simulation moved to W.hydroRelaxStep(); keep this wrapper for diagnostics only
    W.hydroRelaxStep();
    double DEBUG_wsum_after=0.0;
    double DEBUG_E_after   =0.0;
    W.hydro.sumWater( DEBUG_wsum_after, DEBUG_E_after );
    printf( "DEBUG relaxWater[%i] dE %g | %g -> %g \n", frameCount, DEBUG_E_before-DEBUG_E_after,  DEBUG_E_before, DEBUG_E_after );
    if(fabs(DEBUG_wsum_before-DEBUG_wsum_after)>1e-6 ){
        printf( "DEBUG ERROR in hydroRelaxStep water before,after diff %g | %g %g \n", DEBUG_wsum_before-DEBUG_wsum_after, DEBUG_wsum_before, DEBUG_wsum_after );
        exit(0);
    }
    if(frameCount<=200) fprintf( file_relax_debug, " %i %g %g \n", frameCount, DEBUG_E_after,  DEBUG_E_before-DEBUG_E_after );
    if(frameCount==200) fclose(file_relax_debug);
    //SDL_Delay(200);
    //return;

	if(bDrawing){
        Vec2i ind;Vec2d dind;
        //int i = ruler.simplexIndex({mouse_begin_x,mouse_begin_y},ind,dind)>>1;
        ruler.simplexIndexBare({mouse_begin_x,mouse_begin_y},ind);
        ind = ruler.wrap_index(ind);
        int i = ruler.ip2i(ind);
        //printf("i %i (%i,%i) \n",i, ind.x, ind.y );
        W.hydro.ground[i] = randf(0.0,1.0);
	}
}

// ===================== MAIN

LandCraftApp * thisApp;

int main(int argc, char *argv[]){

    // SDL_Init(SDL_INIT_VIDEO);
	// SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
    // SDL_DisplayMode dm;
    // SDL_GetDesktopDisplayMode(0, &dm);
	// int junk;

    // disable stdout buffering
    setbuf(stdout, NULL);
    //Or use the more flexible setvbuf:
    //setvbuf(stdout, NULL, _IONBF, 0); 

    // example: use like : ./spaceCraftEditor -s data/ship_ICF_interceptor_1.lua
    printf( "argc %i \n", argc );
    // Initialize shared world state first so both app and C-API operate on same instance
    W.init("cpp/apps/LandCraft/data");

    SDL_DisplayMode dm = initSDLOGL( 8 );
    int junk;
    thisApp = new LandCraftApp( junk, dm.w-150, dm.h-100 ); 

    // ---- Command-line arguments (e.g., run a Lua script on startup)
    LambdaDict funcs;
#ifdef LUA
    funcs["-lua"] = {1, [&](const char** ss){
        const char* fname = ss[0];
        printf("[LandCraft] CLI: -lua %s\n", fname);
        if(gLua){ Lua::dofile(gLua, fname); }
    }};
#endif
    process_args(argc, argv, funcs, /*bExitOnMiss*/false);

	//thisApp = new LandCraftApp( junk , dm.w-150, dm.h-100 );
	//thisApp = new LandCraftApp( junk , 1000,600 );
	//SDL_SetWindowPosition(thisApp->window, 100, 0 );
	thisApp->zoom  = 3000;
	thisApp->camX0 = 4000;
	thisApp->camY0 = 2000;
	thisApp->loop( 1000000 );
	return 0;
}
















