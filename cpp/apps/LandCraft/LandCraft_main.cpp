
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <unordered_map>
#include <math.h>

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

#include "IO_utils.h"

#include "CommandParser.h"

#include "Economy.h"
#include "Roads.h"

#include "GUI.h"
#include "Plot2D.h"

#include <algorithm>




#include "hydraulics1D.h"


// font rendering:
//  http://www.willusher.io/sdl2%20tutorials/2013/12/18/lesson-6-true-type-fonts-with-sdl_ttf
//  http://stackoverflow.com/questions/28880562/rendering-text-with-sdl2-and-opengl

int   default_font_texture;

void cmapHeight(double g){
    //glColor3f(0.2+0.8*g*g,0.2+0.3*(1-g)*g,0.2);
    //double snow = g*g*g*g;
    //glColor3f(g*0.8+0.2,0.5+0.5*snow,0.2+0.8*snow);
    glColor3f(g,g,g);
    //return g;
}






FILE* file_relax_debug = 0;





class LandCraftApp : public AppSDL2OGL { public:

    Hydraulics1D hydro1d;

    SimplexRuler       ruler;
    Ruler2DFast        square_ruler;
    //TerrainHydraulics  hydraulics;
    HydraulicGrid2D  hydraulics;
    double * ground   = NULL;
    double * water    = NULL;

    Vec2d map_center;
    double maxHeight = 500.0;

    double drawHeight = 0;

    //int doDrain = 0;
    //int       ifaction = 0;
    bool bDrawing = false;
    int terrainViewMode  = 1;

    const int nTraceMax = 256;
    int       nTrace=0;
    int * trace = NULL;

    std::vector<int> river;
    std::vector<int> feeders;

    std::unordered_map<std::string,Commodity*>  commodities;
    std::unordered_map<std::string,Technology*> technologies;


    RoadBuilder roadBuilder;
    std::vector<RoadVehicleType*>  vehicleTypes;
    std::vector<Road*>         roads;
    std::vector<RoadVehicle*>  vehicles;

    CommandParser cmdPars;

    Plot2D riverProfile;
    Plot2D roadProfile;

    int      fontTex = 0;
    //GLuint       itex;
    GUI gui;
    DropDownList* riverList = 0;

    int iBuildStartHex = -1;

    // ==== function declaration
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


    void generateTerrain();
	void terrainColor( int i );
	void drawTerrain( Vec2i i0, Vec2i n, int NX );
	void drawRoad( Road* road );
	void drawRiver( River* road );

	void addRoadStright( Vec2i p1, Vec2i p2 );
	void addRoadPlot(Road* road);
	void addRiverPlot(int iRiver, float scFlow);

	LandCraftApp( int& id, int WIDTH_, int HEIGHT_ );

};

void LandCraftApp::printASCItable( int imin, int imax  ){
    int len = imax-imin;
    char str[len];
    for ( int i=0; i<len; i++ ){
        str[i] = (char)(i+imin);
    }
    printf("%s\n", str );
};


void LandCraftApp::generateTerrain(){
    //hydraulics.genTerrainNoise( 8, 2.0, 1.0,  0.5, 0.8, 45454, {100.0,100.0} );
    //hydraulics.genTerrainNoise( 8, 2.0, 1.0,  0.5, 0.8, rand(), {100.0,100.0} );

    srand(16464);
    hydraulics.ground[0]=0.2;
    bisecNoise( 7, hydraulics.ground, -1.0/256, 1.0/256 );
    hydraulics.initNeighs_6(false);
    //hydraulics.errodeDroples( 1000, 100, +0.1, 0.15, 0.9, {10, 10}, {100, 100} );
    for( int j=0; j<500; j++ ){
        int isz = 25;
        int ix0 = rand()%(hydraulics.n.x-isz);
        int iy0 = rand()%(hydraulics.n.y-isz);
        //printf("%i : %i %i \n", j, ix0, iy0);
        //                         n nStepMax, w,  disolve, sediment, ix0, iy0, ix1, iy1
        //hydraulics.errodeDroples( 400, 500, +0.1, 0.15, 0.9, ix0, iy0, ix0+isz, iy0+isz );
        hydraulics.errodeDroples( 400, 500, +0.1, 0.15, 0.9, {ix0, iy0}, {ix0+isz, iy0+isz} );
        //hydraulics.errodeDroples( 200, 500,   0.02, 0.1,    0.0,   0,0, hydraulics.nx, hydraulics.ny ); // fast
    }
    for(int i=0; i<hydraulics.ntot; i++){ ground[i] *= maxHeight; water[i] = ground[i]; }


    /*
    hydraulics.init_outflow( 500.0 );
    int n = 0;
    for(int ix=0; ix<ruler.na; ix++){ int idx=ix;                     hydraulics.contour2[n]=idx; water[idx]=ground[idx]; n++; };
    for(int ix=0; ix<ruler.na; ix++){ int idx=ix+ruler.ntot-ruler.na; hydraulics.contour2[n]=idx; water[idx]=ground[idx]; n++; };
    for(int iy=0; iy<ruler.nb; iy++){ int idx=iy*    ruler.na;        hydraulics.contour2[n]=idx; water[idx]=ground[idx]; n++; };
    for(int iy=0; iy<ruler.nb; iy++){ int idx=(iy+1)*ruler.na-1;      hydraulics.contour2[n]=idx; water[idx]=ground[idx]; n++; };
    hydraulics.nContour = n;
    hydraulics.isOutflow = true;
    //doDrain = 1;
    */

}

void LandCraftApp::addRoadStright( Vec2i p1, Vec2i p2 ){
    Road* road  = new Road();
    roads.push_back( road );
    roadBuilder.road = road;
    roadBuilder.path.clear();
    //std::list<int> lst;
    //lst.push_back(15454);
    //printf( " %i \n", *(lst.end()) );  //  Returns an iterator to the element following the last element of the container. This element acts as a placeholder; attempting to access it results in undefined behavior.;
    //printf( " %i \n", lst.back() );
    //RoadTile rt =   (RoadTile){(uint16_t)10,(uint16_t)15,(double)154.0};
    //rt.print();
    //roadBuilder.path.push_back( (RoadTile){10,15,0.0} );
    roadBuilder.path.push_back( (RoadTile){p1.x,p1.y,0.0} );
    //roadBuilder.path.push_back ( rt );
    //roadBuilder.path.end()->print();
    //roadBuilder.pushStright ( {55,38}     );
    roadBuilder.pushStright ( p2  );
    for( RoadTile& p : roadBuilder.path ){
        int i     = hydraulics.ip2i( (Vec2i){p.ia,p.ib} );
        //p.height = hydraulics.ground[i];
        //printf( "road (%i,%i)%i  (%i) %g \n", p.ia,p.ib, i, hydraulics.n.x, p.height );
    };
    //exit(0);
    roadBuilder.writeIt();
}

void LandCraftApp::addRiverPlot(int iRiver, float scFlow){

    River* thisRiver = hydraulics.rivers[iRiver];
    DataLine2D * riverH1 = new DataLine2D(thisRiver->path.size());
    riverProfile.lines.push_back( riverH1 );
    riverH1->linspan(0,thisRiver->path.size());
    int ii=0;
    if(scFlow>0){
        for( int i : thisRiver->path ){
            //hydraulics.i2ip(i);
            //riverH1->ys[ii] = hydraulics.water[i] * scFlow;
            riverH1->ys[ii] = thisRiver->flow[ii] * scFlow;
            //riverH1->ys[ii] = thisRiver->flow[ii];
            //printf( "road[%i] h flow %g \n", ii, riverH1->ys[ii] );
            ii++;
        };
    }else{
        for( int i : thisRiver->path ){
            riverH1->ys[ii] = hydraulics.ground[i];
            //printf( "road[%i] h height %g \n", ii, riverH1->ys[ii] );
            ii++;
        };
    }
    //exit(0);
    riverH1->clr = hash_Wang( (iRiver+15464)*5646 );

    //DataLine2D * lDrag = new DataLine2D(nsamp); mainWingL.lines.push_back( lDrag ); lDrag->linspan(-phiRange,phiRange); lDrag->clr = 0xFF0000ff;
    riverProfile.update();
    riverProfile.autoAxes(0.5,0.2);
    riverProfile.render();
}


void LandCraftApp::addRoadPlot(Road* road){
    int n = road->n;

    DataLine2D * roadH1 = new DataLine2D(n);
    roadProfile.lines.push_back( roadH1 );
    roadH1->linspan(0,n);

    DataLine2D * roadH2 = new DataLine2D(n);
    roadProfile.lines.push_back( roadH2 );
    roadH2->linspan(0,n);

    for( int ii=0; ii<n; ii++ ){
        RoadTile& p = road->path[ii];
        int i = hydraulics.ip2i({p.ia,p.ib});
        double h = hydraulics.ground[i];
        double w = hydraulics.water [i];
        roadH1->ys[ii] = w;
        roadH2->ys[ii] = h;
        //roadH2->ys[ii] = h + w;
        //printf( "road[%i] h %g w %g w+h %g p.h %g \n", ii, h,w,h+w, p.height );
    };
    //exit(0);
    roadH2->clr = 0xFF0000ff;
    roadH1->clr = 0xFFff0000;
    //DataLine2D * lDrag = new DataLine2D(nsamp); mainWingL.lines.push_back( lDrag ); lDrag->linspan(-phiRange,phiRange); lDrag->clr = 0xFF0000ff;
}


LandCraftApp::LandCraftApp( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    default_font_texture = makeTexture(  "common_resources/dejvu_sans_mono.bmp" );
    //default_font_texture = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    //itex = makeTexture(  "data/tank.bmp" );
    //itex = makeTexture(  "data/nehe.bmp" );
    printf( "default_font_texture :  %i \n", default_font_texture );




    fontTex     = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );
    GUI_fontTex = fontTex;


    /*
    FILE *ptr_myfile=0;
    char fname[] = "data/FUCK_YOU.txt";
    ptr_myfile=fopen( fname,"w");    printf("DEBUG A.1\n", fname );
    if (!ptr_myfile){ printf("Unable to open file! \n"); exit(-1); }
    fprintf( ptr_myfile, " FUCK YOU BASTARD !!!!! \n" );   printf("DEBUG A.2  >>%s<<\n", fname );
    fclose(ptr_myfile);     printf("DEBUG A.2\n", fname );

    saveBin( "test.bin", 454*sizeof(float), (char*)(new float[454]) );

    exit(0);
    */

    FILE* pFile = fopen("data/Technologies.txt", "r" );
    if (!pFile){ printf("Unable to open file!\n"); exit(0); }
    while( true ){
        Technology* tech = new Technology();
        if( tech->fromFile(pFile) ){
            tech->print();
            technologies[tech->name]=tech;
            //technologies.insert({tech->name,tech});
            //technologies[tech->name]->print();
        }else{
            delete tech;
            break;
        }
    }
    printf("===========\n");
    Factory factory;
    //for( auto it : technologies ){ printf("%s ::::::\n", it.first.c_str() ); it.second->print(); }
    //technologies[ std::string("Charcoal") ]->print();
    factory.setTechnology( technologies[ std::string("Charcoal") ] );
    factory.currentTenchnology->print();
    //factory.stored[ "Wood" ]->ammount = 500.0;
    factory.stored[ "Wood" ] = 500.0;
    printMap_d( factory.stored );
    factory.produce(10.0);
    printMap_d( factory.stored );
     printf("===========\n");
    //exit(0);



    printASCItable( 33, 127  );

    ruler.setSize(128,128);
    ruler.setStep(50);

    map_center = (Vec2d){ruler.na*0.75*ruler.step,ruler.nb*0.5*ruler.step};

    trace = new int [nTraceMax];

    //saveBin( "data/test-1.bin", 454*sizeof(float), (char*)(new float[454]) );

    //exit(-1);

    hydraulics.allocate( {ruler.na,ruler.nb} );
    hydraulics.initNeighs_6(false);
    ground = hydraulics.ground;
    water  = hydraulics.water;
    hydraulics.allocate_outflow();

    bool newMap = false;
    //bool newMap = true;
    if( newMap ){
        generateTerrain();
        //for(int i=0; i<ruler.ntot; i++){ ground[i] = randf(0.0,500.0); };
        saveBin( "data/ground.bin", sizeof(double)*hydraulics.ntot, (char*)ground );
        saveBin( "data/water.bin", sizeof(double)*hydraulics.ntot,  (char*)water  );
    }else{
        loadBin( "data/ground.bin", sizeof(double)*hydraulics.ntot, (char*)ground );
        loadBin( "data/water.bin", sizeof(double)*hydraulics.ntot,  (char*)water  );
    }


    // ======= Test Hydraulic relaxation

    /*
    Vec2i ip0 = {16,18};
    double g0=10.5;
    double w0=14.5;
    int i0=hydraulics.ip2i(ip0);
    hydraulics.ground[i0] = g0;
    hydraulics.water [i0] = w0;
    for(int ing=0; ing<hydraulics.nneigh; ing++){
        Vec2i ip = hydraulics.wrap_index( ip0 + hydraulics.neighs[ing] );
        int j    = hydraulics.ip2i(ip);
        double g = g0 + randf(-5.0,5.0);
        hydraulics.ground[j]=g;
        hydraulics.water [j]=g+randf(0,2.5);
    }
    hydraulics.relaxWatter( ip0 );
    exit(0);
    */

    double wmax = hydraulics.gatherRain( 100.0 ); printf("wmax %f \n",wmax ); // exit(0);
    terrainViewMode = 2;
    hydraulics.findAllRivers( 50.0 );

    std::sort(hydraulics.rivers.begin(), hydraulics.rivers.end(), [](River* a, River* b) { return a->path.size() > b->path.size(); } );

    //DropDownList* riverList = new DropDownList("Rivers",100.0,300.0,200.0,5);
    riverList = new DropDownList("Rivers",20.0,HEIGHT-100.0,200.0,5);
    gui.addPanel( riverList );
    //printf( "FOUND %i RIVERS \n", hydraulics.rivers.size() );
    for(int i=0; i<hydraulics.rivers.size(); i++ ){
        std::string name = "river " + std::to_string(i) +  " ("+ std::to_string( hydraulics.rivers[i]->path.size() )  +")";
        //printf( "%s \n", name.c_str() );
        riverList->addItem( name );
    };
    //exit(0);

    // plot rivers
    riverProfile.init();
    riverProfile.fontTex = fontTex;
    riverProfile.clrGrid = 0xFF404040;
    addRiverPlot(6,-1);
    addRiverPlot(6,0.1);
    //exit(0);
    //addRiverPlot(2);
    //addRiverPlot(3);


    /*
    Road* road  = new Road();
    roads.push_back( road );

    roadBuilder.road = road;

    //std::list<int> lst;
    //lst.push_back(15454);
    //printf( " %i \n", *(lst.end()) );  //  Returns an iterator to the element following the last element of the container. This element acts as a placeholder; attempting to access it results in undefined behavior.;

    //printf( " %i \n", lst.back() );

    //RoadTile rt =   (RoadTile){(uint16_t)10,(uint16_t)15,(double)154.0};
    //rt.print();
    roadBuilder.path.push_back( (RoadTile){10,15,0.0} );
    //roadBuilder.path.push_back ( rt );
    //roadBuilder.path.end()->print();
    roadBuilder.pushStright ( {55,38}     );
    for( RoadTile& p : roadBuilder.path ){
        int i     = hydraulics.ip2i( (Vec2i){p.ia,p.ib} );
        //p.height = hydraulics.ground[i];
        //printf( "road (%i,%i)%i  (%i) %g \n", p.ia,p.ib, i, hydraulics.n.x, p.height );
    };
    //exit(0);
    roadBuilder.writeIt();
    */

    addRoadStright( {10,15}, {55,38} );

    // plot road profile
    roadProfile.init();
    roadProfile.fontTex = fontTex;
    roadProfile.clrGrid = 0xFF404040;
    addRoadPlot(roads[0]);

    roadProfile.update();
    roadProfile.autoAxes(0.5,0.2);
    roadProfile.render();

    //exit(0);

    RoadVehicleType* vehicleType= new RoadVehicleType();
    vehicleTypes.push_back(vehicleType);

    RoadVehicle* vehicle = new RoadVehicle();
    vehicle->road = roads[0];
    vehicle->type = vehicleType;
    vehicles.push_back( vehicle );

    cmdPars.execFile( "data/comands.ini" );

    /*
    for(int i=0; i<100; i++){
        //Vec3d hs = {0.0,1.0,1.0};
        Vec3d hs = {randf(-1.0,1.0),randf(-1.0,1.0),randf(-1.0,1.0)};
        Vec2d p  = {0.2, 0.2};
        double h   = trinagleInterp( toBaricentric(p), hs );
        double hdx = trinagleInterp( toBaricentric({p.x+0.1,p.y }), hs );
        double hdy = trinagleInterp( toBaricentric({p.x,p.y+0.1 }), hs );
        Vec2d  dh  = trinagleDeriv ( hs );
        printf( " %i (%g,%g) (%g,%g) \n", i, (hdx-h)/0.1, (hdy-h)/0.1, dh.x, dh.y );
    }
    exit(0);
    */


    /*
    std::vector<int> feeders;
    River* river = new River();
    river->mouth=NULL;
    int n=hydraulics.trackRiver( hydraulics.sinks[0], 50.0, river->path, feeders );
    //delete river;
    hydraulics.rivers.push_back( river );
    printf( "River length %i \n", n );
    //exit(0);
    */

    //TiledView::init( 6, 6 );
    //tiles    = new int[ nxy ];
    //TiledView::renderAll( -10, -10, 10, 10 );


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


    // Fill random watter
    for(int i=0; i<hydraulics.ntot; i++){
        hydraulics.water[i] = hydraulics.ground[i] + randf(0,20.0);;
    }
    hydraulics.relaxWater( {36, 39});

}

void LandCraftApp::terrainColor( int i ){
    double maxDepth = 50.0;
    double depth;
    double w = water [i];
    double g = ground[i];
    switch( terrainViewMode ){
        case 1:
            depth = clamp( (w-g)/maxDepth, 0.0, 1.0 );
            //w-=g; if(w)depth=clamp( sqrt(w-g)*0.05, 0.0, 1.0 );
            g /= maxHeight;
            glColor3f( (1-depth)*g, (1-depth)*(1-g), depth );
            break;
        case 2:
            //w *= 0.1;
            //w *= log(w*5.0)*0.1;
            w=sqrt(w)*0.05;
            if( hydraulics.known[i] ){ glColor3f( 1, 0, 0 ); }else{ glColor3f( w, w, w ); }
            break;
    }
    //glColor3f( 1.0,1.0,1.0 );
}

void LandCraftApp::drawTerrain( Vec2i i0, Vec2i n, int NX ){
    Vec2f a,b,p;
    a.set( 1.0d, 0.0d           ); //a.mul(scale);
    b.set( 0.5d, 0.86602540378d ); //b.mul(scale);
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

void LandCraftApp::draw(){
    //long tTot = getCPUticks();
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );
    glShadeModel ( GL_SMOOTH );


    /*
    // ---- Hydro 1D
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
    */


    if(!file_relax_debug) file_relax_debug = fopen( "data/relax_debug.log", "w" );
    terrainViewMode = 1;
    double DEBUG_wsum_before=0.0;
    double DEBUG_E_before   =0.0;
    hydraulics.sumWater( DEBUG_wsum_before, DEBUG_E_before );
    hydraulics.relaxWater();
    for(int iy=0; iy<hydraulics.n.y; iy++){
    //for(int iy=0; iy<hydraulics.n.y; iy+=3){
        hydraulics.relaxWaterRasterX( iy, 0, hydraulics.n.x, 1.0 );
        hydraulics.relaxWaterRasterY( iy, 0, hydraulics.n.x, 1.0 );
    }
    double DEBUG_wsum_after=0.0;
    double DEBUG_E_after   =0.0;
    hydraulics.sumWater( DEBUG_wsum_after, DEBUG_E_after );
    printf( "DEBUG relaxWater[%i] dE %g | %g -> %g \n", frameCount, DEBUG_E_before-DEBUG_E_after,  DEBUG_E_before, DEBUG_E_after );
    if(fabs(DEBUG_wsum_before-DEBUG_wsum_after)>1e-6 ){
        printf( "DEBUG ERROR in hydraulics.relaxWater water before,after diff %g | %g %g \n", DEBUG_wsum_before-DEBUG_wsum_after, DEBUG_wsum_before, DEBUG_wsum_after );
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
        ground[i] = randf(0.0,1.0);
	}

    hydraulics.outflow_step();
	//if( doDrain==1 ){ hydraulics.outflow_step(); }else if (doDrain==-1){ hydraulics.inflow_step();  }
	//printf( " nContour %i \n", hydraulics.nContour );
	//if(frameCount>150)exit(0);

	glPushMatrix();
	glScalef(ruler.step,ruler.step,1.0);
	//Draw2D::drawTriaglePatch<cmapHeight>( {0,0}, {128,128}, ruler.na, ground, 0.0, maxHeight );
	//Draw2D::drawTriaglePatch<cmapHeight>( {0,0}, {128,128}, ruler.na, water, 0.0, maxHeight );
	drawTerrain( {0,0}, {128,128}, ruler.na );
	glPopMatrix();

    glColor3f( 1.0, 1.0, 1.0 );
	drawRoad( roads[0] );

	for( RoadVehicle* veh : vehicles ){
        //veh->moveStep( 1.0 );
        veh->move( 2.0 );
        //if( !veh->onWay ){ veh->idir*=-1; veh->onWay=true; }
        if( !veh->onWay ) veh->depart();
	}
	//if(frameCount>100) exit(0);

	for( RoadVehicle* veh : vehicles ){
        Vec2d p;
        //veh->ipath
        RoadTile& tile = veh->road->path[veh->ipath];
        ruler.nodePoint( {tile.ia,tile.ib}, p );
        Draw2D::drawCircle_d( p, 25, 16, false );
	}

    /*
	// test of hexIndex
    glBegin(GL_POINTS);
    for( int ix=0; ix<50; ix++ ){
        for( int iy=0; iy<50; iy++ ){
            double x = mouse_begin_x+ix*2.1;
            double y = mouse_begin_y+iy*2.1;
            int ihex = ruler.hexIndex( {x,y} );
            //glColor3f(  );
            Draw::color_of_hash(ihex+15454);
            glVertex3f( x, y, 100 );
        }
    }
    glEnd();
    */

    // draw sinks
    glColor3f( 0.0, 1.0, 1.0 );
    for(int isink : hydraulics.sinks ){
        Vec2d p;
        ruler.nodePoint( isink,p);
        Draw2D::drawPointCross_d( p, 100.0 );
    }

    // draw droplet trace
    glBegin(GL_LINE_STRIP);
    glColor3f( 1.0, 0.0, 1.0 );
    for(int ii=0; ii<nTrace; ii++){
        Vec2d p;
        ruler.nodePoint( trace[ii],p);
        glVertex3f( p.x, p.y, 100.0 );
    }
    glEnd();

    /*
    // draw river
    glBegin(GL_LINE_STRIP);
    glColor3f( 0.0, 1.0, 1.0 );
    for( i : river ){
        Vec2d p;
        ruler.nodePoint(i,p);
        glVertex3f( p.x, p.y, 100.0 );
    }
    glEnd();

    // draw feeders
    glColor3f( 1.0, 0.0, 1.0 );
    for(int i : feeders ){
        Vec2d p;
        ruler.nodePoint( i,p);
        Draw2D::drawPointCross_d( p, 30.0 );
    }
    */


    /*
    // draw all rivers
    for( River* river: hydraulics.rivers ){
        //printf( " path size: %i \n", river->path.size() );
        Draw::color_of_hash( river->path[0]+54877 );
        drawRiver( river );

        //glBegin(GL_LINE_STRIP);
        //for(int i : river->path ){
        //    Vec2d p;
        //    ruler.nodePoint(i,p);
        //    glVertex3f( p.x, p.y, 100.0 );
        //}
        //glEnd();

    }
    */

    //River* thisRiver = ;
    glLineWidth( 3 );
    glColor3f( 0.0,1.0,0.0 );
    drawRiver( hydraulics.rivers[riverList->iSelected] );
    glLineWidth( 1 );

    //float camMargin = ( camXmax - camXmin )*0.1;
    //float camMargin = 0;
    //TiledView::draw(  camXmin-camMargin, camYmin-camMargin, camXmax+camMargin, camYmax+camMargin  );
    //printf( " camRect  %f %f %f %f \n", camXmin-camMargin, camYmin-camMargin, camXmax+camMargin, camYmax+camMargin );
    //long tComp = getCPUticks();
    //update( );
    //tComp = getCPUticks() - tComp;

    //printf("===== frameCount %i \n", frameCount);

    Draw2D::drawPointCross({mouse_begin_x,mouse_begin_y}, 100.0);

};

/*
int LandCraftApp::tileToList( float x0, float y0, float x1, float y1 ){
	int ilist=glGenLists(1);
	glNewList( ilist, GL_COMPILE );
		terrain.renderRect( x0, y0, x1, y1, 31 );
		//glColor3f(0.9f,0.2f,0.2f); Draw2D::drawRectangle( x0+0.1, y0+0.1, x1-0.1, y1-0.1, false );
	glEndList();
	return ilist;
}
*/

void LandCraftApp::drawHUD(){

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);

    gui.draw();

    glPushMatrix();
    //glTranslatef(WIDTH-300,HEIGHT-200.0,200.0);
    roadProfile.view();
    glTranslatef(100,0.0,0.0);
    riverProfile.view();
    glPopMatrix();
    //glColor3f(1.0,1.0,0.0); mainWingLD.drawVline(phiL);
    //glColor3f(0.0,1.0,1.0); mainWingLD.drawVline(phiR);

}

void LandCraftApp::eventHandling ( const SDL_Event& event  ){
    //printf( "NBodyWorldApp::eventHandling() \n" );
    int ihex;
    double val;
    GUIAbstractPanel* activeGUIPanel = gui.onEvent(mouseX,mouseY,event);
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_n : generateTerrain(); break;
                case SDLK_r:  hydraulics.relaxWater(); break;
                case SDLK_e:
                    ihex = ruler.hexIndex({mouse_begin_x,mouse_begin_y});
                    hydraulics.relaxWater( ruler.i2ip(ihex) );
                    break;
                case SDLK_s :
                    saveBin( "data/ground.bin", sizeof(double)*hydraulics.ntot, (char*)ground );
                    saveBin( "data/water.bin", sizeof(double)*hydraulics.ntot,  (char*)water  );
                    break;
                case SDLK_l :
                    loadBin( "data/ground.bin", sizeof(double)*hydraulics.ntot, (char*)ground );
                    loadBin( "data/water.bin", sizeof(double)*hydraulics.ntot,  (char*)water  );
                    terrainViewMode = 1;
                    break;
                case SDLK_o :
                    for(int i=0; i<hydraulics.ntot; i++){ hydraulics.known[i]=false; }
                    ihex = ruler.hexIndex({mouse_begin_x,mouse_begin_y});
                    hydraulics.contour2[0] = ihex;
                    hydraulics.nContour++;
                    printf( "idrain  %i \n", ihex  );
                    water[ihex]          = ground[ihex];
                    hydraulics.isOutflow = true;
                    terrainViewMode = 1;
                    //doDrain = 1;
                    break;
                case SDLK_i :
                    for(int i=0; i<hydraulics.ntot; i++){ hydraulics.known[i]=false; }
                    ihex = ruler.hexIndex({mouse_begin_x,mouse_begin_y});
                    hydraulics.contour2[0] = ihex;
                    hydraulics.nContour++;
                    printf( "idrain  %i \n", ihex   );
                    water[ihex]          = ground[ihex] + 10.0;
                    hydraulics.isOutflow = false;
                    terrainViewMode = 1;
                    //doDrain = -1;
                    break;
                case SDLK_g :
                    hydraulics.gatherRain( 100.0 );
                    terrainViewMode = 2;
                    //val=0.0; ihex=0;
                    //for( int isink : hydraulics.sinks ){ double w=water[isink]; if(w>val){val=w; ihex=isink; } } // find largest river
                    //hydraulics.trackRiver( ihex, 50.0, river, feeders );
                    hydraulics.findAllRivers( 50.0 );
                    //for(int i=0; i<hydraulics.ntot; i++){ hydraulics.known[i]=false; }
                    //hydraulics.trackRiverRecursive( ihex, 50.0, NULL );
                    break;
                case SDLK_m :
                    terrainViewMode=(terrainViewMode%2)+1;
                    break;
                case SDLK_t :
                    ihex = ruler.hexIndex({mouse_begin_x,mouse_begin_y});
                    nTrace = hydraulics.traceDroplet( {ihex%hydraulics.n.x, ihex/hydraulics.n.x}, nTraceMax, trace );
                    break;
                //case SDLK_0:  formation_view_mode = 0;            printf( "view : default\n" ); break;
                //case SDLK_1:  formation_view_mode = VIEW_INJURY;  printf( "view : injury\n"  ); break;
                //case SDLK_2:  formation_view_mode = VIEW_STAMINA; printf( "view : stamina\n" ); break;
                //case SDLK_3:  formation_view_mode = VIEW_CHARGE;  printf( "view : charge\n"  ); break;
                //case SDLK_4:  formation_view_mode = VIEW_MORAL;   printf( "view : moral\n"   ); break;
                //case   SDLK_n: ifaction++; if(ifaction>=factions.size()) ifaction=0; currentFaction = factions[ifaction]; printf("ifaction %i\n",ifaction); break;
            }

            roadProfile.clear();
            addRoadPlot(roads[0]);
            roadProfile.update();
            roadProfile.autoAxes(0.5,0.2);
            roadProfile.render();

            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:

                    if( activeGUIPanel == 0 ){

                        ihex = ruler.hexIndex( {mouse_begin_x, mouse_begin_y} );
                        drawHeight = ground[ihex];
                        printf( "drawHeight %f \n", drawHeight );
                        //printf( " (%f,%f) %i  \n", mouse_begin_x, mouse_begin_y, ihex );
                        //hydraulics.water[ihex] = 1000.0;
                        //printf( "left button pressed !!!! " );
                        //if( currentFaction != NULL ) currentUnit = currentFaction->getUnitAt( { mouse_begin_x, mouse_begin_y } );
                        //bDrawing=true;

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
                    break;
            }
            break;

        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_RIGHT:
                    if( iBuildStartHex>=0 ){
                        int ihex = ruler.hexIndex( {mouse_begin_x, mouse_begin_y} );
                        for(Road* road : roads){ delete road; }
                        roads.clear();
                        addRoadStright( ruler.i2ip(iBuildStartHex), ruler.i2ip(ihex) );

                        roadProfile.clear();
                        addRoadPlot(roads[0]);
                        roadProfile.update();
                        roadProfile.autoAxes(0.5,0.2);
                        roadProfile.render();
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
        hydraulics.ground[ihex] = drawHeight;
    }
};

// ===================== MAIN

LandCraftApp * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
    SDL_DisplayMode dm;
    SDL_GetDesktopDisplayMode(0, &dm);
	int junk;
	//thisApp = new LandCraftApp( junk , dm.w-150, dm.h-100 );
	thisApp = new LandCraftApp( junk , 1000,600 );
	SDL_SetWindowPosition(thisApp->window, 100, 0 );
	thisApp->zoom  = 3000;
	thisApp->camX0 = 4000;
	thisApp->camY0 = 2000;
	thisApp->loop( 1000000 );
	return 0;
}
















