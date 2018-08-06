
//#define SPEED_TEST

//#include <ctime>
//#include <iostream>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "Draw3D.h"
//#include "testUtils.h"

#include "Body.h"

//#include "SimplexRuler.h"
//#include "Ruler2DFast.h"
//#include "TerrainHydraulics.h"

#include "Terrain25D.h"
#include "Shooter.h"

#include "AeroSurf.h"
#include "AeroCraft.h"
#include "AeroCraftControl.h"
#include "AeroCraftWarrior.h"

#include "AeroControler1.h"

#include "DynamicControl.h"

//#include "AeroCraftControler.h"

//#include "FieldPatch.h"
#include "Solids.h"
//#include "AeroCraftWorld.h"
//#include "AeroCraftGUI.h"

#include "SDL_utils.h"
#include "Draw.h"
#include "Draw3D.h"
#include "AeroDraw.h"
#include "AeroCombatHelpers.h"
#include "AeroTest.h"

#include "GUI.h"
#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

#include "GUI.h"
#include "Plot2D.h"

//#include "AeroCraftDesign.h"

Vec3d opos=(Vec3d){0.0,0.0,0.0};
double traveledDistance=0;

char str[1<<16];

/*
TODO:

- Collision Detection
   - Projectiles-Target
   - Ground
- Use AeroCarftType vs AeroCarft
- Use projectile Burst
- Make Contorler which does not require copy of AeroCraft

*/

// ====================================
//      AeroCraftGUI
// ====================================

class AeroCraftGUI : public AppSDL2OGL_3D { public:
	//AeroCraftWorld * world;
	//Shooter * world;

	bool SimOn = true;

	Shooter            * world  = NULL;
    AeroCraftControler * pilot  = NULL;
    //AeroTester         * tester = NULL;

    //DynamicControl rollControl;
    //double roll;

    const Uint8 *scanKeys;
    Uint32 mouseButtons;

    int      fontTex;
    /*
    GUIPanel   panel;
    MultiPanel mpanel;
    GUITextInput txt;
    GUIAbstractPanel*  focused = NULL;
    */

    GUI gui;
    bool mouseSteer   = false;
    bool autoRetractAirelon    = true;
    bool autoRetractRudder     = true;
    bool autoRetractElevator   = false;

    bool useAutoPilot = false;

	AeroCraftWarrior * myCraft;
	AeroControler1 controler;

    bool staticTest = true;
    // - put to Spline manager ? ... make indepemented AeroCraft Test ?

	//AeroSurfaceDebugRecord leftWingRec,rightWingRec;
	std::vector<AeroSurfaceDebugRecord> dbgRects;
	//AeroSurfaceDebugRecord *leftWingRec,*rightWingRec;

	Plot2D mainWingLD;
	Plot2D mainWingPolar;
	int polarPlotKind=1;

	QuePlot2D historyPlot;
	QuePlot2D wingsTrj;
	QuePlot2D controlTrj;

	static const int nBaloons=100;
	Vec3d baloons[nBaloons];
	int gloBaloon;

	// ==== function declarations

	virtual void camera     ();
	//virtual void cameraHUD();
	virtual void draw   ();
	virtual void drawHUD();

    virtual void eventHandling   ( const SDL_Event& event  );
	virtual void keyStateHandling( const Uint8 *keys );
    virtual void mouseHandling   ( );


    void initPolarPlot();

    void logAeroCraftTrj();
    void controlAeroCraft();
    void drawAerocraftInset();
    //void drawAerocraftDebug();

    void drawControlSurfaceStateHUD();
    void drawPolarPlotHUD();
    void drawAeroCraftTrjHUD();
    void writeAeroCraftStateHUD();
    void drawCompassHUD( const Vec2d& pos2d, const Vec2d& dir2d );

	AeroCraftGUI( int& id, int WIDTH_, int HEIGHT_ );

};

void AeroCraftGUI::camera (){
    glMatrixMode(GL_PROJECTION); //glLoadIdentity();
    camera_FPS ( myCraft->pos, myCraft->rotMat );
    //if(first_person){ camera_FPS ( myCraft->pos, myCraft->rotMat ); }else{ camera_FreeLook( myCraft->pos ); };
    //camera_FwUp( myCraft->pos, myCraft->vel, {0.0,1.0,1.0}, false ); //camera_VelY();
    //camera_FwUp( myCraft->pos, myCraft->vel, {0.0,1.0,1.0}, true  ); //camera_YVel();
    //if(scanKeys[SDL_SCANCODE_LSHIFT]){ camera_FPS ( myCraft->pos, myCraft->rotMat );                    }
    //else                             { camera_FwUp( myCraft->pos, myCraft->vel, {0.0,1.0,1.0}, false ); };
    glMatrixMode(GL_MODELVIEW); glLoadIdentity();
}

void AeroCraftGUI::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);

	//printf( " upT %f t %f v %f | d  %f d %f \n", upTime*1e-3, world->time, myCraft->vel.norm(), traveledDistance, world->time*myCraft->vel.norm() );
	//Mat3d rot;     rot.setT(myCraft->rotMat);
    //Mat3d rot = myCraft->rotMat;

    //controler.bActive = false;
    controlAeroCraft();

	if(SimOn){
        //printf( "y %f vy %f x %f torq=(%f,%f,%f)     \n", rollControl.oy, rollControl.ovy, rollControl.x, myCraft->torq.x, myCraft->torq.y, myCraft->torq.z );
        world->update_world(); // ALL PHYSICS COMPUTATION DONE HERE
        //SimOn = false;
    }
	camera();

    logAeroCraftTrj();
    glColor3f( 1.0,0.0,0.0 );   wingsTrj.drawTrj3DPoints( {0,1,2}, 0.1 );
    glColor3f( 0.0,0.0,1.0 );   wingsTrj.drawTrj3DPoints( {3,4,5}, 0.1 );

	//renderSkyBox(myCraft->pos.x, myCraft->pos.y-1000, myCraft->pos.z, VIEW_DEPTH*0.25 );
	glEnable    (GL_DEPTH_TEST);
	glEnable    (GL_LIGHTING);
	glShadeModel(GL_FLAT);

	//renderAeroCraftVectors( *myCraft, 0.01, 1.0, true );

	glColor3f( 0.7,0.7,0.7 );
	renderAeroCraft( *myCraft, true, 1.0 );
    drawAerocraftInset();

	// ------ Aim & Shoot
	double aimDist = 400.0;
	double aimSize = 5.0;
    Mat3d rot = myCraft->rotMat;
	glColor3f(0.0,1.0,1.0); Draw3D::drawPointCross( myCraft->pos+rot.c*aimDist, aimSize );
	Draw3D::drawLine( myCraft->pos, myCraft->pos+rot.c*aimDist );
    //world->projLifetime=1.0;
    //if( scanKeys[SDL_SCANCODE_LCTRL] && ( frameCount%10 == 0 ) ){
    if( ( mouseButtons & SDL_BUTTON(SDL_BUTTON_LEFT) ) && ( frameCount%10 == 0 ) ){
        myCraft->gun_rot = rot.c;
        world->fireProjectile( myCraft, 700.0, 0 );
        printf("n projectiles: %i \n", world->projectiles.size() );
	}
	glColor3f(1.0,1.0,1.0);
	for(Projectile3D* prj : world->projectiles){
        //Draw3D::drawLine( prj->pos, prj->old_pos );
        Draw3D::drawPointCross( prj->pos, 1.0 );
        //printf( " %f, %f, %f\n", prj->pos.z, prj->pos.y, prj->pos.z);
	}
	for(int i=0; i<nBaloons; i++){
        Draw3D::drawShape( baloons[i], {10.0,0.0,0.0, 0.0,10.0,0.0, 0.0,0.0,10.0} ,gloBaloon);
        Draw3D::drawLine( {baloons[i].x,0.0,baloons[i].z}, baloons[i] );
    }

	//glDisable (GL_LIGHTING);
	glShadeModel(GL_SMOOTH);
	if(world->terrain) glCallList(world->terrain->shape);
	//Draw3D::drawAxis( 1000 );

	/*
	if(useAutoPilot){
        //printf("autoPiloting frame %i\n", frameCount);
        pilot->control(world->dt); return;
    }
    */

};

void AeroCraftGUI::drawHUD(){

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);

    gui.draw();

    drawControlSurfaceStateHUD();
    drawPolarPlotHUD();
    drawAeroCraftTrjHUD();
    writeAeroCraftStateHUD();

    drawCompassHUD( myCraft->pos.xz(), myCraft->rotMat.c.xz() );

    //glColor3f(1.0f,1.0f,1.0f,0.9f);
}


AeroCraftGUI:: AeroCraftGUI( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    printf( " === GUI \n" );
    fontTex     = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );
    GUI_fontTex = fontTex;
    gui.addPanel( new MultiPanel( "MultiPanel_1", 120,5,200,fontSizeDef*4, 4 ) )
        ->bgColor = 0x9090A0;
    SDL_StartTextInput();

    printf( " === world  \n" );

	world = new Shooter();
    world->dt = 0.6*(delay*1.0e-3)/world->perFrame;

    printf( " === Init Environment \n" );

    world->terrain =  prepareTerrain( 128, 2, 100.0, 50 );
    double baloonRande=10000.0;
    for(int i=0; i<nBaloons; i++){ baloons[i]=(Vec3d){randf(-baloonRande,baloonRande),randf(100.0,500.0),randf(-baloonRande,baloonRande) }; };
    gloBaloon=glGenLists(1);
	glNewList( gloBaloon, GL_COMPILE );
        glEnable    ( GL_LIGHTING );
        glShadeModel( GL_FLAT     );
        glColor3f( 0.8f, 0.8f, 0.8f );
        Draw3D::drawPolygons( Solids::Icosahedron_nfaces,        Solids::Icosahedron_ngons,        Solids::Icosahedron_faces,        Solids::Icosahedron_verts        );
	glEndList();

    printf( " === Init AeroCraft \n" );

    //char* fname = "data/AeroCraft1.ini";
    //char* fname = "data/AeroCraftStright1.ini";
    char* fname = "data/AeroCraftMainWing1.ini";
	myCraft     = new AeroCraftWarrior();
	myCraft  ->fromFile(fname);
	//myCraft     = new AeroCraft();   myCraft->fromFile("data/AeroCraft1.ini");

	//myCraft->rotMat.rotate( M_PI/2, {0.0,1.0,0.0} );

	//evalAeroFoceAtRotations( 8, {0.0,0.0,100.0}, {0.0,1.0,0.0}, *myCraft );

    myCraft->pos.y=200.0;
    myCraft->vel.set_mul( myCraft->rotMat.c, 100.0 );
    world->registrWarrior(myCraft);

    dbgRects.resize( myCraft->nPanels );
    for(int i=0; i<myCraft->nPanels; i++){
        myCraft->panels[i].dbgRec=&dbgRects[i];
    };

    printf( " === Init controlers \n" );

    pilot  = new AeroCraftControler();
    pilot->attach( myCraft );
    pilot->rudder  .setSymetricRange(0.2);
    pilot->elevator.setSymetricRange(0.2);

    world->controlers.push_back( &controler );
    controler.setup(false,true,myCraft,new AeroCraft());
    controler.craft_bak->fromFile(fname);
	controler.roll.setup (0,-0.3,+0.3,1.50);
    controler.pitch.setup(0,-0.3,+0.3,0.25);
    controler.yaw.setup  (0,-0.3,+0.3,0.25);

    printf( " === Init Polar \n" );

    initPolarPlot();

    printf( " === Init Trajectories  \n" );

    historyPlot.init( 100, 3 );
    historyPlot.lColors[0] = 0xFFff0000;
    historyPlot.lColors[1] = 0xFF007f00;
    historyPlot.lColors[2] = 0xFFff00ff;
    wingsTrj.init( 100, 6 );
    controlTrj.init(100,3);

    upTime=0;
    SDL_SetRelativeMouseMode( SDL_TRUE );
    VIEW_DEPTH = 10000;
    zoom = 30.0;

    printf(" === SETUP DONE !!! \n");
    //exit(0);
};

void AeroCraftGUI:: eventHandling   ( const SDL_Event& event  ){
    switch( event.type ){
        case SDL_KEYDOWN :
        switch( event.key.keysym.sym ){
            case SDLK_ESCAPE   : SDL_Quit(); exit(1); break;
            case SDLK_KP_PLUS  : zoom/=VIEW_ZOOM_STEP; printf("zoom: %f \n", zoom); break;
            case SDLK_KP_MINUS : zoom*=VIEW_ZOOM_STEP; printf("zoom: %f \n", zoom); break;
            //case SDLK_SPACE    : STOP = !STOP; printf( STOP ? " STOPED\n" : " UNSTOPED\n"); break;
            case SDLK_SPACE    : SimOn = !SimOn; printf( SimOn ? " STOPED\n" : " UNSTOPED\n"); break;

            /*
            case SDLK_LSHIFT   :
                //SDL_WarpMouseInWindow( window, WIDTH/2, HEIGHT/2);
                SDL_WarpMouseInWindow( window, WIDTH/2, HEIGHT*(1-pilot->elevator.getRelativeVal()) );
                mouseHandling();
            */
            case SDLK_u : useAutoPilot = !useAutoPilot; break;
            case SDLK_p : first_person = !first_person; break;
            case SDLK_m : mouseSteer   = !mouseSteer;   break;
            //case SDLK_r : myCraft->qrot.roll( 1.0 ); break;
            //case SDLK_c : pilot->resetSteer( );
        }; break;
        case SDL_QUIT: SDL_Quit(); exit(1); break;

        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_RIGHT:
                    mouseSteer = true;
                    break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_RIGHT:
                    mouseSteer = false;
                    pilot->resetSteer();
                    break;
            }
            break;

    }
};

void AeroCraftGUI:: keyStateHandling( const Uint8 *keys ){
    scanKeys = keys;

    bool RMB = mouseButtons&SDL_BUTTON(SDL_BUTTON_RIGHT);

	if ( keys[ SDL_SCANCODE_DOWN  ] ) { qCamera.pitch( -0.005 ); }
	if ( keys[ SDL_SCANCODE_UP    ] ) { qCamera.pitch(  0.005 ); }
	if ( keys[ SDL_SCANCODE_RIGHT ] ) { qCamera.yaw  (  0.005 ); }
	if ( keys[ SDL_SCANCODE_LEFT  ] ) { qCamera.yaw  ( -0.005 ); }

    if      ( keys[ SDL_SCANCODE_A ] ){ pilot->leftAirelon.inc();   pilot->rightAirelon.dec();   }
	else if ( keys[ SDL_SCANCODE_D ] ){ pilot->leftAirelon.dec();   pilot->rightAirelon.inc();   }


	else if ( keys[ SDL_SCANCODE_X ] ){ double val=0.5+(mouseX)/float(WIDTH);
                                        pilot->leftAirelon .towardRalative( 1.0-val );
                                        pilot->rightAirelon.towardRalative(     val );
                                      }
	else if ( autoRetractAirelon     ){ pilot->leftAirelon.relax(); pilot->rightAirelon.relax(); }

    if      ( keys[ SDL_SCANCODE_W ] ){ pilot->elevator.inc();   }
	else if ( keys[ SDL_SCANCODE_S ] ){ pilot->elevator.dec();   }
	//else if ( keys[ SDL_SCANCODE_Z ] ){ pilot->elevator.towardRelative( mouseY/WIDTH );   }
    else if ( keys[ SDL_SCANCODE_Z ]||keys[ SDL_SCANCODE_LSHIFT ] || RMB ){ pilot->elevator.towardRalative( 0.5+(mouseY)/float(HEIGHT) ); }                                                  //towardRalative
    else if ( autoRetractElevator    ){ pilot->elevator.relax(); }

    if      ( keys[ SDL_SCANCODE_E ] ){ pilot->rudder.inc();   }
	else if ( keys[ SDL_SCANCODE_Q ] ){ pilot->rudder.dec();   }
	else if ( keys[ SDL_SCANCODE_LSHIFT ] || RMB ){ pilot->rudder.towardRalative( 0.5+(mouseX)/float(WIDTH) ); }
    else if ( autoRetractRudder      ){ pilot->rudder.relax(); }


    //if( keys[SDL_SCANCODE_W]||keys[SDL_SCANCODE_S]||keys[SDL_SCANCODE_A]||keys[SDL_SCANCODE_D]||keys[SDL_SCANCODE_E]||keys[SDL_SCANCODE_Q] ){
    if( keys[SDL_SCANCODE_W]||keys[SDL_SCANCODE_S]||keys[SDL_SCANCODE_E]||keys[SDL_SCANCODE_Q] ){
        controler.bActive = false;
        Mat3d rot;
        //rot.setT(myCraft->rotMat);
        rot.set(myCraft->rotMat);
        controler.goalDir = rot.c;
        controler.goalUp  = rot.b;
    }else{
        controler.bActive = true;
    };


};

void AeroCraftGUI:: mouseHandling   ( ){
	// mouse Camera
	int mx,my;
	mouseButtons = SDL_GetMouseState(&mx,&my);
	int dmx = mx - WIDTH/2; 	int dmy = my - HEIGHT/2 ;
	mouseX = dmx;
	mouseY = -dmy;
	SDL_GetRelativeMouseState(&dmx,&dmy);
	qCamera.pitch( 0.005* dmy );
	qCamera.yaw  ( 0.005* dmx );
};


// =========== HELPER FUNCTIONS


void AeroCraftGUI::initPolarPlot(){
    // === mainWingLD
    mainWingLD.init();
    mainWingLD.fontTex = fontTex;
    mainWingLD.clrGrid = 0xFF404040;
    //mainWingLD.clrBg   = 0xFF408080;
    int nsamp = 100;
    double phiRange = M_PI*0.5;
    DataLine2D * lLift = new DataLine2D(nsamp); mainWingLD.lines.push_back( lLift ); lLift->linspan(-phiRange,phiRange); lLift->clr = 0xFFff0000;
    DataLine2D * lDrag = new DataLine2D(nsamp); mainWingLD.lines.push_back( lDrag ); lDrag->linspan(-phiRange,phiRange); lDrag->clr = 0xFF0000ff;
    // === mainWingPolar
    mainWingPolar.init();
    mainWingPolar.fontTex = fontTex;
    mainWingPolar.clrGrid = 0xFF404040;
    DataLine2D * LDpolar = new DataLine2D(nsamp); mainWingPolar.lines.push_back( LDpolar );
    for(int i=0; i<nsamp; i++){
        double phi = lLift->xs[i];
        double ca = cos(phi);
        double sa = sin(phi);
        double CD,CL;
        myCraft->leftAirelon->polarModel(ca,sa, CD, CL);
        lLift->ys[i]   = CL;
        lDrag->ys[i]   = CD;
        LDpolar->xs[i] = CD;
        LDpolar->ys[i] = CL;
    }
    mainWingPolar.update();
    mainWingPolar.autoAxes(0.1,0.5);
    mainWingPolar.render();

    mainWingLD.update();
    mainWingLD.autoAxes(0.5,0.2);
    mainWingLD.render();
}

void AeroCraftGUI::logAeroCraftTrj(){
	if( frameCount%10 == 0 ){
        historyPlot.next( world->time*0.5 );
        historyPlot.set_back( 0, myCraft->vel.norm() *0.02 );
        historyPlot.set_back( 1, myCraft->pos.y      *0.02  );
        historyPlot.set_back( 2, myCraft->vel.y      *0.02 );
    }
    if( frameCount%5 == 0 ){
        controlTrj.next( world->time*0.5 );
        controlTrj.set_back( 0, controler.roll.x   );
        controlTrj.set_back( 1, controler.roll.ovy );
        controlTrj.set_back( 2, controler.roll.oy  );
    }
    if( frameCount%10 == 0 ){
        wingsTrj.next( world->time*0.5 );
        Vec3d gp;
        gp = myCraft->leftAirelon ->getPos();
        wingsTrj.set_back( 0, gp.x );
        wingsTrj.set_back( 1, gp.y );
        wingsTrj.set_back( 2, gp.z );
        gp = myCraft->rightAirelon->getPos();
        wingsTrj.set_back( 3, gp.x );
        wingsTrj.set_back( 4, gp.y );
        wingsTrj.set_back( 5, gp.z );
        //Draw3D::drawPointCross( gp, 0.2 );
    }
    if(frameCount>0){ traveledDistance += (myCraft->pos-opos).norm(); }else{ traveledDistance=0; };
	opos=myCraft->pos;
}

void AeroCraftGUI::controlAeroCraft(){
    //Mat3d rot = myCraft->rotMat;
    controler.goalUp = (Vec3d)( cam.rot.a*mouseX + cam.rot.b*mouseY );
    controler.goalUp.normalize();

    int dmx,dmy;
    SDL_GetRelativeMouseState(&dmx,&dmy);
    //controler.goalDir = (Vec3d)cam.rot.c;
	controler.goalDir.add_mul((Vec3d)cam.rot.a, dmx* 0.003);
	controler.goalDir.add_mul((Vec3d)cam.rot.b, dmy*-0.003);
	controler.goalDir.normalize();
    glColor3f( 1.0,1.0,1.0); Draw3D::drawVecInPos( controler.goalDir*5, myCraft->pos );
}

void AeroCraftGUI::drawAerocraftInset(){
	// ---- Render Aerocraft Inset
	Mat3d rot = myCraft->rotMat;
    glMatrixMode( GL_PROJECTION ); glPushMatrix();
    camera_OrthoInset( {-5,-5}, {30.0,30.0}, {-100.0,100.0}, rot.a, {0.0,1.0,0.0}, true );
    glMatrixMode (GL_MODELVIEW); glPushMatrix();
        glLoadIdentity();
        renderAeroCraft(*myCraft, false, 1.0 );
        renderAeroCraftVectors( *myCraft, 0.01, 1.0, false );
        Draw3D::drawMatInPos( cam.rot, {0.0,0.0,0.0} );
    glMatrixMode( GL_MODELVIEW );  glPopMatrix();
    glMatrixMode( GL_PROJECTION ); glPopMatrix();
}

void AeroCraftGUI::drawControlSurfaceStateHUD(){
    // Control surface state
    glColor3f(1.0,1.0,1.0);
    Vec2d s0={WIDTH*0.5,110.0};
    double sza = 100, szb =10, val;
    val = (pilot->elevator   .getRelativeVal()*2.0-1.0)*sza;
    Draw2D::drawLine_d( {s0.x    ,s0.y+sza},{s0.x    ,s0.y-sza} );
    Draw2D::drawLine_d( {s0.x-szb,s0.y+val},{s0.x+szb,s0.y+val} );
    //s0.y-=szb;
    val =-(pilot->leftAirelon.getRelativeVal()*2.0-1.0)*sza;
    Draw2D::drawLine_d( {s0.x+sza,s0.y    },{s0.x-sza,s0.y    } );
    Draw2D::drawLine_d( {s0.x+val,s0.y-szb},{s0.x+val,s0.y+szb} );
    //s0.y+=szb*2;
    glColor3f(1.0,0.0,1.0);
    val = (pilot->rudder     .getRelativeVal()*2.0-1.0)*sza;
    //Draw2D::drawLine_d( {s0.x+sza,s0.y    },{s0.x-sza,s0.y    } );
    Draw2D::drawLine_d( {s0.x+val,s0.y-szb},{s0.x+val,s0.y+szb} );
    glColor3f(0.0,1.0,0.0); Draw2D::drawPointCross({s0.x+sza*mouseX*2.0/WIDTH,s0.y+sza*mouseY*2.0/HEIGHT},5.0);
}

void AeroCraftGUI::drawPolarPlotHUD(){
    float sc = WIDTH*0.04;
    if(polarPlotKind<=0) return;
    glPushMatrix();
    glTranslatef(WIDTH-300,HEIGHT-200.0,200.0);

    AeroSurfaceDebugRecord* leftWingRec = myCraft->leftAirelon->dbgRec;
    AeroSurfaceDebugRecord* rightWingRec = myCraft->rightAirelon->dbgRec;
    if(polarPlotKind==1){
        glScalef    (sc*4,sc,WIDTH);
        mainWingPolar.view();
        Draw2D::drawLine( {0.0,-1.0}, {0.0,1.0} );
        glColor3f(1.0,1.0,0.0); Draw2D::drawPointCross( {leftWingRec->CD, leftWingRec->CL},  0.05 );
        glColor3f(0.0,1.0,1.0); Draw2D::drawPointCross( {rightWingRec->CD,rightWingRec->CL}, 0.05 );
    }else if( polarPlotKind==2 ){
        glScalef    (sc,sc,WIDTH);
        double phiL = atan2( leftWingRec->sa,  leftWingRec->ca );
        double phiR = atan2( rightWingRec->sa, rightWingRec->ca );
        mainWingLD.view();
        glColor3f(1.0,1.0,0.0); mainWingLD.drawVline(phiL);
        glColor3f(0.0,1.0,1.0); mainWingLD.drawVline(phiR);
    }
    glPopMatrix();
}

void AeroCraftGUI::drawCompassHUD( const Vec2d& pos2d, const Vec2d& dir2d ){
    glPushMatrix();
    glTranslatef(WIDTH-50,50,1.0);
    glScalef(50.0,50.0,1.0);
    Draw2D::drawCircle( {0.0,0.0}, 1.0, 16, false );
    Draw2D::drawLine({-1.0,0.0},{1.0,0.0});
    Draw2D::drawLine({0.0,-1.0},{0.0,1.0});
    Draw2D::drawLine_d({0.0,0.0}, dir2d );
    //Draw::drawText( "x", fontTex, fontSizeDef, 0 );
    //Draw::drawText( "y", fontTex, fontSizeDef, 0 );
    glPopMatrix();
}

void AeroCraftGUI::drawAeroCraftTrjHUD(){
    float sc = WIDTH*0.04;
    glPushMatrix();
        glTranslatef(10.0,HEIGHT*0.5,200.0);
        glScalef    (sc,sc,WIDTH);
        controlTrj.draw( true, true );
        Draw::setRGBA(controlTrj.lColors[0]); Draw::drawText( "x\n",  fontTex, 0.1, 0 );  glTranslatef(0.0,0.2,0.0);
        Draw::setRGBA(controlTrj.lColors[1]); Draw::drawText( "vy\n", fontTex, 0.1, 0 );  glTranslatef(0.0,0.2,0.0);
        Draw::setRGBA(controlTrj.lColors[2]); Draw::drawText( "y\n",  fontTex, 0.1, 0 );  glTranslatef(0.0,0.2,0.0);
    glPopMatrix();
    /*
    glPushMatrix();
        glTranslatef(10.0,HEIGHT*0.5,200.0);
        glScalef    (sc,sc,WIDTH);
        historyPlot.draw( true, true );
        Draw::setRGBA(historyPlot.lColors[0]); Draw::drawText( "vel",        fontTex, 0.1, 0 );  glTranslatef(0.0,0.2,0.0);
        Draw::setRGBA(historyPlot.lColors[1]); Draw::drawText( "attidude\n", fontTex, 0.1, 0 );  glTranslatef(0.0,0.2,0.0);
        Draw::setRGBA(historyPlot.lColors[2]); Draw::drawText( "vVert\n",    fontTex, 0.1, 0 );  glTranslatef(0.0,0.2,0.0);
    glPopMatrix();
    */
}

void AeroCraftGUI::writeAeroCraftStateHUD(){
	Vec3d& pos    = myCraft->pos;
	Vec3d& vel    = myCraft->vel;
    double vtot   = vel.norm();
    double thrust = myCraft->totalThrust.norm();
    double power  = myCraft->getTotalPower();
    Vec3d vdir  = vel*(1.0/vtot);
    char* s=str;
    s+=sprintf(s,"Time[s]: %3.1f(%3.1f)\n", upTime*1e-3, world->time );
    s+=sprintf(s,"altitude[m  ] %g\n",            pos.y  );
    s+=sprintf(s,"speed   [m/s] %g [km/h] %g \n", vtot, vtot*3.6 );
    s+=sprintf(s,"climb   [m/s] %g           \n", vel.y  );
    s+=sprintf(s,"thrust  [N  ] %g [kg]   %g \n", thrust, thrust/9.81   );
    s+=sprintf(s,"thrust0 [N  ] %g P[kW]  %g\n", power/vtot, power*1e-3 );
    s+=sprintf(s,"pos  (%g,%g,%g) \n", pos .x, pos .y, pos .z  );
    s+=sprintf(s,"vel  (%g,%g,%g) \n", vel .x, vel .y, vel .z  );
    s+=sprintf(s,"vdir (%g,%g,%g) \n", vdir.x, vdir.y, vdir.z );
    //s+=sprintf(s,"vdir (%g,%g,%g) \n", vdir.x, vdir.y, vdir.z );
    glPushMatrix();
    glColor4f(1.0f,1.0f,1.0f,0.9f);
    glTranslatef( 10.0, HEIGHT-20, 0.0 );
    Draw::drawText( str, fontTex, fontSizeDef, {40,40} );
    glPopMatrix();
    //Draw::drawText( str, fontTex, 10, 0 );
}


// =========== MAIN

AeroCraftGUI * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	SDL_ShowCursor(SDL_DISABLE);
	SDL_DisplayMode dm;
    SDL_GetDesktopDisplayMode(0, &dm);
	int junk;
	thisApp = new AeroCraftGUI( junk , dm.w-150, dm.h-100 );
	SDL_SetWindowPosition(thisApp->window, 100, 0 );
	thisApp->loop( 1000000 );
	return 0;
}


