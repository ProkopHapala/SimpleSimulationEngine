

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



// ====================================
//      AeroCraftGUI
// ====================================

class AeroCraftGUI : public AppSDL2OGL_3D { public:
	//AeroCraftWorld * world;
	//Shooter * world;

	bool SimOn = true;

	Shooter            * world  = NULL;
    AeroCraftControler * pilot  = NULL;
    AeroTester         * tester = NULL;

    //DynamicControl rollControl;
    //double roll;

    const Uint8 *scanKeys;
    Uint32 mouseButtons;

    int      fontTex;
    GUIPanel   panel;
    MultiPanel mpanel;
    GUITextInput txt;

    GUIAbstractPanel*  focused = NULL;

    bool mouseSteer   = false;
    bool autoRetractAirelon    = true;
    bool autoRetractRudder     = true;
    bool autoRetractElevator   = false;

    bool useAutoPilot = false;

	AeroCraftWarrior * myCraft;
	AeroControler1 controler;

    bool staticTest = true;
    // - put to Spline manager ? ... make indepemented AeroCraft Test ?

	AeroSurfaceDebugRecord leftWingRec,rightWingRec;
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

	if(staticTest) return;

	if(frameCount>0){ traveledDistance += (myCraft->pos-opos).norm(); }else{ traveledDistance=0; };
	opos=myCraft->pos;

	//printf( " upT %f t %f v %f | d  %f d %f \n", upTime*1e-3, world->time, myCraft->vel.norm(), traveledDistance, world->time*myCraft->vel.norm() );

	Mat3d rot;     rot.    setT(myCraft->rotMat);
    //Mat3d camMatT; camMatT.setT(camMat);

    controler.goalUp = ( camMat.a*mouseX + camMat.b*mouseY ); controler.goalUp.normalize();

    int dmx,dmy;
    SDL_GetRelativeMouseState(&dmx,&dmy);
    //SDL_WarpMouse(0, 0);
    //SDL_WarpMouseInWindow( window, WIDTH*0.5, HEIGHT*0.5 );

    //printf( "dmYX %i %i \n", dmx, dmy );
	controler.goalDir.add_mul(camMat.a, dmx*0.01);
	controler.goalDir.add_mul(camMat.b, dmy*-0.01);
	controler.goalDir.normalize();

    //glColor3f( 1.0,1.0,1.0); Draw3D::drawVecInPos( controler.goalRoll*5, myCraft->pos );
    glColor3f( 1.0,1.0,1.0); Draw3D::drawVecInPos( controler.goalDir*5, myCraft->pos );

    //Draw3D::drawMatInPos( camMat, myCraft->pos );
    //Draw3D::drawMatInPos( camMatT, myCraft->pos );
    //Draw3D::drawMatInPos( rot   , myCraft->pos );



	if(SimOn){
        //printf( "y %f vy %f x %f torq=(%f,%f,%f)     \n", rollControl.oy, rollControl.ovy, rollControl.x, myCraft->torq.x, myCraft->torq.y, myCraft->torq.z );
        world->update_world(); // ALL PHYSICS COMPUTATION DONE HERE
        //SimOn = false;
    }
	camera();



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
    glColor3f( 1.0,0.0,0.0 );
    //wingsTrj.drawTrj3D( {0,1,2} );
    wingsTrj.drawTrj3DPoints( {0,1,2}, 0.1 );
    glColor3f( 0.0,0.0,1.0 );
    //wingsTrj.drawTrj3D( {3,4,5} );
    wingsTrj.drawTrj3DPoints( {3,4,5}, 0.1 );


	//renderSkyBox(myCraft->pos.x, myCraft->pos.y-1000, myCraft->pos.z, VIEW_DEPTH*0.25 );
	glEnable(GL_DEPTH_TEST);

	glEnable    (GL_LIGHTING);
	glShadeModel(GL_FLAT);

	//world->myCraft->render();
	renderAeroCraft( *myCraft, true );

	double aimDist = 400.0;
	double aimSize = 5.0;

	glColor3f(0.0,1.0,1.0); Draw3D::drawPointCross( myCraft->pos+rot.c*aimDist, aimSize );
	Draw3D::drawLine( myCraft->pos, myCraft->pos+rot.c*aimDist );

	// ---- Render Aerocraft Inset
    glMatrixMode( GL_PROJECTION ); glPushMatrix();
    camera_OrthoInset( {-5,-5}, {30.0,30.0}, {-100.0,100.0}, rot.a, {0.0,1.0,0.0}, true );
    glMatrixMode (GL_MODELVIEW); glPushMatrix();
        glLoadIdentity();
        renderAeroCraft(*myCraft, false);
        Draw3D::drawMatInPos( camMat, {0.0,0.0,0.0} );
    glMatrixMode( GL_MODELVIEW );  glPopMatrix();
    glMatrixMode( GL_PROJECTION ); glPopMatrix();
    //Draw3D::drawMatInPos( camMat, myCraft->pos );
    //Draw3D::drawMatInPos( myCraft->rotMat, myCraft->pos );
    Draw3D::drawMatInPos( rot, myCraft->pos );


    // ----  draw AreroCraftDebug
	glColor3f(1.0,0.0,0.0);
	double fsc = 0.001;
	double vsc = 1.1;
	Draw3D::drawVecInPos( leftWingRec .force*fsc, leftWingRec .gdpos+myCraft->pos );
	Draw3D::drawVecInPos( rightWingRec.force*fsc, rightWingRec.gdpos+myCraft->pos );
	glColor3f(0.0,0.0,1.0);
    Draw3D::drawVecInPos( leftWingRec .uair*vsc, leftWingRec .gdpos+myCraft->pos );
	Draw3D::drawVecInPos( rightWingRec.uair*vsc, rightWingRec.gdpos+myCraft->pos );
	//glColor3f(1.0,0.0,1.0);
	//Draw3D::drawVecInPos( (leftWingRec.force+rightWingRec.force)*fsc, rightWingRec.gdpos+myCraft->pos );
	Draw3D::drawVecInPos( (Vec3d){0.0,world->gravity,0.0} *(myCraft->mass*fsc), myCraft->pos );
	glColor3f(1.0,0.0,1.0); Draw3D::drawVecInPos( myCraft->force*fsc, myCraft->pos );
	glColor3f(0.0,0.0,1.0); Draw3D::drawVecInPos( myCraft->vel*vsc,   myCraft->pos );

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

	if(useAutoPilot){
        //printf("autoPiloting frame %i\n", frameCount);
        pilot->control(world->dt); return;
    }

    if(mouseSteer){
        if (first_person){
            //double dpitch=mouseY*0.005;
            //double dyaw  =mouseX*0.002;
            //double droll =0.2*dyaw;
            //pilot->resetSteer( );
            //myCraft->steerTo(droll, dpitch , dyaw);
        }else{
            Mat3d matCam;
            qCamera.toMatrix_T( matCam );
            Draw3D::drawMatInPos(matCam, myCraft->pos);
            //pilot->steerToDir( matCam.c );
        }
    }

};

void AeroCraftGUI::drawHUD(){

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
	//panel .tryRender();  panel.draw();
	//mpanel.tryRender(); mpanel.draw();
	//if(focused) Draw2D::drawRectangle(focused->xmin,focused->ymin,focused->xmax,focused->ymax,false);

	/*
	if(first_person){
        glPushMatrix();
            float sc = WIDTH*0.01;
            glTranslatef(sc*8,sc*8,100.0);
            glScalef(sc,sc,1.0);
            renderAeroCraft(*myCraft, false);
        glPopMatrix();
    }
    */

    //printf( " (%i,%i) (%i,%i) \n", WIDTH,HEIGHT, mouseX,mouseY  );
    //glColor3f(1.0,1.0,1.0); Draw2D::drawLine( {WIDTH*0.5,HEIGHT*0.5}, {WIDTH*0.5+mouseX*2,HEIGHT*0.5+mouseY*2} );



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


    //printf( "time sec %li \n", time(0)-appTime0 );
    //printf( "time sec %g sec %li \n", (clock()-appTime0)/clock_per_second, (time(0)-appSec0) );
    //printf( "time sec %g \n", SDL_GetTicks()*1e-3 );



    //glPopMatrix();

	//if(staticTest){ drawStaticTest2D( *tester, fontTex, WIDTH, HEIGHT ); return; }

	char str[256];
	//sprintf(str, "speed %3.3f attitude %4.3f glideRatio %3.3f \0",world->myCraft->vel.norm(), world->myCraft->pos.y,   -sqrt(sq(world->myCraft->vel.x)+sq(world->myCraft->vel.z))/world->myCraft->vel.y );
	double vtot   = myCraft->vel.norm();
	//double thrust = world->myCraft->propelers[0].getThrust( vtot );
	double thrust = myCraft->totalThrust.norm();
	sprintf(str, "T=%3.1fs(%3.1f) attitude %4.3f speed %3.3f vVert %3.3f tgAlfa %3.3f thrust %3.3f \0", upTime*1e-3, world->time, myCraft->pos.y, vtot, myCraft->vel.y, myCraft->vel.y/vtot, thrust );
	glColor4f(1.0f,1.0f,1.0f,0.9f); Draw::drawText( str, fontTex, 10, 0 );

	//if(first_person){
    //    glColor4f(1.0f,1.0f,1.0f,0.9f); Draw2D::drawPointCross({mouseX+WIDTH*0.5,mouseY+HEIGHT*0.5},5.0);
    //}

    float sc = WIDTH*0.04;
    glPushMatrix();
    if(polarPlotKind==1){
        glTranslatef(200.0,HEIGHT-150.0,200.0);
        glScalef    (sc*4,sc,WIDTH);
        mainWingPolar.view();
        Draw2D::drawLine( {0.0,-1.0}, {0.0,1.0} );
        glColor3f(1.0,1.0,0.0); Draw2D::drawPointCross( {leftWingRec.CD, leftWingRec.CL},  0.05 );
        glColor3f(0.0,1.0,1.0); Draw2D::drawPointCross( {rightWingRec.CD,rightWingRec.CL}, 0.05 );
    }else if( polarPlotKind==2 ){
        double phiL = atan2( leftWingRec.sa,  leftWingRec.ca );
        double phiR = atan2( rightWingRec.sa, rightWingRec.ca );
        glTranslatef(200.0,HEIGHT-150.0,200.0);
        glScalef    (sc,sc,WIDTH);
        mainWingLD.view();
        glColor3f(1.0,1.0,0.0); mainWingLD.drawVline(phiL);
        glColor3f(0.0,1.0,1.0); mainWingLD.drawVline(phiR);
    }
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
    glTranslatef( 100.0, 100.0, 0.0 );
    glScalef(100.0,100.0,0.0);
    Draw2D::drawVecInPos_d( {cos(roll),sin(roll)}, {0.0,0.0} );
    Draw2D::drawLine_d({-1.0, 0.0},{ 1.0, 0.0});
    Draw2D::drawLine_d({ 0.0,-1.0},{ 0.0, 1.0});
    glPopMatrix();
    */

    //glColor3f(1.0f,1.0f,1.0f,0.9f);
}

AeroCraftGUI:: AeroCraftGUI( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    //long tick_per_sec = CLOCKS_PER_SEC;
    //printf( "CLOCKS_PER_SEC %e \n", (double)tick_per_sec );
    //exit(0);

    VIEW_DEPTH = 10000;
    zoom = 30.0;

    printf( " === GUI \n" );

    //fontTex = makeTexture( "common_resources/dejvu_sans_mono.bmp" );
    //fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA.bmp" );
    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    //fontTex = makeTexture( "common_resources/dejvu_sans_mono_Alpha.bmp" );
    panel.init( 5,5,105,35,  fontTex );
    panel.caption   = "rotation [Rad]"; panel.vmin = -3.14159265359; panel.vmax = 3.14159265359;
    mpanel.initMulti( 120,5,200,120, fontTex , 4 );
    mpanel.caption="MultiPanel_1";
    txt.inputText = "insert number using =+-*/";
    SDL_StartTextInput ();

    printf( " === world  \n" );

	world = new Shooter();
    //world->perFrame = 1;
    //world->dt       = 0.005d;
    //world->dt       = 0.002d;
    world->dt = 0.6*(delay*1.0e-3)/world->perFrame;

    printf( " === Environment \n" );
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

    printf( " === aerocraft \n" );

    //char* fname = "data/AeroCraft1.ini";
    //char* fname = "data/AeroCraftStright1.ini";
    char* fname = "data/AeroCraftMainWing1.ini";
	myCraft     = new AeroCraftWarrior();   myCraft  ->fromFile(fname);

	//myCraft     = new AeroCraft();   myCraft->fromFile("data/AeroCraft1.ini");
    myCraft->pos.y=200.0;
    myCraft->vel.set_mul( myCraft->rotMat.c, 100.0 );
    world->registrWarrior(myCraft);

    myCraft->leftAirelon->dbgRec  = &leftWingRec;
    myCraft->rightAirelon->dbgRec = &rightWingRec;


    world->controlers.push_back( &controler );
    controler.craft     = myCraft;
    controler.craft_bak = new AeroCraft();  controler.craft_bak->fromFile(fname);
	//myCraft->controler   = &controler;
	controler.goalDir = myCraft->rotMat.c;
	controler.bUp   = false;
    controler.roll.y0 =  0; // M_PI*0.5; //sqrt(0.5);
    //controler.rollControl.dxdt_max =  1.0;
    //controler.rollControl.dydx     =  0.2; //5.5;
    //controler.rollControl.T        =  10.0;
    controler.roll.xmin = -0.3;
    controler.roll.xmax =  0.3;
    controler.roll.K    =  1.5;

    controler.bDir      = true;
    //controler.roll.y0   =  0;
    controler.pitch.xmin = -0.3;
    controler.pitch.xmax =  0.3;
    controler.pitch.K    =  0.25;
    //controler.roll.y0   =  0;
    controler.yaw.xmin = -0.3;
    controler.yaw.xmax =  0.3;
    controler.yaw.K    =  0.25;


    printf( " === autoPilot1 \n" );

    pilot  = new AeroCraftControler();
    pilot->attach( myCraft );
    //pilot->craft0=myCraft_bak;

    pilot->rudder  .setSymetricRange(0.2);
    pilot->elevator.setSymetricRange(0.2);

    printf( " === tester \n" );

    /*
    tester      = new AeroTester();
    tester->autoPilot = pilot;
    tester->craft     = myCraft;
    tester->gravityG  = world->gravity;
    //tester->reallocateTrj(int n);
    */

    printf( " === tester \n" );
    first_person = true;

    staticTest=false;
    //if( staticTest ) tester->doStaticTesting( 500, 0.01, 300.0, 5.0 );

    // Polar Plotting
    mainWingLD.init();
    mainWingLD.fontTex = fontTex;
    mainWingLD.clrGrid = 0xFF404040;
    //mainWingLD.clrBg   = 0xFF408080;
    int nsamp = 100;
    double phiRange = M_PI*0.5;
    DataLine2D * lLift = new DataLine2D(nsamp); mainWingLD.lines.push_back( lLift ); lLift->linspan(-phiRange,phiRange); lLift->clr = 0xFFff0000;
    DataLine2D * lDrag = new DataLine2D(nsamp); mainWingLD.lines.push_back( lDrag ); lDrag->linspan(-phiRange,phiRange); lDrag->clr = 0xFF0000ff;

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
    mainWingLD.update();
    mainWingLD.autoAxes(0.5,0.2);
    mainWingLD.render();

    mainWingPolar.update();
    mainWingPolar.autoAxes(0.1,0.5);
    mainWingPolar.render();

    historyPlot.init( 100, 3 );
    historyPlot.lColors[0] = 0xFFff0000;
    historyPlot.lColors[1] = 0xFF007f00;
    historyPlot.lColors[2] = 0xFFff00ff;

    wingsTrj.init( 100, 6 );

    controlTrj.init(100,3);

    upTime=0;

    SDL_SetRelativeMouseMode( SDL_TRUE );

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

            case SDLK_LSHIFT   :
                //SDL_WarpMouseInWindow( window, WIDTH/2, HEIGHT/2);
                SDL_WarpMouseInWindow( window, WIDTH/2, HEIGHT*(1-pilot->elevator.getRelativeVal()) );
                mouseHandling();

            case SDLK_u : useAutoPilot = !useAutoPilot; break;
            case SDLK_p : first_person = !first_person; break;
            case SDLK_m : mouseSteer   = !mouseSteer;   break;


            case SDLK_r : myCraft->qrot.roll( 1.0 ); break;

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
        Mat3d rot; rot.setT(myCraft->rotMat);
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
    /*
	controler.goalDir.add_mul(camMat.a, dmx*0.2);
	controler.goalDir.add_mul(camMat.b, dmy*0.2);
	controler.goalDir.normalize();
	*/
};

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


