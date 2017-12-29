

//#define SPEED_TEST

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

//#include "AeroCraftControler.h"

#include "FieldPatch.h"
#include "Solids.h"
//#include "AeroCraftWorld.h"
//#include "AeroCraftGUI.h"

#include "SDL_utils.h"
#include "Draw.h"
#include "Draw3D.h"
#include "AeroDraw.h"
#include "AeroTest.h"

#include "GUI.h"
#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

#include "GUI.h"
#include "Plot2D.h"

#include "AeroCraftDesign.h"

// ===============================
// ===== GLOBAL CONSTAMNTS
// ===============================

Shooter            * world  = NULL;
AeroCraftControler * pilot  = NULL;
AeroTester         * tester = NULL;

//AeroTester         tester     ;
//AeroCraftControler autoPilot1 ;
//Shooter            world      ;

//SDL_Surface *screen;
SDL_Event event;

int frameCount=0;
double tickSum=0;
int    stepSum=0;

bool loopEnd = false;
bool STOP    = false;

// ===============================
// ===== Free Functions
// ===============================

//void rotateTo( int pivot, Mat3d& rot, const Mat3d& rot0, const Vec3d& xhat, const Vec3d& yhat, double step ){

Terrain25D * prepareTerrain( int nsz, int nsub, double step, double hmax ){
    Terrain25D_bicubic * terrain = new Terrain25D_bicubic();
    terrain->ruler.setup( (Vec2d){nsz*0.5,nsz*0.5}*-step, (Vec2d){step,step} );
    terrain->allocate( {nsz,nsz} );
    terrain->makeRandom( 0.0, hmax );

    terrain->shape = glGenLists(1);
    glNewList( terrain->shape , GL_COMPILE );
    //int na=100,nb=100;

    int na = (terrain->ruler.n.a - 3)*nsub;
    int nb = (terrain->ruler.n.b - 3)*nsub;
    float da=terrain->ruler.step.a/float(nsub);
    float db=terrain->ruler.step.b/float(nsub);
    float x0=terrain->ruler.pos0.x;
    float y0=terrain->ruler.pos0.y;

    glEnable(GL_LIGHTING);
    glColor3f (0.5f,0.5f,0.5f);
    glNormal3f(0.0f,1.0f,0.0f);

    float * oldvals = new float[na*3];
    for(int ia=0; ia<na; ia++){
        glBegin(GL_TRIANGLE_STRIP);
        for(int ib=0; ib<nb; ib++){
            int i3 = 3*ib;
            Vec2d dv1,dv2;
            Vec2d p1; p1.set( (ia  )*da+x0, ib*db+y0 );
            Vec2d p2; p2.set( (ia+1)*da+x0, ib*db+y0 );
            float v1,v2;
            if( ia == 0 ){
                v1 = (float)terrain->eval( p1, dv1 );
            }else{
                v1 = oldvals[i3]; dv1.x=oldvals[i3+1]; dv1.y=oldvals[i3+2];
            }
            v2 = (float)terrain->eval( p2, dv2 );
            oldvals[i3] = v2; oldvals[i3+1] = dv2.x; oldvals[i3+2] = dv2.y;
            glNormal3f(-dv1.x,1.0,-dv1.y); glVertex3f( (float)p1.x,  v1, (float)p1.y );
            glNormal3f(-dv2.x,1.0,-dv2.y); glVertex3f( (float)p2.x,  v2, (float)p2.y );

            //glColor3f(v1,0.5,-v1); glVertex3f( (float)p1.x,  v1, (float)p1.y );
            //glColor3f(v2,0.5,-v2); glVertex3f( (float)p2.x,  v2, (float)p2.y );
            //printf( " %i (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f)\n", p1.x, p1.y, v1 ,  p2.x, p2.y, v2  );
        }
        glEnd();
    }

    // Normals
    /*
    glBegin(GL_LINES);
    for(int ia=0; ia<na; ia++){
        for(int ib=0; ib<nb; ib++){
            int i3 = 3*ib;
            Vec2d p,dv; p.set( ia*da+x0, ib*db+y0 );
            double v = (float)terrain->eval( p, dv );
            glVertex3f( (float)p.x,         v, (float)p.y );
            glVertex3f( (float)(p.x-dv.x),  v+1.0, (float)(p.y-dv.y) );
        }

    }
    glEnd();
    */
    glEndList();
    delete [] oldvals;
    return terrain;
}


// ====================================
//      AeroCraftGUI
// ====================================

class AeroCraftGUI : public AppSDL2OGL_3D { public:
	//AeroCraftWorld * world;
	//Shooter * world;

    int      fontTex;
    GUIPanel   panel;
    MultiPanel mpanel;
    GUITextInput txt;

    GUIAbstractPanel*  focused = NULL;

    bool mouseSteer   = false;
    bool autoRetractAirelon    = true;
    bool autoRetractRudder     = true;
    bool autoRetractElevator   = false;
    float  ElevatorRate = 0.01;
	float  AirelonRate  = 0.002;
	float  RudderRate   = 0.01;
    double AirelonRetractRate  = AirelonRate;
    double RudderRetractRate   = RudderRate;
    double ElevatorRetractRate = ElevatorRate;

    bool useAutoPilot = false;

	int perFrame = 10;
	//double dt = 0.001;

	AeroCraftWarrior * myCraft;
	//AeroCraft * myCraft;
	AeroCraft * myCraft_bak;

	//FieldPatch fieldPatch;
	//int buildings_shape = -1;
	//int terrain_shape   = -1;

    bool staticTest = true;
    // - put to Spline manager ? ... make indepemented AeroCraft Test ?

	//int fontTex_DEBUG;

	AeroSurfaceDebugRecord leftWingRec,rightWingRec;
	Plot2D mainWingLD;
	Plot2D mainWingPolar;
	int polarPlotKind=1;

	QuePlot2D historyPlot;
	QuePlot2D wingsTrj;

	static const int nBaloons=100.0;
	Vec3d baloons[nBaloons];
	int gloBaloon;

	// ==== function declarations

	//void reallocateTrj(int n);
	//void resetSteer( );
	//void steerToDir( const Vec3d& dir );

	//void renderSkyBox( float x0, float y0, float z0 );
	//void drawStaticTest2D();

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
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    //float fov = VIEW_ZOOM_DEFAULT/zoom;
    //glFrustum( -ASPECT_RATIO, ASPECT_RATIO, -1, 1, 1*fov, VIEW_DEPTH*fov );
    float camDist = 100.0;
    glFrustum( -ASPECT_RATIO, ASPECT_RATIO, -1, 1, camDist/zoom, VIEW_DEPTH );

    Mat3d camMat;
    Vec3f camPos;
    convert(myCraft->pos, camPos );

    if(first_person){
        // third person camera attached to aero-craft
        camMat.setT( myCraft->rotMat );
        //glTranslatef ( -camPos.x, -camPos.y, -camPos.z );
    }else{
        // third person camera attached to aero-craft
        qCamera.toMatrix( camMat );
        camMat.T();
	}
	float glMat[16];
	Draw3D::toGLMatCam( { 0.0f, 0.0f, 0.0f}, camMat, glMat );
	glMultMatrixf( glMat );
    glTranslatef ( -camPos.x+camMat.cx*camDist, -camPos.y+camMat.cy*camDist, -camPos.z+camMat.cz*camDist );
    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity();

}


void AeroCraftGUI::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);


	if(staticTest) return;

	world->update_world(); // ALL PHYSICS COMPUTATION DONE HERE
	camera ();

	if( frameCount%10 == 0 ){
        historyPlot.next( world->time*0.5 );
        historyPlot.set_back( 0, myCraft->vel.norm() *0.02 );
        historyPlot.set_back( 1, myCraft->pos.y      *0.02  );
        historyPlot.set_back( 2, myCraft->vel.y      *0.02 );
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


    Mat3d rot; rot.setT(myCraft->rotMat);
    glMatrixMode( GL_PROJECTION );
    glPushMatrix();
        glLoadIdentity();
        glOrtho( -ASPECT_RATIO*5.0, ASPECT_RATIO*30.0, -5.0,30.0,  -100.0, 100.0);
        Mat3d camMat;
        //camMat.setOne();
        //camMat.fromDirUp({0.0,1.0,0.0},myCraft->rotMat.a);
        // camera y-axis than along wing
        camMat.b.set( 0.0,1.0,0.0 );
        camMat.c.set( rot.a );
        camMat.c.makeOrtho( camMat.b );
        camMat.c.normalize();
        camMat.a.set_cross(camMat.b,camMat.c);
        // camera along wing than y-axis
        /*
        camMat.c.set( rot.a );
        camMat.b.set( 0.0,1.0,0.0 );
        camMat.b.makeOrtho( camMat.c );
        camMat.b.normalize();
        camMat.a.set_cross(camMat.b,camMat.c);
        */
        float glMat[16];

        Draw3D::toGLMatCam( { 0.0f, 0.0f, 0.0f}, camMat, glMat );
        //Draw3D::toGLMat( { 0.0f, 0.0f, 0.0f}, camMat, glMat );
        glMultMatrixf( glMat );
    glMatrixMode (GL_MODELVIEW);
    glPushMatrix();
        glLoadIdentity();
        renderAeroCraft(*myCraft, false);

        Draw3D::drawMatInPos( camMat, {0.0,0.0,0.0} );
    glMatrixMode( GL_MODELVIEW );
    glPopMatrix();
    glMatrixMode( GL_PROJECTION );
    glPopMatrix();
    //Draw3D::drawMatInPos( camMat, myCraft->pos );
    //Draw3D::drawMatInPos( myCraft->rotMat, myCraft->pos );
    Draw3D::drawMatInPos( rot, myCraft->pos );


	glColor3f(1.0,0.0,0.0);
	double fsc = 0.001;
	double vsc = 1.1;
	Draw3D::drawVecInPos( leftWingRec.force*fsc,  leftWingRec.gdpos+myCraft->pos );
	Draw3D::drawVecInPos( rightWingRec.force*fsc, rightWingRec.gdpos+myCraft->pos );
	//glColor3f(1.0,0.0,1.0);
	//Draw3D::drawVecInPos( (leftWingRec.force+rightWingRec.force)*fsc, rightWingRec.gdpos+myCraft->pos );
	Draw3D::drawVecInPos( (Vec3d){0.0,world->gravity,0.0} *(myCraft->mass*fsc), myCraft->pos );
	glColor3f(1.0,0.0,1.0); Draw3D::drawVecInPos( myCraft->force*fsc, myCraft->pos );
	glColor3f(0.0,0.0,1.0); Draw3D::drawVecInPos( myCraft->vel*vsc,   myCraft->pos );

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
            double dpitch=mouseY*0.005;
            double dyaw  =mouseX*0.002;
            double droll =0.2*dyaw;
            pilot->resetSteer( );
            myCraft->steerTo(droll, dpitch , dyaw);
        }else{
            Mat3d matCam;
            qCamera.toMatrix_T( matCam );
            Draw3D::drawMatInPos(matCam, myCraft->pos);
            pilot->steerToDir( matCam.c );
        }
    }

};

void AeroCraftGUI::drawHUD(){

    glDisable(GL_DEPTH_TEST);
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

	if(staticTest){ drawStaticTest2D( *tester, fontTex, WIDTH, HEIGHT ); return; }

	char str[256];
	//sprintf(str, "speed %3.3f attitude %4.3f glideRatio %3.3f \0",world->myCraft->vel.norm(), world->myCraft->pos.y,   -sqrt(sq(world->myCraft->vel.x)+sq(world->myCraft->vel.z))/world->myCraft->vel.y );
	double vtot   = myCraft->vel.norm();
	//double thrust = world->myCraft->propelers[0].getThrust( vtot );
	double thrust = myCraft->totalThrust.norm();
	sprintf(str, "attitude %4.3f speed %3.3f vVert %3.3f tgAlfa %3.3f thrust %3.3f \0", myCraft->pos.y, vtot, myCraft->vel.y, myCraft->vel.y/vtot, thrust );
	glColor4f(1.0f,1.0f,1.0f,0.9f); Draw::drawText( str, fontTex, 10, 0 );

	if(first_person){
        glColor4f(1.0f,1.0f,1.0f,0.9f); Draw2D::drawPointCross({mouseX+WIDTH*0.5,mouseY+HEIGHT*0.5},5.0);
    }

    float sc = WIDTH*0.04;
    glPushMatrix();
    if(polarPlotKind==1){
        glTranslatef(200.0,HEIGHT-150.0,200.0);
        glScalef    (sc*4,sc,WIDTH);
        mainWingPolar.view();
        Draw2D::drawLine( {0.0,-1.0}, {0.0,1.0} );
        glColor3f(1.0,1.0,0.0); Draw2D::drawPointCross( {leftWingRec.CD, leftWingRec.CL}, 0.05 );
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

    glPushMatrix();
        glTranslatef(10.0,HEIGHT*0.5,200.0);
        glScalef    (sc,sc,WIDTH);
        historyPlot.draw( true, true );
        Draw::setRGBA(historyPlot.lColors[0]); Draw::drawText( "vel",        fontTex, 0.1, 0 );  glTranslatef(0.0,0.2,0.0);
        Draw::setRGBA(historyPlot.lColors[1]); Draw::drawText( "attidude\n", fontTex, 0.1, 0 );  glTranslatef(0.0,0.2,0.0);
        Draw::setRGBA(historyPlot.lColors[2]); Draw::drawText( "vVert\n",    fontTex, 0.1, 0 );  glTranslatef(0.0,0.2,0.0);
    glPopMatrix();

    //glColor3f(1.0f,1.0f,1.0f,0.9f);
}

AeroCraftGUI:: AeroCraftGUI( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

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
    world->perFrame = 1;
    world->dt       = 0.005d;

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
    char* fname = "data/AeroCraftStright1.ini";
	myCraft_bak = new AeroCraft();          myCraft_bak->fromFile(fname);
	myCraft     = new AeroCraftWarrior();   myCraft    ->fromFile(fname);
	//myCraft     = new AeroCraft();   myCraft->fromFile("data/AeroCraft1.ini");
    myCraft->pos.y=200.0;
    myCraft->vel.set_mul( myCraft->rotMat.c, 100.0 );
    world->registrWarrior(myCraft);

    myCraft->leftAirelon->dbgRec  = &leftWingRec;
    myCraft->rightAirelon->dbgRec = &rightWingRec;

    printf( " === autoPilot1 \n" );

    pilot  = new AeroCraftControler();
    pilot->attach( myCraft );
    pilot->craft0=myCraft_bak;

    pilot->rudder  .setSymetricRange(0.2);
    pilot->elevator.setSymetricRange(0.2);

    printf( " === tester \n" );

    tester      = new AeroTester();
    tester->autoPilot = pilot;
    tester->craft     = myCraft;
    tester->gravityG  = world->gravity;
    //tester->reallocateTrj(int n);

    first_person = true;

    staticTest=false;
    if( staticTest ) tester->doStaticTesting( 500, 0.01, 300.0, 5.0 );

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

};

void AeroCraftGUI:: eventHandling   ( const SDL_Event& event  ){
    switch( event.type ){
        case SDL_KEYDOWN :
        switch( event.key.keysym.sym ){
            case SDLK_ESCAPE   : SDL_Quit(); exit(1); break;
            case SDLK_KP_PLUS  : zoom/=VIEW_ZOOM_STEP; printf("zoom: %f \n", zoom); break;
            case SDLK_KP_MINUS : zoom*=VIEW_ZOOM_STEP; printf("zoom: %f \n", zoom); break;
            case SDLK_SPACE    : STOP = !STOP; printf( STOP ? " STOPED\n" : " UNSTOPED\n"); break;

            case SDLK_u : useAutoPilot = !useAutoPilot;    break;
            case SDLK_p : first_person = !first_person; break;
            case SDLK_m : mouseSteer   = !mouseSteer;   break;
            case SDLK_c : pilot->resetSteer( );
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

	if ( keys[ SDL_SCANCODE_DOWN  ] ) { qCamera.pitch( -0.005 ); }
	if ( keys[ SDL_SCANCODE_UP    ] ) { qCamera.pitch(  0.005 ); }
	if ( keys[ SDL_SCANCODE_RIGHT ] ) { qCamera.yaw  (  0.005 ); }
	if ( keys[ SDL_SCANCODE_LEFT  ] ) { qCamera.yaw  ( -0.005 ); }

    if      ( keys[ SDL_SCANCODE_A ] ){ pilot->leftAirelon.inc();   pilot->rightAirelon.dec();   }
	else if ( keys[ SDL_SCANCODE_D ] ){ pilot->leftAirelon.dec();   pilot->rightAirelon.inc();   }
	else if ( autoRetractAirelon     ){ pilot->leftAirelon.relax(); pilot->rightAirelon.relax(); }

    if      ( keys[ SDL_SCANCODE_W ] ){ pilot->elevator.inc();   }
	else if ( keys[ SDL_SCANCODE_S ] ){ pilot->elevator.dec();   }
    else if ( autoRetractElevator    ){ pilot->elevator.relax(); }

    if      ( keys[ SDL_SCANCODE_E ] ){ pilot->rudder.inc();   }
	else if ( keys[ SDL_SCANCODE_Q ] ){ pilot->rudder.dec();   }
    else if ( autoRetractRudder      ){ pilot->rudder.relax(); }


};

void AeroCraftGUI:: mouseHandling   ( ){
	// mouse Camera
	int mx,my;
	SDL_GetMouseState(&mx,&my);
	int dmx = mx - WIDTH/2; 	int dmy = my - HEIGHT/2 ;
	mouseX = dmx;
	mouseY = -dmy;

	SDL_GetRelativeMouseState(&dmx,&dmy);
	qCamera.pitch( 0.005* dmy );
	qCamera.yaw  ( 0.005* dmx );
};

AeroCraftGUI * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	SDL_SetRelativeMouseMode( SDL_TRUE );
	SDL_DisplayMode dm;
    SDL_GetDesktopDisplayMode(0, &dm);
	int junk;
	thisApp = new AeroCraftGUI( junk , dm.w-150, dm.h-100 );
	SDL_SetWindowPosition(thisApp->window, 100, 0 );
	thisApp->loop( 1000000 );
	return 0;
}


