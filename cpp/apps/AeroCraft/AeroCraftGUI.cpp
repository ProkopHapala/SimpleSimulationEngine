
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "SDL_utils.h"
#include "Draw.h"
#include "Draw2D.h"
#include "Draw3D.h"
//#include "drawDrawUtils.h"

#include "AeroCraftGUI.h" // THE HEADER

void AeroCraftGUI:: camera (){
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    float fov = VIEW_ZOOM_DEFAULT/zoom;
    glFrustum( -ASPECT_RATIO, ASPECT_RATIO, -1, 1, 1*fov, VIEW_DEPTH*fov );
    Mat3d camMat;
    Vec3f camPos;
    convert(world->myCraft->pos, camPos );
    float camDist = 10.0;
    if(first_person){
        // third person camera attached to aero-craft
        camMat.setT( world->myCraft->rotMat );
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

void AeroCraftGUI:: renderSkyBox( float x0, float y0, float z0){
	glDepthMask(0);
	glDisable (GL_LIGHTING);
	glDisable (GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);
	float skysz = VIEW_DEPTH*0.25;
	float R0=0.1,G0=0.1,B0=0.5;
	float R1=0.7,G1=0.8,B1=0.8;
	glBegin(GL_QUADS);
		glColor3f( R0, G0, B0 );  glVertex3f( -skysz+x0, skysz+y0, -skysz+z0 );
		glColor3f( R0, G0, B0 );  glVertex3f( -skysz+x0, skysz+y0, +skysz+z0 );
		glColor3f( R0, G0, B0 );  glVertex3f( +skysz+x0, skysz+y0, +skysz+z0 );
		glColor3f( R0, G0, B0 );  glVertex3f( +skysz+x0, skysz+y0, -skysz+z0 );

		glColor3f( R1, G1, B1 );  glVertex3f( -skysz+x0,     0+y0, -skysz+z0 );
		glColor3f( R1, G1, B1 );  glVertex3f( -skysz+x0,     0+y0, +skysz+z0 );
		glColor3f( R0, G0, B0 );  glVertex3f( -skysz+x0, skysz+y0, +skysz+z0 );
		glColor3f( R0, G0, B0 );  glVertex3f( -skysz+x0, skysz+y0, -skysz+z0 );

		glColor3f( R1, G1, B1 );  glVertex3f( +skysz+x0,     0+y0, -skysz+z0 );
		glColor3f( R1, G1, B1 );  glVertex3f( +skysz+x0,     0+y0, +skysz+z0 );
		glColor3f( R0, G0, B0 );  glVertex3f( +skysz+x0, skysz+y0, +skysz+z0 );
		glColor3f( R0, G0, B0 );  glVertex3f( +skysz+x0, skysz+y0, -skysz+z0 );

		glColor3f( R1, G1, B1 );  glVertex3f( -skysz+x0,     0+y0, -skysz+z0 );
		glColor3f( R1, G1, B1 );  glVertex3f( +skysz+x0,     0+y0, -skysz+z0 );
		glColor3f( R0, G0, B0 );  glVertex3f( +skysz+x0, skysz+y0, -skysz+z0 );
		glColor3f( R0, G0, B0 );  glVertex3f( -skysz+x0, skysz+y0, -skysz+z0 );

		glColor3f( R1, G1, B1 );  glVertex3f( -skysz+x0,     0+y0, +skysz+z0 );
		glColor3f( R1, G1, B1 );  glVertex3f( +skysz+x0,     0+y0, +skysz+z0 );
		glColor3f( R0, G0, B0 );  glVertex3f( +skysz+x0, skysz+y0, +skysz+z0 );
		glColor3f( R0, G0, B0 );  glVertex3f( -skysz+x0, skysz+y0, +skysz+z0 );
	glEnd();
	glDepthMask(1);
}


void AeroCraftGUI:: draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	if(world->staticTest) return;

	world->update(); // ALL PHYSICS COMPUTATION DONE HERE
	camera ();

	renderSkyBox(world->myCraft->pos.x, world->myCraft->pos.y-1000, world->myCraft->pos.z );
	glEnable(GL_DEPTH_TEST);

	glEnable    (GL_LIGHTING);
	glShadeModel(GL_FLAT);

	world->myCraft->render();

	//glDisable (GL_LIGHTING);
	glShadeModel(GL_SMOOTH);

	if ( world->buildings_shape >0 ) glCallList( world->buildings_shape );
	if ( world->terrain_shape >0)  {
		glCallList( world->terrain_shape   );
	} else {
	 	// terrain
		float groundsz = VIEW_DEPTH;
		glBegin(GL_QUADS);
			glColor3f( 0.3, 0.6, 0.1 );
			glNormal3f(0,1,0);
			glVertex3f( -groundsz, 0, -groundsz );
			glVertex3f( +groundsz, 0, -groundsz );
			glVertex3f( +groundsz, 0, +groundsz );
			glVertex3f( -groundsz, 0, +groundsz );
		glEnd();
	};

	//Draw3D::drawAxis( 1000 );

	/*
	glColor4f(1.0f,1.0f,1.0f,0.9f);
	char str[256];
	sprintf(str, "speed %3.3f\0",world->myCraft->vel.norm());
	Draw3D::drawText(str, world->myCraft->pos, fontTex, 0.2, 0, 0 );
	//Draw3D::drawText( "AHOJ!\0", world->myCraft->pos, fontTex, 0.2, 0, 0 );
	//Draw3D::drawText( "AHOJ!\0", {0.0,0.0,0.0}, fontTex, 0.5, 0, 0 );
	*/

	if(autoPilot){
        //printf("autoPiloting frame %i\n", frameCount);
        world->autoPilot1.control(world->dt); return;
    }

    if(mouseSteer){
        if (first_person){
            double dpitch=mouseY*0.005;
            double dyaw  =mouseX*0.002;
            double droll =0.2*dyaw;
            world->resetSteer( );
            world->myCraft->steerTo(droll, dpitch , dyaw);
        }else{
            Mat3d matCam;
            qCamera.toMatrix_T( matCam );
            Draw3D::drawMatInPos(matCam, world->myCraft->pos);
            world->steerToDir( matCam.c );
        }
    }

};

void AeroCraftGUI::drawHUD(){
	//panel .tryRender();  panel.draw();
	//mpanel.tryRender(); mpanel.draw();
	//if(focused) Draw2D::drawRectangle(focused->xmin,focused->ymin,focused->xmax,focused->ymax,false);

	if(world->staticTest){ drawStaticTest2D(); return; }

	char str[256];
	//sprintf(str, "speed %3.3f attitude %4.3f glideRatio %3.3f \0",world->myCraft->vel.norm(), world->myCraft->pos.y,   -sqrt(sq(world->myCraft->vel.x)+sq(world->myCraft->vel.z))/world->myCraft->vel.y );
	double vtot   = world->myCraft->vel.norm();
	//double thrust = world->myCraft->propelers[0].getThrust( vtot );
	double thrust = world->myCraft->totalThrust.norm();
	sprintf(str, "attitude %4.3f speed %3.3f vVert %3.3f tgAlfa %3.3f thrust %3.3f \0", world->myCraft->pos.y, vtot, world->myCraft->vel.y, world->myCraft->vel.y/vtot, thrust );
	glColor4f(1.0f,1.0f,1.0f,0.9f); Draw::drawText( str, fontTex, 10, 0,0 );

	if(first_person){ glColor4f(1.0f,1.0f,1.0f,0.9f); Draw2D::drawPointCross({mouseX+WIDTH*0.5,mouseY+HEIGHT*0.5},5.0); }

	/*
	int npol = 101;
	double phi0 =0.0;
	double dphi =3.14/(npol-1);
    for( int i=0; i<npol; i++){
        double phi = phi0 + dphi*i;
        double sa = sin( phi );
        double ca = cos( phi );
        double CD,CL;
        world->myCraft->panels[0].polarModel( ca, sa, CD, CL );
        printf("%i: %3.3f (%3.3f,%3.3f) (%3.3f,%3.3f) \n", i, phi, ca, sa, CD, CL);
    };
    exit(0);
    */

}


void AeroCraftGUI::drawStaticTest2D(){
    char str[256];
    /*
        Vec3d pos   = world->trjPos  [i];
        Vec3d vel   = world->trjVel  [i];
        Vec3d Force = world->trjForce[i];
        Vec3d Fw    = world->trjFw   [i];
        Vec3d Up    = world->trjUp   [i];
    */

    int n = world->ntrj;
     double t;
    Vec3d pos,vel,force,Up,Fw;
    // --- attitude
    double attitude;
    //glPushMatrix();
    //glScalef(1.0,0.5,1.0);
    glBegin(GL_LINE_STRIP);
    glColor3f(0.0f,0.0f,0.0f);
    for(int i=0; i<n; i++){
        t    = world->trjT    [i];
        pos = world->trjPos[i];
        glVertex3f(t,pos.y, 1.0);
    }
    glEnd();
    glColor4f(1.0f,1.0f,1.0f,0.9f);
    t        = world->trjT  [0]; attitude = world->trjPos[0].y;
    sprintf(str, "attitude=%4.3f m\0", attitude );
    Draw2D::drawText( str, {t,attitude},0.0,fontTex, 10, 0,0 );
    t        = world->trjT[n-1]; attitude = world->trjPos[n-1].y;
    sprintf(str, "attitude=%4.3f m\0", attitude );
    Draw2D::drawText( str, {t,attitude},0.0,fontTex, 10, 0,0 );
    //glPopMatrix();

    // --- pos.zy
    glBegin(GL_LINE_STRIP);
    glColor3f(0.0f,0.8f,0.8f);
    for(int i=0; i<n; i++){
        pos = world->trjPos[i];
        glVertex3f(pos.z*0.02+WIDTH*0.5,pos.y*0.02+HEIGHT*0.5, 1.0);
    }
    glEnd();

    // --- pos.xy
    glBegin(GL_LINE_STRIP);
    glColor3f(0.8f,0.8f,0.0f);
    for(int i=0; i<n; i++){
        pos = world->trjPos[i];
        glVertex3f(pos.x*0.02+WIDTH*0.5,pos.y*0.02+HEIGHT*0.5, 1.0);
    }
    glEnd();

    // --- pos.xz
    int iturn_old=0;
    int iturn    =0;
    double xold  =0.0;
    double xmin=1e+8,xmax=-1e+8,zmin=1e+8,zmax=-1e+8;
    glBegin(GL_LINE_STRIP);
    glColor3f(0.8f,0.8f,0.8f);
    for(int i=0; i<n; i++){
        pos = world->trjPos[i];
        vel = world->trjVel[i];
        if(pos.x<xmin) xmin=pos.x;
        if(pos.x>xmax) xmax=pos.x;
        if(pos.z<zmin) zmin=pos.z;
        if(pos.z>zmax) zmax=pos.z;
        if(vel.x*xold <0 ){ iturn_old=iturn; iturn=i; } xold=vel.x;
        glVertex3f(pos.x*0.02+WIDTH*0.5,pos.z*0.02+HEIGHT*0.5, 1.0);
    }
    glEnd();
    //sprintf(str, "xspan=%4.3f[m] zspan=%4.3f[m] T=%4.3f [s]\0", xmax-xmin, zmax-zmin, world->trjT[iturn]-world->trjT[iturn_old] );
    sprintf(str, "xspan=%4.3f[m] T=%4.3f [s]\0", xmax-xmin, world->trjT[iturn]-world->trjT[iturn_old] );
	glColor4f(1.0f,1.0f,1.0f,0.9f); Draw2D::drawText( str, {WIDTH*0.5,HEIGHT*0.5},0.0,fontTex, 10, 0,0 );
	pos=world->trjPos[iturn    ]; Draw2D::drawPointCross({pos.x*0.02+WIDTH*0.5,pos.z*0.02+HEIGHT*0.5},10);
	pos=world->trjPos[iturn_old]; Draw2D::drawPointCross({pos.x*0.02+WIDTH*0.5,pos.z*0.02+HEIGHT*0.5},10);


    // --- v.y
    float d0 = HEIGHT*0.5f;
    glBegin(GL_LINE_STRIP);
    glColor3f(1.0f,0.0f,1.0f);
    for(int i=0; i<n; i++){
        t     = world->trjT    [i];
        vel   = world->trjVel  [i];
        glVertex3f(t,vel.y+d0, 1.0);
    }
    glEnd();
    glBegin(GL_LINES);
        glVertex3f(0,d0, 1.0); glVertex3f(WIDTH,d0, 1.0);
    glEnd();


    // --- speed
    glBegin(GL_LINE_STRIP);
    glColor3f(0.0f,0.0f,1.0f);
    double speed;
    for(int i=0; i<n; i++){
        t     = world->trjT    [i];
        vel   = world->trjVel  [i];
        speed = vel.norm();
        glVertex3f(t,speed, 1.0);
    }
    glEnd();
    sprintf(str, "v=%4.3fm/s(%4.3fkm/h)\0", speed, speed*3.6 );
	glColor4f(1.0f,1.0f,1.0f,0.9f); Draw2D::drawText( str, {t,speed},0.0,fontTex, 10, 0,0 );


    //exit(0);
}

//void AeroCraftGUI:: drawHUD(){};

//AeroCraftGUI:: AeroCraftGUI( int& id, int WIDTH_, int HEIGHT_ ) : ScreenSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {
AeroCraftGUI:: AeroCraftGUI( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {
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
    //panel.nChars = 6;
};

