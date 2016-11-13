
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "SDL_utils.h"
#include "Draw.h"
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
    float camDist = 5.0;
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

void AeroCraftGUI:: renderSkyBox(){
	glDepthMask(0);
	glDisable (GL_LIGHTING);
	glShadeModel(GL_SMOOTH);
	float skysz = VIEW_DEPTH/2;
	glBegin(GL_QUADS);
		glColor3f( 0.1, 0.1, 0.5 );  glVertex3f( -skysz, skysz, -skysz );
		glColor3f( 0.1, 0.1, 0.5 );  glVertex3f( -skysz, skysz, +skysz );
		glColor3f( 0.1, 0.1, 0.5 );  glVertex3f( +skysz, skysz, +skysz );
		glColor3f( 0.1, 0.1, 0.5 );  glVertex3f( +skysz, skysz, -skysz );

		glColor3f( 0.3, 0.5, 0.5 );  glVertex3f( -skysz,     0, -skysz );
		glColor3f( 0.3, 0.5, 0.5 );  glVertex3f( -skysz,     0, +skysz );
		glColor3f( 0.1, 0.1, 0.5 );  glVertex3f( -skysz, skysz, +skysz );
		glColor3f( 0.1, 0.1, 0.5 );  glVertex3f( -skysz, skysz, -skysz );

		glColor3f( 0.3, 0.5, 0.5 );  glVertex3f( +skysz,     0, -skysz );
		glColor3f( 0.3, 0.5, 0.5 );  glVertex3f( +skysz,     0, +skysz );
		glColor3f( 0.1, 0.1, 0.5 );  glVertex3f( +skysz, skysz, +skysz );
		glColor3f( 0.1, 0.1, 0.5 );  glVertex3f( +skysz, skysz, -skysz );

		glColor3f( 0.3, 0.5, 0.5 );  glVertex3f( -skysz,     0, -skysz );
		glColor3f( 0.3, 0.5, 0.5 );  glVertex3f( +skysz,     0, -skysz );
		glColor3f( 0.1, 0.1, 0.5 );  glVertex3f( +skysz, skysz, -skysz );
		glColor3f( 0.1, 0.1, 0.5 );  glVertex3f( -skysz, skysz, -skysz );

		glColor3f( 0.3, 0.5, 0.5 );  glVertex3f( -skysz,     0, +skysz );
		glColor3f( 0.3, 0.5, 0.5 );  glVertex3f( +skysz,     0, +skysz );
		glColor3f( 0.1, 0.1, 0.5 );  glVertex3f( +skysz, skysz, +skysz );
		glColor3f( 0.1, 0.1, 0.5 );  glVertex3f( -skysz, skysz, +skysz );
	glEnd();
	glDepthMask(1);
}


void AeroCraftGUI:: draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	renderSkyBox();

	glEnable(GL_DEPTH_TEST);

	world->update(); // ALL PHYSICS COMPUTATION DONE HERE

	camera ();

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

	Draw3D::drawAxis( 1000 );

	glColor4f(1.0f,1.0f,1.0f,0.9f);
	char str[256];
	sprintf(str, "speed %3.3f\0",world->myCraft->vel.norm());
	Draw3D::drawText(str, world->myCraft->pos, fontTex, 0.2, 0, 0 );
	//Draw3D::drawText( "AHOJ!\0", world->myCraft->pos, fontTex, 0.2, 0, 0 );
	//Draw3D::drawText( "AHOJ!\0", {0.0,0.0,0.0}, fontTex, 0.5, 0, 0 );

};

void AeroCraftGUI::drawHUD(){
	//panel .tryRender();  panel.draw();
	//mpanel.tryRender(); mpanel.draw();
	//if(focused) Draw2D::drawRectangle(focused->xmin,focused->ymin,focused->xmax,focused->ymax,false);

	char str[256];
	//sprintf(str, "speed %3.3f attitude %4.3f glideRatio %3.3f \0",world->myCraft->vel.norm(), world->myCraft->pos.y,   -sqrt(sq(world->myCraft->vel.x)+sq(world->myCraft->vel.z))/world->myCraft->vel.y );
	double vtot = world->myCraft->vel.norm();
	sprintf(str, "attitude %4.3f speed %3.3f vVert %3.3f tgAlfa %3.3f \0", world->myCraft->pos.y, vtot, world->myCraft->vel.y, world->myCraft->vel.y/vtot );
	Draw::drawText( str, fontTex, 10, 0,0 );
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

