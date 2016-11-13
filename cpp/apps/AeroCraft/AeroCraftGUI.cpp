
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "SDL_utils.h"
#include "Draw.h"
#include "Draw3D.h"
//#include "drawDrawUtils.h"

#include "AeroCraftGUI.h" // THE HEADER

void AeroCraftGUI:: camera (){

    // third person camera attached to aero-craft
    float camDist = 5.0;
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
	float fov = VIEW_ZOOM_DEFAULT/zoom;
	glFrustum( -ASPECT_RATIO, ASPECT_RATIO, -1, 1, 1*fov, VIEW_DEPTH*fov );
	Mat3f camMat;
	Vec3f camPos;
	convert(world->myCraft->pos, camPos );
    qCamera.toMatrix( camMat );
    camMat.T();
	float glMat[16];
	Draw3D::toGLMatCam( { 0.0f, 0.0f, 0.0f}, camMat, glMat );
	glMultMatrixf( glMat );
	//glTranslatef ( -camPos.x, -camPos.y, -camPos.z );
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
    printf( "AeroCraftGUI draw\n" );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	renderSkyBox();
	world->update(); // ALL PHYSICS COMPUTATION DONE HERE

	/*
	printQuat( myCraft.qrot );  printf("qrot\n");
	printVec( myCraft.force );  printf("force\n");
	printVec( myCraft.torq );   printf("torq\n");
	printVec( myCraft.L );      printf("L\n");
	printVec( myCraft.omega );  printf("omega\n");
	printf("invI\n");
	printMat(  myCraft.invI );
	printf("rotMat\n");
	printMat(  myCraft.rotMat  );
	*/

	camera ();

	glEnable    (GL_LIGHTING);
	glShadeModel(GL_FLAT);
	glEnable(GL_DEPTH_TEST);
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


	//glDisable ( GL_LIGHTING );
	//glColor4f(1.0f,1.0f,1.0f,0.5f);
	//Draw3D::drawText( "AHOJ!\0", world->myCraft->pos, default_font_texture, 0.5, 0, 0 );
	//Draw3D::drawText( "AHOJ!\0", {0.0,0.0,0.0}, default_font_texture, 0.5, 0, 0 );

};

//void AeroCraftGUI:: drawHUD(){};

AeroCraftGUI:: AeroCraftGUI( int& id, int WIDTH_, int HEIGHT_ ) : ScreenSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {
    //default_font_texture = makeTexture( "common_resource/dejvu_sans_mono.bmp" );
    default_font_texture = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
};

