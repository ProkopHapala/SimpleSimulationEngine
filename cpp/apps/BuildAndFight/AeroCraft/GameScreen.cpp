
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "drawMath.h"
//#include "drawDrawUtils.h"

#include "GameScreen.h" // THE HEADER

void GameScreen:: camera (){
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();

/*
	float zoomo = zoom*10;
	glOrtho      ( -zoomo*ASPECT_RATIO, zoomo*ASPECT_RATIO, -zoomo, zoomo, -VIEW_DEPTH, +VIEW_DEPTH );	
	Mat3d camMat;
	qmouse.toMatrix(camMat);
	float glMat[16];
	toGLMatCam( {0,0,0}, camMat, glMat);
	glMultMatrixf( glMat );
*/

/*
	glFrustum ( -zoom*ASPECT_RATIO, zoom*ASPECT_RATIO, -zoom, zoom, 0.5, VIEW_DEPTH);
	Quat4d qcam;
	qcam.setQmul( qmouse, myCraft.qrot );
	//qcam.set(qmouse);
	//qcam.set( myCraft.qrot );
	Mat3d camMat;
	qcam.toMatrix(camMat);
	float glMat[16];
	toGLMatCam( {0,0,0}, camMat, glMat);
	glMultMatrixf( glMat );
*/

	glTranslatef( (float)-world->myCraft->pos.x, (float)-world->myCraft->pos.y, (float)-world->myCraft->pos.z );

	glMatrixMode (GL_MODELVIEW);

}

void GameScreen:: renderSkyBox(){
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

void GameScreen:: draw(){
	glClearColor( 0.9, 0.9, 0.9, 0.0);                      
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );   

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT ); 

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

	drawAxis( 5 );

	//glFinish();
	//SDL_GL_SwapBuffers();
 
};

//void GameScreen:: drawHUD(){};

GameScreen:: GameScreen( int& id, int WIDTH_, int HEIGHT_ ) : ScreenSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

};

