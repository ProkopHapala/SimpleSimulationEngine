
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "ScreenSDL2OGL_3D.h"

/*

	// TODO 

void ScreenSDL2OGL_3D::mouse_camera (float x, float y){
	if (mouse_spinning){
		trackball ( qCameraOld,
			(2.0f*mouse_begin_x-WIDTH)/WIDTH,  (HEIGHT-2.0f*mouse_begin_y)/HEIGHT, 
			(2.0f*x-WIDTH            )/WIDTH,  (HEIGHT-2.0f*y            )/HEIGHT 
		);
		add_quats ( qCameraOld, qCamera, qCamera );
		mouse_begin_x = x; mouse_begin_y = y; 
	}
}

void ScreenSDL2OGL_3D::getCameraDirections( ){
	float mat[4][4];
	glGetFloatv (GL_MODELVIEW_MATRIX, &mat[0][0]);
	camRight.set( mat[0][0], mat[1][0], mat[2][0] );
	camUp   .set( mat[0][1], mat[1][1], mat[2][1] );
	camDir  .set( mat[0][2], mat[1][2], mat[2][2] );
	camDir.mul( -1 ); // for some reason it is inverted
}

*/

void ScreenSDL2OGL_3D::camera(){
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
	glOrtho ( -zoom*ASPECT_RATIO, zoom*ASPECT_RATIO, -zoom, zoom, -VIEW_DEPTH, +VIEW_DEPTH );
	//glOrtho ( -zoom, zoom, -zoom, zoom, -VIEW_DEPTH, +VIEW_DEPTH );
	//glOrtho ( -zoom, zoom, -zoom*ASPECT_RATIO, zoom*ASPECT_RATIO, -VIEW_DEPTH, +VIEW_DEPTH );
	glMatrixMode (GL_MODELVIEW);
	//float camMatrix[4][4];
	//build_rotmatrix (camMatrix, qCamera );
	//glLoadMatrixf(&camMatrix[0][0]);
}

ScreenSDL2OGL_3D::ScreenSDL2OGL_3D( int& id, int WIDTH_, int HEIGHT_ ) : ScreenSDL2OGL( id, WIDTH_, HEIGHT_ ) { }
