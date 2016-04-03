
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "Draw3D.h"

#include "AppSDL2OGL_3D.h" // THE HEADER

void AppSDL2OGL_3D::camera(){

    float camMatrix[16];
    qCamera.toMatrix_unitary( camMat );
    //first_person = true;
    //perspective  = true;
	if(first_person){
        glMatrixMode (GL_PROJECTION);
        glLoadIdentity();
        if( perspective ){
            //glFrustum( left, right, bottom, top, near, far );
            //glFrustum( -zoom*ASPECT_RATIO, zoom*ASPECT_RATIO, -zoom, zoom, 5, 20 );
            float fov = VIEW_ZOOM_DEFAULT/zoom;
            glFrustum( -ASPECT_RATIO, ASPECT_RATIO, -1, 1, 1*fov, VIEW_DEPTH*fov );
        }else{
            glOrtho  ( -zoom*ASPECT_RATIO, zoom*ASPECT_RATIO, -zoom, zoom, -VIEW_DEPTH, +VIEW_DEPTH );
        }
        Draw3D::toGLMatCam( {0.0f,0.0f,0.0f}, camMat, camMatrix );
        glMultMatrixf( camMatrix );
        glMatrixMode (GL_MODELVIEW);
        glLoadIdentity();
	}else{
        glMatrixMode( GL_PROJECTION );
        glLoadIdentity();
        if( perspective ){
            //glFrustum( left, right, bottom, top, near, far );
            //glFrustum( -zoom*ASPECT_RATIO, zoom*ASPECT_RATIO, -zoom, zoom, 5, 20 );
            float fov = VIEW_ZOOM_DEFAULT/zoom;
            glFrustum( -ASPECT_RATIO, ASPECT_RATIO, -1, 1, 1*fov, VIEW_DEPTH*fov );
        }else{
            glOrtho ( -zoom*ASPECT_RATIO, zoom*ASPECT_RATIO, -zoom, zoom, -VIEW_DEPTH, +VIEW_DEPTH );
        }
        glMatrixMode (GL_MODELVIEW);
        Draw3D::toGLMatCam( {0.0f,0.0f,0.0f}, camMat, camMatrix );
        glLoadMatrixf(camMatrix);
	}
	//glMatrixMode (GL_MODELVIEW);
}

void AppSDL2OGL_3D::draw   (){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glEnable    ( GL_LIGHTING );
	glShadeModel( GL_FLAT     );

	Draw3D::drawBox       ( -1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 1.0f, 0.8f, 0.8f, 0.8f );

	glShadeModel( GL_SMOOTH     );
	Draw3D::drawSphere_oct( 5, 1.0f, {3.0d,3.0d,3.0d} );

	glDisable ( GL_LIGHTING );
	Draw3D::drawAxis ( 3.0f );

};

void AppSDL2OGL_3D::eventHandling ( const SDL_Event& event  ){
    switch( event.type ){
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //mouseSpinning = true;
					//SDL_GetMouseState( &spinning_start.x, &spinning_start.y );
                break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
					//mouseSpinning = false;
					//int mx,my;
					//SDL_GetMouseState( &mx, &my );
					//qCamera.fromTrackballQ( spinning_start.x, spinning_start.y, mx, my );
                    break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}

void AppSDL2OGL_3D::keyStateHandling( const Uint8 *keys ){
    if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.dyaw  (  0.01 ); }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.dyaw  ( -0.01 ); }
	if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.dpitch(  0.01 ); }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.dpitch( -0.01 ); }
/*
    if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.dyaw  (  0.01 ); printf( "yaw   qCamera (%3.3f,%3.3f,%3.3f,%3.3f) \n", qCamera.x, qCamera.y, qCamera.z, qCamera.w ); }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.dyaw  ( -0.01 ); printf( "yaw   qCamera (%3.3f,%3.3f,%3.3f,%3.3f) \n", qCamera.x, qCamera.y, qCamera.z, qCamera.w ); }
	if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.dpitch(  0.01 ); printf( "pitch qCamera (%3.3f,%3.3f,%3.3f,%3.3f) \n", qCamera.x, qCamera.y, qCamera.z, qCamera.w ); }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.dpitch( -0.01 ); printf( "pitch qCamera (%3.3f,%3.3f,%3.3f,%3.3f) \n", qCamera.x, qCamera.y, qCamera.z, qCamera.w ); }
*/
/*
    if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.yaw  (  0.01 ); printf( "yaw   qCamera (%3.3f,%3.3f,%3.3f,%3.3f) \n", qCamera.x, qCamera.y, qCamera.z, qCamera.w ); }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.yaw  ( -0.01 ); printf( "yaw   qCamera (%3.3f,%3.3f,%3.3f,%3.3f) \n", qCamera.x, qCamera.y, qCamera.z, qCamera.w ); }
	if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.pitch(  0.01 ); printf( "pitch qCamera (%3.3f,%3.3f,%3.3f,%3.3f) \n", qCamera.x, qCamera.y, qCamera.z, qCamera.w ); }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.pitch( -0.01 ); printf( "pitch qCamera (%3.3f,%3.3f,%3.3f,%3.3f) \n", qCamera.x, qCamera.y, qCamera.z, qCamera.w ); }
*/

};

/*
void AppSDL2OGL_3D::mouse_camera (float x, float y){
	if (mouse_spinning){
		trackball ( qCameraOld,
			(2.0f*mouse_begin_x-WIDTH)/WIDTH,  (HEIGHT-2.0f*mouse_begin_y)/HEIGHT,
			(2.0f*x-WIDTH            )/WIDTH,  (HEIGHT-2.0f*y            )/HEIGHT
		);
		add_quats ( qCameraOld, qCamera, qCamera );
		mouse_begin_x = x; mouse_begin_y = y;
	}
}

void AppSDL2OGL_3D::getCameraDirections( ){
	float mat[4][4];
	glGetFloatv (GL_MODELVIEW_MATRIX, &mat[0][0]);
	camRight.set( mat[0][0], mat[1][0], mat[2][0] );
	camUp   .set( mat[0][1], mat[1][1], mat[2][1] );
	camDir  .set( mat[0][2], mat[1][2], mat[2][2] );
	camDir.mul( -1 ); // for some reason it is inverted
}
*/

void AppSDL2OGL_3D::mouseHandling( ){
    int mx,my;
    //SDL_GetMouseState( &mouseX, &mouseY );
    SDL_GetRelativeMouseState( &mx, &my);
    //printf( " %i %i \n", mx,my );
    Quat4d q; q.fromTrackball( 0, 0, mx*0.001, my*0.001 );
    //qCamera.qmul( q );
    qCamera.qmul_T( q );
}

AppSDL2OGL_3D::AppSDL2OGL_3D( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {
	qCamera.setOne();
	qCamera.toMatrix_unitary( camMat );
}

