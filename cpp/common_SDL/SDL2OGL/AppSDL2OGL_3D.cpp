
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
        //Draw3D::toGLMatCam( camPos, camMat, camMatrix ); // this does not work properly
        //camMat.a.mul(-1.0);
        //camMat.b.mul(-1.0);
        //camMat.c.mul(-1.0);
        Draw3D::toGLMatCam( {0.0d,0.0d,0.0d}, camMat, camMatrix );
        glMultMatrixf( camMatrix );
        glTranslatef ( -camPos.x, -camPos.y, -camPos.z );
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

        Draw3D::toGLMatCam( camPos*-1.0, camMat, camMatrix );
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
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_o:  perspective   = !perspective; break;
                case SDLK_p:  first_person  = !first_person ;   break;
            }
            break;
        /*
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
            */
    };
    AppSDL2OGL::eventHandling( event );
}

void AppSDL2OGL_3D::keyStateHandling( const Uint8 *keys ){

    if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.dyaw  (  keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.dyaw  ( -keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.dpitch(  keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.dpitch( -keyRotSpeed ); }

/*
    if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.yaw  (  keyRotSpeed ); qCamera.normalize(); }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.yaw  ( -keyRotSpeed ); qCamera.normalize(); }
	if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.pitch(  keyRotSpeed ); qCamera.normalize(); }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.pitch( -keyRotSpeed ); qCamera.normalize(); }
 */

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


void AppSDL2OGL_3D::mouseHandling( ){
    int mx,my;
    SDL_GetMouseState( &mouseX, &mouseY );   mouseY=HEIGHT-mouseY;
    mouse_begin_x = (2*mouseX-WIDTH )*zoom/HEIGHT;
    mouse_begin_y = (2*mouseY-HEIGHT)*zoom/HEIGHT;
    Uint32 buttons = SDL_GetRelativeMouseState( &mx, &my);
    //printf( " %i %i \n", mx,my );
    if ( buttons & SDL_BUTTON(SDL_BUTTON_RIGHT)) {
        Quat4d q; q.fromTrackball( 0, 0, -mx*mouseRotSpeed, my*mouseRotSpeed );
        qCamera.qmul_T( q );
    }
    //qCamera.qmul( q );
}

void AppSDL2OGL_3D::drawCrosshair( float sz ){
    glBegin( GL_LINES );
    float whalf = WIDTH *0.5;
    float hhalf = HEIGHT*0.5;
    glVertex3f( whalf-10,hhalf, 0 ); glVertex3f( whalf+10,hhalf, 0 );
    glVertex3f( whalf,hhalf-10, 0 ); glVertex3f( whalf,hhalf+10, 0 );
    glEnd();
}

AppSDL2OGL_3D::AppSDL2OGL_3D( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {
	qCamera.setOne();
	qCamera.toMatrix_unitary( camMat );
	camPos.set(0.0d);
}

