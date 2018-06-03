
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "Draw3D.h"

#include "AppSDL2OGL_3D.h" // THE HEADER

//  camera_FPS( pos, rotMat ){

void AppSDL2OGL_3D::camera_FPS( const Vec3d& pos, const Mat3d& rotMat ){
    //glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    glFrustum( -ASPECT_RATIO, ASPECT_RATIO, -1, 1, camDist/zoom, VIEW_DEPTH );
    //Mat3d camMat;
    Vec3f camPos;
    convert( pos, camPos );
    camMat.setT( rotMat );
	float glMat[16];
	Draw3D::toGLMatCam( { 0.0f, 0.0f, 0.0f}, camMat, glMat );
	glMultMatrixf( glMat );
    glTranslatef ( -camPos.x+camMat.cx*camDist, -camPos.y+camMat.cy*camDist, -camPos.z+camMat.cz*camDist );
};

// camera( pos, dir, Up )
void AppSDL2OGL_3D::camera_FwUp( const Vec3d& pos, const Vec3d& fw, const Vec3d& up, bool upDominant ){
    //glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    glFrustum( -ASPECT_RATIO, ASPECT_RATIO, -1, 1, camDist/zoom, VIEW_DEPTH );
    //Mat3d camMat;
    Vec3f camPos;
    convert( pos, camPos );
    camMat.b = up;
    camMat.c = fw;
    if( upDominant ){
        camMat.b.normalize();
        camMat.c.makeOrtho( camMat.b );
        camMat.c.normalize();
    }else{
        camMat.c.normalize();
        camMat.b.makeOrtho( camMat.c );
        camMat.b.normalize();
    }
    camMat.a.set_cross(camMat.b,camMat.c);
	float glMat[16];
	Draw3D::toGLMatCam( { 0.0f, 0.0f, 0.0f}, camMat, glMat );
	glMultMatrixf( glMat );
    glTranslatef ( -camPos.x+camMat.cx*camDist, -camPos.y+camMat.cy*camDist, -camPos.z+camMat.cz*camDist );
};

void AppSDL2OGL_3D::camera_FreeLook( const Vec3d& pos ){
    //glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    glFrustum( -ASPECT_RATIO, ASPECT_RATIO, -1, 1, camDist/zoom, VIEW_DEPTH );
    //Mat3d camMat;
    Vec3f camPos;
    convert( pos, camPos );
    qCamera.toMatrix( camMat );
    camMat.T();
	float glMat[16];
	Draw3D::toGLMatCam( { 0.0f, 0.0f, 0.0f}, camMat, glMat );
	glMultMatrixf( glMat );
    glTranslatef ( -camPos.x+camMat.cx*camDist, -camPos.y+camMat.cy*camDist, -camPos.z+camMat.cz*camDist );
};

void AppSDL2OGL_3D::camera_OrthoInset( const Vec2d& p1, const Vec2d& p2, const Vec2d& zrange, const Vec3d& fw, const Vec3d& up, bool upDominant ){
    //glMatrixMode( GL_PROJECTION ); glPushMatrix();
    glLoadIdentity();
    //glOrtho( -ASPECT_RATIO*5.0, ASPECT_RATIO*30.0, -5.0, 30.0,  -100.0, 100.0);
    //printf( "--- %f %f  %f %f  %f %f \n", -ASPECT_RATIO*5.0, ASPECT_RATIO*30.0, -5.0, 30.0,  -100.0, 100.0  );
    //printf( "    %f %f  %f %f  %f %f \n", ASPECT_RATIO*p1.x, ASPECT_RATIO*p2.x, p1.y, p2.y,   zrange.a, zrange.b );
    glOrtho( ASPECT_RATIO*p1.x, ASPECT_RATIO*p2.x, p1.y, p2.y, zrange.a, zrange.b );
    //Mat3d camMat;
    camMat.b = up;
    camMat.c = fw;
    if( upDominant ){
        camMat.b.normalize();
        camMat.c.makeOrtho( camMat.b );
        camMat.c.normalize();
    }else{
        camMat.c.normalize();
        camMat.b.makeOrtho( camMat.c );
        camMat.b.normalize();
    }
    camMat.a.set_cross(camMat.b,camMat.c);
    float glMat[16];
    Draw3D::toGLMatCam( {0.0f,0.0f,0.0f}, camMat, glMat );
    //Draw3D::toGLMat( { 0.0f, 0.0f, 0.0f}, camMat, glMat );
    glMultMatrixf( glMat );
    //glMatrixMode (GL_MODELVIEW);
}

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
        //glMatrixMode (GL_MODELVIEW);
        //glTranslatef ( -camPos.x, -camPos.y, -camPos.z );
        //Draw3D::toGLMatCam( camPos*-1.0, camMat, camMatrix );
        //Mat3d camMatT; camMatT.setT(camMat);
        //Draw3D::toGLMat( {0.0,0.0,0.0}, camMatT, camMatrix );
        Draw3D::toGLMatCam( {0.0,0.0,0.0}, camMat, camMatrix );
        glTranslatef ( -camPos.x, -camPos.y, -camPos.z );
        //glLoadMatrixf(camMatrix);
        //glLMultMatrixf(camMatrix);
        glMultMatrixf(camMatrix);
        glMatrixMode (GL_MODELVIEW);
        glTranslatef ( -camPos.x, -camPos.y, -camPos.z );

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
	Draw3D::drawSphere_oct( 5, 1.0f, (Vec3f){3.0,3.0,3.0} );

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
    SDL_GetMouseState( &mouseX, &mouseY );   mouseY=HEIGHT-mouseY;
    mouse_begin_x = (2*mouseX-WIDTH )*zoom/HEIGHT;
    mouse_begin_y = (2*mouseY-HEIGHT)*zoom/HEIGHT;
    int mx,my; Uint32 buttons = SDL_GetRelativeMouseState( &mx, &my);
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
	GLbyte* s;
	// http://stackoverflow.com/questions/40444046/c-how-to-detect-graphics-card-model
	printf( "GL_VENDOR  : %s \n", glGetString(GL_VENDOR)  );
	printf( "GL_VERSION : %s \n", glGetString(GL_VERSION) );
}

