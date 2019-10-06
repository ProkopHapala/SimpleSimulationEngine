
#ifndef  AppSDL2OGL_3D_h
#define  AppSDL2OGL_3D_h

#include "Vec2.h"
#include "quaternion.h"

#include "Camera.h"
#include "cameraOGL.h"

#include "AppSDL2OGL.h"

//#include "Camera.h"
//#include "cameraOGL.h"

class AppSDL2OGL_3D : public AppSDL2OGL{ public:
	bool mouseSpinning = false;
	//Quat4f qCamera;
	//Mat3f  camMat;
	float  mouseRotSpeed   = 0.001;
	float  keyRotSpeed     = 0.01;
    float  cameraMoveSpeed = 0.2f;
    Quat4f qCamera;
	//Mat3d  camMat;
	//Vec3d  camPos;
	Camera cam;

	float camDist = 50.0;
	Vec2i spinning_start;

	bool perspective  = false;
	bool first_person = false;

// ============ function declarations

	//virtual void quit(       );
	//virtual void loop( int n );
	//virtual void inputHanding();
	virtual void keyStateHandling( const Uint8 *keys       );
	virtual void eventHandling   ( const SDL_Event& event  );
	virtual void mouseHandling   (                         );

	virtual void draw     ();
	virtual void camera   ();

	void camera_FPS       ( const Vec3d& pos, const Mat3d& rotMat );
	void camera_FwUp      ( const Vec3d& pos, const Vec3d& fw, const Vec3d& up, bool upDominant );
	void camera_FreeLook  ( const Vec3d& pos );
	void camera_OrthoInset( const Vec2d& p1, const Vec2d& p2, const Vec2d& zrange, const Vec3d& fw, const Vec3d& up, bool upDominant );

	//void orthoCamera      ( );
	//void perspectiveCamera( );

	void drawCrosshair( float sz );

	AppSDL2OGL_3D( int& id, int WIDTH_, int HEIGHT_ );

};

#endif
