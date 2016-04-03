
#ifndef  AppSDL2OGL_3D_h
#define  AppSDL2OGL_3D_h

#include "Vec2.h"
#include "quaternion.h"


#include "AppSDL2OGL.h"

class AppSDL2OGL_3D : public AppSDL2OGL{
	public:
	bool mouseSpinning = false;
	//Quat4f qCamera;
	//Mat3f  camMat;
    Quat4d qCamera;
	Mat3d  camMat;

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

	//void orthoCamera      ( );
	//void perspectiveCamera( );

	AppSDL2OGL_3D( int& id, int WIDTH_, int HEIGHT_ );

};

#endif
