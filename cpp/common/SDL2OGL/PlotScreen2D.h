
#ifndef  PlotScreen2D_h
#define  PlotScreen2D_h

#include "Vec3.h"
#include "Vec2.h"
#include "Mat3.h"
#include "quaternion.h"
#include "geom2D.h"

#include "ScreenSDL2OGL.h"

class PlotScreen2D : public ScreenSDL2OGL{
	public:

    Rect2d defaultExtend;

	// ==== function declarations

	//virtual void camera      ();

	PlotScreen2D ( int& id, int WIDTH_, int HEIGHT_ );

};

#endif







