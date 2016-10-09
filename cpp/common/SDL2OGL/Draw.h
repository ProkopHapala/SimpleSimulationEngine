
#ifndef  Draw_h
#define  Draw_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"

namespace Draw{

	void setRGB ( uint32_t i );
	void setRGBA( uint32_t i );

	void color_of_hash( int i  );

    void billboardCam( );
    void drawText    ( const char * str, int itex, float sz, int istart, int iend );





}; // namespace Draw


#endif

