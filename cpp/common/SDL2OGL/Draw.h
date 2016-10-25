
#ifndef  Draw_h
#define  Draw_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"

namespace Draw{

    constexpr int     ncolors = 5;
    static uint32_t   colors_rainbow[ncolors] = { 0xFF000000, 0xFFFF0000, 0xFF00FF00, 0xFF00FFFF, 0xFFFFFFFF };

	void setRGB ( uint32_t i );
	void setRGBA( uint32_t i );


	void colorScale( double d, int ncol=ncolors, const uint32_t * colors=&colors_rainbow[0] );

	void color_of_hash( int i  );

    void billboardCam( );
    void drawText    ( const char * str, int itex, float sz, int istart, int iend );


}; // namespace Draw


#endif

