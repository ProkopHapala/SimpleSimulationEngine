
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
    void billboardCamProj( );
    void drawText    ( const char * str, int itex, float sz, int istart, int iend );


    inline void printGLmat( float * glMat ){
    	printf( "  %f %f %f %f  \n", glMat[ 0], glMat[ 1], glMat[ 2], glMat[ 3] );
        printf( "  %f %f %f %f  \n", glMat[ 4], glMat[ 5], glMat[ 6], glMat[ 7] );
        printf( "  %f %f %f %f  \n", glMat[ 8], glMat[ 9], glMat[10], glMat[11] );
        printf( "  %f %f %f %f  \n", glMat[12], glMat[13], glMat[14], glMat[15] );
    }


}; // namespace Draw


#endif

