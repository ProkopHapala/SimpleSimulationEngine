
#ifndef  Draw_h
#define  Draw_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#define  fontSizeDef 7



namespace Draw{
    constexpr int     ncolors = 5;
    static uint32_t   colors_rainbow[ncolors] = { 0xFF000000, 0xFFFF0000, 0xFF00FF00, 0xFF00FFFF, 0xFFFFFFFF };

	void setRGB ( uint32_t i );
	void setRGBA( uint32_t i );

	void      colorScale( double d, int ncol=ncolors, const uint32_t * colors=&colors_rainbow[0] );
    uint32_t icolorScale( double d, int ncol=ncolors, const uint32_t * colors=&colors_rainbow[0] );

	void color_of_hash( int i  );

    void billboardCam( );
    void billboardCamProj( );
    void drawText ( const char * str, int itex, float sz, int iend         );
    void drawText ( const char * str, int itex, float sz, Vec2i block_size );

    //GLuint makeTexture( char * fname );
    //GLuint makeTexture( int nx, int ny, float * data );

    inline void printGLmat( float * glMat ){
    	printf( "  %f %f %f %f  \n", glMat[ 0], glMat[ 1], glMat[ 2], glMat[ 3] );
        printf( "  %f %f %f %f  \n", glMat[ 4], glMat[ 5], glMat[ 6], glMat[ 7] );
        printf( "  %f %f %f %f  \n", glMat[ 8], glMat[ 9], glMat[10], glMat[11] );
        printf( "  %f %f %f %f  \n", glMat[12], glMat[13], glMat[14], glMat[15] );
    }

    inline uint32_t float2RGBA(float f){
        uint8_t c = (uint8_t)(f*255.0);
        return 0xFF000000|(c<<16)|(c<<8)|c;
    }

    template<uint32_t _float2RGBA_(float f)>
    GLuint makeTexture( int nx, int ny, float * data ){
        GLuint itex=0;
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        glGenTextures  ( 1, &itex );
        glBindTexture  ( GL_TEXTURE_2D, itex );
        //glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, nx, ny, 0, GL_RED, GL_FLOAT, data);
        //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        int ntot = nx*ny;
        uint32_t * data_ = new uint32_t[ntot];
        for(int i=0;i<ntot;i++){ data_[i] = _float2RGBA_(data[i]); }
        //for(int i=0;i<ntot;i++){ data_[i] = _float2RGBA_(0.5); }
        //glTexImage2D( GL_TEXTURE_2D, 0, GL_RGB,  surf->w,  surf->h, 0, GL_BGR,       GL_UNSIGNED_BYTE, surf->pixels );
        //glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, nx, ny, 0, 0x8000,  GL_UNSIGNED_BYTE, data );
        glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, nx,  ny, 0, GL_RGBA,  GL_UNSIGNED_BYTE, data_ );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
        delete[] data_;
        return itex;
    };

}; // namespace Draw


#endif

