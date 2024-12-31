
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


    inline void setRGB( uint32_t i ){
        constexpr float inv255 = 1.0f/255.0f;
        //glColor3f( (i&0xFF)*inv255, ((i>>8)&0xFF)*inv255, ((i>>16)&0xFF)*inv255 );
        glColor3f( ((i>>16)&0xFF)*inv255, ((i>>8)&0xFF)*inv255, (i&0xFF)*inv255  );
    };

    inline void setRGB( uint32_t i, Vec3f& color ){
        constexpr float inv255 = 1.0f/255.0f;
        //glColor3f( (i&0xFF)*inv255, ((i>>8)&0xFF)*inv255, ((i>>16)&0xFF)*inv255 );
        color.set( ((i>>16)&0xFF)*inv255, ((i>>8)&0xFF)*inv255, (i&0xFF)*inv255  );
    };

    inline void setRGBA( uint32_t i ){
        constexpr float inv255 = 1.0f/255.0f;
        glColor4f( (i&0xFF)*inv255, ((i>>8)&0xFF)*inv255, ((i>>16)&0xFF)*inv255, ((i>>24)&0xFF)*inv255 );
    };

    /*
    void Draw::setRGBA( uint32_t i, Quat4f& color ){
        constexpr float inv255 = 1.0f/255.0f;
        color.set( (i&0xFF)*inv255, ((i>>8)&0xFF)*inv255, ((i>>16)&0xFF)*inv255, ((i>>24)&0xFF)*inv255 );
    };
    */

    inline void color_of_hash( int i ){
        //constexpr float inv255 = 1.0f/255.0f;
        int h = hash_Wang( i );
        Draw::setRGB( h );
        //glColor3f( (h&0xFF)*inv255, ((h>>8)&0xFF)*inv255, ((h>>16)&0xFF)*inv255 );
    };

    inline void color_of_hash( int i, Vec3f& color ){
        int h = hash_Wang( i );
        Draw::setRGB( h, color );
    };

    static int  fontTex = 0;
    static char tmp_str[1024];

    constexpr int     ncolors = 5;
    //static uint32_t   colors_rainbow[ncolors] = { 0xFF000000, 0xFFFF0000, 0xFF00FF00, 0xFF00FFFF, 0xFFFFFFFF };
    //static uint32_t   colors_rainbow[ncolors]   = { 0xFF000000, 0xFF0000FF, 0xFF00FF00, 0xFF00FFFF, 0xFFFFFFFF };
    //static uint32_t   colors_rainbow[ncolors]   = { 0xFF000000, 0xFFFF0000, 0xFF0000FF, 0xFF00FFFF, 0xFFFFFFFF };
    //static uint32_t   colors_rainbow[ncolors]   = { 0xFF000000, 0xFFFF0000, 0xFFFF00FF, 0xFF00FFFF, 0xFFFFFFFF };
    static uint32_t   colors_rainbow[ncolors]   = { 0xFF000000, 0xFFFF0000, 0xFF8000FF, 0xFF00FFFF, 0xFFFFFFFF };
    static uint32_t   colors_RWB    [ncolors]   = { 0xFFFF0000, 0xFFFFFF00, 0xFFFFFFFF, 0xFF00FFFF, 0xFF0000FF };

    //void setRGB ( uint32_t i );
    //void setRGBA( uint32_t i );
    //void setRGB ( uint32_t i, Vec3f& color );
    //void setRGBA( uint32_t i, Vec3f& color );

    void      colorScale( double d, int ncol=ncolors, const uint32_t * colors=&colors_rainbow[0] );
    uint32_t icolorScale( double d, int ncol=ncolors, const uint32_t * colors=&colors_rainbow[0] );

    void color_of_hash( int i  );
    void color_of_hash( int i, Vec3f& c );

    void billboardCam( );
    void billboardCamProj( float scale=200.0 );
    void drawText ( const char * str, int itex, float sz, int iend=0       );
    void drawText ( const char * str, int itex, float sz, Vec2i block_size );

    //GLuint makeTexture( char * fname );
    //GLuint makeTexture( int nx, int ny, float * data );

    inline int list(int ogl=0){ if(ogl)glDeleteLists(ogl,1); ogl=glGenLists(1); glNewList(ogl,GL_COMPILE); return ogl; };

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

