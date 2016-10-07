
#ifndef  SDL_utils_h
#define  SDL_utils_h

#include <stdlib.h>
#include <stdio.h>

GLuint makeTexture( char * fname ){
    GLuint itex;
    //SDL_Surface * surf = IMG_Load( fname );
    SDL_Surface * surf = SDL_LoadBMP( fname );
    if ( surf ){
        glGenTextures  ( 1, &itex );
        glBindTexture  ( GL_TEXTURE_2D, itex );
        glTexImage2D   ( GL_TEXTURE_2D, 0, 3, surf->w,  surf->h, 0, GL_BGR, GL_UNSIGNED_BYTE, surf->pixels );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
    }else{
        printf( "cannot load %s\n", fname  );
    }
    if ( surf ) SDL_FreeSurface( surf );
    //glGenTextures( 1, &itex );
    return itex;
};

#endif
