
#ifndef  SDLplot_h
#define  SDLplot_h

#include "functions.h"

// ==== globals

class SDLplot{ public:

    int SCREEN_WIDTH  = 512;
    int SCREEN_HEIGHT = 512;

    int IY0 = SCREEN_HEIGHT/2;
    int IX0 = IY0;

    double ZOOM  = 30.0;
    double IZOOM = 1.0d/ZOOM;

    SDL_Renderer*	render			= NULL;
    SDL_Window*		window        	= NULL;
    SDL_Surface*	screenSurface 	= NULL;
    SDL_Surface*	tempSurf;
    SDL_Texture*	tempTex;

    SDL_Rect SrcR;
    SDL_Rect DestR;

    // ==== inline functions

    inline void setZoom( double zoom_ ){ ZOOM=zoom_; IZOOM=1.0d/ZOOM; }

    inline double i2x( int ix ){ return (ix - IX0)*IZOOM;  };
    inline double i2y( int iy ){ return (iy - IY0)*IZOOM;  };

    inline int x2i( double  x ){ return x*ZOOM + IX0;     };
    inline int y2i( double  y ){ return y*ZOOM + IY0;     };

    // ==== functions

    void init( char * title, int SCREEN_WIDTH_, int SCREEN_HEIGHT_ ){
        SCREEN_WIDTH  = SCREEN_WIDTH_;
        SCREEN_HEIGHT = SCREEN_HEIGHT_;
        int IY0       = SCREEN_HEIGHT/2;
        int IX0       = IY0;
    	window   = SDL_CreateWindow( title, SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN );
        render   = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC );

        SDL_UpdateWindowSurface( window );
        tempSurf        = SDL_CreateRGBSurface(0,SCREEN_WIDTH,SCREEN_HEIGHT,32,0,0,0,0 );
        SrcR = {0,0,tempSurf->w,tempSurf->h};
        DestR = SrcR;
    }

    void plotHline( double y ){
        SDL_RenderDrawLine( render, 0, y2i(-y), SCREEN_WIDTH, y2i(-y) );
    }

    void plotVline( double x ){
        SDL_RenderDrawLine( render, x2i(x), 0, x2i(x), SCREEN_HEIGHT );
    }

    void plotAxes( ){
        plotHline( 0 );
        plotVline( 0 );
    }

    void plotFunc( int n, double * xs, double * ys, double yscale ){
        int oix=x2i(  xs[0]        );
        int oiy=y2i( -ys[0]*yscale );
        for (int i=1; i<n; i++){
            int ix = x2i(  xs[i]        );
            int iy = y2i( -ys[i]*yscale );
            SDL_RenderDrawLine( render, oix, oiy, ix, iy );
            //printf( " %i  :  %f %f %i %i \n", i, xs[i], ys[i] , ix, iy  );
            oix = ix;
            oiy = iy;
        }
    }

    void setPixelsFunc2i( SDL_Surface* surface, int ix0, int iy0, int ix1, int iy1, Func2i func ){
        SDL_LockSurface( surface  );
        int nx = ix1-ix0;
        int ny = iy1-iy0;
        for (int iy=iy0; iy<iy1; iy++){
            Uint32 * pixel  = ((Uint32*)surface->pixels) + iy*surface->w + ix0;
            for (int ix=ix0; ix<ix1; ix++){
                int f = func( ix, iy );
                Uint32 color = f | 0xFF000000;
                *pixel = color;
                pixel++;
            }
        }
        SDL_UnlockSurface( surface  );
    }

    void setPixelsFunction( SDL_Surface* surface, int ix0, int iy0, int ix1, int iy1, Function2d func ){
        SDL_LockSurface( surface  );
        int nx = ix1-ix0;
        int ny = iy1-iy0;
        for (int iy=iy0; iy<iy1; iy++){
            Uint32 * pixel  = ((Uint32*)surface->pixels) + iy*surface->w + ix0;
            double y = i2y( iy );
            for (int ix=ix0; ix<ix1; ix++){
                double x = i2x( ix );
                double f = func( x, y );
                Uint8  fc    = ( (Uint8)(255.0*f) ) & 0xFF;
                Uint32 color = (fc<<16)|(fc<<8)|fc | 0xFF000000;
                *pixel = color;
                pixel++;
            }
        }
        SDL_UnlockSurface( surface  );
    }

    void setPixelsFunctionClamped( SDL_Surface* surface, int ix0, int iy0, int ix1, int iy1, Function2d func, float vmin, float vmax ){
        SDL_LockSurface( surface  );
        int nx = ix1-ix0;
        int ny = iy1-iy0;
        float fsc = 1.0f/(vmax-vmin);
        for (int iy=iy0; iy<iy1; iy++){
            Uint32 * pixel  = ((Uint32*)surface->pixels) + iy*surface->w + ix0;
            double y = i2y( iy );
            for (int ix=ix0; ix<ix1; ix++){
                double x = i2x( ix );
                double f = (func( x, y )-vmin)*fsc;
                if(f<0.0f) f=0.0;
                if(f>1.0f) f=1.0;
                Uint8  fc    = ( (Uint8)(255.0*f) ) & 0xFF;
                Uint32 color = (fc<<16)|(fc<<8)|fc | 0xFF000000;
                *pixel = color;
                pixel++;
            }
        }
        SDL_UnlockSurface( surface  );
    }

};

#endif
