
#include <SDL2/SDL_opengl.h>

#include "Draw.h"  // THE HEADER

void Draw::setRGB( uint32_t i ){
	constexpr float inv255 = 1.0f/255.0f;
	//glColor3f( (i&0xFF)*inv255, ((i>>8)&0xFF)*inv255, ((i>>16)&0xFF)*inv255 );
	glColor3f( ((i>>16)&0xFF)*inv255, ((i>>8)&0xFF)*inv255, (i&0xFF)*inv255  );
};

void Draw::setRGBA( uint32_t i ){
	constexpr float inv255 = 1.0f/255.0f;
	glColor4f( (i&0xFF)*inv255, ((i>>8)&0xFF)*inv255, ((i>>16)&0xFF)*inv255, ((i>>24)&0xFF)*inv255 );
};

void Draw::color_of_hash( int i ){
	constexpr float inv255 = 1.0f/255.0f;
	int h = hash_Wang( i );
	glColor3f( (h&0xFF)*inv255, ((h>>8)&0xFF)*inv255, ((h>>16)&0xFF)*inv255 );
};

void Draw::colorScale( double d, int ncol, const uint32_t * colors ){
    constexpr float inv255 = 1.0f/255.0f;
    d*=(ncol-1);
    int icol = (int)d;
    d-=icol; double md = 1-d;
    uint32_t clr1=colors[icol];
    uint32_t clr2=colors[icol];
    glColor3f(
        ( d*(  colors[icol]    &0xFF) + md*( colors[icol]     &0xFF ))*inv255,
        ( d*((colors[icol]>>8 )&0xFF) + md*((colors[icol]>>8 )&0xFF ))*inv255,
        ( d*((colors[icol]>>16)&0xFF) + md*((colors[icol]>>16)&0xFF ))*inv255
    );
};

/*
void Draw::setColorInt32( uint32_t clr ) {
    constexpr float i255 = 1/255.0f;
    uint8_t b = ( ( clr       ) & 0xFF );
    uint8_t g = ( ( clr >> 8  ) & 0xFF );
    uint8_t r = ( ( clr >> 16 ) & 0xFF );
    uint8_t a = (   clr >> 24          );
    glColor4f( i255*r, i255*g, i255*b, i255*a );
    //printf( " r %i g %i b %i a %i     %f %f %f %f  \n", r, g, b, a,  i255*r, i255*g, i255*b, i255*a   );
};
*/

void Draw::billboardCam( ){
    float glMat[16];
    //glMatrixMode(GL_MODELVIEW);
    glGetFloatv(GL_MODELVIEW_MATRIX , glMat);
    glMat[0 ] = 1;   glMat[1 ] = 0;   glMat[2 ] = 0;
    glMat[4 ] = 0;   glMat[5 ] = 1;   glMat[6 ] = 0;
    glMat[8 ] = 0;   glMat[9 ] = 0;   glMat[10] = 1;
    glLoadMatrixf(glMat);
};

void Draw::billboardCamProj( ){
    float glCam  [16];
    float glModel[16];
    glGetFloatv (GL_MODELVIEW_MATRIX,  glModel);
    glGetFloatv (GL_PROJECTION_MATRIX, glCam);

    Mat3f mat;
    mat.a.set(glCam[0],glCam[1],glCam[2]);       mat.a.mul(1/mat.a.norm2());
    mat.b.set(glCam[4],glCam[5],glCam[6]);       mat.b.mul(1/mat.b.norm2());
    mat.c.set(glCam[8],glCam[9],glCam[10]);      mat.c.mul(1/mat.c.norm2());

    glModel[0 ] = mat.a.x;   glModel[1 ] = mat.b.x;   glModel[2 ] = mat.c.x;
    glModel[4 ] = mat.a.y;   glModel[5 ] = mat.b.y;   glModel[6 ] = mat.c.y;
    glModel[8 ] = mat.a.z;   glModel[9 ] = mat.b.z;   glModel[10] = mat.c.z;

    //glModel[0 ] = glCam[0];   glModel[1 ] = glCam[4];   glModel[2 ] = glCam[8];
    //glModel[4 ] = glCam[1];   glModel[5 ] = glCam[5];   glModel[6 ] = glCam[9];
    //glModel[8 ] = glCam[2];   glModel[9 ] = glCam[6];   glModel[10] = glCam[10];

    glLoadMatrixf(glModel);
};

void Draw::drawText( const char * str, int itex, float sz, int istart, int iend ){
    const int nchars = 95;
    float persprite = 1.0f/nchars;
    glEnable     ( GL_TEXTURE_2D );
    glBindTexture( GL_TEXTURE_2D, itex );
    glEnable(GL_BLEND);
    glEnable(GL_ALPHA_TEST);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    glBegin(GL_QUADS);
    int terminator = 0xFFFF;
    if(iend<=0) { terminator=-iend; iend=256; };
    for(int i=istart; i<iend; i++){
        if  (str[i]==terminator) break;
        int isprite = str[i] - 33;
        float offset  = isprite*persprite+(persprite*0.57);
        float xi = i*sz;
        glTexCoord2f( offset          , 1.0f ); glVertex3f( xi   ,    0, 0.0f );
        glTexCoord2f( offset+persprite, 1.0f ); glVertex3f( xi+sz,    0, 0.0f );
        glTexCoord2f( offset+persprite, 0.0f ); glVertex3f( xi+sz, sz*2, 0.0f );
        glTexCoord2f( offset          , 0.0f ); glVertex3f( xi   , sz*2, 0.0f );
    }
    glEnd();
    glDisable  ( GL_BLEND );
    glDisable  ( GL_ALPHA_TEST );
    glDisable  ( GL_TEXTURE_2D );
    glBlendFunc( GL_ONE, GL_ZERO );
};


/*
GLuint Draw::makeTexture( char * fname ){

    //SDL_Surface * surf = IMG_Load( fname );
    SDL_Surface * surf = SDL_LoadBMP( fname );
    if ( surf ){
        GLuint itex=0;
        glGenTextures  ( 1, &itex );
        glBindTexture  ( GL_TEXTURE_2D, itex );
        //if      (surf->format->BytesPerPixel == 1) { glTexImage2D( GL_TEXTURE_2D, 0, 1,  surf->w,  surf->h, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, surf->pixels ); }
        //if      (surf->format->BytesPerPixel == 1) { glTexImage2D( GL_TEXTURE_2D, 0, 1,   surf->w,  surf->h, 0, GL_INTENSITY, GL_UNSIGNED_BYTE, surf->pixels ); }
        //if      (surf->format->BytesPerPixel == 1) { glTexImage2D( GL_TEXTURE_2D, 0, 1,       surf->w,  surf->h, 0, GL_ALPHA,     GL_UNSIGNED_BYTE, surf->pixels ); }
        if      (surf->format->BytesPerPixel == 1) {
            glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
            glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
            glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
            glTexImage2D( GL_TEXTURE_2D, 0, 1,       surf->w,  surf->h, 0, GL_RED,       GL_UNSIGNED_BYTE, surf->pixels );
        }
        else if (surf->format->BytesPerPixel == 3) { glTexImage2D( GL_TEXTURE_2D, 0, GL_RGB,  surf->w,  surf->h, 0, GL_BGR,       GL_UNSIGNED_BYTE, surf->pixels ); }
        else if (surf->format->BytesPerPixel == 4) { glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, surf->w,  surf->h, 0, 0x8000,       GL_UNSIGNED_BYTE, surf->pixels ); }
        else return 0;
        printf( "surface->format->Rmask : %i itex %i \n", surf->format->BytesPerPixel, itex  ) ;// surface->format->Rmask/ == 0x000000ff;

        //glTexImage2D   ( GL_TEXTURE_2D, 0, 3, surf->w,  surf->h, 0, GL_BGR, GL_UNSIGNED_BYTE, surf->pixels );
        //glTexImage2D   ( GL_TEXTURE_2D, 0, 3, surf->w,  surf->h, 0, GL_BGRA, GL_UNSIGNED_BYTE, surf->pixels );
        //glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, surf->w,  surf->h, 0, GL_RGBA, GL_UNSIGNED_BYTE, surf->pixels );
        //glTexImage2D(GL_TEXTURE_2D, 0, GL_BGRA, surf->w,  surf->h, 0, GL_BGRA, GL_UNSIGNED_BYTE, surf->pixels );
        //glTexImage2D(GL_TEXTURE_2D, 0,  GL_RGBA, surf->w,  surf->h, 0, GL_RGBA, GL_UNSIGNED_BYTE, surf->pixels );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
        glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
        SDL_FreeSurface( surf );
        return itex;
    }else{
        printf( "cannot load %s\n", fname  );
    }
    return 0;
};
*/

/*
GLuint Draw::makeTexture( int nx, int ny, float * data ){

    GLuint itex=0;
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glGenTextures  ( 1, &itex );
    glBindTexture  ( GL_TEXTURE_2D, itex );

    //glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, nx, ny, 0, GL_RED, GL_FLOAT, data);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    int ntot = nx*ny;
    uint32_t * data_ = new uint32_t[ntot];
    for(int i=0;i<ntot;i++){
        //data_[i] = (int)(255*data[i]);
        //data_[i] = (int)(255*data[i]);
        //data_[i] = (int)(255*data[i]);
        //data_[i] = (int)(255*data[i]);
        uint8_t R = 0xFF; uint8_t G = 0xFF; uint8_t B = 0xFF; uint8_t A = 0xFF;
        data_[i] = (A<<24)|(B<<16)|(G<<8)|R;
    }


    //glTexImage2D( GL_TEXTURE_2D, 0, GL_RGB,  surf->w,  surf->h, 0, GL_BGR,       GL_UNSIGNED_BYTE, surf->pixels );
    //glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, nx, ny, 0, 0x8000,  GL_UNSIGNED_BYTE, data );
    glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, nx,  ny, 0, GL_RGBA,  GL_UNSIGNED_BYTE, data_ );

    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );

    delete[] data_;
    return itex;
};
*/



