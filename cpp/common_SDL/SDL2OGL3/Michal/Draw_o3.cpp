
#include "Renderer.h"
#include "Vec2.h"
#include "Vec3.h"

#include "MeshBuilder.h"
#include "TextRenderer.h"
#include "Draw_o3.h"  // THE HEADER

void Draw::colorScale( double d, int ncol, const uint32_t * colors ){
    constexpr float inv255 = 1.0f/255.0f;
    //double d_bak = d;
    _clamp(d,0.000001,0.999999);
    d*=(ncol-1);
    int icol = (int)d;
    // if(icol<0) Draw::setRGB( colors[0] );
    // else if(icol>=ncol-1){ uint32_t i=colors[ncol-1];  opengl1renderer.color3f( ((i)&0xFF)*inv255, ((i>>8)&0xFF)*inv255, ((i>>16)&0xFF)*inv255  ); return; }
    // else if(icol<      0){ uint32_t i=colors[    0 ];  opengl1renderer.color3f( ((i)&0xFF)*inv255, ((i>>8)&0xFF)*inv255, ((i>>16)&0xFF)*inv255  ); return; }
    //if(icol>=ncol-1){ printf( "ERROR: Draw::colorScale() icol(%i)>=ncol(%i)-1 d=%g  d_bak=%g \n", icol, ncol, d, d_bak ); }
    d-=icol; double md = 1-d; // linear interpolation coefficients
    //printf( "Draw::colorScale() d,md %g %g icol/ncol %i/%i \n", d, md, icol, ncol );    
    
    if( (icol<0)||((icol+1)>ncol)){
        printf( "Draw::colorScale() icol/ncol %i/%i d,md %g %g \n", icol, ncol, d, md );   
        exit(0);
    }
    uint32_t clr1=colors[icol  ];
    uint32_t clr2=colors[icol+1];
    opengl1renderer.color3f(
        ( d*( clr2     &0xFF) + md*( clr1     &0xFF ))*inv255,
        ( d*((clr2>>8 )&0xFF) + md*((clr1>>8 )&0xFF ))*inv255,
        ( d*((clr2>>16)&0xFF) + md*((clr1>>16)&0xFF ))*inv255
    );
    // From setRGB( uint32_t i ){
    //opengl1renderer.color3f( ((i>>16)&0xFF)*inv255, ((i>>8)&0xFF)*inv255, (i&0xFF)*inv255  );
};

uint32_t Draw::icolorScale( double d, int ncol, const uint32_t * colors ){
    constexpr float inv255 = 1.0f/255.0f;
    d*=(ncol-1);
    int icol = (int)d;
    d-=icol; double md = 1-d;
    uint32_t clr1=colors[icol  ];
    uint32_t clr2=colors[icol+1];
    return 0xFF000000
        | (((uint32_t)( d*( clr2     &0xFF) + md*( clr1     &0xFF )))    )
        | (((uint32_t)( d*((clr2>>8 )&0xFF) + md*((clr1>>8 )&0xFF )))<<8 )
        | (((uint32_t)( d*((clr2>>16)&0xFF) + md*((clr1>>16)&0xFF )))<<16);
};

/*
void Draw::setColorInt32( uint32_t clr ) {
    constexpr float i255 = 1/255.0f;
    uint8_t b = ( ( clr       ) & 0xFF );
    uint8_t g = ( ( clr >> 8  ) & 0xFF );
    uint8_t r = ( ( clr >> 16 ) & 0xFF );
    uint8_t a = (   clr >> 24          );
    opengl1renderer.color4f( i255*r, i255*g, i255*b, i255*a );
    //printf( " r %i g %i b %i a %i     %f %f %f %f  \n", r, g, b, a,  i255*r, i255*g, i255*b, i255*a   );
};
*/

void Draw::billboardCam( ){
    float glMat[16];
    //opengl1renderer.matrixMode(GL_MODELVIEW);
    opengl1renderer.getFloatv(GL_MODELVIEW_MATRIX , glMat);
    glMat[0 ] = 1;   glMat[1 ] = 0;   glMat[2 ] = 0;
    glMat[4 ] = 0;   glMat[5 ] = 1;   glMat[6 ] = 0;
    glMat[8 ] = 0;   glMat[9 ] = 0;   glMat[10] = 1;
    opengl1renderer.loadMatrixf(glMat);
};

void Draw::billboardCamProj( float scale_ ){
    //printf( "billboardCamProj(%g) \n", scale );
    float glCam  [16];
    float glModel[16];
    opengl1renderer.getFloatv (GL_MODELVIEW_MATRIX,  glModel );
    opengl1renderer.getFloatv (GL_PROJECTION_MATRIX, glCam   );
    Mat3f mat;
    mat.a.set(glCam[0],glCam[1],glCam[2]);       //mat.a.mul(1/mat.a.norm2());
    mat.b.set(glCam[4],glCam[5],glCam[6]);       //mat.b.mul(1/mat.b.norm2());
    mat.c.set(glCam[8],glCam[9],glCam[10]);      //mat.c.mul(1/mat.c.norm2());
    //float scale = 1/( scale_ * ( mat.a.norm2() + mat.b.norm2() + mat.c.norm2() ) );
    float scale = 1/( mat.a.norm2() + mat.b.norm2() + mat.c.norm2() );
    mat.a.mul(scale); mat.b.mul(scale);mat.c.mul(scale);
    glModel[0 ] = mat.a.x;   glModel[1 ] = mat.b.x;   glModel[2 ] = mat.c.x;
    glModel[4 ] = mat.a.y;   glModel[5 ] = mat.b.y;   glModel[6 ] = mat.c.y;
    glModel[8 ] = mat.a.z;   glModel[9 ] = mat.b.z;   glModel[10] = mat.c.z;
    opengl1renderer.loadMatrixf(glModel);
};

/*
// ---------- Backup
void Draw::billboardCamProj( ){
    printf( "billboardCamProj() \n" );
    float glCam  [16];
    float glModel[16];
    opengl1renderer.getFloatv (GL_MODELVIEW_MATRIX,  glModel);
    opengl1renderer.getFloatv (GL_PROJECTION_MATRIX, glCam);
    //opengl1renderer.matrixMode(GL_MODELVIEW);

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

    opengl1renderer.loadMatrixf(glModel);
};
*/

TextRenderer text;
void Draw::drawText( const char * str, Vec3f pos, float sz, int iend, Vec3f color ){
    text.set(str);
    text.draw2D(pos, sz, color);
}

void Draw::drawText( const char * str, Vec3f pos, float sz, Vec2i block_size, Vec3f color ){
    text.set(str, block_size.x, block_size.y);
    text.draw2D(pos, sz, color);
};
