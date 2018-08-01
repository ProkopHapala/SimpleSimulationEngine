
#include <stdlib.h>
#include <stdio.h>

#include <vector>
#include <GL/glew.h>
//#define GL_GLEXT_PROTOTYPES
//#include <GL/gl.h>
//#include <SDL2/SDL.h>

#include <fastmath.h>
#include <Vec2.h>
#include <Vec3.h>
#include <Mat3.h>
#include <quaternion.h>
#include <raytrace.h>
//#include <Body.h>

#include "Shader.h"
#include "GLObject.h"
#include "DrawOGL3.h"
#include "SceneOGL3.h"
#include "ScreenSDL2OGL3.h"
#include "AppSDL2OGL3.h"

#include "Mesh.h"
#include "Solids.h"
#include "GLfunctions.h"
#include "GLobjects.h"
#include "GLObject.h"
#include "GL3Utils.h"
#include "Shader.h"

// =============== Global variables

// http://www.real-time-volume-graphics.org/?page_id=28
// http://moddb.wikia.com/wiki/OpenGL:Tutorials:3D_Textures

int frameCount = 0;
double lastTime = 0.0;

// =============== Functions

/*
GLMesh* makeQuad3D( Vec2f p0, Vec2f p1, Vec2f u0, Vec2f u1 ){
    GLfloat verts[]  = { p0.x,p0.y,0.0, p1.x,p1.y,0.0, p0.x,p1.y,0.0,  p0.x,p0.y,0.0, p1.x,p1.y,0.0, p1.x,p0.y,0.0 };
    //Vec2f vUVs   [] = { u0.x,u0.y,     u1.x,u1.y,     u0.x,u1.y,       u0.x,u0.y,     u1.x,u1.y,     u1.x,u0.y     };
    GLfloat colors[] = { u0.x,u0.y,0.0, u1.x,u1.y,0.0, u0.x,u1.y,0.0,  u0.x,u0.y,0.0, u1.x,u1.y,0.0, u1.x,u0.y,0.0 };
    GLMesh* glquad = new GLMesh();
    //glquad->init( 6, 0, NULL, verts, NULL, NULL, vUVs );
    glquad->init( 6, 0, NULL, verts, NULL, colors, NULL );
    return glquad;
}
*/

uint8_t* makeImage3DRGBA( int nx, int ny, int nz ){
    int i = 0;
    uint8_t * buff = new uint8_t[nx*ny*nz*4];
    float dx=1.0/nx; float dy=1.0/ny; float dz=1.0/nz;
    for(int iz=0; iz<nz; iz++){
        for(int iy=0; iy<ny; iy++){
            for(int ix=0; ix<nx; ix++){
                int ioff=i*4;
                uint8_t b = ix^iy^iz;
                //uint8_t b = ix*iz + iy*iz;
                /*
                buff[ioff+0]=(int)255/(1+(ix-nx/2)*(ix-nx/2)*0.1) + b;
                buff[ioff+1]=(int)255/(1+(iy-ny/2)*(iy-ny/2)*0.1) + b;
                buff[ioff+2]=(int)255/(1+(iz-nz/2)*(iz-nz/2)*0.1) + b;
                buff[ioff+3]=255;
                */
                float x = ix*dx-0.5; float y = iy*dy-0.5; float z = iz*dz-0.5;
                int r=rand();
                buff[ioff+0] = (r    )&0xff;
                buff[ioff+1] = (r>>8 )&0xff;
                buff[ioff+2] = (r>>16)&0xff;
                //buff[ioff+3] =(int)255*exp( -10*(x*x + y*y + z*z) );
                float f1 = sin(x*10)-10*y;
                float f2 = sin(z);
                buff[ioff+3] =(int)255*( 1/(1.0+20*f1*f1+20*f2*f2) );
                i++;
            }
        }
    }
    return buff;
}

uint8_t* makeImage3DDistance( int nx, int ny, int nz, int sz, uint8_t transhold, uint8_t* inBuff ){
    uint8_t * outBuff = new uint8_t[nx*ny*nz];
    float dx=1.0/nx; float dy=1.0/ny; float dz=1.0/nz;
    for(int iz=0; iz<nz; iz++){
        int jzmin = iz-sz; if(jzmin<0 ) jzmin=0;
        int jzmax = iz+sz; if(jzmax>nz) jzmax=nz;
        for(int iy=0; iy<ny; iy++){
            int jymin = iy-sz; if(jymin<0 ) jymin=0;
            int jymax = iy+sz; if(jymax>ny) jymax=ny;
            for(int ix=0; ix<nx; ix++){
                int jxmin = ix-sz; if(jxmin<0 ) jxmin=0;
                int jxmax = ix+sz; if(jxmax>nx) jxmax=nx;
                int ddmin = sz*sz;
                for(int jz=jzmin; jz<jzmax; jz++){
                    for(int jy=jymin; jy<jymax; jy++){
                        for(int jx=jxmin; jx<jxmax; jx++){
                            int j = jx + nx*(jy + ny*jz);
                            if( inBuff[j] > transhold ){
                                int dx = jx-ix;
                                int dy = jy-iy;
                                int dz = jz-iz;
                                int dd = dx*dx + dy*dy + dz*dz;
                                if(dd<ddmin){
                                    ddmin = dd;
                                }
                            }
                        }
                    }
                }
                int i = ix + nx*(iy + ny*iz);
                //printf( " (%i,%i,%i) %g   ", ix,iy,iz, sqrt(ddmin) );
                outBuff[i] = 255*( sqrt(ddmin)/nz );
            }
        }
    }
    return outBuff;
}

uint8_t* makeImage3D( int nx, int ny, int nz ){
    int i = 0;
    uint8_t * buff = new uint8_t[nx*ny*nz];
    float dx=1.0/nx; float dy=1.0/ny; float dz=1.0/nz;
    for(int iz=0; iz<nz; iz++){
        for(int iy=0; iy<ny; iy++){
            for(int ix=0; ix<nx; ix++){
                float x = ix*dx-0.5; float y = iy*dy-0.5; float z = iz*dz-0.5;
                //float x = ix*dx; float y = iy*dy; float z = iz*dz;
                float damp = exp( -20.0*(x*x + y*y + z*z) )*(x*x+y*y)*5.0;
                //damp*=damp;
                float noise = (0.1+randf());
                if((ix==0)||(ix==nx-1)) damp=0;
                if((iy==0)||(iy==ny-1)) damp=0;
                if((iz==0)||(iz==nz-1)) damp=0;
                buff[i] =(int)255*( damp*noise );
                //buff[i] = ix^iy^iz;

                i++;
            }
        }
    }
    return buff;
}

uint8_t* makeImage3DWhiteNoise( int nx, int ny, int nz ){
    int i = 0;
    uint8_t * buff = new uint8_t[nx*ny*nz];
    float dx=1.0/nx; float dy=1.0/ny; float dz=1.0/nz;
    for(int iz=0; iz<nz; iz++){
        for(int iy=0; iy<ny; iy++){
            for(int ix=0; ix<nx; ix++){
                buff[i] = rand()&0xFF;
                i++;
            }
        }
    }
    return buff;
}

class TestAppScreenOGL3: public AppSDL2OGL3, public SceneOGL3 { public:
    //virtual void draw(){}
    //Shader *sh1;

    Shader *sh1,*sh2;
    GLMesh *msh1,*msh2;
    GLuint  tx3D_1,tx3D_2;

    //int npoints=0;
    //Vec3f* points=NULL;

    TestAppScreenOGL3():AppSDL2OGL3(800,600),SceneOGL3(){
        // FIXME: second window does not work :((((
        //SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
        //screens.push_back( new ScreenSDL2OGL3( 800, 600) );
        for( ScreenSDL2OGL3* screen: screens ) screen->scenes.push_back( this );

        // from here: http://moddb.wikia.com/wiki/OpenGL:Tutorials:3D_Textures
        glEnable(GL_TEXTURE_3D);
        int NX=32,NY=32,NZ=32;
        //int NX=64,NY=64,NZ=64;
        uint8_t* buff     = makeImage3D( NX, NY, NZ );
        uint8_t* buffDist = makeImage3DDistance( NX, NY, NZ, 16, 16, buff );

        glGenTextures(1, &tx3D_1);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_3D, tx3D_1);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        //glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        //glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP);
        //glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA8, NX, NY, NZ, 0, GL_RGBA, GL_UNSIGNED_BYTE, buff);
        //glTexImage3D(GL_TEXTURE_3D,0,GL_INTENSITY,NX,NY,NZ,0,GL_LUMINANCE,GL_UNSIGNED_BYTE,buff);
        //glTexImage3D(GL_TEXTURE_3D,0,GL_RED,NX,NY,NZ,0,GL_RED,GL_UNSIGNED_BYTE,buff);
        glTexImage3D(GL_TEXTURE_3D,0,GL_RED,NX,NY,NZ,0,GL_RED,GL_UNSIGNED_BYTE,buffDist);
        //glBindTexture(GL_TEXTURE_3D, tx3D_1);

        buff = makeImage3DWhiteNoise(NX,NY,NZ);
        glGenTextures(1, &tx3D_2);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_3D, tx3D_2);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_REPEAT);
        glTexImage3D(GL_TEXTURE_3D,0,GL_RED,NX,NY,NZ,0,GL_RED,GL_UNSIGNED_BYTE,buff);

        delete [] buffDist;
        delete [] buff;

        sh1=new Shader();
        //sh1->init( "common_resources/shaders/const3D.glslv",   "common_resources/shaders/const3D.glslf"   );
        //sh1->init( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/color3D.glslf"   );
        //sh1->init( "common_resources/shaders/rayMarch3DTexture.glslv",   "common_resources/shaders/rayMarch3DTexture.glslf"   );
        sh1->init( "common_resources/shaders/rayMarch3DTexture.glslv","common_resources/shaders/rayMarchCloudTexture3D.glslf"   );
        sh1->getDefaultUniformLocation();

        sh2=new Shader();
        sh2->init( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/cut3DTexture.glslf"   );
        sh2->getDefaultUniformLocation();

        
        //msh1 = new GLMesh(); msh1->init_wireframe( Solids::Octahedron );
        //msh1 = new GLMesh(); msh1->init_wireframe( Solids::Octahedron );
        //msh1 = hardTriangles2mesh( Solids::Octahedron );
        //msh1 = hardTriangles2mesh( Solids::Icosahedron );
        msh1 = hardTriangles2mesh( Solids::Cube );

        msh2 = makeQuad3D( {0.0f,0.0f}, {1.0f,1.0f}, {0.0f,0.0f}, {1.0f,1.0f} );

        Camera& cam = screens[0]->cam;
        cam.zmin = 1.0; cam.zmax = 1000.0; cam.zoom = 20.00f;
        cam.aspect = screens[0]->HEIGHT/(float)screens[0]->WIDTH;
        screen->camLookAt=new Vec3f(); *screen->camLookAt = (Vec3f){0.0,0.0,0.0};
        printf( "SETUP DONE \n" );
    }

    virtual void draw( Camera& cam ){
        //glClearColor(1.0, 1.0, 1.0, 1.0);
        glClearColor(0.0, 0.0, 1.0, 1.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );
        glEnable(GL_DEPTH_TEST);


        glEnable(GL_BLEND);
        glEnable(GL_CULL_FACE);
        glCullFace(GL_FRONT); // This is hack since we have probably inverted normals
        //glBlendEquation( GL_FUNC_ADD )
        glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);


        Shader* sh=sh1;
        sh->use();
        Mat3f mrot; mrot.setOne();
        sh->set_modelMat( (GLfloat*)&mrot );
        sh->set_modelPos( (const GLfloat[]){0.0f,0.0f,0.0f} );
        setCamera( *sh1, cam );

        glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, tx3D_1);
        glActiveTexture(GL_TEXTURE2); glBindTexture(GL_TEXTURE_2D, tx3D_2);
        glBindTexture(GL_TEXTURE_3D, tx3D_1);
        glUniform1i(sh->getUloc("texture_1"),     0  );
        glUniform1i(sh->getUloc("texture_noise"), 1  );
        glUniform1f(sh->getUloc("txScale"  ), 0.5);
        glUniform3f(sh->getUloc("txOffset"), 0.0,0.0,frameCount*0.01);

        msh1->draw(GL_TRIANGLES);
        msh1->draw();

        double fjunk;

        sh=sh2;
        sh->use();
        mrot.setOne();
        sh->set_modelMat( (GLfloat*)&mrot );
        sh->set_modelPos( (const GLfloat[]){0.0f,0.0f,0.0f} );
        setCamera( *sh1, cam );
        glUniform3f(sh->getUloc("txOffset"), 0.0,0.0, modf( frameCount/500.0, &fjunk )  );

        msh2->draw();

    };


};


// ================== main

TestAppScreenOGL3 * app;

int main(int argc, char *argv[]){
    app = new TestAppScreenOGL3( );
    app->loop( 1000000 );
    app->quit();
    return 0;
}




