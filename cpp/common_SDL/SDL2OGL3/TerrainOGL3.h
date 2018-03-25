
#ifndef  TerrainOGL3_h
#define  TerrainOGL3_h

#include "Vec2.h"
#include "Vec3.h"
#include "GLfunctions.h"
#include "Shader.h"
#include "GLObject.h"
#include "GLobjects.h"
#include "GLUtils.h"
#include "CMesh.h"
#include "GL3Utils.h"
#include "fastmath.h"

class TerrainOGL3{ public:
    Shader sh;
    //GLMesh
    Vec2i nSamp;
    double dist0;

    Vec3f pos      = (Vec3f){ 0.0,   0.0,    0.0 };
    Vec3f mapScale = (Vec3f){ 0.002, 0.002, 20.0 };
    Vec2f uv0      = (Vec2f){ 0.0,   0.0 };
    Vec2f viewMin  = (Vec2f){ 1.0,   0.0 };
    Vec2f viewMax  = (Vec2f){ 0.0,   1.0 };

    int nVertDrawn=0;
    int nDrawCalls=0;

    GLuint txHeight;
    GLuint stripUVs;
    //GLuint hexUVs;  // central hexagon - todo later
    //GLuint hexInds
    struct { GLuint
        uv0,
        pa,pb,
        mapScale,
        txHeight,
        lightColor,
        diffuseColor,
        ambientColor,
        specularColor,
        lightPos;
    } ulocs;

    void makeStrip( int n, float f ){
        float z = 1.0;
        Vec2f* uvs = new Vec2f[n*2];
        for(int i=0; i<n; i++){
            //float z_ = z-z0;
            //uvs[2*i].set(0.0,z_); uvs[2*i+1].set(1.0,z_);
            uvs[2*i].set(0.0,z); uvs[2*i+1].set(1.0,z);
            //printf( "%i (%f,%f) (%f,%f)\n", i, uvs[2*i].x,uvs[2*i].y,  uvs[2*i+1].x,uvs[2*i+1].y );
            z*=f;
        }
        printf( "z_max: %f ", z );
        newArrayBuffer( stripUVs, n*2*2*sizeof(GLfloat), (GLfloat*)uvs, GL_STATIC_DRAW );
        delete [] uvs;
        //exit(0);
    }

    /*
    void init( Vec2i nSamp_, float dist0_, Vec2i nHeighs, float* height_map ){
        nSamp = nSamp_;
        dist0 = dist0_;
        sh.init( "common_resources/shaders/terrain_strip.glslv", "common_resources/shaders/terrain_world.glslf" );
        //sh.init( "common_resources/shaders/terrain_strip.glslv", "common_resources/shaders/color3D.glslf" );
        sh.use();
        sh.getDefaultUniformLocation();
        ulocs.pa=sh.getUloc("pa"); ulocs.pb=sh.getUloc("pb"); ulocs.uv0=sh.getUloc("uv0"); ulocs.mapScale=sh.getUloc("mapScale"); ulocs.txHeight=sh.getUloc("txHeight");
        printf( "ulocs : %i %i %i %i %i \n", ulocs.pa, ulocs.pb, ulocs.uv0, ulocs.mapScale, txHeight );
        newTexture2D( txHeight, nHeighs.x, nHeighs.y, height_map, GL_RED, GL_FLOAT );
        //makeStrip( nSamp.b, dist0, 1.0 + 2.0/nSamp.a );
        makeStrip( nSamp.b, 1.0 + 2.0/nSamp.a );
    }
    */

    void fromHeightsDerivs_byte( int nx, int ny, float* height_map, float hsc ){
        uint8_t * hxyz = new uint8_t[nx*ny*4];
        for( int iy=0; iy<ny; iy++ ){
            for( int ix=0; ix<nx; ix++ ){
                int i  = iy*nx + ix;
                int i4 = i<<2;
                Vec3d nv;
                nv.x = hsc*(height_map[wrap_index2d_fast(ix+1,iy  ,nx,ny)]-height_map[wrap_index2d_fast(ix-1,iy  ,nx,ny)]);
                nv.y = hsc*(height_map[wrap_index2d_fast(ix  ,iy+1,nx,ny)]-height_map[wrap_index2d_fast(ix  ,iy-1,nx,ny)]);
                nv.z = 1;
                nv.normalize();
                //printf( " %i %i   |  %f %f %f  %f \n", ix, iy, nv.x, nv.y, nv.z, height_map[ i ] );
                hxyz[ i4   ] = (uint8_t)( (nv.x+1.0)*128 ); // a
                hxyz[ i4+1 ] = (uint8_t)( (nv.y+1.0)*128 );
                hxyz[ i4+2 ] = (uint8_t)( (nv.z+1.0)*128 );
                hxyz[ i4+3 ] = (int)( height_map[ i ]*255 );
            }
        }
        newTexture2D( txHeight, nx, ny, hxyz, GL_RGBA, GL_UNSIGNED_BYTE );
        delete [] hxyz;
    }

    void fromHxy_byte( int nx, int ny, float* height_map, float hsc ){
        uint8_t * hxy = new uint8_t[nx*ny*3];
        for( int iy=0; iy<ny; iy++ ){
            for( int ix=0; ix<nx; ix++ ){
                int i  = iy*nx + ix;
                int i3 = i*3;
                float dx = hsc*(height_map[wrap_index2d_fast(ix+1,iy  ,nx,ny)]-height_map[wrap_index2d_fast(ix-1,iy  ,nx,ny)]);
                float dy = hsc*(height_map[wrap_index2d_fast(ix  ,iy+1,nx,ny)]-height_map[wrap_index2d_fast(ix  ,iy-1,nx,ny)]);
                //printf( " %i %i   |  %f %f %f  %f \n", ix, iy, nv.x, nv.y, nv.z, height_map[ i ] );
                hxy[ i3   ] = (uint8_t)( _clamp(dx,-1.0,1.0)*127 + 128 ); // a
                hxy[ i3+1 ] = (uint8_t)( _clamp(dy,-1.0,1.0)*127 + 128 );
                hxy[ i3+2 ] = (uint8_t)( _clamp(height_map[ i ], 0.0,1.0 )*255 );
            }
        }
        newTexture2D( txHeight, nx, ny, hxy, GL_RGB, GL_UNSIGNED_BYTE );
        delete [] hxy;
    }

    void fromHeightsDerivs_float( int nx, int ny, float* height_map, float hsc ){
        float * hxyz = new float[nx*ny*4];
        for( int iy=0; iy<ny; iy++ ){
            for( int ix=0; ix<nx; ix++ ){
                int i  = iy*nx + ix;
                int i4 = i<<2;
                Vec3d nv;
                nv.x = hsc*(height_map[wrap_index2d_fast(ix+1,iy  ,nx,ny)]-height_map[wrap_index2d_fast(ix-1,iy  ,nx,ny)]);
                nv.y = hsc*(height_map[wrap_index2d_fast(ix  ,iy+1,nx,ny)]-height_map[wrap_index2d_fast(ix  ,iy-1,nx,ny)]);
                nv.z = 1;
                nv.normalize();
                //printf( " %i %i   |  %f %f %f  %f \n", ix, iy, nv.x, nv.y, nv.z, height_map[ i ] );
                hxyz[ i4   ] = nv.x;
                hxyz[ i4+1 ] = nv.y;
                hxyz[ i4+2 ] = nv.z;
                hxyz[ i4+3 ] = height_map[ i ];
            }
        }
        newTexture2D( txHeight, nx, ny, hxyz, GL_RGBA, GL_FLOAT ); // this does not seem to work
        delete [] hxyz;
    }

    void init( Vec2i nSamp_, float dist0_, Vec2i nHeighs, float* height_map ){
        nSamp = nSamp_;
        dist0 = dist0_;
        sh.init( "common_resources/shaders/terrain_strip.glslv", "common_resources/shaders/terrain_world.glslf" );
        //sh.init( "common_resources/shaders/terrain_strip.glslv", "common_resources/shaders/color3D.glslf" );
        sh.use();
        sh.getDefaultUniformLocation();
        ulocs.pa=sh.getUloc("pa"); ulocs.pb=sh.getUloc("pb"); ulocs.uv0=sh.getUloc("uv0"); ulocs.mapScale=sh.getUloc("mapScale"); ulocs.txHeight=sh.getUloc("txHeight");
        ulocs.lightColor    = sh.getUloc("lightColor"   );
        ulocs.diffuseColor  = sh.getUloc("diffuseColor" );
        ulocs.ambientColor  = sh.getUloc("ambientColor" );
        ulocs.specularColor = sh.getUloc("specularColor");
        ulocs.lightPos      = sh.getUloc("lightPos"     );

        printf( "ulocs : %i %i %i %i %i \n", ulocs.pa, ulocs.pb, ulocs.uv0, ulocs.mapScale, txHeight );
        glUniform3f(ulocs.lightColor ,    0.50, 0.45, 0.40 );
        glUniform3f(ulocs.diffuseColor,   1.00, 1.00, 1.00 );
        glUniform3f(ulocs.ambientColor ,  0.20, 0.25, 0.30 );
        glUniform3f(ulocs.specularColor , 2.00, 2.00, 2.00 );
        glUniform3f(ulocs.lightPos,        0.0,+1000.0,0.0 );

        //fromHeightsDerivs_byte ( nHeighs.x, nHeighs.y, height_map, 5.0 );
        fromHxy_byte ( nHeighs.x, nHeighs.y, height_map, 100.0 );
        //fromHeightsDerivs_float( nHeighs.x, nHeighs.y, height_map, 5.0 );
        //newTexture2D( txHeight, nHeighs.x, nHeighs.y, height_map, GL_RED, GL_FLOAT );


        //makeStrip( nSamp.b, dist0, 1.0 + 2.0/nSamp.a );
        makeStrip( nSamp.b, 1.0 + 2.0/nSamp.a );
    }

    void setViewRange( Vec2f camDir, float camTg ){
        Vec2f rot; rot.x=sqrt(1.0/(1.0+camTg*camTg)); rot.y=rot.x*camTg;
        viewMin.set_udiv_cmplx(camDir,rot);
        viewMax.set_mul_cmplx (camDir,rot);
        //printf(" (%f,%f), (%f,%f), (%f,%f) (%f,%f) \n", camDir.x,camDir.y,  viewMin.x,viewMin.y,   viewMax.x,viewMax.y,    rot.x,rot.y );
        //exit(0);
    }

    void drawStrips( Vec2f pa, Vec2f pb ){
        Vec2f op,dp,p;
        dp = (pb-pa)*(1.0/nSamp.a);
        p=pa+dp; op = pa;
        for(int i=0; i<nSamp.a; i++ ){
            if( p.isBetweenRotations(viewMin,viewMax) ){
                nDrawCalls++;
                nVertDrawn+=nSamp.b*2;
                glUniform2f( ulocs.pa, op.x, op.y );
                glUniform2f( ulocs.pb,  p.x,  p.y );
                glDrawArrays( GL_TRIANGLE_STRIP, 0, nSamp.b );
                //printf( "%i (%f,%f) (%f,%f)\n", i, op.x, op.y,  p.x, p.y  );
            }
            op = p; p.add(dp);
        }
    }

    void draw(){
        nDrawCalls=0;
        nVertDrawn=0;
        //sh.use();
        glUniform2f ( ulocs.uv0,      uv0.x, uv0.y );
        glUniform3fv( ulocs.mapScale,1, (GLfloat*)&mapScale );
        sh.set_modelPos( (GLfloat*)&pos );

        Vec2f drot; drot.fromAngle( M_PI/3.0 );
        Vec2f p = (Vec2f){dist0,0.0};
        bindVertexAttribPointer( 0, stripUVs, 2, GL_FLOAT, GL_FALSE );
        //bindTexture( 0, txHeight, ulocs.txHeight );
        bindTexture( 0, txHeight, ulocs.txHeight  );
        for(int i=0; i<6; i++){
            Vec2f p_; p_.set_mul_cmplx(p,drot);
            //printf( " === %i (%f,%f) (%f,%f)\n", i, p.x, p.y,  p_.x, p_.y  );
            drawStrips( p, p_ );
            p = p_;
        }
        //exit(0);
    }

};

#endif
