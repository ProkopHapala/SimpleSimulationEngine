
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

class TerrainOGL3{ public:
    Shader sh;
    //GLMesh
    Vec2i nSamp;
    double dist0;

    Vec3f pos      = (Vec3f){ 0.0,  0.0,  0.0 };
    Vec3f mapScale = (Vec3f){0.002, 0.002, 20.0 };
    Vec2f uv0      = (Vec2f){ 0.0,  0.0  };

    Vec2f viewMin = (Vec2f){1.0,0.0};
    Vec2f viewMax = (Vec2f){0.0,1.0};


    GLuint txHeight;
    GLuint stripUVs;
    //GLuint hexUVs;  // central hexagon - todo later
    //GLuint hexInds
    struct { GLuint uv0,pa,pb,mapScale,txHeight; } ulocs;

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

    void setViewRange( Vec2f camDir, float camTg ){
        Vec2f rot; rot.x=sqrt(1.0/(1.0+camTg*camTg)); rot.y=rot.x*camTg;
        viewMin.set_udiv_cmplx(camDir,rot);
        viewMax.set_mul_cmplx (camDir,rot);
        printf(" (%f,%f), (%f,%f), (%f,%f) (%f,%f) \n", camDir.x,camDir.y,  viewMin.x,viewMin.y,   viewMax.x,viewMax.y,    rot.x,rot.y );
        //exit(0);
    }

    void drawStrips( Vec2f pa, Vec2f pb ){
        Vec2f op,dp,p;
        dp = (pb-pa)*(1.0/nSamp.a);
        p=pa+dp; op = pa;
        for(int i=0; i<nSamp.a; i++ ){
            if( p.isBetweenRotations(viewMin,viewMax) ){
                glUniform2f( ulocs.pa, op.x, op.y );
                glUniform2f( ulocs.pb,  p.x, p.y );
                glDrawArrays( GL_TRIANGLE_STRIP, 0, nSamp.b );
                //printf( "%i (%f,%f) (%f,%f)\n", i, op.x, op.y,  p.x, p.y  );
            }
            op = p; p.add(dp);
        }
    }

    void draw(){
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
