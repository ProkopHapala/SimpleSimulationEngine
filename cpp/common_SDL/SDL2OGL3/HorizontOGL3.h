
#ifndef  HorizontOGL3_h
#define  HorizontOGL3_h

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

//
/*
Motivation:

- pre-render supersampled far-horizont with distance (depth), and map it on cylindrical mesh arund camera at large distance (e.g. 1000 m)
  assuming limited velocity of observe, pralax chnages (i.e. change of view angle) of distant object are very slow, therefore background sceen does not have to be updated often
- In case we need recalculate paralax, it ti best interpolate 3-4 textures captured at corners of simplex (triangle,tetrahedron) arround current position. when observer exit simplex we need to recalculate only one point

*/

class HorizontOGL3{ public:
    Shader sh;
    Vec2i nSamp;
    double dist0;

    Vec3f pos      = (Vec3f){ 0.0,   0.0,    0.0 };
    Vec3f mapScale = (Vec3f){ 0.002, 0.002, 20.0 };
    Vec2f uv0      = (Vec2f){ 0.0,   0.0 };
    Vec2f viewMin  = (Vec2f){ 1.0,   0.0 };
    Vec2f viewMax  = (Vec2f){ 0.0,   1.0 };

    int nVertDrawn=0;
    int nDrawCalls=0;

    GLuint tx;
    GLuint stripUVs;
    GLuint stripVerts;

    struct { GLuint
        //uv0,
        //pa,pb,
        //mapScale,
        tx
        //lightColor,
        //diffuseColor,
        //ambientColor,
        //specularColor,
        //lightPos;
        ;
    } ulocs;

    void makeStrip( int n, float ymin, float ymax){
        Vec2d ph = Vec2dX;
        float du = 1.0/n;
        Vec2d dph; dph.fromAngle(2*M_PI/n);
        Vec3f* verts = new Vec3f[n*2];
        Vec2f* uvs   = new Vec2f[n*2];
        for(int i=0; i<n*2; i+=2){
            verts[i].set(ph.x,ymin,ph.y);
            verts[i].set(ph.x,ymax,ph.y);
            uvs  [i].set(du*i,0.0);
            uvs  [i].set(du*i,1.0);
            ph.mul_cmplx(dph);
        }
        newArrayBuffer( stripVerts, n*2*3*sizeof(GLfloat), (GLfloat*)verts, GL_STATIC_DRAW );
        newArrayBuffer( stripUVs,   n*2*2*sizeof(GLfloat), (GLfloat*)uvs,   GL_STATIC_DRAW );
        delete [] uvs;
        delete [] verts;
    }


    void init( Vec2i nSamp_, float dist0_, Vec2i nHeighs, float* height_map ){
        nSamp = nSamp_;
        dist0 = dist0_;
        //sh.init( "common_resources/shaders/terrain_strip.glslv", "common_resources/shaders/terrain_world.glslf" );
        sh.init( "common_resources/shaders/texture3D.glslv",   "common_resources/shaders/texture.glslf"  );
        //sh.init( "common_resources/shaders/terrain_strip.glslv", "common_resources/shaders/color3D.glslf" );
        sh.use();
        sh.getDefaultUniformLocation();

        /*
        ulocs.pa=sh.getUloc("pa"); ulocs.pb=sh.getUloc("pb"); ulocs.uv0=sh.getUloc("uv0"); ulocs.mapScale=sh.getUloc("mapScale"); ulocs.txHeight=sh.getUloc("txHeight");
        ulocs.lightColor    = sh.getUloc("lightColor"   );
        ulocs.diffuseColor  = sh.getUloc("diffuseColor" );
        ulocs.ambientColor  = sh.getUloc("ambientColor" );
        ulocs.specularColor = sh.getUloc("specularColor");
        ulocs.lightPos      = sh.getUloc("lightPos"     );
       // printf( "ulocs : %i %i %i %i %i \n", ulocs.pa, ulocs.pb, ulocs.uv0, ulocs.mapScale, txHeight );
        glUniform3f(ulocs.lightColor ,    0.50, 0.45, 0.40 );
        glUniform3f(ulocs.diffuseColor,   1.00, 1.00, 1.00 );
        glUniform3f(ulocs.ambientColor ,  0.20, 0.25, 0.30 );
        glUniform3f(ulocs.specularColor , 2.00, 2.00, 2.00 );
        glUniform3f(ulocs.lightPos,        0.0,+1000.0,0.0 );
        */

        //fromHeightsDerivs_byte ( nHeighs.x, nHeighs.y, height_map, 5.0 );
        //fromHxy_byte ( nHeighs.x, nHeighs.y, height_map, 100.0 );
        //fromHeightsDerivs_float( nHeighs.x, nHeighs.y, height_map, 5.0 );
        //newTexture2D( txHeight, nHeighs.x, nHeighs.y, height_map, GL_RED, GL_FLOAT );


        //makeStrip( nSamp.b, dist0, 1.0 + 2.0/nSamp.a );
        makeStrip( nSamp.b,  0.0, 1.0 );
    }

/*
    void setViewRange( Vec2f camDir, float camTg ){
        Vec2f rot; rot.x=sqrt(1.0/(1.0+camTg*camTg)); rot.y=rot.x*camTg;
        viewMin.set_udiv_cmplx(camDir,rot);
        viewMax.set_mul_cmplx (camDir,rot);
        //printf(" (%f,%f), (%f,%f), (%f,%f) (%f,%f) \n", camDir.x,camDir.y,  viewMin.x,viewMin.y,   viewMax.x,viewMax.y,    rot.x,rot.y );
        //exit(0);
    }
*/

    void draw(){
        nDrawCalls=0;
        nVertDrawn=0;
        //sh.use();
        //glUniform2f ( ulocs.uv0,      uv0.x, uv0.y );
        //glUniform3fv( ulocs.mapScale,1, (GLfloat*)&mapScale );
        sh.set_modelPos( (GLfloat*)&pos );

        Vec2f drot; drot.fromAngle( M_PI/3.0 );
        Vec2f p = (Vec2f){dist0,0.0};
        bindVertexAttribPointer( 0, stripVerts, 3, GL_FLOAT, GL_FALSE );
        bindVertexAttribPointer( 1, stripUVs,   2, GL_FLOAT, GL_FALSE );
        //bindTexture( 0, txHeight, ulocs.txHeight );
        bindTexture( 0, tx, ulocs.tx  );
        glDrawArrays( GL_TRIANGLE_STRIP, 0, nSamp.x*2);
        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(1);
        /*
        for(int i=0; i<6; i++){
            Vec2f p_; p_.set_mul_cmplx(p,drot);
            //printf( " === %i (%f,%f) (%f,%f)\n", i, p.x, p.y,  p_.x, p_.y  );
            drawStrips( p, p_ );
            p = p_;
            #ifdef _DEBUG_VIEW_
            DEBUG_mesh->addLine( pos, {p.x,100.0,p.y},{1.0,0.0,1.0});
            #endif // _DEBUG_VIEW_
        }
        */
        //exit(0);
    }

};

#endif
