
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
#include "DrawOGL3.h"
#include "fastmath.h"

/*

TODO:
 - Procedural terrain shader (Perilin/Simplex Noise)
 - Multi respolution terrain texture - Main features by standard texture, smaller features by random noise texture
 - Color map except height map
 - Cubic terrain texture interpolation
 - Terrain shading using pixel-shader ray

*/

int heightTextureFromHeightsDerivs_byte( int nx, int ny, float* height_map, float hsc ){
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
    GLuint txHeight = 0;
    newTexture2D( txHeight, nx, ny, hxyz, GL_RGBA, GL_UNSIGNED_BYTE );
    delete [] hxyz;
    return txHeight;
}

int heightTextureFromHxy_byte( int nx, int ny, float* height_map, float hsc ){
    uint8_t * hxy = new uint8_t[nx*ny*3];
    float inv255 = 1.0/255.0;


    for( int iy=0; iy<ny; iy++ ){
        for( int ix=0; ix<nx; ix++ ){
            int i  = iy*nx + ix;
            int i3 = i*3;
            float dx = hsc*(height_map[wrap_index2d_fast(ix+1,iy  ,nx,ny)]-height_map[wrap_index2d_fast(ix-1,iy  ,nx,ny)]);
            float dy = hsc*(height_map[wrap_index2d_fast(ix  ,iy+1,nx,ny)]-height_map[wrap_index2d_fast(ix  ,iy-1,nx,ny)]);
            //printf( " (%i,%i)  (%f,%f) %f %\n", ix, iy, dx, dy, height_map[i] );
            //hxy[ i3   ] = (uint8_t)( _clamp(dx,-1.0,1.0)*127 + 128 ); // a
            //hxy[ i3+1 ] = (uint8_t)( _clamp(dy,-1.0,1.0)*127 + 128 );
            //hxy[ i3+2 ] = (uint8_t)( _clamp(height_map[ i ], 0.0,1.0 )*255 );

            hxy[ i3   ] = (uint8_t)( dx*127 + 128 ); // a
            hxy[ i3+1 ] = (uint8_t)( dy*127 + 128 );
            hxy[ i3+2 ] = (uint8_t)( height_map[ i ]*255 );

            //printf( " (%i,%i)  (%i,%i,%i) \n", ix,iy, hxy[i3], hxy[i3+1], hxy[i3+2] );

            //printf( " (%i,%i) (%f,%f)  (%f,%f) \n", ix,iy, dx, dy, 2.0*(hxy[i3]*inv255-0.5), 2.0*(hxy[i3+1]*inv255-0.5)  );

            //hxy[ i3   ] = (uint8_t)( _clamp(dy,-1.0,1.0)*127 + 128 ); // a
            //hxy[ i3+1 ] = (uint8_t)( _clamp(dx,-1.0,1.0)*127 + 128 );

            //hxy[ i3+2 ] = (uint8_t)( _clamp(dx,-1.0,1.0)*127 + 128 ); // a
            //hxy[ i3+1 ] = (uint8_t)( _clamp(dy,-1.0,1.0)*127 + 128 );
            //hxy[ i3+0 ] = (uint8_t)( _clamp(height_map[ i ], 0.0,1.0 )*255 );
        }
    }
    GLuint txHeight = 0;
    newTexture2D( txHeight, nx, ny, hxy, GL_RGB, GL_UNSIGNED_BYTE );

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);


    delete [] hxy;
    return txHeight;
}

int heightTextureFromHeightsDerivs_float( int nx, int ny, float* height_map, float hsc ){
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
    GLuint txHeight = 0;
    newTexture2D( txHeight, nx, ny, hxyz, GL_RGBA, GL_FLOAT ); // this does not seem to work
    delete [] hxyz;
    return txHeight;
}

// ===========================
// ====  TerrainOGL3Prototype
// ===========================

class TerrainOGL3Prototype{ public:
    Shader sh;
    Vec2i nSamp;

    Vec2f viewMin  = (Vec2f){ 1.0,   0.0 };
    Vec2f viewMax  = (Vec2f){ 0.0,   1.0 };

    double dist0;
    Vec3f pos      = (Vec3f){ 0.0,   0.0,    0.0 };
    Vec3f mapScale = (Vec3f){ 0.002, 0.002, 20.0 };
    Vec2f uv0      = (Vec2f){ 0.0,   0.0 };
    Vec2f txStep   = (Vec2f){ 1.0,   1.0 };
    float derivScale = 1.0;

    float flexibility = 0.0;

    GLuint txHeight;

    struct { GLuint
        uv0,
        mapScale,
        derivScale,
        txHeight,
        txStep,
        flexibility,
        lightColor,
        diffuseColor,
        ambientColor,
        specularColor,
        lightPos;
    } ulocs0;

    void getUlocs(){
        //ulocs.pa=sh.getUloc("pa");
        //ulocs.pb=sh.getUloc("pb");
        ulocs0.uv0           = sh.getUloc("uv0"          );
        ulocs0.mapScale      = sh.getUloc("mapScale"     );
        ulocs0.derivScale    = sh.getUloc("derivScale"   );
        ulocs0.txStep        = sh.getUloc("txStep"       );
        ulocs0.txHeight      = sh.getUloc("txHeight"     );
        ulocs0.flexibility   = sh.getUloc("flexibility"  );
        ulocs0.lightColor    = sh.getUloc("lightColor"   );
        ulocs0.diffuseColor  = sh.getUloc("diffuseColor" );
        ulocs0.ambientColor  = sh.getUloc("ambientColor" );
        ulocs0.specularColor = sh.getUloc("specularColor");
        ulocs0.lightPos      = sh.getUloc("lightPos"     );
    }

    void setDefaultColors(){
        glUniform3f(ulocs0.lightColor ,    0.50, 0.45, 0.40 );
        glUniform3f(ulocs0.diffuseColor,   1.00, 1.00, 1.00 );
        glUniform3f(ulocs0.ambientColor ,  0.10, 0.15, 0.20 );
        glUniform3f(ulocs0.specularColor , 1.00, 1.00, 1.00 );
        glUniform3f(ulocs0.lightPos,        0.0,+1000.0,0.0 );
    }

    void setViewRange( Vec2f camDir, float camTg ){
        Vec2f rot; rot.x=sqrt(1.0/(1.0+camTg*camTg)); rot.y=rot.x*camTg;
        viewMin.set_udiv_cmplx(camDir,rot);
        viewMax.set_mul_cmplx (camDir,rot);
        //printf(" (%f,%f), (%f,%f), (%f,%f) (%f,%f) \n", camDir.x,camDir.y,  viewMin.x,viewMin.y,   viewMax.x,viewMax.y,    rot.x,rot.y );
        //exit(0);
    }

    virtual void draw() = 0;

    virtual void draw(const Camera& cam){
        //cam.lookAt( Vec3dZ, 20.0 );
        sh.use();
        setCamera( sh, cam );
        draw();
    }

};

// ===========================
// ====  TerrainOGL3
// ===========================

class TerrainOGL3 : public TerrainOGL3Prototype { public:

    int nVertDrawn=0;
    int nDrawCalls=0;

    //GLuint txHeight;
    GLuint stripUVs;
    //GLuint hexUVs;  // central hexagon - todo later
    //GLuint hexInds
    struct { GLuint
        uv0,
        pa,pb;
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

    void init( Vec2i nSamp_, float dist0_, Vec2i nHeighs, float* height_map ){
        nSamp = nSamp_;
        dist0 = dist0_;
        sh.init( "common_resources/shaders/terrain_strip.glslv", "common_resources/shaders/terrain_world.glslf" );
        //sh.init( "common_resources/shaders/terrain_strip.glslv", "common_resources/shaders/color3D.glslf" );
        sh.use();
        sh.getDefaultUniformLocation();
        TerrainOGL3Prototype::getUlocs();

        ulocs.pa=sh.getUloc("pa");
        ulocs.pb=sh.getUloc("pb");

        printf( "ulocs : %i %i %i %i %i \n", ulocs.pa, ulocs.pb, ulocs0.uv0, ulocs0.mapScale, txHeight );
        setDefaultColors();

        //txHeight = heightTextureFromHeightsDerivs_byte ( nHeighs.x, nHeighs.y, height_map, 5.0 );
        txHeight = heightTextureFromHxy_byte ( nHeighs.x, nHeighs.y, height_map, 100.0 );
        //txHeight = heightTextureFromHeightsDerivs_float( nHeighs.x, nHeighs.y, height_map, 5.0 );
        //newTexture2D( txHeight, nHeighs.x, nHeighs.y, height_map, GL_RED, GL_FLOAT );

        //makeStrip( nSamp.b, dist0, 1.0 + 2.0/nSamp.a );
        makeStrip( nSamp.b, 1.0 + 2.0/nSamp.a );
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

    virtual void draw() override {
        nDrawCalls=0;
        nVertDrawn=0;
        //sh.use();
        glUniform2f ( ulocs0.uv0,        uv0.x, uv0.y );
        glUniform3fv( ulocs0.mapScale,1, (GLfloat*)&mapScale );
        glUniform1f ( ulocs0.derivScale, derivScale );
        glUniform1f ( ulocs0.flexibility, flexibility );
        glUniform2fv( ulocs0.txStep,  1,  (GLfloat*)&txStep );
        sh.set_modelPos( (GLfloat*)&pos );

        Vec2f drot; drot.fromAngle( M_PI/3.0 );
        Vec2f p = (Vec2f){dist0,0.0};
        bindVertexAttribPointer( 0, stripUVs, 2, GL_FLOAT, GL_FALSE );
        //bindTexture( 0, txHeight, ulocs.txHeight );
        bindTexture( 0, txHeight, ulocs0.txHeight  );
        for(int i=0; i<6; i++){
            Vec2f p_; p_.set_mul_cmplx(p,drot);
            //printf( " === %i (%f,%f) (%f,%f)\n", i, p.x, p.y,  p_.x, p_.y  );
            drawStrips( p, p_ );
            p = p_;
            #ifdef _DEBUG_VIEW_
            DEBUG_mesh->addLine( pos, {p.x,100.0,p.y},{1.0,0.0,1.0});
            #endif // _DEBUG_VIEW_
        }

        //exit(0);
    }

};

// ===========================
// ====  TerrainOGL3_patch
// ===========================

class TerrainOGL3_patch : public TerrainOGL3Prototype { public:

    //float sc  = 100.0;
    Vec2f dirHex  = (Vec2f){ 0.5,0.86602540378};
    Vec2f dirHex2 = (Vec2f){-0.5,0.86602540378};


    GLMesh *mesh=0;
    //GLuint hexUVs;  // central hexagon - todo later
    //GLuint hexInds
    struct { GLuint
        //p00,p01,p10,p11;
        ps;
    } ulocs;

    void init( Vec2i nSamp_, float dist0_, Vec2i nHeighs, float* height_map, float dhsc, bool wire ){
        nSamp = nSamp_;
        dist0 = dist0_;
        sh.init( "common_resources/shaders/terrain_patch.glslv", "common_resources/shaders/terrain_world.glslf" );
        //sh.init( "common_resources/shaders/terrain_strip.glslv", "common_resources/shaders/color3D.glslf" );

        //mesh = glHalfHexGrid( nSamp, true );
        mesh = glHalfHexGrid( nSamp, wire );

        sh.use();
        sh.getDefaultUniformLocation();

        TerrainOGL3Prototype::getUlocs();

        //ulocs.p00=sh.getUloc("p00");
        //ulocs.p01=sh.getUloc("p01");
        //ulocs.p10=sh.getUloc("p10");
        //ulocs.p11=sh.getUloc("p11");
        ulocs.ps=sh.getUloc("ps");

        //printf( "ulocs : ps(%i,%i,%i,%i) %i %i %i \n", ulocs.p11, ulocs.p01, ulocs.p10, ulocs.p11,  ulocs0.uv0, ulocs0.mapScale, txHeight );
        setDefaultColors();
        //txHeight = heightTextureFromHeightsDerivs_byte ( nHeighs.x, nHeighs.y, height_map, 5.0 );
        txHeight = heightTextureFromHxy_byte ( nHeighs.x, nHeighs.y, height_map, dhsc );
        //txHeight = heightTextureFromHeightsDerivs_float( nHeighs.x, nHeighs.y, height_map, 5.0 );
        //newTexture2D( txHeight, nHeighs.x, nHeighs.y, height_map, GL_RED, GL_FLOAT );
    }

    virtual void draw() override {
        //sh.use();
        glUniform2f ( ulocs0.uv0,         uv0.x, uv0.y );
        glUniform3fv( ulocs0.mapScale,1, (GLfloat*)&mapScale );
        glUniform1f ( ulocs0.derivScale, derivScale );
        glUniform1f ( ulocs0.flexibility, flexibility );
        glUniform2fv( ulocs0.txStep,1,   (GLfloat*)&txStep );
        sh.set_modelPos( (GLfloat*)&pos );

        //Vec2f drot; drot.fromAngle( M_PI/3.0 );
        //Vec2f p = (Vec2f){dist0,0.0};
        //bindVertexAttribPointer( 0, stripUVs, 2, GL_FLOAT, GL_FALSE );
        //bindTexture( 0, txHeight, ulocs.txHeight );
        bindTexture( 0, txHeight, ulocs0.txHeight  );


        //Vec2f rot  = Vec2fX;

        // TODO
        float fc  = nSamp.y/nSamp.x;
        float mfc = 1-fc;
        Vec2f drot = dirHex;
        Vec2f d    = dirHex*dist0;
        Vec2f d2   = d*(1+fc);


        Vec2f ps[4] = { -d.x, d.y,  d.x, d.y,   -d2.x, d2.y,  d2.x, d2.y  };

        // TODO : we should better use glUniform2fv();
        int nind = 0;
        if      ( mesh->draw_mode == GL_TRIANGLES ) nind = nSamp.x*(nSamp.x*2+nSamp.x  )*3;
        else if ( mesh->draw_mode == GL_LINES     ) nind = nSamp.x*(nSamp.x*2+nSamp.x+1)*3 + nSamp.x*2;
        int narg = mesh->preDraw();

        glBindBuffer  ( GL_ELEMENT_ARRAY_BUFFER, mesh->inds );



        for(int i=0; i<6; i++){
            for(int j=0; j<4; j++ ){  ps[j].mul_cmplx(drot); };

            float side1 = viewMin.cross(ps[0]);
            float side2 = viewMax.cross(ps[1]);
            if( (side1<0)||(side2>0) ) continue;

            //float side1 = viewMin.cross(ps[1]);
            //float side2 = viewMax.cross(ps[0]);
            //if( (side1>0)||(side2<0) ) continue;  // frustrum inside block
            //if( (side1>0)&&(side2>0) ) continue;
            //glUniform2f( ulocs.p00,  ps[0].x, ps[0].y );
            //glUniform2f( ulocs.p01,  ps[1].x, ps[1].y );
            //glUniform2f( ulocs.p10,  ps[2].x, ps[2].y );
            //glUniform2f( ulocs.p11,  ps[3].x, ps[3].y );
            glUniform2fv( ulocs.ps, 4, (const GLfloat*)ps );
            //drawElements( mesh->draw_mode, mesh->inds, nSamp.x*nfac );
            glDrawElements( mesh->draw_mode, mesh->nInds, GL_UNSIGNED_INT, (void*)0 );
            //mesh->draw();
            //rot.mul_cmplx(drot);
        }

        //glUniform2fv( ulocs.ps, 4, (const GLfloat*){ -d.x,d.y, d.x,d.y, -d2.x,d.y*mfc, d2.x,d.y*mfc } );
        //glUniform2f( ulocs.p00, -d.x,  d.y     );
        //glUniform2f( ulocs.p01,  d.x,  d.y     );
        //glUniform2f( ulocs.p10, -d2.x, d.y*mfc );
        //glUniform2f( ulocs.p11,  d2.x, d.y*mfc );
        ps[0].set(-d.x,d.y); ps[1].set(d.x,d.y); ps[2].set(-d2.x,d.y*mfc); ps[3].set(d2.x, d.y*mfc);
        glUniform2fv( ulocs.ps, 4, (const GLfloat*)ps );
        //mesh->draw();
        //drawElements( mesh->draw_mode, mesh->inds, nSamp.x*nfac );
        glDrawElements( mesh->draw_mode, nind, GL_UNSIGNED_INT, (void*)0 );

        //glUniform2f( ulocs.p00, -d.x,  -d.y     );
        //glUniform2f( ulocs.p01,  d.x,  -d.y     );
        //glUniform2f( ulocs.p10, -d2.x, -d.y*mfc );
        //glUniform2f( ulocs.p11,  d2.x, -d.y*mfc );
        ps[0].set(-d.x,-d.y); ps[1].set(d.x,-d.y); ps[2].set(-d2.x,-d.y*mfc); ps[3].set(d2.x, -d.y*mfc);
        glUniform2fv( ulocs.ps, 4, (const GLfloat*)ps );
        //drawElements( mesh->draw_mode, mesh->inds, nSamp.x*nfac );
        glDrawElements( mesh->draw_mode, nind, GL_UNSIGNED_INT, (void*)0 );

        mesh->postDraw(narg);

    }

    virtual void draw(const Camera& cam) override{
        //float cab = dirHex.dot(dirHex2);
        float f = dist0/nSamp.x;
        Vec2f p2  = cam.pos.xz();
        //float ca  = p2.dot(dirHex )/f; ca=round(ca);
        //float cb  = p2.dot(dirHex2)/f; cb=round(cb);
        Vec2f cs; cs.fromLinearSolution( dirHex, dirHex2, p2*(1/f) );
        cs.a=round(cs.a); cs.b=round(cs.b);
        Vec2f p2_ = ( dirHex*cs.a + dirHex2*cs.b )*f;
        //Vec2f p2_ = dirHex*cs.a + dirHex2*cs.b;
        //printf( " (%f,%f) -> (%f,%f)   |  (%f,%f) \n",  p2.x, p2.y, p2_.x, p2_.y,   cs.a, cs.b );
        pos.x=p2_.x;
        pos.z=p2_.y;

        Vec2f camDir = (Vec2f){cam.rot.c.x,cam.rot.c.z};
        camDir.normalize();
        setViewRange( camDir, cam.getTgX() );
        //printf( "camDir (%f,%f) vieMin (%f,%f)\n", camDir.x, camDir.y,    viewMin.x, viewMin.y  );

        //TerrainOGL3Prototype::draw(cam);

        sh.use();
        setCamera( sh, cam );
        draw();
    }

};
















// ===========================
// ====  TerrainOGL3_patch
// ===========================

class TerrainOGL3_normals : public TerrainOGL3Prototype { public:

    GLMesh *mesh=0;
    //GLuint hexUVs;  // central hexagon - todo later
    //GLuint hexInds
    struct { GLuint
        //p00,p01,p10,p11;
        ps;
    } ulocs;

    void init( Vec2i nSamp_, float dist0_, Vec2i nHeighs, float* height_map, float dhsc ){
        nSamp = nSamp_;
        dist0 = dist0_;
        //sh.init( "common_resources/shaders/terrain_normals.glslv", "common_resources/shaders/terrain_world.glslf" );
        //sh.init( "common_resources/shaders/terrain_normals.glslv", "common_resources/shaders/const3D.glslf" );
        sh.init( "common_resources/shaders/terrain_normals.glslv", "common_resources/shaders/color3D.glslf" );
        //sh.init( "common_resources/shaders/terrain_strip.glslv", "common_resources/shaders/color3D.glslf" );
        //mesh = glHalfHexGrid( nSamp, true );
        mesh = glNormalGrid( nSamp );
        sh.use();
        sh.getDefaultUniformLocation();
        TerrainOGL3Prototype::getUlocs();
        ulocs.ps=sh.getUloc("ps");
        setDefaultColors();
        txHeight = heightTextureFromHxy_byte ( nHeighs.x, nHeighs.y, height_map, dhsc );
    }

    virtual void draw() override {
        //sh.use();
        glUniform2f ( ulocs0.uv0,         uv0.x, uv0.y );
        glUniform3fv( ulocs0.mapScale,1, (GLfloat*)&mapScale );
        glUniform1f ( ulocs0.derivScale, derivScale );
        glUniform1f ( ulocs0.flexibility, flexibility );
        glUniform2fv( ulocs0.txStep,1,   (GLfloat*)&txStep );
        sh.set_modelPos( (GLfloat*)&pos );

        bindTexture( 0, txHeight, ulocs0.txHeight  );
        Vec2f ps[4];
        ps[0].set(-dist0,-dist0);
        ps[1].set(+dist0,-dist0);
        ps[2].set(-dist0,dist0);
        ps[3].set(+dist0,dist0);
        glUniform2fv( ulocs.ps, 4, (const GLfloat*)ps );

        mesh->draw();

    }

    virtual void draw(const Camera& cam) override{
        pos.x=cam.pos.x;
        pos.z=cam.pos.z;
        sh.use();
        setCamera( sh, cam );
        draw();
    }

};

#endif
