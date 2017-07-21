
#ifndef  GL3Utils_h
#define  GL3Utils_h

//#include <GL/glew.h>

#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>
#include "Vec3.h"
#include "GLObject.h"
#include "GLUtils.h"

#include "CMesh.h"

//GLObject * qaudPatchHard( int n, Vec2d p0, Vec2d da, Vec2d db, Vec3d (vertFunc)(Vec2d) ){


GLObject * makeOgl_flat( const CMesh& mesh ){
    GLObject * ogl = new GLObject();
    ogl->setup( countVerts( mesh.nfaces, mesh.ngons ) );
    hardFace( mesh.nfaces, mesh.ngons, mesh.faces, mesh.verts, ogl->buffs[0].cbuff, ogl->buffs[1].cbuff );
    ogl->init();
    return ogl;
}

template<typename Func>
GLObject * qaudPatchHard( int n, Vec2d p0, Vec2d da, Vec2d db, Func vertFunc ){
    //GLObject * ogl = new GLObject( 3*2*n*n );
    GLObject * ogl = new GLObject();
    ogl->setup( 3*2*n*n );
    Vec3f * verti = (Vec3f*)ogl->buffs[0].cbuff;
    Vec3f * normi = (Vec3f*)ogl->buffs[1].cbuff;
    double d = 1.0d/(n);
    Vec3d p00,p01,p10,p11,nv;
    for(int i=0; i<n; i++){
        Vec2d uv  = p0 + da*i;
        p00 = vertFunc( uv );
        p01 = vertFunc( uv+da );
        for(int j=0; j<n; j++){
            uv.add(db);
            p10 = vertFunc( uv    );
            p11 = vertFunc( uv+da );
            nv.set_cross(p10-p00,p01-p00); nv.normalize();
            convert(p00,verti[0]); convert(p01,verti[1]); convert(p10,verti[2]);
            convert(nv ,normi[0]); convert(nv ,normi[1]); convert(nv ,normi[2]);
            nv.set_cross(p01-p11,p10-p11); nv.normalize();
            convert(p11,verti[3]); convert(p01,verti[4]); convert(p10,verti[5]);
            convert(nv ,normi[3]); convert(nv ,normi[4]); convert(nv ,normi[5]);
            p00=p10; p01=p11;
            verti+=6;
            normi+=6;
        }
    }
    ogl->init();
    return ogl;
}

template<typename Func>
GLObject * makeNVerts( int n, Func vertFunc ){
    //GLObject * ogl = new GLObject( 3*2*n*n );
    GLObject * ogl = new GLObject();
    ogl->setup( n );
    Vec3f * verti = (Vec3f*)ogl->buffs[0].cbuff;
    Vec3f * normi = (Vec3f*)ogl->buffs[1].cbuff;
    for(int i=0; i<n; i++){
        Vec3d p,nv;
        vertFunc( i, p, nv );
        convert(p,verti[i]); convert(nv,normi[i]);
    }
    ogl->init();
    return ogl;
}


template<typename Func>
GLObject * makeNTris( int n, Func vertFunc ){
    //GLObject * ogl = new GLObject( 3*2*n*n );
    GLObject * ogl = new GLObject();
    ogl->setup( n*3 );
    Vec3f * verti = (Vec3f*)ogl->buffs[0].cbuff;
    Vec3f * normi = (Vec3f*)ogl->buffs[1].cbuff;
    int ii = 0;
    for(int i=0; i<n; i++){
        Vec3d p,nv;
        vertFunc( i, 0, p, nv ); convert(p,verti[ii]); convert(nv,normi[ii]); ii++;
        vertFunc( i, 1, p, nv ); convert(p,verti[ii]); convert(nv,normi[ii]); ii++;
        vertFunc( i, 2, p, nv ); convert(p,verti[ii]); convert(nv,normi[ii]); ii++;
    }
    ogl->init();
    return ogl;
}


#endif
