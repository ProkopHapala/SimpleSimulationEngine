
#ifndef  GL3Utils_h
#define  GL3Utils_h

#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>
#include "Vec3.h"
#include "GLObject.h"
#include "GLUtils.h"

GLObject * qaudPatchHard( int n, Vec2d p0, Vec2d da, Vec2d db, Vec3d (vertFunc)(Vec2d) ){

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

/*
void triPatchHard( int n, Vec2d p0, Vec2d da, Vec2d db, Vec3d (vertFunc)(Vec2d), GLfloat* verts, GLfloat* normals ){
    Vec3f * verti = (Vec3f*)verts;
    Vec3f * normi = (Vec3f*)normals;
    double d = 1.0d/(n);
    int nVert = 3*n*n;
    for(int i=0; i<n; i++){
        uv.x    = p0 + da*i;
        Vec3d oa = vertFunc( uv );
        Vec3d ob = vertFunc( uv );
        for(int j=0; j<n-i; j++){
            uv.add(da);
            Vec3d a  = vertFunc( uv+da );
            Vec3d n; n.set_cross(a-oa,b-oa);
            n.normalize();
            convert(a,verti[0]); convert(b,verti[1]); convert(a,verti[2]);
            convert(n,normi[0]); convert(n,normi[1]); convert(n,normi[2]);
            if( j<i ){
                Vec3d b  = vertFunc( uv+db );
                Vec3d n; n.set_cross(a-oa,b-oa);
                n.normalize();
                convert(a,verti[0]); convert(b,verti[1]); convert(a,verti[2]);
                convert(n,normi[0]); convert(n,normi[1]); convert(n,normi[2]);
                oa=a; ob=b;
                verti+=3;
                normi+=3;
            }
        }
        facei += ni;
    }
}
*/


#endif
