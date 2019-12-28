
#ifndef  GL3Utils_h
#define  GL3Utils_h

#include <GL/glew.h>
//#define GL_GLEXT_PROTOTYPES
//#include <GL/gl.h>

#include "Vec2.h"
#include "Vec3.h"
#include "GLObject.h"
#include "GLobjects.h"
#include "GLUtils.h"
#include "CMesh.h"
#include "Camera.h"

//GLObject * qaudPatchHard( int n, Vec2d p0, Vec2d da, Vec2d db, Vec3d (vertFunc)(Vec2d) ){

GLMesh* cam2mesh( Camera& cam ){
    GLMesh* msh = new GLMesh();
    msh->draw_mode = GL_LINES;

    printf( "zoom %g aspect %g zmin %g zmax %g\n", cam.zoom, cam.aspect, cam.zmin, cam.zmax  );

    //float tgx = 1.0/(cam.zoom*cam.aspect);
    //float tgy = 1.0/(cam.zoom);
    float tgx  = cam.getTgX();
    float tgy  = cam.getTgY();
    float zmin = cam.zmin;
    float zmax = cam.zmax;
    //float tgx = 0.5;
    //float tgy = 0.5;
    //float zmin = 2.0;
    //float zmax = 6.0;
    float xzmin = tgx*zmin;
    float xzmax = tgx*zmax;
    float yzmin = tgy*zmin;
    float yzmax = tgy*zmax;
    float verts[] = {
        -xzmin,-yzmin,zmin,
        -xzmax,-yzmax,zmax,
        -xzmin,+yzmin,zmin,
        -xzmax,+yzmax,zmax,
        +xzmin,-yzmin,zmin,
        +xzmax,-yzmax,zmax,
        +xzmin,+yzmin,zmin,
        +xzmax,+yzmax,zmax
    };
    int edges[] = { 0,1, 0,2, 0,4,  1,3,1,5, 2,3, 2,6,    7,5, 7,6, 7,3,  5,4, 6,4    };
    msh->init( 8, 12*2, edges, verts, NULL, NULL, NULL );
    return msh;
}

GLMesh* makeQuad3D( Vec2f p0, Vec2f p1, Vec2f u0, Vec2f u1 ){
    //GLfloat verts = new GLfloat[3*2*3];
    //GLfloat vUVs  = new GLfloat[3*2*2];
    //delete [] verts;
    //GLfloat verts[] = { p0.x,p0.y,0.0, p1.x,p1.y,0.0, p0.x,p1.y,0.0,   p0.x,p0.y,0.0, p1.x,p1.y,0.0, p1.x,p0.y,0.0 };
    //Vec2f vUVs   [] = { u0.x,u0.y,     u1.x,u1.y,     u0.x,u1.y,       u0.x,u0.y,     u1.x,u1.y,     u1.x,u0.y     };
    //GLfloat verts[] = { p0.x,p0.y,0.0, p1.x,p1.y,0.0, p0.x,p1.y,0.0,   p0.x,p0.y,0.0, p1.x,p0.y,0.0,  p1.x,p1.y,0.0  };
    //Vec2f vUVs   [] = { u0.x,u0.y,     u1.x,u1.y,     u0.x,u1.y,       u0.x,u0.y,     u1.x,u0.y,      u1.x,u1.y      };
    GLfloat verts[] = { p0.x,p0.y,0.0, p0.x,p1.y,0.0, p1.x,p1.y,0.0,    p0.x,p0.y,0.0, p1.x,p1.y,0.0, p1.x,p0.y,0.0 };
    Vec2f vUVs   [] = { u0.x,u0.y,     u0.x,u1.y,     u1.x,u1.y,        u0.x,u0.y,     u1.x,u1.y,     u1.x,u0.y     };
    GLMesh* glquad = new GLMesh();
    glquad->init( 6, 0, NULL, verts, NULL, NULL, vUVs );
    return glquad;
}

template<typename PixFunc>
GLuint makeTestTextureRGBA( Vec2i sz, PixFunc pixFunc ){
    double dx = 1.0d/sz.x;
    double dy = 1.0d/sz.y;
    uint32_t * c_img1 = new uint32_t[sz.x*sz.y];
    for( int iy=0; iy<sz.y; iy++ ){
        for( int ix=0; ix<sz.x; ix++ ){
            uint8_t r,g,b,a;
            pixFunc( (Vec2i){ix,iy}, r,g,b,a );
            //r=ix; g=ix^iy; b=iy; a=255;
            c_img1[ iy*sz.x + ix ] = (a<<24) | (b<<16) | (g<<8) | (r);
        }
    }
    GLuint texID;
    newTexture2D( texID, sz.x, sz.y, c_img1, GL_RGBA, GL_UNSIGNED_BYTE );
    return texID;
}

GLuint makeTestTextureRGBA_1( Vec2i sz ){
    return makeTestTextureRGBA( sz, [](Vec2i uv, uint8_t& r,uint8_t& g,uint8_t& b,uint8_t& a ){ r=uv.x; g=uv.x^uv.y; b=uv.y; a=255; } );
}


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
GLMesh * qaudPatchSmooth( Vec2i n, Vec2f p0, Vec2f da, Vec2f db, Func vertFunc, void**cbuffs=NULL ){
    int nVerts = n.a*n.b;
    int nTris  =(n.a-1)*(n.b-1)*2;
    Vec3f * vpos = new Vec3f[nVerts];
    Vec3f * vnor = new Vec3f[nVerts];
    Vec3i * tris = new Vec3i[nTris ];
    int iv = 0;
    int ii = 0;
    //printf( "(%i,%i) (%i,%i) (%g,%g)(%g,%g)(%g,%g)\n", n.a,n.b, nVerts, nInds , p0.x,p0.y, da.x,da.y, db.x,db.y );
    for(int ia=0; ia<n.a; ia++){
        for(int ib=0; ib<n.b; ib++){
            Vec2f p = p0 + da*ia + db*ib;
            vertFunc( p, vpos[iv], vnor[iv] );
            //printf( "%i %i %i (%g,%g) (%g,%g,%g) \n", ia, ib, iv, p.x, p.y, vpos[iv].x, vpos[iv].y, vpos[iv].z );
            if(ia<(n.a-1)&&(ib<(n.b-1))){ tris[ii].set(iv,iv+1,iv+n.b); ii++; tris[ii].set(iv+n.b,iv+1,iv+n.b+1); ii++; }
            iv++;
        }
    }
    GLMesh * glmesh = new GLMesh();
    glmesh->init( nVerts, nTris*3, tris, vpos, vnor, NULL, NULL );
    if(cbuffs){ cbuffs[0]=tris; cbuffs[1]=vpos; cbuffs[2]=vnor; }
    else{ delete [] vpos; delete [] vnor; delete [] tris; }
    return glmesh;
}



template<typename Func>
GLMesh * qaudPatchSmooth( Vec2i n, Vec2f span, int pattern, Func vertFunc, void**cbuffs=NULL ){
    int nVerts=0, nTris=0;
    switch(pattern){
        case 0: nTris=(n.a-1)*(n.b-1)*2; nVerts=n.a*n.b;       break;
        case 1: nTris=(n.a-1)*(n.b-1)*4; nVerts=n.a*(2*n.b-1); break;
    }
    Vec3f * vpos = new Vec3f[nVerts];
    Vec3f * vnor = new Vec3f[nVerts];
    Vec3i * tris = new Vec3i[nTris ];
    int iq = 0; int iv = 0; int ii = 0;
    float da = span.a/(n.a-1);
    float db = span.b/(n.b-1);
    //printf( "(%i,%i) (%i,%i) (%g,%g)(%g,%g)(%g,%g)\n", n.a,n.b, nVerts, nInds , p0.x,p0.y, da.x,da.y, db.x,db.y );
    for(int ia=0; ia<n.a; ia++){
        for(int ib=0; ib<n.b; ib++){
            Vec2f p = { da*ia, db*ib };
            vertFunc( p, vpos[iv], vnor[iv] ); iv++;
            if(ia<(n.a-1)&&(ib<(n.b-1))){
                switch(pattern){
                    case 0: tris[ii].set(iv-1,iv,iv+n.b-1); ii++; tris[ii].set(iv+n.b-1,iv,iv+n.b); ii++; break;
                    case 1:
                        p.add(0.5f*da,0.5f*db);
                        vertFunc( p, vpos[iv], vnor[iv] );
                        int iv2 = (n.b<<1)-1 + iv;
                        //if(ia<(n.a-2)){ iv2+=iv; }else{ iv2+=(iv>>1); }
                        if(ia==(n.a-2)){ iv2-=ib; if(ib==(n.b-2)) iv2--; };
                        tris[ii].set(iv,iv -1,iv +1); ii++;
                        tris[ii].set(iv,iv -1,iv2-1); ii++;
                        tris[ii].set(iv,iv2+1,iv2-1); ii++;
                        tris[ii].set(iv,iv2+1,iv +1); ii++;
                        iv++;
                        break;
                }
            }
            //printf( "%i %i %i (%g,%g) (%g,%g,%g) \n", ia, ib, iv, p.x, p.y, vpos[iv].x, vpos[iv].y, vpos[iv].z );
            iq++;
        }
    }
    GLMesh * glmesh = new GLMesh();
    glmesh->init( nVerts, nTris*3, tris, vpos, vnor, NULL, NULL );
    if(cbuffs){ cbuffs[0]=tris; cbuffs[1]=vpos; cbuffs[2]=vnor; }
    else{ delete [] vpos; delete [] vnor; delete [] tris; }
    return glmesh;
}

GLMesh * qaudPatchUV( Vec2i n, void**cbuffs=NULL ){
    int nVerts = n.a*n.b;
    int nTris  =(n.a-1)*(n.b-1)*2;
    Vec2f * vUVs = new Vec2f[nVerts];
    Vec3i * tris = new Vec3i[nTris ];
    int iv = 0;
    int ii = 0;
    float da = 1.0f/(n.a-1);
    float db = 1.0f/(n.b-1);
    //printf( "(%i,%i) (%i,%i) (%g,%g)(%g,%g)(%g,%g)\n", n.a,n.b, nVerts, nInds , p0.x,p0.y, da.x,da.y, db.x,db.y );
    for(int ia=0; ia<n.a; ia++){
        for(int ib=0; ib<n.b; ib++){
            vUVs[iv] = { da*ia, db*ib };
            //printf( "%i %i %i (%g,%g) (%g,%g,%g) \n", ia, ib, iv, p.x, p.y, vpos[iv].x, vpos[iv].y, vpos[iv].z );
            if(ia<(n.a-1)&&(ib<(n.b-1))){ tris[ii].set(iv,iv+1,iv+n.b); ii++; tris[ii].set(iv+n.b,iv+1,iv+n.b+1); ii++; }
            iv++;
        }
    }
    GLMesh * glmesh = new GLMesh();
    glmesh->init( nVerts, nTris*3, tris, NULL, NULL, NULL, vUVs );
    if(cbuffs){ cbuffs[0]=tris; cbuffs[1]=vUVs; }
    else{ delete [] tris; delete [] vUVs; }
    return glmesh;
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

GLMesh* glNormalGrid( Vec2i ns ){
    int nverts = ns.x*ns.y*2;
    //int nquads = ns.x*ns.y;
    float dx = 1.0/ns.x;
    float dy = 1.0/ns.y;
    Vec3f* verts = new Vec3f[nverts  ];
    //Vec3i* inds  = new Vec3i[nquads*2];
    int i=0;
    for(int iy=0; iy<ns.y; iy++){
        for(int ix=0; ix<ns.x; ix++){
            verts[i].set( ix*dx, iy*dy, 0.0 ); i++;
            verts[i].set( ix*dx, iy*dy, 1.0 ); i++;
        }
    }
    //printf( "i %i \n", i, nquads );
    GLMesh* mesh=new GLMesh();
    mesh->draw_mode = GL_LINES;
    mesh->init( nverts, 0, NULL, (GLfloat*)verts, NULL, NULL, NULL);
    delete [] verts;
    return mesh;
}

GLMesh* glQuadGrid( Vec2i ns, bool wire ){
    int nverts = (ns.x+1)*(ns.y+1);
    //int nquads = ns.x*ns.y;
    float dx = 1.0/ns.x;
    float dy = 1.0/ns.y;
    Vec3f* verts = new Vec3f[nverts  ];
    //Vec3i* inds  = new Vec3i[nquads*2];
    int i=0;
    for(int iy=0; iy<=ns.y; iy++){
        for(int ix=0; ix<=ns.x; ix++){
            verts[i].set( ix*dx, iy*dy, 0.0 );
            //printf( "%i %g %g %g\n", i, verts[i].x, verts[i].y, verts[i].z);
            i++;
        }
    }

    // indexes
    i=0;
    GLMesh* mesh=new GLMesh();
    int ninds=0;
    GLuint* inds = 0;
    if(wire){ // Lines
        mesh->draw_mode = GL_LINES;
        ninds = (ns.x*ns.y*3 + ns.x + ns.y)*2;
        //ninds = ( ns.x + ns.y)*2;
        //ninds = ns.x*2;
        //ninds = ns.y*2;
        inds  = new GLuint[ninds];
        Vec2i* edges = (Vec2i*)inds;

        for(int iy=0; iy<ns.y; iy++){
            int i1 = iy*(ns.x+1);
            int i2 = i1+(ns.x+1);
            //printf( "i %i nx %i i1 %i i2 %i \n", iy, nx, i1, i2 );
            for(int ix=0; ix<ns.x; ix++){
                edges[i].set( i1+ix, i1+ix+1 ); i++;
                edges[i].set( i1+ix, i2+ix   ); i++;
                edges[i].set( i1+ix, i2+ix+1 ); i++;
            }
        }
        int i0=ns.y*(ns.x+1);
        for(int ix=0; ix<ns.x; ix++){ edges[i].set( i0+ix, i0+ix+1 ); i++; };
        for(int iy=0; iy<ns.y; iy++){ edges[i].set( iy*(ns.x+1)+ns.x, (iy+1)*(ns.x+1)+ns.x ); i++; };
        //printf( "ninds_ %i ninds %i \n",      ninds_,   ninds );
    }else{    // Triangles
        mesh->draw_mode = GL_TRIANGLES;
        ninds = ns.x*ns.y*2*3;
        inds  = new GLuint[ninds];
        Vec3i* tris = (Vec3i*)inds;
        for(int iy=0; iy<ns.y; iy++){
            int i1 = iy*(ns.x+1);
            int i2 = i1+(ns.x+1);
            for(int ix=0; ix<ns.x; ix++){
                tris[i  ].set( i1+ix, i2+ix,   i2+ix+1 );
                tris[i+1].set( i1+ix, i1+ix+1, i2+ix+1 );
                i+=2;
            }
        }
    }

    //printf( "i %i \n", i, nquads );
    mesh->init( nverts, ninds, inds, (GLfloat*)verts, NULL, NULL, NULL);
    delete [] verts;
    delete [] inds;
    return mesh;
}


GLMesh* glTriangleGrid( int n, bool wire ){
    int nverts = (n+1)*n;
    //int nquads = n*n;
    float d = 1.0/n;
    Vec3f* verts = new Vec3f[nverts];
    //Vec3i* inds  = new Vec3i[nquads];
    int i=0;
    for(int iy=0; iy<n+1; iy++){
        for(int ix=0; ix<=(n-iy); ix++){
            verts[i].set( ix*d, iy*d, 0.0 );
            //printf( "%i %g %g %g\n", i, verts[i].x, verts[i].y, verts[i].z);
            i++;
        }
        //printf( "iy %i i %i \n", iy, i );
    }
    // indexes
    i=0;
    GLMesh* mesh=new GLMesh();
    int ninds=0;
    GLuint* inds = 0;
    if(wire){ // Lines
        mesh->draw_mode = GL_LINES;
        ninds = n*(n+1)*3;
        inds  = new GLuint[ninds];
        Vec2i* edges = (Vec2i*)inds;
        int i1 = 0;
        for(int iy=0; iy<n; iy++){
            int nx = (n-iy);
            int i2 = i1+(nx+1);
            //printf( "i %i nx %i i1 %i i2 %i \n", iy, nx, i1, i2 );
            for(int ix=0; ix<nx; ix++){
                edges[i].set( i1+ix,   i1+ix+1 ); i++;
                edges[i].set( i1+ix+1, i2+ix   ); i++;
                edges[i].set( i2+ix,   i1+ix   ); i++;
            }
            i1 = i2;
        }
        //printf( "ninds_ %i ninds %i \n",      ninds_,   ninds );
    }else{    // Triangles
        mesh->draw_mode = GL_TRIANGLES;
        ninds = n*n*3;
        inds  = new GLuint[ninds];
        Vec3i* tris = (Vec3i*)inds;
        int i1 = 0;
        for(int iy=0; iy<n; iy++){
            int nx = (n-iy);
            int i2 = i1+(nx+1);
            //printf( "i %i nx %i i1 %i i2 %i \n", iy, nx, i1, i2 );
            for(int ix=0; ix<nx; ix++){
                tris[i].set( i1+ix, i1+ix+1, i2+ix ); i++;
                if(ix==(nx-1))continue;
                tris[i].set( i1+ix+1, i2+ix, i2+ix+1 ); i++;
            }
            i1 = i2;
        }
    }

    //printf( "i %i \n", i, nquads );
    mesh->init( nverts, ninds, inds, (GLfloat*)verts, NULL, NULL, NULL);
    delete [] verts;
    delete [] inds;
    return mesh;
}

GLMesh* glHalfHexGrid( Vec2i ns, bool wire ){
    int nverts = (ns.y+1)*(ns.x+1+ns.y);
    float dy = 1.0/ns.y;
    Vec3f* verts = new Vec3f[nverts];
    int i=0;
    for(int iy=0; iy<ns.y+1; iy++){
        int nx = ns.x + iy;
        float dx = 1.0/nx;
        for(int ix=0; ix<(nx+1); ix++){
            verts[i].set( ix*dx, iy*dy, 0.0 );
            //printf( "%i %g %g %g\n", i, verts[i].x, verts[i].y, verts[i].z);
            i++;
        }
        //printf( "iy %i i %i \n", iy, i );
    }
    // indexes
    i=0;
    GLMesh* mesh=new GLMesh();
    int ninds=0;
    GLuint* inds = 0;
    if(wire){ // Lines
        mesh->draw_mode = GL_LINES;
        ninds = ns.y*(ns.x*2+ns.y+1)*3 + ns.x*2;
        inds  = new GLuint[ninds];
        Vec2i* edges = (Vec2i*)inds;
        int i1 = 0;
        //int ninds_=0;
        for(int ix=0; ix<ns.x; ix++){ edges[i].set( ix, ix+1 ); i++; };
        for(int iy=0; iy<ns.y; iy++){
            int nx = ns.x + iy;
            int i2 = i1+(nx+1);
            for(int ix=0; ix<=nx; ix++){
                edges[i].set( i1+ix, i2+ix   ); i++;
                edges[i].set( i1+ix, i2+ix+1 ); i++;
                edges[i].set( i2+ix, i2+ix+1 ); i++;
            }
            i1 = i2;
        }
        //printf( "ninds_ %i ninds %i \n",      ninds_,   ninds );
    }else{    // Triangles
        mesh->draw_mode = GL_TRIANGLES;
        ninds = ns.y*(ns.x*2+ns.y)*3;
        inds  = new GLuint[ninds];
        Vec3i* tris = (Vec3i*)inds;
        int i1 = 0;
        for(int iy=0; iy<ns.y; iy++){
            int nx = ns.x + iy;
            int i2 = i1+(nx+1);
            for(int ix=0; ix<=nx; ix++){
                tris[i].set( i1+ix, i2+ix, i2+ix+1 ); i++;
                if(ix==(nx))continue;
                tris[i].set( i1+ix, i1+ix+1, i2+ix+1 ); i++;
            }
            i1 = i2;
        }
    }
    //printf( "i %i \n", i, nquads );

    mesh->init( nverts, ninds, inds, (GLfloat*)verts, NULL, NULL, NULL);
    delete [] verts;
    delete [] inds;
    return mesh;
}


#endif
