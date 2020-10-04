
#ifndef  DrawOGL3_h
#define  DrawOGL3_h

#include <vector>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"

#include <GL/glew.h>
#include "Shader.h"
//#include "GLObject.h"
#include "GLobjects.h"

#include "CMesh.h"
#include "UVfuncs.h"

Vec3f (*UVfunc)(Vec2f p);

/*
inline void push3f( std::vector<float>& vs, const Vec3f& v ){ vs.push_back(v.x); vs.push_back(v.y); vs.push_back(v.z); }
inline void push2f( std::vector<float>& vs, const Vec2f& v ){ vs.push_back(v.x); vs.push_back(v.y); }
*/
inline void push3i( std::vector<int>&   vs, const Vec3i& v ){ vs.push_back(v.x); vs.push_back(v.y); vs.push_back(v.z); }
inline void push2i( std::vector<int>&   vs, const Vec2i& v ){ vs.push_back(v.x); vs.push_back(v.y); }

GLMesh* vecs2mesh( int n, Vec3f* ps, Vec3f* ds, float sc ){
    //int n    = vpos.size();
    Vec3f* vs =new Vec3f[n*2];
    for(int i=0; i<n; i++){ int i2=i*2; vs[i2]=ps[i]; vs[i2+1]=ps[i]+ds[i]*sc; };
    GLMesh* mesh=new GLMesh();
    mesh->init( n*2, 0, NULL, (GLfloat*)vs, NULL, NULL, NULL);
    mesh->draw_mode = GL_LINES;
    delete [] vs;
    return mesh;
}

GLMesh* hardTriangles2mesh( const CMesh& msh ){
    Vec3f *vs = new Vec3f[msh.ntri*3];
    Vec3f *ns = new Vec3f[msh.ntri*3];
    for(int i=0; i<msh.ntri; i++){
        int i3=i*3;
        Vec3i iv= msh.tris[i];
        Vec3f a,b,c,nor;
        convert( msh.verts[iv.a], a );
        convert( msh.verts[iv.b], b );
        convert( msh.verts[iv.c], c );
        nor.set_cross(a-c,b-c); nor.normalize();
        vs[i3  ]=a; ns[i3  ]=nor;
        vs[i3+1]=b; ns[i3+1]=nor;
        vs[i3+2]=c; ns[i3+2]=nor;
    };
    GLMesh* glmesh=new GLMesh();
    glmesh->init( msh.ntri*3, 0, NULL, vs, ns, NULL, NULL );
    glmesh->draw_mode = GL_TRIANGLES;
    delete [] vs;
    delete [] ns;
    return glmesh;
}

//hardMesh( int n, Vec3i* tris, Vec3f* ps ){}


/*
class LineMeshBuilder{ public:
    std::vector<Vec3f> vpos;
    std::vector<Vec3f> vpos;
}
*/


GLMesh* polyLineMesh( int n, float* ps, GLMesh* mesh=NULL ){
    if(mesh==0) mesh=new GLMesh();
    //Vec3f *pnor=NULL;
    //if( bnor ) pnor=&vnor[0];
    mesh->init( n, 0, 0, ps, 0, NULL, NULL );
    mesh->draw_mode = GL_LINE_STRIP;
    return mesh;
}









class GLMeshBuilder{ public:
    bool bnor = true;
    bool bUVs = true;
    std::vector<Vec2i> subs; //=std::vector<Vec2i>({(Vec2i){0,0}});
    std::vector<Vec3f> vpos;
    std::vector<Vec3f> vnor;
    //std::vector<Vec3f> vcol;
    std::vector<Vec2f> vUVs;
    std::vector<int>   inds;
    GLenum draw_mode = GL_TRIANGLES;

    void clear     (){ vpos.clear(); vnor.clear(); vUVs.clear(); inds.clear(); }
    GLMesh* makeGLmesh(){
        GLMesh* mesh=new GLMesh();
        Vec3f *pnor=NULL;Vec2f *pUV=NULL;
        if( bnor ) pnor=&vnor[0];
        if( bUVs ) pUV =&vUVs[0];
        mesh->init( vpos.size(), inds.size(), &inds[0],&vpos[0],pnor, NULL, pUV);
        mesh->draw_mode = draw_mode;
        return mesh;
    }

    GLMesh* makeLineMesh(){
        GLMesh* mesh=new GLMesh();
        Vec3f *pnor=NULL;
        if( bnor ) pnor=&vnor[0];
        mesh->init( vpos.size(), 0, 0, &vpos[0], pnor, NULL, NULL );
        mesh->draw_mode = GL_LINES;
        return mesh;
    }

    void newSub(){ subs.push_back({vpos.size(),inds.size()}); }
    inline Vec2i subVertRange(int i){ int i0=0; if(i>0)i0=subs[i-1].a; return {i0,subs[i].a}; }
    inline Vec2i subIndRange (int i){ int i0=0; if(i>0)i0=subs[i-1].b; return {i0,subs[i].b}; }

    void move ( Vec2i iv, Vec3f shift ){ for(int i=iv.a; i<iv.b; i++){ vpos[i].add(shift); } }
    void scale( Vec2i iv, Vec3f sc    ){
        Vec3f invSc = {1/sc.x,1/sc.y,1/sc.z};
        for(int i=iv.a; i<iv.b; i++){
            vpos[i].mul(sc);
            if(bnor) vnor[i].mul(invSc); vnor[i].normalize();
        }
    }
    void rotate( Vec2i iv, Vec3f p0, Vec3f p1, float angle ){
        Vec3f uax=p1-p0; uax.normalize();
        Vec2f cs; cs.fromAngle(angle);
        for(int i=iv.a; i<iv.b; i++){
            Vec3f v = vpos[i]-p0;  v.rotate_csa(cs.a,cs.b,uax); vpos[i]=v+p0;
            if(bnor) vnor[i].rotate_csa(cs.a,cs.b,uax);
        }
    }

    void applyMatrix( Vec2i iv, Mat3f mat ){
        for(int i=iv.a; i<iv.b; i++){
            mat.dot_to(vpos[i], vpos[i]);
            if(bnor) mat.dot_to(vnor[i],vnor[i]);
        }
    }

    void applyMatrixT( Vec2i iv, Mat3f mat ){
        for(int i=iv.a; i<iv.b; i++){
            mat.dot_to_T(vpos[i], vpos[i]);
            if(bnor) mat.dot_to_T(vnor[i],vnor[i]);
        }
    }

    void duplicateSub( int i ){
        Vec2i ivs = subVertRange(i);
        Vec2i iis = subIndRange (i);
        int di    = vpos.size() - ivs.a;
        for(int i=iis.a; i<iis.b; i++){
            inds.push_back( inds[i]+di );
        }
        for(int i=ivs.a; i<ivs.b; i++){
            vpos         .push_back( vpos[i] );
            if(bnor) vnor.push_back( vnor[i] );
            if(bUVs) vUVs.push_back( vUVs[i] );
        }
        newSub();
    }

    void addLine( Vec3f p1, Vec3f p2, Vec3f c ){
        vpos.push_back(p1);  vnor.push_back(c);
        vpos.push_back(p2);  vnor.push_back(c);
    };

    void addLine( Vec3d p1, Vec3d p2, Vec3f c ){
        Vec3f p;
        convert(p1,p); vpos.push_back(p);  vnor.push_back(c);
        convert(p2,p); vpos.push_back(p);  vnor.push_back(c);
    };

    void addPointCross( Vec3f p, float d, Vec3f c ){
        addLine( p+(Vec3f){d,0,0}, p+(Vec3f){-d, 0, 0}, c );
        addLine( p+(Vec3f){0,d,0}, p+(Vec3f){ 0,-d, 0}, c );
        addLine( p+(Vec3f){0,0,d}, p+(Vec3f){ 0, 0,-d}, c );
    };

    void addArrow( Vec3f p1, Vec3f p2, float d, Vec3f c ){

    };

    void addLines( int n, int * inds, Vec3f* verts, Vec3f c ){
        for(int i=0; i<n; i++){
            int i2=i*2;
            addLine( verts[inds[i2]], verts[inds[i2+1]], c );
        }
    }


    void addLineStrip( int n, Vec3f* verts, Vec3f c ){
        for(int i=0; i<n; i++){
            int i2=i*2;
            addLine( verts[inds[i2]], verts[inds[i2+1]], c );
        }
    }

    void moveSub     ( int i, Vec3f shift ){ move ( subVertRange(i), shift); }
    void scaleSub    ( int i, Vec3f sc    ){ scale( subVertRange(i), sc   ); }
    void rotateSub   ( int i, Vec3f p0, Vec3f p1, float angle ){ rotate( subVertRange(i), p0, p1, angle ); }

    GLMesh* normals2GLmesh( float sc ){ return vecs2mesh( vpos.size(), &vpos[0], &vnor[0], sc ); }
}; // class GLMeshBuilder



//   ======  UV Draw Templates



template<typename UVfunc>
void UVFunc2smooth( Vec2i n, Vec2f UVmin, Vec2f UVmax, float voff, UVfunc func, GLMeshBuilder& mesh ){
    mesh.draw_mode = GL_TRIANGLES;
    Vec2f duv = UVmax-UVmin; duv.mul( {1.0f/n.a,1.0f/n.b} );
    //int i0=mesh.vpos.size()/3;
    for(int ia=0;ia<=n.a;ia++){
        Vec2f uv = { UVmin.a+duv.a*ia, UVmin.b+voff*duv.b*ia };
        for(int ib=0;ib<=n.b;ib++){
            int i=mesh.vpos.size();
            Vec3f p = func(uv);
            mesh.vpos.push_back( p );
            if( (ia>0) && (ib>0) ){
                int nb = n.b+1;
                push3i( mesh.inds, {i-nb,i-1,i     } );
                push3i( mesh.inds, {i-nb,i-1,i-nb-1} );
            }
            if(mesh.bUVs) mesh.vUVs.push_back( uv );
            if(mesh.bnor){
                /*
                float h=0.01;
                Vec2f o;
                Vec3f nor,da,db;
                o=uv; if(ia>0  )o.a+=h; da.set(func(o));
                o=uv; if(ia<n.a)o.a-=h; da.sub(func(o));
                o=uv; if(ib>0  )o.b+=h; db.set(func(o));
                o=uv; if(ib<n.b)o.b-=h; db.sub(func(o));
                nor.set_cross(db,da); nor.normalize();
                */
                mesh.vnor.push_back( getUVFuncNormal(uv,0.01,func) );
                //push3f(mesh.vnor,p);
            };
            uv.b+=duv.b;
        }
    }
    mesh.newSub();
}

template<typename UVfunc>
void UVFunc2wire( Vec2i n, Vec2f UVmin, Vec2f UVmax, float voff, UVfunc func, GLMeshBuilder& mesh ){
    mesh.draw_mode = GL_LINES;
    Vec2f duv = UVmax-UVmin; duv.mul( {1.0f/n.a,1.0f/n.b} );
    //int i0=mesh.vpos.size()/3;
    for(int ia=0;ia<=n.a;ia++){
        Vec2f uv = { UVmin.a+duv.a*ia, UVmin.b+voff*duv.b*ia };
        for(int ib=0;ib<=n.b;ib++){
            int i=mesh.vpos.size();
            Vec3f p = func(uv);
            mesh.vpos.push_back( p );
            int nb = n.b+1;
            if (ia<n.a){ push2i( mesh.inds, {i,i+nb } ); }
            if (ib<n.b){ push2i( mesh.inds, {i,i+1  } ); }
            //if(mesh.bUVs) push(mesh.vUVs, uv );
            uv.b+=duv.b;
        }
    }
    mesh.newSub();
}


template<typename UVfunc>
void drawExtrudedWireUVFunc( Vec2i n, float thick, Vec2f UVmin, Vec2f UVmax, float voff, UVfunc func, GLMeshBuilder& mesh ){
    mesh.draw_mode = GL_TRIANGLES;

    float eps = 0.001;
    Vec2f duv = UVmax-UVmin; duv.mul( {1.0f/n.a,1.0f/n.b} );
    //int i0=mesh.vpos.size()/3;
    for(int ia=0;ia<=n.a;ia++){
        Vec2f uv = { UVmin.a+duv.a*ia, UVmin.b+voff*duv.b*ia };
        for(int ib=0;ib<=n.b;ib++){
            int i=mesh.vpos.size();
            Vec3f p = func(uv);
            Vec3f nr = getUVFuncNormal(uv,eps,func);
            mesh.vpos.push_back( p          );
            mesh.vpos.push_back( p+nr*thick );
            if( (ia>0) && (ib>0) ){
                int nb = n.b+1;
                push3i( mesh.inds, {i-nb,i-1,i     } );
                push3i( mesh.inds, {i-nb,i-1,i-nb-1} );
            }
            if(mesh.bnor){
                Vec3f nrl;
                nrl.set_cross(func(uv+(Vec2f){0.0,duv.b*0.5})-p, nr); nrl.normalize();
                mesh.vnor.push_back( nrl );
                mesh.vnor.push_back( nrl );
                //push3f(mesh.vnor,p);
            };
            //int nb = n.b+1;
            //if (ia<n.a){ push2i( mesh.inds, {i,i+nb } ); }
            //if (ib<n.b){ push2i( mesh.inds, {i,i+1  } ); }
            //if(mesh.bUVs) push(mesh.vUVs, uv );
            uv.b+=duv.b;
        }
    }
    mesh.newSub();


    /*
    Vec2f uv;
    Vec3f p,nr,nrl;
    for(int ia=0;ia<=n.a;ia++){
        // strips
        glBegin(GL_TRIANGLE_STRIP);
        for(int ib=0;ib<=n.b;ib++){
            uv.a = UVmin.a+duv.a*ia;
            uv.b = UVmin.b+duv.b*ib + voff*duv.b*ia;
            p  = func(uv);
            nr = getUVFuncNormal(uv,eps,func);
            nrl.set_cross(func(uv+(Vec2f){0.0,duv.b*0.5})-p, nr); nrl.normalize();
            glNormal3f(nrl.x,nrl.y,nrl.z);
            glVertex3f(p.x,p.y,p.z);
            p.add_mul(nr,thick);
            glVertex3f(p.x,p.y,p.z);
        }
        glEnd();
        if(ia<n.a){
            glBegin(GL_TRIANGLE_STRIP);
            for(int ib=0;ib<=n.b;ib++){
                uv.a = UVmin.a+duv.a*ia;
                uv.b = UVmin.b+duv.b*ib + voff*duv.b*ia;

                p  = func(uv);
                nr = getUVFuncNormal(uv,eps,func);
                nrl.set_cross(func(uv+(Vec2f){0.0,duv.b*0.5})-p, nr); nrl.normalize();
                glNormal3f(nrl.x,nrl.y,nrl.z);
                glVertex3f(p.x,p.y,p.z);
                p.add_mul(nr,thick);
                glVertex3f(p.x,p.y,p.z);

                uv.a += duv.a;
                uv.b += voff*duv.b;

                p  = func(uv);
                nr = getUVFuncNormal(uv,eps,func);
                nrl.set_cross(func(uv+(Vec2f){0.0,duv.b*0.5})-p, nr); nrl.normalize();
                glNormal3f(nrl.x,nrl.y,nrl.z);
                glVertex3f(p.x,p.y,p.z);
                p.add_mul(nr,thick);
                glVertex3f(p.x,p.y,p.z);

            }
            glEnd();
        }
    }
    */
}

// ================== UV FUNCS

void Cone2Mesh( Vec2i n, Vec2f UVmin, Vec2f UVmax, float R1, float R2, float L, float voff, bool wire, GLMeshBuilder& mesh ){
    auto uvfunc = [&](Vec2f uv){return ConeUVfunc(uv,R1,R2,L);};
    if(wire){ UVFunc2wire  ( n, UVmin, UVmax, voff, uvfunc, mesh ); }
    else    { UVFunc2smooth( n, UVmin, UVmax, voff, uvfunc, mesh ); }
}

void Sphere2Mesh( Vec2i n, Vec2f UVmin, Vec2f UVmax, float R, float voff, bool wire, GLMeshBuilder& mesh ){
    auto uvfunc = [&](Vec2f uv){return SphereUVfunc(uv,R);};
    if(wire){ UVFunc2wire  ( n, UVmin, UVmax, voff, uvfunc, mesh ); }
    else    { UVFunc2smooth( n, UVmin, UVmax, voff, uvfunc, mesh ); }
}

void Torus2Mesh( Vec2i n, Vec2f UVmin, Vec2f UVmax, float r, float R, float voff, bool wire, GLMeshBuilder& mesh ){
    auto uvfunc = [&](Vec2f uv){return TorusUVfunc(uv,r,R);};
    if(wire){ UVFunc2wire  ( n, UVmin, UVmax, voff, uvfunc, mesh ); }
    else    { UVFunc2smooth( n, UVmin, UVmax, voff, uvfunc, mesh ); }
}

void Teardrop2Mesh( Vec2i n, Vec2f UVmin, Vec2f UVmax, float R1, float R2, float L, float voff, bool wire, GLMeshBuilder& mesh ){
    auto uvfunc = [&](Vec2f uv){return TeardropUVfunc(uv,R1,R2,L);};
    if(wire){ UVFunc2wire  ( n, UVmin, UVmax, voff, uvfunc, mesh ); }
    else    { UVFunc2smooth( n, UVmin, UVmax, voff, uvfunc, mesh ); }
}

void NACASegment2Mesh( Vec2i n, Vec2f UVmin, Vec2f UVmax, float *coefs1, float *coefs2, float L, float voff, bool wire, GLMeshBuilder& mesh ){
    auto uvfunc = [&](Vec2f uv){return NACA4digitUVfunc(uv,coefs1,coefs2,L);};
    if(wire){ UVFunc2wire  ( n, UVmin, UVmax, voff, uvfunc, mesh ); }
    else    { UVFunc2smooth( n, UVmin, UVmax, voff, uvfunc, mesh ); }
}




inline void HarmonicTube2Mesh( Vec2i n, Vec2f UVmin, Vec2f UVmax, float R1, float R2, float L, float voff, float freq, float amp, bool wire, GLMeshBuilder& mesh ){
    auto uvfunc = [&](Vec2f uv){return HarmonicTubeUVfunc(uv,R1,R2,L,freq,amp);};
    if(wire){ UVFunc2wire( n, UVmin, UVmax, voff, uvfunc, mesh ); }
    else    { UVFunc2smooth( n, UVmin, UVmax,voff, uvfunc, mesh ); }
}

inline void Parabola2Mesh( Vec2i n, Vec2f UVmin, Vec2f UVmax, float R, float L, float voff, bool wire, GLMeshBuilder& mesh ){
    float K = L/(R*R);
    UVmin.a*=R; UVmax.a*=R;
    //printf( "drawUV_Parabola: R %f L %f K %f \n", R, L, K  );
    auto uvfunc = [&](Vec2f uv){return ParabolaUVfunc(uv,K);};
    if(wire){ UVFunc2wire( n, UVmin, UVmax,voff, uvfunc, mesh ); }
    else    { UVFunc2smooth( n, UVmin, UVmax,voff, uvfunc, mesh ); }
}

/*
inline void Parabola_ExtrudedWire( Vec2i n,Vec2f UVmin, Vec2f UVmax, float R, float L, float voff, double thick, GLMeshBuilder& mesh ){
    float K = L/(R*R);
    UVmin.a*=R; UVmax.a*=R;
    auto uvfunc = [&](Vec2f uv){return ParabolaUVfunc(uv,K);};
    drawExtrudedWireUVFunc( n,thick, UVmin, UVmax,voff, uvfunc );
}
*/

inline void Hyperbola2Mesh( Vec2i n, Vec2f UVmin, Vec2f UVmax, float r, float R, float L, float voff, bool wire, GLMeshBuilder& mesh ){
    //printf( "drawUV_Hyperbola: r %f R %f L %f \n", r, R, L  );
    if(r>0){
        float K = R/L;
        UVmin.a*=L; UVmax.a*=L;
        auto uvfunc = [&](Vec2f uv){return HyperbolaLUVfunc(uv,r,K);};
        if(wire){ UVFunc2wire( n, UVmin, UVmax,voff, uvfunc, mesh ); }
        else    { UVFunc2smooth( n, UVmin, UVmax,voff, uvfunc, mesh ); }
    }else{
        r=-r;
        float K = L/R;
        UVmin.a*=R; UVmax.a*=R;
        auto uvfunc = [&](Vec2f uv){return HyperbolaRUVfunc(uv,r,K);};
        if(wire){ UVFunc2wire( n, UVmin, UVmax,voff, uvfunc, mesh ); }
        else    { UVFunc2smooth( n, UVmin, UVmax,voff, uvfunc, mesh ); }
    }
}

#endif
