
#ifndef  MeshBuilder_h
#define  MeshBuilder_h

#include <vector>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"

//#include "CMesh.h"    // ToDo : possily we should derive or ouput to CMesh
//#include "Mesh.h"       // ToDo : we should connect with Mesh.h ( Genral polygon-mesh with usefull editation tools )
#include "UVfuncs.h"

// ===== MeshBuilder =====
//   Abstraction over rendering polygonn meshes
//   This should allow same algorithm to be used to create mesh which can be later easily transformed into
//    * OpenCL 1.2
//    * OpenGL 3./4.
//    * .obj file format
//    * .mesh generators for finite-element calculations ( RayTracing, Radiosity, Scattering, soft-bondy mechanis, linear flow, viscoelasticity ... )
//     There should be tools for both creation of mesh (e.g. UV function, constructive solid geometry (CSG), polygon edity, besier patches ...  as well as eddition
// ==========================

inline void push3i( std::vector<int>&   vs, const Vec3i& v ){ vs.push_back(v.x); vs.push_back(v.y); vs.push_back(v.z); }
inline void push2i( std::vector<int>&   vs, const Vec2i& v ){ vs.push_back(v.x); vs.push_back(v.y); }


class MeshBuilder{ public:
    bool bnor = true;
    bool bUVs = true;
    std::vector<Vec3i>  subs;   //=std::vector<Vec2i>({(Vec2i){0,0}});
    //std::vector<int>   modes;  // mode for each sub
    std::vector<Vec3f>  vpos;
    std::vector<Vec3f>  vnor;
    //std::vector<Vec3f> vcol;
    std::vector<Vec2f> vUVs;
    std::vector<int>   inds;
    GLenum draw_mode = GL_TRIANGLES;
    Vec3f penColor;

    void clear(){ vpos.clear(); vnor.clear(); vUVs.clear(); inds.clear(); }

/*
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
*/

    void   newSub( GLenum draw_mode = GL_TRIANGLES ){ subs.push_back({vpos.size(),inds.size(), draw_mode}); }
    inline Vec2i subVertRange(int i){ int i0=0; if(i<0)i=subs.size()+i; i0=subs[i-1].a; return {i0,subs[i].a}; }
    inline Vec2i subIndRange (int i){ int i0=0; if(i<0)i=subs.size()+i; i0=subs[i-1].b; return {i0,subs[i].b}; }

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

    void addLine( Vec3f p1, Vec3f p2 ){
        vpos.push_back(p1);  vnor.push_back(penColor);
        vpos.push_back(p2);  vnor.push_back(penColor);
    };

    void addLine( Vec3d p1, Vec3d p2 ){
        vpos.push_back((Vec3f)p1);  vnor.push_back(penColor);
        vpos.push_back((Vec3f)p2);  vnor.push_back(penColor);
    };

    void addPointCross( Vec3f p, float d ){
        addLine( p+(Vec3f){d,0,0}, p+(Vec3f){-d, 0, 0} );
        addLine( p+(Vec3f){0,d,0}, p+(Vec3f){ 0,-d, 0} );
        addLine( p+(Vec3f){0,0,d}, p+(Vec3f){ 0, 0,-d} );
    };

    void addArrow( Vec3f p1, Vec3f p2, float d ){

    };

    void addLines( int n, int * inds, Vec3f* verts ){
        for(int i=0; i<n; i++){
            int i2=i*2;
            addLine( verts[inds[i2]], verts[inds[i2+1]] );
        }
    }

    void addTriangle( const Vec3f& a, const Vec3f& b, const Vec3f& c ){
        Vec3f nor; nor.set_cross( a-c, b-c ); nor.normalize();
        vpos.push_back( a ); vnor.push_back( nor );
        vpos.push_back( b ); vnor.push_back( nor );
        vpos.push_back( c ); vnor.push_back( nor );
    }

    void addQuad( const Vec3f& p00, const Vec3f& p01, const Vec3f& p10, const Vec3f& p11 ){
        addTriangle( p00, p01, p10 );
        addTriangle( p11, p10, p01 );
    }

    void addTube4( Vec3f p1, Vec3f p2, Vec3f up, float r1, float r2 ){
        up.normalize();
        Vec3f lf; lf.set_cross( p2-p1, up ); lf.normalize();
        addQuad( p1+up*r1, p1+lf*r1, p2+up*r2, p2+lf*r2 );
        addQuad( p1+lf*r1, p1-up*r1, p2+lf*r2, p2-up*r2 );
        addQuad( p1-up*r1, p1-lf*r1, p2-up*r2, p2-lf*r2 );
        addQuad( p1-lf*r1, p1+up*r1, p2-lf*r2, p2+up*r2 );
    };

    int addCircleAxis( int n, const Vec3f& pos, const Vec3f& v0, const Vec3f& uaxis, float R ){
        float dphi = 2*M_PI/n;
        float dca  = cos( dphi );
        float dsa  = sin( dphi );
        int nvert=0;
        Vec3f v; v.set(v0);
        //glBegin( GL_LINE_LOOP );
        for( int i=0; i<n; i++ ){
            //glVertex3f( pos.x+v.x*R, pos.y+v.y*R, pos.z+v.z*R ); nvert++;
            //printf( " drawCircleAxis %i (%3.3f,%3.3f,%3.3f) \n", i, v.x, v.y, v.z );
            Vec3f ov = v;
            v.rotate_csa( dca, dsa, uaxis );
            addLine( ov, v );
        }
        //glEnd();
        return nvert;
    }

    void moveSub     ( int i, Vec3f shift ){ move ( subVertRange(i), shift); }
    void scaleSub    ( int i, Vec3f sc    ){ scale( subVertRange(i), sc   ); }
    void rotateSub   ( int i, Vec3f p0, Vec3f p1, float angle ){ rotate( subVertRange(i), p0, p1, angle ); }



    int  addCapsula( Vec3f p0, Vec3f p1, float r1, float r2, float theta1, float theta2, float dTheta, int nPhi, bool capped ){
        int nvert=0;
        Vec3f ax   = p1-p0;  float L = ax.normalize();
        Vec3f up,left;       ax.getSomeOrtho(up,left);
        Vec2f cph=Vec2fX, dph;
        dph.fromAngle( 2*M_PI/nPhi );
        // Cylinder
        Vec2f cth,dth;
        float dr = (r2-r1);
        float cv = sqrt(L*L+dr*dr);
        cth.set( L/cv, -dr/cv );
        //glBegin(GL_TRIANGLE_STRIP);

        Vec3f opa=p0,opb=p1;
        for(int iph=0; iph<(nPhi+1); iph++){
            Vec3f pa = p0 + left*(cph.x*r1) + up*(cph.y*r1);
            Vec3f pb = p1 + left*(cph.x*r2) + up*(cph.y*r2);
            Vec3f na = (left*(cph.x) + up*(cph.y)*cth.x + ax*cth.y)*-1.0;
            Vec3f nb = (left*(cph.x) + up*(cph.y)*cth.x + ax*cth.y)*-1.0;
            //glNormal3f(na.x,na.y,na.z); glVertex3f(pa.x,pa.y,pa.z);
            //glNormal3f(nb.x,nb.y,nb.z); glVertex3f(pb.x,pb.y,pb.z);
            addQuad( opa, opb, pa, pb );
            opa=pa; opb=pb;
            cph.mul_cmplx(dph);
            nvert+=2;
        }
        glEnd();

        float DTh,h;
        int nTheta;
        // Spherical Cap
        cph=Vec2fX;
        cth.set( L/cv, -dr/cv );
        //dth.fromSin(v1/r1);
        DTh = (-theta1 - asin(cth.y));
        nTheta = (int)(fabs(DTh)/dTheta);
        dth.fromAngle( DTh/nTheta );
        //printf( " cth (%f,%f)  dth (%f,%f) \n", cth.x, cth.y,  dth.x, dth.y );
        r1/=cth.x;
        h  =-cth.y*r1;
        // Left
        opa=p0;opb=p1;
        for(int ith=0; ith<(nTheta+1); ith++){
            Vec2f cth_ = Vec2f::mul_cmplx(cth,dth);
            //glBegin(GL_TRIANGLE_STRIP);
            //glBegin(GL_LINES);
            for(int iph=0; iph<(nPhi+1); iph++){
                Vec3f pa = p0 + (left*(cph.x*r1) + up*(cph.y*r1))*cth.x  + ax*(h+cth.y*r1);
                Vec3f pb = p0 + (left*(cph.x*r1) + up*(cph.y*r1))*cth_.x + ax*(h+cth_.y*r1);
                Vec3f na = (left*(cph.x) + up*(cph.y)*cth.x  + ax*cth.y)*1.0;
                Vec3f nb = (left*(cph.x) + up*(cph.y)*cth_.x + ax*cth_.y)*1.0;
                //glNormal3f(na.x,na.y,na.z); glVertex3f(pa.x,pa.y,pa.z);
                //glNormal3f(nb.x,nb.y,nb.z); glVertex3f(pb.x,pb.y,pb.z);
                addQuad( opa, opb, pa, pb );
                opa=pa; opb=pb;
                nvert+=2;
                cph.mul_cmplx(dph);
            }
            //glEnd();
            //printf( "%i cth (%f,%f)  cth_ (%f,%f) \n", ith, cth.x, cth.y,  cth_.x, cth_.y );
            cth=cth_;
        }
        //return 0;
        cph=Vec2fX;
        cth.set( L/cv, -dr/cv );
        //cth = Vec2fX;
        //cth.set( dr/cv, L/cv);
        //dth.fromAngle( asin(v2/r2)/nTheta );
        DTh    = (theta2-asin(cth.y));
        nTheta = (int)(fabs(DTh)/dTheta);
        dth.fromAngle(DTh/nTheta );
        r2/= cth.x;
        h  =-cth.y*r2;
        // Right
        for(int ith=0; ith<(nTheta+1); ith++){
            Vec2f cth_ = Vec2f::mul_cmplx(cth,dth);
            glBegin(GL_TRIANGLE_STRIP);
            for(int iph=0; iph<(nPhi+1); iph++){
                Vec3f pa = p1 + (left*(cph.x*r2) + up*(cph.y*r2))*cth.x  + ax*(h+cth.y*r2);
                Vec3f pb = p1 + (left*(cph.x*r2) + up*(cph.y*r2))*cth_.x + ax*(h+cth_.y*r2);
                Vec3f na = (left*(cph.x) + up*(cph.y)*cth.x  + ax*cth.y)*-1.0;
                Vec3f nb = (left*(cph.x) + up*(cph.y)*cth_.x + ax*cth_.y)*-1.0;
                //glNormal3f(na.x,na.y,na.z); glVertex3f(pa.x,pa.y,pa.z);
                //glNormal3f(nb.x,nb.y,nb.z); glVertex3f(pb.x,pb.y,pb.z);
                addQuad( opa, opb, pa, pb );
                opa=pa; opb=pb;
                nvert+=2;
                cph.mul_cmplx(dph);
            }
            glEnd();
            cth=cth_;
        }
        return nvert;
    }
    //GLMesh* normals2GLmesh( float sc ){ return vecs2mesh( vpos.size(), &vpos[0], &vnor[0], sc ); }




    void write_obj( char* fname ){
        int nsubs  = subs.size();
        Vec3i osub = subs[0];

        FILE * pFile;
        pFile = fopen (fname,"w");

        int nobj =0;
        int nvert=0;
        int nnor =0;
        for( int i=1; i<nsubs; i++ ){
            Vec3i sub = subs[i];
            int mode  = osub.z;
            if      (mode == GL_TRIANGLES ){      // un-indexed triangles

                printf(         "o OBJ_TRIANGLES.%i  [ %i ... %i ] \n", nobj, osub.x, sub.x );
                fprintf( pFile, "o OBJ_TRIANGLES.%i \n", nobj ); nobj++;
                int iii = 0;
                for(int j=osub.x; j<sub.x; j++){
                    Vec3f vp = vpos[j];
                    Vec3f vn = vnor[j];
                    fprintf(pFile, "v   %f %f %f\n", vp.x, vp.y, vp.z ); nvert++;
                    fprintf(pFile, "vn  %f %f %f\n", vn.x, vn.y, vn.z ); nnor ++;
                    if(iii%3==2) fprintf( pFile, "f %i//%i %i//%i %i//%i \n", nvert-2,nnor-2,  nvert-1,nnor-1,   nvert,nnor );
                    iii++;
                }

            }else if(mode == GL_TRIANGLE_STRIP ) {  // Indexed Triangles
                printf(         "o OBJ_TRIANGLE_STRIP.%i  [ %i ... %i ] \n", nobj, osub.x, sub.x );
                fprintf( pFile, "o OBJ_TRIANGLE_STRIP.%i \n", nobj ); nobj++;
                /*
                //  ToDo - This does not work for some reason => brute force polygonization
                int iv0 = nvert-osub.x;
                int in0 = nnor -osub.x;
                for(int j=osub.x; j<sub.x; j++){
                    Vec3f vp = vpos[j];
                    Vec3f vn = vnor[j];
                    fprintf(pFile, "v   %f %f %f\n", vp.x, vp.y, vp.z ); nvert++;
                    fprintf(pFile, "vn  %f %f %f\n", vn.x, vn.y, vn.z ); nnor ++;
                }
                int iii = 0;
                for(int j=osub.y; j<sub.y; j+=3 ){
                    int i0=inds[j  ];
                    int i1=inds[j+1];
                    int i2=inds[j+2];
                    fprintf( pFile, "f %i//%i %i//%i %i//%i \n",    i0+iv0, i0+in0,       i1+iv0, i1+in0,     i2+iv0, i2+in0 );
                    iii++;
                }
                */
                int iii = 0;
                for(int j=osub.y; j<sub.y; j++ ){
                    int ii = inds[j  ];
                    Vec3f vp = vpos[ii];
                    Vec3f vn = vnor[ii];
                    fprintf(pFile, "vn  %f %f %f\n", vn.x, vn.y, vn.z ); nnor ++;
                    fprintf(pFile, "v   %f %f %f\n", vp.x, vp.y, vp.z ); nvert++;
                    if(iii%3==2) fprintf( pFile, "f %i//%i %i//%i %i//%i \n", nvert-2,nnor-2,  nvert-1,nnor-1,   nvert,nnor );
                    iii++;
                }

            }else if(mode == GL_LINES ) {
                /*
                printf(         "o OBJ_LINES.%i  [ %i ... %i ] \n", nobj, osub.x, sub.x );
                fprintf( pFile, "o OBJ_LINES.%i \n", nobj ); nobj++;
                int iii = 0;
                for(int j=osub.x; j<sub.x; j++){
                    //Vec3f vnor = mesh.vnor[j];
                    Vec3f vp = vpos[j];
                    fprintf(pFile, "v   %f %f %f\n", vp.x, vp.y, vp.z ); nvert++;
                    //if(iii%2==1) printf( "l %i %i \n", nvert-1, nvert );
                    if(iii%2==1) printf( "f %i %i \n", nvert-1, nvert );
                    iii++;
                }
                */
            }
            osub=sub;
        }
        fclose(pFile);
    }

}; // class MeshBuilder



//   ======  UV Draw Templates



template<typename UVfunc>
void UVFunc2smooth( Vec2i n, Vec2f UVmin, Vec2f UVmax, float voff, UVfunc func, MeshBuilder& mesh ){
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
                mesh.vnor.push_back( getUVFuncNormal(uv,0.01,func) );
            };
            uv.b+=duv.b;
        }
    }
    mesh.newSub();
}

template<typename UVfunc>
void UVFunc2wire( Vec2i n, Vec2f UVmin, Vec2f UVmax, float voff, UVfunc func, MeshBuilder& mesh ){
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
void ExtrudedWireUVFunc( Vec2i n, float thick, Vec2f UVmin, Vec2f UVmax, float voff, UVfunc func, MeshBuilder& mesh ){
    //mesh.draw_mode = GL_TRIANGLES;
    //mesh.draw_mode = GL_TRIANGLE_STRIP;

    float eps = 0.001;
    Vec2f duv = UVmax-UVmin; duv.mul( {1.0f/n.a,1.0f/n.b} );
    //int i0=mesh.vpos.size()/3;

    for(int ia=0;ia<=n.a;ia++){
        Vec2f uv = { UVmin.a+duv.a*ia, UVmin.b+voff*duv.b*ia };
        for(int ib=0;ib<=n.b;ib++){
            int i=mesh.vpos.size();
            Vec3f p  = func(uv);
            Vec3f nr = getUVFuncNormal(uv,eps,func);
            mesh.vpos.push_back( p          );
            mesh.vpos.push_back( p+nr*thick );
            if( ib>0 ){
                //int nb = n.b+1;
                //push3i( mesh.inds, {i-nb,i-1,i     } );
                //push3i( mesh.inds, {i-nb,i-1,i-nb-1} );
                push3i( mesh.inds, {i-2,i-1,i } );
                push3i( mesh.inds, {i-1,i  ,i+1 } );
            }
            if( ia>0 ){
                int nb = n.b*2+1;
                push3i( mesh.inds, {i-nb-1,i-nb,i} );
                push3i( mesh.inds, {i-nb,i,i+1} );

                push3i( mesh.inds, {i-nb+1,i-nb+2,i} );
                push3i( mesh.inds, {i-nb+2,i,i+1} );
            }
            if(mesh.bnor){
                Vec3f nrl;
                nrl.set_cross(func(uv+(Vec2f){0.0,duv.b*0.5})-p, nr); nrl.normalize();
                mesh.vnor.push_back( nrl );
                mesh.vnor.push_back( nrl );

                // ToDo : Normals are not correct
                mesh.vnor.push_back( nrl );
                mesh.vnor.push_back( nrl );

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

}

// ================== UV FUNCS

void Cone2Mesh( Vec2i n, Vec2f UVmin, Vec2f UVmax, float R1, float R2, float L, float voff, bool wire, MeshBuilder& mesh ){
    auto uvfunc = [&](Vec2f uv){return ConeUVfunc(uv,R1,R2,L);};
    if(wire){ UVFunc2wire  ( n, UVmin, UVmax, voff, uvfunc, mesh ); }
    else    { UVFunc2smooth( n, UVmin, UVmax, voff, uvfunc, mesh ); }
}

void Sphere2Mesh( Vec2i n, Vec2f UVmin, Vec2f UVmax, float R, float voff, bool wire, MeshBuilder& mesh ){
    auto uvfunc = [&](Vec2f uv){return SphereUVfunc(uv,R);};
    if(wire){ UVFunc2wire  ( n, UVmin, UVmax, voff, uvfunc, mesh ); }
    else    { UVFunc2smooth( n, UVmin, UVmax, voff, uvfunc, mesh ); }
}

void Torus2Mesh( Vec2i n, Vec2f UVmin, Vec2f UVmax, float r, float R, float voff, bool wire, MeshBuilder& mesh ){
    auto uvfunc = [&](Vec2f uv){return TorusUVfunc(uv,r,R);};
    if(wire){ UVFunc2wire  ( n, UVmin, UVmax, voff, uvfunc, mesh ); }
    else    { UVFunc2smooth( n, UVmin, UVmax, voff, uvfunc, mesh ); }
}

void Teardrop2Mesh( Vec2i n, Vec2f UVmin, Vec2f UVmax, float R1, float R2, float L, float voff, bool wire, MeshBuilder& mesh ){
    auto uvfunc = [&](Vec2f uv){return TeardropUVfunc(uv,R1,R2,L);};
    if(wire){ UVFunc2wire  ( n, UVmin, UVmax, voff, uvfunc, mesh ); }
    else    { UVFunc2smooth( n, UVmin, UVmax, voff, uvfunc, mesh ); }
}

void NACASegment2Mesh( Vec2i n, Vec2f UVmin, Vec2f UVmax, float *coefs1, float *coefs2, float L, float voff, bool wire, MeshBuilder& mesh ){
    auto uvfunc = [&](Vec2f uv){return NACA4digitUVfunc(uv,coefs1,coefs2,L);};
    if(wire){ UVFunc2wire  ( n, UVmin, UVmax, voff, uvfunc, mesh ); }
    else    { UVFunc2smooth( n, UVmin, UVmax, voff, uvfunc, mesh ); }
}




inline void HarmonicTube2Mesh( Vec2i n, Vec2f UVmin, Vec2f UVmax, float R1, float R2, float L, float voff, float freq, float amp, bool wire, MeshBuilder& mesh ){
    auto uvfunc = [&](Vec2f uv){return HarmonicTubeUVfunc(uv,R1,R2,L,freq,amp);};
    if(wire){ UVFunc2wire( n, UVmin, UVmax, voff, uvfunc, mesh ); }
    else    { UVFunc2smooth( n, UVmin, UVmax,voff, uvfunc, mesh ); }
}

inline void Parabola2Mesh( Vec2i n, Vec2f UVmin, Vec2f UVmax, float R, float L, float voff, bool wire, MeshBuilder& mesh ){
    float K = L/(R*R);
    UVmin.a*=R; UVmax.a*=R;
    //printf( "drawUV_Parabola: R %f L %f K %f \n", R, L, K  );
    auto uvfunc = [&](Vec2f uv){return ParabolaUVfunc(uv,K);};
    if(wire){ UVFunc2wire( n, UVmin, UVmax,voff, uvfunc, mesh ); }
    else    { UVFunc2smooth( n, UVmin, UVmax,voff, uvfunc, mesh ); }
}

inline void Parabola2Mesh_ExtrudedWire( Vec2i n,Vec2f UVmin, Vec2f UVmax, float R, float L, float voff, double thick, MeshBuilder& mesh ){
    float K = L/(R*R);
    UVmin.a*=R; UVmax.a*=R;
    auto uvfunc = [&](Vec2f uv){return ParabolaUVfunc(uv,K);};
    ExtrudedWireUVFunc( n,thick, UVmin, UVmax,voff, uvfunc, mesh );
}

/*
inline void Parabola_ExtrudedWire( Vec2i n,Vec2f UVmin, Vec2f UVmax, float R, float L, float voff, double thick, MeshBuilder& mesh ){
    float K = L/(R*R);
    UVmin.a*=R; UVmax.a*=R;
    auto uvfunc = [&](Vec2f uv){return ParabolaUVfunc(uv,K);};
    drawExtrudedWireUVFunc( n,thick, UVmin, UVmax,voff, uvfunc );
}
*/

inline void Hyperbola2Mesh( Vec2i n, Vec2f UVmin, Vec2f UVmax, float r, float R, float L, float voff, bool wire, MeshBuilder& mesh ){
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


