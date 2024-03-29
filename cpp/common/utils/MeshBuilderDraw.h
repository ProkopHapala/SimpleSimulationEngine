
#ifndef  MeshBuilderDraw_h
#define  MeshBuilderDraw_h

#include <vector>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"

//#include "CMesh.h"    // ToDo : possily we should derive or ouput to CMesh
//#include "Mesh.h"     // ToDo : we should connect with Mesh.h ( Genral polygon-mesh with usefull editation tools )
#include "UVfuncs.h"
#include "MeshBuilder.h"

namespace Mesh{

//   ======  UV Draw Templates

template<typename UVfunc>
void UVFunc2smooth( Vec2i n, Vec2f UVmin, Vec2f UVmax, float voff, UVfunc func, Builder& mesh ){
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
void UVFunc2wire( Vec2i n, Vec2f UVmin, Vec2f UVmax, float voff, UVfunc func, Builder& mesh ){
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
void ExtrudedWireUVFunc( Vec2i n, float thick, Vec2f UVmin, Vec2f UVmax, float voff, UVfunc func, Builder& mesh ){
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

void Cone2Mesh( Vec2i n, Vec2f UVmin, Vec2f UVmax, float R1, float R2, float L, float voff, bool wire, Builder& mesh ){
    auto uvfunc = [&](Vec2f uv){return ConeUVfunc(uv,R1,R2,L);};
    if(wire){ UVFunc2wire  ( n, UVmin, UVmax, voff, uvfunc, mesh ); }
    else    { UVFunc2smooth( n, UVmin, UVmax, voff, uvfunc, mesh ); }
}

void Sphere2Mesh( Vec2i n, Vec2f UVmin, Vec2f UVmax, float R, float voff, bool wire, Builder& mesh ){
    auto uvfunc = [&](Vec2f uv){return SphereUVfunc(uv,R);};
    if(wire){ UVFunc2wire  ( n, UVmin, UVmax, voff, uvfunc, mesh ); }
    else    { UVFunc2smooth( n, UVmin, UVmax, voff, uvfunc, mesh ); }
}

void Torus2Mesh( Vec2i n, Vec2f UVmin, Vec2f UVmax, float r, float R, float voff, bool wire, Builder& mesh ){
    auto uvfunc = [&](Vec2f uv){return TorusUVfunc(uv,r,R);};
    if(wire){ UVFunc2wire  ( n, UVmin, UVmax, voff, uvfunc, mesh ); }
    else    { UVFunc2smooth( n, UVmin, UVmax, voff, uvfunc, mesh ); }
}

void Teardrop2Mesh( Vec2i n, Vec2f UVmin, Vec2f UVmax, float R1, float R2, float L, float voff, bool wire, Builder& mesh ){
    auto uvfunc = [&](Vec2f uv){return TeardropUVfunc(uv,R1,R2,L);};
    if(wire){ UVFunc2wire  ( n, UVmin, UVmax, voff, uvfunc, mesh ); }
    else    { UVFunc2smooth( n, UVmin, UVmax, voff, uvfunc, mesh ); }
}

void NACASegment2Mesh( Vec2i n, Vec2f UVmin, Vec2f UVmax, float *coefs1, float *coefs2, float L, float voff, bool wire, Builder& mesh ){
    auto uvfunc = [&](Vec2f uv){return NACA4digitUVfunc(uv,coefs1,coefs2,L);};
    if(wire){ UVFunc2wire  ( n, UVmin, UVmax, voff, uvfunc, mesh ); }
    else    { UVFunc2smooth( n, UVmin, UVmax, voff, uvfunc, mesh ); }
}

inline void HarmonicTube2Mesh( Vec2i n, Vec2f UVmin, Vec2f UVmax, float R1, float R2, float L, float voff, float freq, float amp, bool wire, Builder& mesh ){
    auto uvfunc = [&](Vec2f uv){return HarmonicTubeUVfunc(uv,R1,R2,L,freq,amp);};
    if(wire){ UVFunc2wire( n, UVmin, UVmax, voff, uvfunc, mesh ); }
    else    { UVFunc2smooth( n, UVmin, UVmax,voff, uvfunc, mesh ); }
}

inline void Parabola2Mesh( Vec2i n, Vec2f UVmin, Vec2f UVmax, float R, float L, float voff, bool wire, Builder& mesh ){
    float K = L/(R*R);
    UVmin.a*=R; UVmax.a*=R;
    //printf( "drawUV_Parabola: R %f L %f K %f \n", R, L, K  );
    auto uvfunc = [&](Vec2f uv){return ParabolaUVfunc(uv,K);};
    if(wire){ UVFunc2wire( n, UVmin, UVmax,voff, uvfunc, mesh ); }
    else    { UVFunc2smooth( n, UVmin, UVmax,voff, uvfunc, mesh ); }
}

inline void Parabola2Mesh_ExtrudedWire( Vec2i n,Vec2f UVmin, Vec2f UVmax, float R, float L, float voff, double thick, Builder& mesh ){
    float K = L/(R*R);
    UVmin.a*=R; UVmax.a*=R;
    auto uvfunc = [&](Vec2f uv){return ParabolaUVfunc(uv,K);};
    ExtrudedWireUVFunc( n,thick, UVmin, UVmax,voff, uvfunc, mesh );
}

/*
inline void Parabola_ExtrudedWire( Vec2i n,Vec2f UVmin, Vec2f UVmax, float R, float L, float voff, double thick, Builder& mesh ){
    float K = L/(R*R);
    UVmin.a*=R; UVmax.a*=R;
    auto uvfunc = [&](Vec2f uv){return ParabolaUVfunc(uv,K);};
    drawExtrudedWireUVFunc( n,thick, UVmin, UVmax,voff, uvfunc );
}
*/

inline void Hyperbola2Mesh( Vec2i n, Vec2f UVmin, Vec2f UVmax, float r, float R, float L, float voff, bool wire, Builder& mesh ){
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

//   ======  Capsula

int drawCapsula( Builder& mesh, Vec3f p0, Vec3f p1, float r1, float r2, float theta1, float theta2, float dTheta, int nPhi, bool capped ){
    printf( "Mesh::drawCapsula() p0(%g,%g,%g), p1(%g,%g,%g), rs(%g,%g) theta(%g,%g) dTheta=%g nPhi=%i capped=%i\n", p0.x,p0.y,p0.z, p1.x,p1.y,p1.z, r1,r2, theta1,theta2, dTheta, nPhi, capped );
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
        mesh.addQuad( opa, opb, pa, pb );
        opa=pa; opb=pb;
        cph.mul_cmplx(dph);
        nvert+=2;
    }
    //glEnd();

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
            mesh.addQuad( opa, opb, pa, pb );
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
        //glBegin(GL_TRIANGLE_STRIP);
        for(int iph=0; iph<(nPhi+1); iph++){
            Vec3f pa = p1 + (left*(cph.x*r2) + up*(cph.y*r2))*cth.x  + ax*(h+cth.y*r2);
            Vec3f pb = p1 + (left*(cph.x*r2) + up*(cph.y*r2))*cth_.x + ax*(h+cth_.y*r2);
            Vec3f na = (left*(cph.x) + up*(cph.y)*cth.x  + ax*cth.y)*-1.0;
            Vec3f nb = (left*(cph.x) + up*(cph.y)*cth_.x + ax*cth_.y)*-1.0;
            //glNormal3f(na.x,na.y,na.z); glVertex3f(pa.x,pa.y,pa.z);
            //glNormal3f(nb.x,nb.y,nb.z); glVertex3f(pb.x,pb.y,pb.z);
            mesh.addQuad( opa, opb, pa, pb );
            opa=pa; opb=pb;
            nvert+=2;
            cph.mul_cmplx(dph);
        }
        //glEnd();
        cth=cth_;
    }
    return nvert;
}
//GLMesh* normals2GLmesh( float sc ){ return vecs2mesh( vpos.size(), &vpos[0], &vnor[0], sc ); }

// ============== ToDo Sphere
// ============== Solids

} // namespace Mesh

#endif


