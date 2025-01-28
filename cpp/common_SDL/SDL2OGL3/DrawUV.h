#ifndef DrawUV_h
#define DrawUV_h

#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "MeshBuilder2.h"
#include "geometry/UVfuncs.h"

namespace Mesh {

// template<typename UVfunc> Vec3f getUVFuncNormal(Vec2f uv, float h, UVfunc func) { 
//     Vec2f o; Vec3f nor,da,db; 
//     o=uv; o.a+=h; da.set(func(o)); 
//     o=uv; o.a-=h; da.sub(func(o)); 
//     o=uv; o.b+=h; db.set(func(o)); 
//     o=uv; o.b-=h; db.sub(func(o)); 
//     nor.set_cross(db,da); nor.normalize(); 
//     return nor; 
// }

template<typename UVfunc> void UVFunc2smooth(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float voff, UVfunc func) {
    Vec2f duv = UVmax-UVmin; duv.mul({1.0f/n.a,1.0f/n.b});
    std::vector<int> verts((n.a+1)*(n.b+1)); int iv = 0;
    for(int ia=0; ia<=n.a; ia++) { 
        Vec2f uv = {UVmin.a+duv.a*ia, UVmin.b+voff*duv.b*ia}; 
        for(int ib=0; ib<=n.b; ib++) { 
            Vec3d p,nor; 
            convert(func(uv), p); 
            convert(getUVFuncNormal(uv,0.01,func), nor); 
            verts[iv++] = builder.vert(p, nor, Vec2d{uv.x,uv.y}); 
            uv.b += duv.b; 
        }
    }
    for(int ia=0; ia<n.a; ia++) for(int ib=0; ib<n.b; ib++) { 
        int i=ia*(n.b+1)+ib; 
        builder.tri(verts[i], verts[i+1], verts[i+n.b+1]); 
        builder.tri(verts[i+1], verts[i+n.b+2], verts[i+n.b+1]); 
    }
    builder.chunk({builder.tris.size()-2*n.a*n.b, 2*n.a*n.b, -1, (int)Builder2::ChunkType::face});
}

template<typename UVfunc> void UVFunc2wire(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float voff, UVfunc func) {
    Vec2f duv = UVmax-UVmin; duv.mul({1.0f/n.a,1.0f/n.b});
    //std::vector<int> verts((n.a+1)*(n.b+1)); int iv = 0;
    int iv = builder.verts.size();
    for(int ia=0; ia<=n.a; ia++) { 
        Vec2f uv = {UVmin.a+duv.a*ia, UVmin.b+voff*duv.b*ia};
        for(int ib=0; ib<=n.b; ib++) { 
            Vec3d p; 
            convert(func(uv), p); 
            //verts[iv] = 
            builder.vert(p); 
            if(ia<n.a) builder.edge( iv, iv+n.b+1 ); 
            if(ib<n.b) builder.edge( iv, iv+1     ); 
            iv++; 
            uv.b += duv.b; 
        }
    }
    builder.chunk({builder.edges.size()-n.a*(n.b+1)-n.b*(n.a+1), n.a*(n.b+1)+n.b*(n.a+1), -1, (int)Builder2::ChunkType::edgestrip});
}

template<typename UVfunc> void drawExtrudedWireUVFunc(Builder2& builder, Vec2i n, float thick, Vec2f UVmin, Vec2f UVmax, float voff, UVfunc func) {
    Vec2f duv = UVmax-UVmin; duv.mul({1.0f/n.a,1.0f/n.b}); float eps = 0.001;
    std::vector<int> verts((n.a+1)*(n.b+1)*2); int iv = 0;
    for(int ia=0; ia<=n.a; ia++) { 
        Vec2f uv = {UVmin.a+duv.a*ia, UVmin.b+voff*duv.b*ia};
        for(int ib=0; ib<=n.b; ib++) { 
            Vec3d p,nr; 
            convert(func(uv), p); 
            convert(getUVFuncNormal(uv,eps,func), nr);
            verts[iv*2  ] = builder.vert(p); 
            verts[iv*2+1] = builder.vert(p+nr*thick);
            if(ia>0 && ib>0) { 
                int i=iv*2, i_prev=(iv-1)*2, i_up=(iv-(n.b+1))*2;
                builder.tri(verts[i_prev], verts[i  ], verts[i_prev+1]); 
                builder.tri(verts[i     ], verts[i+1], verts[i_prev+1]);
                builder.tri(verts[i_up  ], verts[i  ], verts[i_up  +1]); 
                builder.tri(verts[i     ], verts[i+1], verts[i_up  +1]);
            }
            iv++; 
            uv.b += duv.b;
        }
    }
    builder.chunk({builder.tris.size()-4*n.a*n.b, 4*n.a*n.b, -1, (int)Builder2::ChunkType::face});
}

void Cone2Mesh(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float R1, float R2, float L, float voff, bool wire) { 
    auto uvfunc = [&](Vec2f uv){return ConeUVfunc(uv,R1,R2,L);}; 
    if(wire) UVFunc2wire(builder, n,UVmin,UVmax,voff,uvfunc); 
    else     UVFunc2smooth(builder, n,UVmin,UVmax,voff,uvfunc); 
}

void Sphere2Mesh(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float R, float voff, bool wire) { 
    auto uvfunc = [&](Vec2f uv){return SphereUVfunc(uv,R);}; 
    if(wire) UVFunc2wire(builder, n,UVmin,UVmax,voff,uvfunc); 
    else     UVFunc2smooth(builder, n,UVmin,UVmax,voff,uvfunc); 
}

void Torus2Mesh(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float r, float R, float voff, bool wire) { 
    auto uvfunc = [&](Vec2f uv){return TorusUVfunc(uv,r,R);}; 
    if(wire) UVFunc2wire(builder, n,UVmin,UVmax,voff,uvfunc); 
    else     UVFunc2smooth(builder, n,UVmin,UVmax,voff,uvfunc); 
}

void Teardrop2Mesh(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float R1, float R2, float L, float voff, bool wire) { 
    auto uvfunc = [&](Vec2f uv){return TeardropUVfunc(uv,R1,R2,L);}; 
    if(wire) UVFunc2wire(builder, n,UVmin,UVmax,voff,uvfunc); 
    else     UVFunc2smooth(builder, n,UVmin,UVmax,voff,uvfunc); 
}

void NACASegment2Mesh(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float *coefs1, float *coefs2, float L, float voff, bool wire) { 
    auto uvfunc = [&](Vec2f uv){return NACA4digitUVfunc(uv,coefs1,coefs2,L);}; 
    if(wire) UVFunc2wire(builder, n,UVmin,UVmax,voff,uvfunc); 
    else     UVFunc2smooth(builder, n,UVmin,UVmax,voff,uvfunc); 
}

void HarmonicTube2Mesh(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float R1, float R2, float L, float voff, float freq, float amp, bool wire) { 
    auto uvfunc = [&](Vec2f uv){return HarmonicTubeUVfunc(uv,R1,R2,L,freq,amp);}; 
    if(wire) UVFunc2wire(builder, n,UVmin,UVmax,voff,uvfunc); 
    else     UVFunc2smooth(builder, n,UVmin,UVmax,voff,uvfunc); 
}

void Parabola2Mesh(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float R, float L, float voff, bool wire) { 
    float K = L/(R*R); 
    UVmin.a*=R; UVmax.a*=R; 
    auto uvfunc = [&](Vec2f uv){return ParabolaUVfunc(uv,K);}; 
    if(wire) UVFunc2wire(builder, n,UVmin,UVmax,voff,uvfunc); 
    else     UVFunc2smooth(builder, n,UVmin,UVmax,voff,uvfunc); 
}

void Hyperbola2Mesh(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float r, float R, float L, float voff, bool wire) {
    if(r>0){ 
        float K = R/L; 
        UVmin.a*=L; UVmax.a*=L; 
        auto uvfunc = [&](Vec2f uv){return HyperbolaLUVfunc(uv,r,K);}; 
        if(wire) UVFunc2wire(builder, n,UVmin,UVmax,voff,uvfunc); 
        else     UVFunc2smooth(builder, n,UVmin,UVmax,voff,uvfunc); 
    }else{ 
        r=-r; 
        float K = L/R; 
        UVmin.a*=R; UVmax.a*=R; 
        auto uvfunc = [&](Vec2f uv){return HyperbolaRUVfunc(uv,r,K);}; 
        if(wire) UVFunc2wire(builder, n,UVmin,UVmax,voff,uvfunc); 
        else     UVFunc2smooth(builder, n,UVmin,UVmax,voff,uvfunc); 
    }
}

void Cone_ExtrudedWire(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float R1, float R2, float L, float voff, float thick) { 
    auto uvfunc = [&](Vec2f uv){return ConeUVfunc(uv,R1,R2,L);}; 
    drawExtrudedWireUVFunc(builder, n,thick,UVmin,UVmax,voff,uvfunc); 
}

void Sphere_ExtrudedWire(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float R, float voff, float thick) { 
    auto uvfunc = [&](Vec2f uv){return SphereUVfunc(uv,R);}; 
    drawExtrudedWireUVFunc(builder, n,thick,UVmin,UVmax,voff,uvfunc); 
}

void Torus_ExtrudedWire(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float r, float R, float voff, float thick) { 
    auto uvfunc = [&](Vec2f uv){return TorusUVfunc(uv,r,R);}; 
    drawExtrudedWireUVFunc(builder, n,thick,UVmin,UVmax,voff,uvfunc); 
}

void Teardrop_ExtrudedWire(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float R1, float R2, float L, float voff, float thick) { 
    auto uvfunc = [&](Vec2f uv){return TeardropUVfunc(uv,R1,R2,L);}; 
    drawExtrudedWireUVFunc(builder, n,thick,UVmin,UVmax,voff,uvfunc); 
}

void NACASegment_ExtrudedWire(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float *coefs1, float *coefs2, float L, float voff, float thick) { 
    auto uvfunc = [&](Vec2f uv){return NACA4digitUVfunc(uv,coefs1,coefs2,L);}; 
    drawExtrudedWireUVFunc(builder, n,thick,UVmin,UVmax,voff,uvfunc); 
}

void HarmonicTube_ExtrudedWire(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float R1, float R2, float L, float voff, float freq, float amp, float thick) { 
    auto uvfunc = [&](Vec2f uv){return HarmonicTubeUVfunc(uv,R1,R2,L,freq,amp);}; 
    drawExtrudedWireUVFunc(builder, n,thick,UVmin,UVmax,voff,uvfunc); 
}

void Parabola_ExtrudedWire(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float R, float L, float voff, float thick) { 
    float K = L/(R*R); UVmin.a*=R; UVmax.a*=R; 
    auto uvfunc = [&](Vec2f uv){return ParabolaUVfunc(uv,K);}; 
    drawExtrudedWireUVFunc(builder, n,thick,UVmin,UVmax,voff,uvfunc); 
}

void Hyperbola_ExtrudedWire(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float r, float R, float L, float voff, float thick) {
    if(r>0){ 
        float K = R/L; UVmin.a*=L; UVmax.a*=L; 
        auto uvfunc = [&](Vec2f uv){return HyperbolaLUVfunc(uv,r,K);}; 
        drawExtrudedWireUVFunc(builder, n,thick,UVmin,UVmax,voff,uvfunc); 
    }else{ 
        r=-r; float K = L/R; UVmin.a*=R; UVmax.a*=R; 
        auto uvfunc = [&](Vec2f uv){return HyperbolaRUVfunc(uv,r,K);}; 
        drawExtrudedWireUVFunc(builder, n,thick,UVmin,UVmax,voff,uvfunc); 
    }
}

void ConeFan(Builder2& builder, int n, float r, const Vec3f& base, const Vec3f& tip) {
    Vec3f a,b,c,c_hat; c.set_sub(tip,base); c_hat.set_mul(c,1/c.norm()); c_hat.getSomeOrtho(a,b); a.normalize(); b.normalize();
    float alfa=2*M_PI/n; Vec2f rot; rot.set(1.0f,0.0f); Vec2f drot; drot.set(cos(alfa),sin(alfa));
    Vec3f q=c; q.add_mul(a,-r); float pnab=c_hat.dot(q)/q.norm(), pnc=sqrt(1-pnab*pnab);
    int i0=builder.verts.size(); builder.vert(Vec3d{tip.x,tip.y,tip.z});
    for(int i=0; i<=n; i++) {
        Vec3f p; p.set(rot.x*a.x+rot.y*b.x, rot.x*a.y+rot.y*b.y, rot.x*a.z+rot.y*b.z);
        Vec3f pn; pn.set(pnab*p.x+pnc*c_hat.x, pnab*p.y+pnc*c_hat.y, pnab*p.z+pnc*c_hat.z);
        builder.vert(Vec3d{base.x+r*p.x,base.y+r*p.y,base.z+r*p.z}, Vec3d{pn.x,pn.y,pn.z});
        if(i>0) builder.tri(i0,i0+i,i0+i+1);
        rot.mul_cmplx(drot);
    }
    builder.chunk({builder.tris.size()-n, n, -1, (int)Builder2::ChunkType::face});
}

void CylinderStrip(Builder2& builder, int n, float r1, float r2, const Vec3f& base, const Vec3f& tip) {
    Vec3f a,b,c,c_hat; c.set_sub(tip,base); c_hat.set_mul(c,1/c.norm()); c_hat.getSomeOrtho(a,b); a.normalize(); b.normalize();
    float alfa=2*M_PI/n; Vec2f rot; rot.set(1.0f,0.0f); Vec2f drot; drot.set(cos(alfa),sin(alfa));
    Vec3f q=c; q.add_mul(a,-(r1-r2)); float pnab=c_hat.dot(q)/q.norm(), pnc=sqrt(1-pnab*pnab);
    int i0=builder.verts.size();
    for(int i=0; i<=n; i++) {
        Vec3f p; p.set(rot.x*a.x+rot.y*b.x, rot.x*a.y+rot.y*b.y, rot.x*a.z+rot.y*b.z);
        Vec3f pn; pn.set(pnab*p.x+pnc*c_hat.x, pnab*p.y+pnc*c_hat.y, pnab*p.z+pnc*c_hat.z);
        builder.vert(Vec3d{base.x+r1*p.x,base.y+r1*p.y,base.z+r1*p.z}, Vec3d{pn.x,pn.y,pn.z});
        builder.vert(Vec3d{tip.x+r2*p.x,tip.y+r2*p.y,tip.z+r2*p.z}, Vec3d{pn.x,pn.y,pn.z});
        if(i>0) { builder.tri(i0+(i-1)*2,i0+i*2,i0+(i-1)*2+1); builder.tri(i0+i*2,i0+i*2+1,i0+(i-1)*2+1); }
        rot.mul_cmplx(drot);
    }
    builder.chunk({builder.tris.size()-2*n, 2*n, -1, (int)Builder2::ChunkType::face});
}

void CylinderStrip_wire(Builder2& builder, int n, float r1, float r2, const Vec3f& base, const Vec3f& tip) {
    Vec3f a,b,c,c_hat; c.set_sub(tip,base); c_hat.set_mul(c,1/c.norm()); c_hat.getSomeOrtho(a,b); a.normalize(); b.normalize();
    float alfa=2*M_PI/n; Vec2f rot; rot.set(1.0f,0.0f); Vec2f drot; drot.set(cos(alfa),sin(alfa));
    int i0=builder.verts.size();
    for(int i=0; i<n; i++) {
        Vec3f p; p.set(rot.x*a.x+rot.y*b.x, rot.x*a.y+rot.y*b.y, rot.x*a.z+rot.y*b.z);
        builder.vert(Vec3d{base.x+r1*p.x,base.y+r1*p.y,base.z+r1*p.z});
        builder.vert(Vec3d{tip.x+r2*p.x,tip.y+r2*p.y,tip.z+r2*p.z});
        if(i>0) { builder.edge(i0+(i-1)*2,i0+i*2); builder.edge(i0+(i-1)*2+1,i0+i*2+1); }
        builder.edge(i0+i*2,i0+i*2+1);
        rot.mul_cmplx(drot);
    }
    builder.edge(i0+(n-1)*2,i0); builder.edge(i0+(n-1)*2+1,i0+1);
    builder.chunk({builder.edges.size()-3*n, 3*n, -1, (int)Builder2::ChunkType::edgestrip});
}

void SphereTriangle_wire(Builder2& builder, int n, float r, const Vec3f& pos, const Vec3f& a, const Vec3f& b, const Vec3f& c) {
    float d=1.0f/n; Vec3f da(a-c),db(b-c); da.mul(d); db.mul(d);
    int i0=builder.verts.size();
    for(int ia=0; ia<=n; ia++) {
        Vec3f p0=c; p0.add_mul(da,ia);
        for(int ib=0; ib<=n-ia; ib++) {
            Vec3f p=p0; p.add_mul(db,ib); p.mul(r/p.norm());
            builder.vert(Vec3d{pos.x+p.x,pos.y+p.y,pos.z+p.z});
            if(ia<n && ib<n-ia) {
                int i1=i0+ia*(2*n-ia+1)/2+ib, i2=i0+(ia+1)*(2*n-ia)/2+ib;
                builder.edge(i1,i2); if(ib<n-ia-1) builder.edge(i2,i2+1);
            }
        }
    }
    builder.chunk({builder.edges.size()-2*n*n, 2*n*n, -1, (int)Builder2::ChunkType::edgestrip});
}

void SphereTriangle(Builder2& builder, int n, float r, const Vec3f& pos, const Vec3f& a, const Vec3f& b, const Vec3f& c) {
    float d=1.0f/n; Vec3f da(a-c),db(b-c); da.mul(d); db.mul(d);
    int i0=builder.verts.size();
    for(int ia=0; ia<=n; ia++) {
        Vec3f p0=c; p0.add_mul(da,ia);
        for(int ib=0; ib<=n-ia; ib++) {
            Vec3f p=p0; p.add_mul(db,ib); p.mul(r/p.norm());
            builder.vert(Vec3d{pos.x+p.x,pos.y+p.y,pos.z+p.z}, Vec3d{p.x/r,p.y/r,p.z/r});
            if(ia<n && ib<n-ia) {
                int i1=i0+ia*(2*n-ia+1)/2+ib, i2=i0+(ia+1)*(2*n-ia)/2+ib;
                builder.tri(i1,i2,i1+1); if(ib<n-ia-1) builder.tri(i2,i2+1,i1+1);
            }
        }
    }
    builder.chunk({builder.tris.size()-n*n, n*n, -1, (int)Builder2::ChunkType::face});
}

void Sphere_oct(Builder2& builder, int n, float r, const Vec3f& pos, bool wire) {
    Vec3f a,b,c; 
    a.set(1.0f,0.0f,0.0f); b.set(0.0f,1.0f,0.0f); c.set(0.0f,0.0f,1.0f);
    Vec3f na=a; na.mul(-1.0f);
    Vec3f nb=b; nb.mul(-1.0f);
    Vec3f nc=c; nc.mul(-1.0f);
    if(wire) {
        SphereTriangle_wire(builder,n,r,pos, a, b, c); SphereTriangle_wire(builder,n,r,pos,na, b, c);
        SphereTriangle_wire(builder,n,r,pos, a,nb, c); SphereTriangle_wire(builder,n,r,pos,na,nb, c);
        SphereTriangle_wire(builder,n,r,pos, a, b,nc); SphereTriangle_wire(builder,n,r,pos,na, b,nc);
        SphereTriangle_wire(builder,n,r,pos, a,nb,nc); SphereTriangle_wire(builder,n,r,pos,na,nb,nc);
    } else {
        SphereTriangle(builder,n,r,pos, a, b, c); SphereTriangle(builder,n,r,pos,na, b, c);
        SphereTriangle(builder,n,r,pos, a,nb, c); SphereTriangle(builder,n,r,pos,na,nb, c);
        SphereTriangle(builder,n,r,pos, a, b,nc); SphereTriangle(builder,n,r,pos,na, b,nc);
        SphereTriangle(builder,n,r,pos, a,nb,nc); SphereTriangle(builder,n,r,pos,na,nb,nc);
    }
}

void Capsula(Builder2& builder, Vec3f p0, Vec3f p1, float r1, float r2, float theta1, float theta2, float dTheta, int nPhi, bool capped) {
    Vec3f ax=p1-p0; float L=ax.normalize();
    Vec3f up,lf; ax.getSomeOrtho(up,lf);
    float dPhi=2*M_PI/nPhi, R=(r1+r2)*0.5f, dR=(r2-r1)*0.5f;
    int nTheta=(theta2-theta1)/dTheta;
    int i0=builder.verts.size();
    for(int iTheta=0; iTheta<=nTheta; iTheta++) {
        float cT=cos(theta1+iTheta*dTheta), sT=sin(theta1+iTheta*dTheta);
        float r=R+dR*cT, z=L*0.5f*sT;
        for(int iPhi=0; iPhi<=nPhi; iPhi++) {
            float cP=cos(iPhi*dPhi), sP=sin(iPhi*dPhi);
            Vec3f d=up*cP+lf*sP;
            Vec3f p=p0+ax*((L+z)*0.5f)+d*r;
            Vec3f n=d+ax*sT; n.normalize();
            builder.vert(Vec3d{p.x,p.y,p.z}, Vec3d{n.x,n.y,n.z});
            if(iTheta<nTheta && iPhi<nPhi) {
                int i=i0+iTheta*(nPhi+1)+iPhi;
                builder.tri(i,i+1,i+nPhi+1); builder.tri(i+1,i+nPhi+2,i+nPhi+1);
            }
        }
    }
    if(capped) {
        int ic1=builder.verts.size(), ic2=ic1+1;
        builder.vert(Vec3d{p0.x,p0.y,p0.z}, Vec3d{-ax.x,-ax.y,-ax.z});
        builder.vert(Vec3d{p1.x,p1.y,p1.z}, Vec3d{ax.x,ax.y,ax.z});
        for(int iPhi=0; iPhi<nPhi; iPhi++) {
            builder.tri(i0+iPhi,i0+iPhi+1,ic1);
            builder.tri(i0+(nTheta)*(nPhi+1)+iPhi,ic2,i0+(nTheta)*(nPhi+1)+iPhi+1);
        }
    }
    builder.chunk({builder.tris.size()-(2*nTheta*nPhi+(capped?2*nPhi:0)), 2*nTheta*nPhi+(capped?2*nPhi:0), -1, (int)Builder2::ChunkType::face});
}

} // namespace Mesh

#endif
