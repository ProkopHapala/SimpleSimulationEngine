#ifndef  SphereSampling_h
#define  SphereSampling_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"

namespace SphereSampling{

static const double oct_edge_planes[] = {
// Equatorial edges [0..10]
 0.2628655560595669, 0.5257311121191337,-0.8090169943749476,
-0.6881909602355869, 0.5257311121191336, 0.5000000000000000,
 0.85065080835204,   0.5257311121191337, 0.0,
-0.6881909602355869, 0.5257311121191336,-0.5000000000000000,
 0.262865556059567,  0.5257311121191337, 0.8090169943749475,
 0.2628655560595667, 0.5257311121191337,-0.8090169943749476,
-0.6881909602355868, 0.5257311121191336, 0.5000000000000000,
 0.85065080835204,   0.5257311121191336, 0.0,
-0.688190960235587,  0.5257311121191337,-0.5000000000000000,
 0.2628655560595671, 0.5257311121191337, 0.8090169943749475,
// Upper Cap edges
 0.42532540417602,   0.85065080835204,  0.3090169943749475,
-0.1624598481164531, 0.85065080835204,  0.5,
-0.5257311121191336, 0.85065080835204,  0.0,
-0.1624598481164533, 0.85065080835204, -0.5,
 0.42532540417602,   0.85065080835204, -0.3090169943749476,
 // Bottom Cap edges
-0.5257311121191336, 0.85065080835204,  0.0,
-0.1624598481164533, 0.85065080835204, -0.5,
 0.4253254041760200, 0.85065080835204, -0.3090169943749476,
 0.4253254041760201, 0.85065080835204,  0.3090169943749473,
-0.1624598481164531, 0.85065080835204,  0.5,

};

static const double oct_polar_verts[] = {
 // upper
-0.8944271909999159,0.4472135954999579,0.0,
 -0.2763932022500211,0.4472135954999579,-0.8506508083520399,
  0.7236067977499788,0.4472135954999579,-0.5257311121191338,
  0.7236067977499789,0.4472135954999579,0.5257311121191336,
-0.276393202250021,0.4472135954999579,0.85065080835204,
// lower band
-0.7236067977499788,-0.4472135954999579,0.5257311121191337,
-0.723606797749979,-0.4472135954999579,-0.5257311121191335,
 0.2763932022500208,-0.4472135954999579,-0.85065080835204,
 0.8944271909999159,-0.4472135954999579,0,
 0.2763932022500211,-0.4472135954999579,0.8506508083520399,
 // caps
 0.0, 1.0,0.0,
 0.0,-1.0,0.0
};


static const int oct_fac2verts[] = {
0,1,6,10,
1,2,7,10,
2,3,8,10,
3,4,9,10,
4,0,5,10,
5,6,0,11,
6,7,1,11,
7,8,2,11,
8,9,3,11,
9,5,4,11,
};

static const int oct_fac2neigh[] = {
0,1,6,10,
1,2,7,10,
2,3,8,10,
3,4,9,10,
4,0,5,10,
5,6,0,11,
6,7,1,11,
7,8,2,11,
8,9,3,11,
9,5,4,11,
};

void sampleOctahedron( Vec3d p, uint8_t& iface, double& a, double& b ){
    //iface = (p.x>0) | ( ((uint8_t)(p.y>0))<<1) | ( ((uint8_t)(p.z>0) )<<2);
    iface = ( ((uint8_t)(p.x>0))<<2) | ( ((uint8_t)(p.y>0))<<1) | ( ((uint8_t)(p.z>0) ));
    double renorm = 1.0d/( fabs(p.x) + fabs(p.y) + fabs(p.z) );
    a = fabs(p.x);
    b = fabs(p.y);
}


void getIcosaFace( Vec3d p, uint8_t& iface ){
//  http://www.kjmaclean.com/Geometry/Icosahedron.html
//  https://en.wikipedia.org/wiki/Regular_icosahedron
//  rc = sqrt(10+2*sqrt(5))/4 = sqrt(phi^2 + 1)/2 = 0.95105651629 * a
//  plate height = phi/(2*sqrt(phi^2+1)).a = 0.42532540417.a  = 0.4472135955 rc
//  top height   = 1/sqrt(phi^2+1).a       = 0.52573111211.a  = 0.5527864045 rc
    double phi10 = (atan2( p.z, p.x )+ M_PI) * 1.5915494309;  //  1.5915494309 = 10/(2*pi)
    int iphi     = (int)phi10;
    int ioff,i;
    Vec3d& d1=((Vec3d*)oct_edge_planes)[iphi];
    if( d1.dot(p)>0 ){ ioff=0; i=iphi/2; }else{ ioff=5; i=(iphi+1)/2;  if(i>=5) i=0; };
    int i2 = i+1; if(i2>=5) i2=0;
    iface=i+ioff;
    Vec3d& d2=((Vec3d*)oct_edge_planes)[iface+10];
    iface<<=1;
    if( d2.dot(p) > 0 ){ // uper half of rect
        iface++;
    }
}

void sampleIcosa2quads( Vec3d p, uint8_t& iface, double& a, double& b ){
//  http://www.kjmaclean.com/Geometry/Icosahedron.html
//  https://en.wikipedia.org/wiki/Regular_icosahedron
//  rc = sqrt(10+2*sqrt(5))/4 = sqrt(phi^2 + 1)/2 = 0.95105651629 * a
//  plate height = phi/(2*sqrt(phi^2+1)).a = 0.42532540417.a  = 0.4472135955 rc
//  top height   = 1/sqrt(phi^2+1).a       = 0.52573111211.a  = 0.5527864045 rc
    double phi10 = (atan2( p.z, p.x )+ M_PI) * ( 0.15915494309 * 10 );
    int iphi     = (int)phi10;
    double dphi  = phi10 - iphi;
    // 0.4472135955;
    const double hbound = 0.4472135954999579;
    //Vec3d *a,*b,*c;
    int ioff,i;
    Vec3d& d1=((Vec3d*)oct_edge_planes)[iphi];
    if( d1.dot(p) > 0 ){ ioff=0; i=iphi/2; }else{ ioff=5; i=(iphi+1)/2;  if(i>=5) i=0; };
    int i2 = i+1; if(i2>=5) i2=0;
    iface=i+ioff;
    Vec3d& va=((Vec3d*)oct_polar_verts)[iface  ];
    Vec3d& vb=((Vec3d*)oct_polar_verts)[i2+ioff];
    Vec3d& d2=((Vec3d*)oct_edge_planes)[iface+10];
    iface<<=1;
    int ic;
    bool bTop = d2.dot(p) > 0;
    if(bTop){ // uper half of rect
        iface++;
        if(ioff){ ic=i; }else{ic=10;};
    }else{
        if(ioff){ic=11;}else{ ic=i+5+1; if(ic>=10) ic=5; };
    }
    Vec3d& vc = ((Vec3d*)oct_polar_verts)[ic];
    Vec3d cs; cs.fromLinearSolution( va, vb, vc,  p );
    cs.mul(1/(cs.a+cs.b+cs.c));
    if(bTop){
        a=1-cs.b;
        b=1-cs.a;
    }else{
        a = cs.a;
        b = cs.b;
    }
    /*
    if(debugPlot ){
        //Draw3D::drawTriangle(va,vb,vc);
        glColor3f(0.0,0.0,0.0); Draw3D::drawLine(p,vc);
        glColor3f(1.0,0.0,0.0); Draw3D::drawLine(va,vc);
        glColor3f(0.0,0.0,1.0); Draw3D::drawLine(vb,vc);
    }
    */
}


void sampleTri( Vec3d p, Vec3i tri, Vec3d* verts, Vec3d& c ){
    c.a = p.dot( verts[tri.a] );
    c.b = p.dot( verts[tri.b] );
    c.c = p.dot( verts[tri.c] );
    c.mul( 1/(c.a + c.b + c.c) );
}


template<typename T>
void diTri2cartes( T fa, T fb, const Vec3T<T>& a, const Vec3T<T>& b, const Vec3T<T>& c, const Vec3T<T>& d,  Vec3T<T>& p ){
    T fd = (fa+fb)*0.5;
    if( fa>fb ){ T f=(fa-fb); p = d*(fd-0.5*f) + c*(1-fd-0.5*f) + a*f; }
    else       { T f=(fb-fa); p = d*(fd-0.5*f) + c*(1-fd-0.5*f) + b*f; }
};

void icosa2cartes( Vec2i ns, int iface, float fa, float fb, Vec3d& p ){
    //Vec2i n = {nsamp,nsamp};
    int nab = ns.a*ns.b;
    Quat4i* f2v  = (Quat4i*)oct_fac2verts;
    Vec3d*  vs   = (Vec3d*)oct_polar_verts;
    Quat4i& iv   = f2v[iface];
    diTri2cartes<double>( fa, fb, vs[iv.z], vs[iv.w], vs[iv.x], vs[iv.y], p );
};

};

#endif
