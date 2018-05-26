
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "fastmath.h"
#include "Vec2.h"
#include "Solids.h"

#include "Draw.h"
#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"

#include "SDL_utils.h"
#include "GLUtils.h"

// ======================  TestApp
const int nsamp = 64*4;
//const int nsamp = 8;
const int npix  = 8*nsamp*nsamp;
uint32_t pix   [npix];


bool debugPlot = false;
Vec3d curPos = (Vec3d){5.0,0.0,0.0};

const double oct_edge_planes[] = {
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

const double oct_polar_verts[] = {
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

void sampleOctahedron( Vec3d p, uint8_t& iface, double& a, double& b ){
    //iface = (p.x>0) | ( ((uint8_t)(p.y>0))<<1) | ( ((uint8_t)(p.z>0) )<<2);
    iface = ( ((uint8_t)(p.x>0))<<2) | ( ((uint8_t)(p.y>0))<<1) | ( ((uint8_t)(p.z>0) ));
    double renorm = 1.0d/( fabs(p.x) + fabs(p.y) + fabs(p.z) );
    a = fabs(p.x);
    b = fabs(p.y);
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
    if(d2.dot(p) > 0){ // uper half of rect
        iface++;
        if(ioff){ ic=i; }else{ic=10;};
    }else{
        if(ioff){ic=11;}else{ ic=i+5+1; if(ic>=10) ic=5; };
    }
    Vec3d& vc = ((Vec3d*)oct_polar_verts)[ic];
    Vec3d pvc = p-vc;
    const double renorm = 0.95105651629; // 1/|va-vc|
    a = pvc.dot(va-vc)*renorm;
    b = pvc.dot(vb-vc)*renorm;

    if(debugPlot ){
        //Draw3D::drawTriangle(va,vb,vc);
        glColor3f(0.0,0.0,0.0); Draw3D::drawLine(p,vc);
        glColor3f(1.0,0.0,0.0); Draw3D::drawLine(va,vc);
        glColor3f(0.0,0.0,1.0); Draw3D::drawLine(vb,vc);
    }
}




void sampleTri( Vec3d p, Vec3i tri, Vec3d* verts, Vec3d& c ){
    c.a = p.dot( verts[tri.a] );
    c.b = p.dot( verts[tri.b] );
    c.c = p.dot( verts[tri.c] );
    c.mul( 1/(c.a + c.b + c.c) );
}

class TestAppSolids : public AppSDL2OGL_3D {
	public:

	int point_cloud;
	int shape;

	char str[2048];
	int  fontTex;

	// ---- function declarations

	virtual void draw   ();

	TestAppSolids( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppSolids::TestAppSolids( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
/*
    //normals =
	shape=glGenLists(1);
	glNewList( shape, GL_COMPILE );
        printf( " Solids::nTetrahedron_tris %i \n", Solids::tetrahedron.nTris );
        glEnable    ( GL_LIGHTING );
        glShadeModel( GL_FLAT     );
        glColor3f( 0.8, 0.8, 0.8 );
        Draw3D::drawTriangles( Solids::tetrahedron.nTris, (int*)(&Solids::tetrahedron.tris[0][0]), Solids::tetrahedron.verts );
	glEndList();
*/


/*
// TEST BY PYTHON
https://www.tutorialspoint.com/execute_python_online.php
import numpy as np
mat = np.array([
  [1.,0.0,-0.1],
  [1.,-1.0,0.0],
  [1.0,1.0,1.0],
])
print np.linalg.solve(mat, [0.5,0.9,0.5])
print np.linalg.solve(mat.transpose(), [0.5,0.9,0.5])
*/

    Vec3d va,vb,vc,p,o;
    va.set(1.,0.0,-0.1);
    vb.set(1.,-1.0,0.0);
    vc.set(1.0,1.0,1.0);
    p .set(0.5,0.9,0.5);
    o.fromLinearSolution(va,vb,vc,p);
    print(o);


    //exit(0);


    for( int i=0; i<npix; i++ ){ pix[i]=0; }

    debugPlot = false;
	shape=glGenLists(1);
	glNewList( shape, GL_COMPILE );
        //printf( " Solids::nTetrahedron_tris %i \n", Solids::nTetrahedron_tris );

        glPushMatrix();
        glDisable ( GL_LIGHTING );
        glColor3f( 1.0f, 0.0f, 1.0f );
        /*
        Draw3D::drawLines    ( Solids::Tetrahedron_nedges, Solids::Tetrahedron_edges, Solids::Tetrahedron_verts );
        glTranslatef(  2.0f, 0.0f, 0.0f );
        Draw3D::drawLines    ( Solids::Octahedron_nedges, Solids::Octahedron_edges, Solids::Octahedron_verts );
        glTranslatef( -4.0f, 0.0f, 0.0f );
        Draw3D::drawLines    ( Solids::Cube_nedges, Solids::Cube_edges, Solids::Cube_verts );
        */
        //Draw3D::drawLines    ( Solids::RhombicDodecahedron_nedges, Solids::RhombicDodecahedron_edges, Solids::RhombicDodecahedron_verts );
        //Draw3D::drawLines    ( Solids::Icosahedron_nedges, (int*)Solids::Icosahedron_edges, Solids::Icosahedron_verts );
        glPopMatrix();

        glPushMatrix();
        glEnable    ( GL_LIGHTING );
        glShadeModel( GL_FLAT     );
        glColor3f( 0.8f, 0.8f, 0.8f );
        //Draw3D::drawTriangles( Solids::Tetrahedron_ntris,  Solids::Tetrahedron_tris,  Solids::Tetrahedron_verts );

        //Draw3D::drawPolygons( Solids::Tetrahedron_nfaces, Solids::Tetrahedron_ngons, Solids::Tetrahedron_faces, Solids::Tetrahedron_verts );
        //glTranslatef(  2.0f, 0.0f, 0.0f );
        //Draw3D::drawPolygons( Solids::Octahedron_nfaces,  Solids::Octahedron_ngons,  Solids::Octahedron_faces,  Solids::Octahedron_verts  );
        //glTranslatef(  -4.0f, 0.0f, 0.0f );
        //Draw3D::drawPolygons( Solids::Cube_nfaces,        Solids::Cube_ngons,        Solids::Cube_faces,        Solids::Cube_verts        );

        //Draw3D::drawPolygons( Solids::RhombicDodecahedron_nfaces,        Solids::RhombicDodecahedron_ngons,        Solids::RhombicDodecahedron_faces,        Solids::RhombicDodecahedron_verts        );

        glScalef(0.5,0.5,0.5);
        //glRotatef( 180.0, 0,1,0 );
        //glRotatef( -90.0, 0,1,0 );
        glRotatef( 31.7174744, 0,0,1 ); // arctan (phi-1) in degreers
        //Draw3D::drawPolygons( Solids::Icosahedron_nfaces,        Solids::Icosahedron_ngons,        Solids::Icosahedron_faces,        Solids::Icosahedron_verts        );
        glDisable( GL_LIGHTING ); glColor3f(1.0,1.0,1.0); Draw3D::drawLines    ( Solids::Icosahedron_nedges, (int*)Solids::Icosahedron_edges, Solids::Icosahedron_verts );
        glRotatef( -31.7174744, 0,0,1 ); // arctan (phi-1) in degreers
        // Vec3i ivs = Solids::Octahedron_tris[iface];


        //double sc = 1.6180339887498948482;
        //double sc = 1.70;
        double sc = 1.90211303259;
        glColor3f(0.0,0.0,0.0);
        glDisable( GL_LIGHTING );
        double dphi = (M_PI*2)/10;

        double tgt = 0.5;
        double ct = 1/sqrt(1+tgt*tgt);
        double st = ct*tgt;
        for( int i=0; i<10; i++ ){
            double sign = 2*(i&1)-1;
            double phi = dphi*i;
            double ca = cos(phi);
            double sa = sin(phi);
            Vec3d a = (Vec3d){ ct*cos(phi),        st*sign  ,ct*sin(phi)};
            Vec3d b = (Vec3d){ ct*cos(phi+dphi),  -st*sign  ,ct*sin(phi+dphi)};
            //printf("|a| %f |b| %f \n", a.norm(), b.norm() );
            //printf( "%.16g,%.16g,%.16g,\n", a.x, a.y, a.z );
            printf("|a| %.16g |b| %.16g |a-b| %.16g \n", a.norm(), b.norm(), (a-b).norm() );
            //printf( "%.16g,%.16g,%.16g,\n", c.x, c.y, c.z );
            Vec3d c; c.set_cross(a+b,a-b);
            //Draw3D::drawLine( a*sc, b*sc);
            //Draw3D::drawVecInPos( c, (a+b)*sc*0.5 );
            c.normalize();
            //printf( "%.16g,%.16g,%.16g,\n", c.x, c.y, c.z );
        }

        for( int i=0; i<5; i++ ){
            double phi = 2*dphi*i  + dphi*2;
            double ca = cos(phi);
            double sa = sin(phi);
            Vec3d a = (Vec3d){ ct*cos(phi),        -st  ,ct*sin(phi)     };
            Vec3d b = (Vec3d){ ct*cos(phi+2*dphi), -st  ,ct*sin(phi+2*dphi)};
            //printf("|a| %f |b| %f \n", a.norm(), b.norm() );
            //printf( "%.16g,%.16g,%.16g,\n", a.x, a.y, a.z );
            //printf( "%.16g,%.16g,%.16g,\n", c.x, c.y, c.z );

            Vec3d c; c.set_cross(a+b,a-b);
            Draw3D::drawLine( a*sc, b*sc);
            Draw3D::drawVecInPos( c, (a+b)*sc*0.5 );
            c.normalize();
            printf( "%.16g,%.16g,%.16g,\n", c.x, c.y, c.z );
        }


        glDisable( GL_LIGHTING );
        glPointSize(1);
        for( int i=0; i<Solids::Octahedron_ntris; i++ ){
            Vec3i ivs = Solids::Octahedron_tris[i];
            //p   = ( Solids::Octahedron_verts[ivs.a] + Solids::Octahedron_verts[ivs.b] + Solids::Octahedron_verts[ivs.c] ) * 0.3;
            //sprintf( str, "%i (%i,%i,%i)", i, (int)(p.x>0), (int)(p.y>0), (int)(p.z>0) );
            //sprintf( str, "%i_%i%i%i", i, (int)(p.x>0), (int)(p.y>0), (int)(p.z>0) );
            //Draw3D::drawText(str, p, fontTex, 0.03, 0);

            uint32_t * pixF = pix + nsamp*nsamp*i;

            glBegin(GL_POINTS);
            Vec3d& a = Solids::Octahedron_verts[ivs.a];
            Vec3d& b = Solids::Octahedron_verts[ivs.b];
            Vec3d& c = Solids::Octahedron_verts[ivs.c];
            float step = 1.0/nsamp;
            for( int ia = 0; ia<nsamp; ia++ ){
                float ca = ia*step;
                for( int ib = 0; ib<(nsamp-ia); ib++ ){
                    float cb = ib*step;
                    float cc = 1-ca-cb;
                    Vec3d p = (a*ca + b*cb + c*cc);
                    p.normalize();
                    //Draw::setRGB ( pixF[ ia*nsamp + ib ] );

                    double fa,fb;
                    uint8_t iface;
                    //sampleIcosa( p, iface, fa, fb );
                    sampleIcosa2quads( p, iface, fa, fb );
                    //Draw::color_of_hash( (iface+1545)*151 );

                    glColor3f(fa,fb,1-fa-fb);

                    glVertex3f( sc*(float)p.x, sc*(float)p.y, sc*(float)p.z );
                    //printf( "%i: %i %i %i | (%f,%f,%f) \n", i, ia, ib, iface, p.x, p.y, p.z );
                }
            }
            glEnd();

        }


        glPopMatrix();

	glEndList();

	zoom = 3.0;

}

void TestAppSolids::draw   (){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glEnable(GL_DEPTH_TEST);
	glCallList( shape );


	glColor3f( 0.0f,1.0f,1.0f );

	curPos = camMat.c*-3.0;
	//curPos = (Vec3d){0.2,0.0,1.0};

	Draw3D::drawPointCross( curPos, 0.1 );



	uint8_t iface; double a,b;
	sampleOctahedron( curPos, iface, a, b );

	//printf( "iface %i | (%f,%f,%f) \n", iface, p.x, p.y, p.z );

	Vec3i ivs = Solids::Octahedron_tris[iface];
	Vec3d c;
	sampleTri( curPos, ivs, Solids::Octahedron_verts, c );

    if ( SDL_GetMouseState(NULL, NULL) & SDL_BUTTON(SDL_BUTTON_LEFT) ) {
        //SDL_Log("Mouse Button 1 (left) is pressed.");
        uint32_t * pixF = pix + nsamp*nsamp*iface;
        //float step = 1.0/nsamp;
        int ia = (int)(c.a*nsamp);
        int ib = (int)(c.b*nsamp);
        int i = ia*nsamp + ib;
        pixF[i] = 0xffffffff;
    }



    debugPlot = true;
    sampleIcosa2quads( curPos, iface, a, b );
    sprintf( str, "(%.3f,%.3f)", a,b );   Draw3D::drawText(str,curPos, fontTex, 0.02, 0);
    debugPlot = false;

	double fsc = 1.1;
	Vec3d p_;
	p_.set_lincomb( c.a,c.b,c.c,  Solids::Octahedron_verts[ivs.a], Solids::Octahedron_verts[ivs.b], Solids::Octahedron_verts[ivs.c] );
	p_.mul(fsc);
	Draw3D::drawLine( Solids::Octahedron_verts[ivs.a]*fsc, p_ );
	Draw3D::drawLine( Solids::Octahedron_verts[ivs.b]*fsc, p_ );
	Draw3D::drawLine( Solids::Octahedron_verts[ivs.c]*fsc, p_ );


    for( int i=0; i<5; i++ ){
        Vec3d p;
        glColor3f(1.0,0.0,0.0);
        p=  ((Vec3d*)oct_polar_verts)[i];
        sprintf( str, "%i", i );   Draw3D::drawText(str, p, fontTex, 0.02, 0);
        p =  ((Vec3d*)oct_polar_verts)[i+5];
        sprintf( str, "%i", i+5 ); Draw3D::drawText(str, p, fontTex, 0.01, 0);
    }


	glColor3f( 1.0f,0.0f,1.0f );
	fsc = 1.1;
    //Draw3D::drawTriangle( Solids::Octahedron_verts[ivs.a]*fsc, Solids::Octahedron_verts[ivs.b]*fsc, Solids::Octahedron_verts[ivs.c]*fsc );

    glDisable ( GL_LIGHTING );

    glColor3f( 1.0f,0.0f,0.0f );


    //exit(0);

    /*
    glColor3f( 0.0f,1.0f,0.0f );
    for( int i=0; i<Solids::Octahedron_nverts; i++ ){
        sprintf( str, "%i", i );
        Draw3D::drawText(str, Solids::Octahedron_verts[i], fontTex, 0.03, 0);
    }
    */

    //if( SDL_MOUSEBUTTONDOWN ){ printf("mouse down %i \n", SDL_MOUSEBUTTONDOWN  ); };



	Draw3D::drawAxis ( 3.0f );

};

// ===================== MAIN

TestAppSolids * testApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestAppSolids( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}
















