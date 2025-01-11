
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Solids.h"
#include "Noise.h"
#include "SphereSampling.h"

#include "Draw.h"
#include "Draw3D.h"
#include "DrawSphereMap.h"
#include "AppSDL2OGL_3D.h"

#include "SDL_utils.h"
#include "GLUtils.h"

using namespace SphereSampling;

// ======================  TestApp
//const int nsamp = 8;
const int nsamp = 32;
//const int nsamp = 8;
const int npix  = 10*nsamp*nsamp;
uint32_t pix   [npix];

float heights[npix];

Vec3f samplePs[npix];

const int nCrater = 50;
Vec3d  craterPos[nCrater];
double craterSz[nCrater];

bool debugPlot = false;
Vec3d curPos = (Vec3d){5.0,0.0,0.0};


char str[2048];
//int  fontTex;

class TestAppSphereSampling : public AppSDL2OGL_3D {
	public:

	int point_cloud;
	int shape;

	// ---- function declarations

	virtual void draw   ();

	TestAppSphereSampling( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppSphereSampling::TestAppSphereSampling( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    //glEnable(GL_AUTO_NORMALS);

    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );


    /*
    for( int i=0; i<npix; i++ ){
        pix[i]=0;


        Vec3d p;
        icosa2cartes( int iface, float fa, float fb, p );

        heights[i] = randf();
    }
    */

    srand(34646);
    for(int i=0; i<nCrater; i++){
        //craterPos[i].setHomogenousSphericalSample( randf(), randf() );
        craterPos[i].fromRandomSphereSample();
        //printf( "%i %f %f %f \n", craterPos[i].x, craterPos[i].y, craterPos[i].z);
        craterSz[i] = randf()+0.4;
    }

    Vec2d dab = (Vec2d){ 1.0/nsamp, 1.0/nsamp };
    for(int iface=0; iface<10; iface++ ){
        for(int ia=0; ia<nsamp; ia++ ){
            for(int ib=0; ib<nsamp; ib++ ){
                Vec3d p;
                icosa2cartes( (Vec2i){nsamp,nsamp}, iface, ia*dab.a, ib*dab.b, p );
                int i = iface*nsamp*nsamp + ia*nsamp + ib;

                //heights[i] = getSphericalHarmonicRand( {5,4}, p );

                p.normalize();
                double h;
                h = Noise::getQuadruRand( p, 35456 )*2*0;

                h += Noise::getCraterHeight( p, nCrater, 1.0, craterPos, craterSz )*0.3;
                heights[i] = h; //+ sin( h*4 );

            }
        }
    };

    debugPlot = false;
	shape=glGenLists(1);
	glNewList( shape, GL_COMPILE );
        //printf( " Solids::nTetrahedron_tris %i \n", Solids::nTetrahedron_tris );

        for(int i=0; i<nCrater; i++){
            Draw3D::drawPointCross( craterPos[i], craterSz[i] );
        }

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

        /*
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
        */

        /*

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
        */


        /*
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
        */


        glScalef(sc,sc,sc);
        //glEnable(GL_LIGHTING);
        glDisable(GL_LIGHTING);
        glShadeModel( GL_SMOOTH     );
        glColor3f(1.0,1.0,1.0);
        Vec3d* vs = ((Vec3d*)oct_polar_verts);
        //drawDiTri( (Vec2i){nsamp,nsamp},  (Vec3f)vs[6], (Vec3f)vs[5], (Vec3f)vs[0], (Vec3f)vs[11], height );

        //drawDiTri_Inner( (Vec2i){nsamp,nsamp},  (Vec3f)vs[6], (Vec3f)vs[5], (Vec3f)vs[0], (Vec3f)vs[11], height );
        drawIcosaMap( (Vec2i){nsamp,nsamp}, heights, 0.1 );

        glPopMatrix();

	glEndList();

	zoom = 3.0;

}

void TestAppSphereSampling::draw   (){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glEnable(GL_DEPTH_TEST);
	glCallList( shape );


	glColor3f( 0.0f,1.0f,1.0f );

	curPos = (Vec3d)(cam.rot.c*-3.0);
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
        sprintf( str, "%i", i+5 ); Draw3D::drawText(str, p, fontTex, 0.02, 0);
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


    /*
    Quat4i* f2v = (Quat4i*)oct_fac2verts;
    Quat4i* f2n = (Quat4i*)oct_fac2neigh;
    Vec3d* vs = (Vec3d*)oct_polar_verts;
    Quat4i& iv = f2v[0];
    drawDiTri( {nsamp,nsamp}, (Vec3f)vs[iv.z], (Vec3f)vs[iv.w], (Vec3f)vs[iv.x], (Vec3f)vs[iv.y],  heights  );

    drawDiTri_seam( nsamp, nsamp, (Vec3f)vs[iv.x], (Vec3f)vs[iv.z], (Vec3f)vs[iv.y], heights, heights, {1,nsamp*(nsamp-1)}, {1,0} );
    */

	Draw3D::drawAxis ( 3.0f );

};

// ===================== MAIN

TestAppSphereSampling * testApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestAppSphereSampling( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}
















