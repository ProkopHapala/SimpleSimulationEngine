
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


const double oct_edge_planes[] = { // there should be normals
0.2542288679069233,0.5684730304826415,-0.7824360014318398,
-0.6655798171017323,0.5684730304826415,0.4835720429064384,
0.8227018983896183,0.5684730304826415,-4.801936448283155e-17,
-0.6655798171017325,0.5684730304826415,-0.4835720429064383,
0.2542288679069233,0.5684730304826416,0.7824360014318397,
0.2542288679069232,0.5684730304826415,-0.7824360014318398,
-0.6655798171017323,0.5684730304826415,0.4835720429064385,
0.8227018983896183,0.5684730304826415,-1.440580934484946e-16,
-0.6655798171017325,0.5684730304826415,-0.4835720429064382,
0.2542288679069234,0.5684730304826416,0.7824360014318397
};



void sampleOctahedron( Vec3d p, uint8_t& iface, double& a, double& b ){
    //iface = (p.x>0) | ( ((uint8_t)(p.y>0))<<1) | ( ((uint8_t)(p.z>0) )<<2);
    iface = ( ((uint8_t)(p.x>0))<<2) | ( ((uint8_t)(p.y>0))<<1) | ( ((uint8_t)(p.z>0) ));
    double renorm = 1.0d/( fabs(p.x) + fabs(p.y) + fabs(p.z) );
    a = fabs(p.x);
    b = fabs(p.y);
}


void sampleIcosa( Vec3d p, uint8_t& iface, double& a, double& b ){
//  http://www.kjmaclean.com/Geometry/Icosahedron.html
//  https://en.wikipedia.org/wiki/Regular_icosahedron
//  rc = sqrt(10+2*sqrt(5))/4 = sqrt(phi^2 + 1)/2 = 0.95105651629 * a
//  plate height = phi/(2*sqrt(phi^2+1)).a = 0.42532540417.a  = 0.4472135955 rc
//  top height   = 1/sqrt(phi^2+1).a       = 0.52573111211.a  = 0.5527864045 rc
    double phi10 = (atan2( p.z, p.x ) + M_PI) * 0.15915494309 * 10 ;
    int iphi     = (int)phi10;
    double dphi  = phi10 - iphi;
    double hbound = 0.4472135955;

    Vec3d line = (Vec3d){0.58778525229,0.4472135955,0.80901699437};

    //double slope  =
    if( p.y>hbound ){
        iface = iphi/2;
    }else if(p.y<-hbound){
        iphi--; if(iphi<0){iface=5+4;}else{ iface = 5+iphi/2; }
    }else{
        if( iphi&1 ) {
            if( p.y > ( (dphi-0.5)*hbound*2) ){ iface = 10+iphi/2; }else{ iface = 15+ iphi/2+1;  if(iface>=20) iface=15; };
        }else{
            if( p.y > ( (0.5-dphi)*hbound*2) ){ iface = 10+iphi/2; }else{ iface = 15+ iphi/2; };
        }
    };
    //iface = iphi;
    //if      ( p.y>hbound ){ iface = 0; }
    //else if ( p.y<-hbound ){ iface = 1; }
    //else                  { iface = 2; };
}

void sampleIcosa2quads( Vec3d p, uint8_t& iface, double& a, double& b ){
    double phi10 = (atan2( p.z, p.x ) + M_PI) * ( 0.15915494309 * 10 );
    int iphi     = (int)phi10;
    double dphi  = phi10 - iphi;
    double hbound = 0.4472135955;
    //double slope  =
    if( p.y>hbound ){
        iface = iphi/2;
    }else if(p.y<-hbound){
        iface = 5+(iphi+1)/2;
        if(iface>=10) iface=5;
    }else{
        Vec3d& dir=((Vec3d*)oct_edge_planes)[iphi];
        if( dir.dot(p) > 0 ){ iface = iphi/2; }else{ iface = 5+ (iphi+1)/2;  if(iface>=10) iface=5; };

        /*
        double ycut = (2*dphi-1.0)*hbound;
        if( iphi&1 ) {
            if( p.y >  ycut ){ iface = iphi/2; }else{ iface = 5+ (iphi+1)/2;  if(iface>=10) iface=5; };
        }else{
            if( p.y > -ycut ){ iface = iphi/2; }else{ iface = 5+ iphi/2; };
        }
        */
    };
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

    for( int i=0; i<npix; i++ ){ pix[i]=0; }

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
        Draw3D::drawPolygons( Solids::Icosahedron_nfaces,        Solids::Icosahedron_ngons,        Solids::Icosahedron_faces,        Solids::Icosahedron_verts        );
        glRotatef( -31.7174744, 0,0,1 ); // arctan (phi-1) in degreers
        // Vec3i ivs = Solids::Octahedron_tris[iface];


        double sc = 1.80;
        glColor3f(0.0,0.0,0.0);
        glDisable( GL_LIGHTING );
        double dphi = (M_PI*2)/10;
        for( int i=0; i<10; i++ ){
            double sign = 2*(i&1)-1;
            double phi = dphi*i;
            double ca = cos(phi);
            double sa = sin(phi);
            Vec3d a = (Vec3d){ cos(phi),       0.4472135955*sign  ,sin(phi)};
            Vec3d b = (Vec3d){ cos(phi+dphi), -0.4472135955*sign  ,sin(phi+dphi)};

            Vec3d c; c.set_cross(a+b,a-b);
            Draw3D::drawLine( a*sc, b*sc);
            Draw3D::drawVecInPos( c, a+b );
            c.normalize();
            printf( "%.16g,%.16g,%.16g,\n", c.x, c.y, c.z );
        }



        glDisable( GL_LIGHTING );
        glPointSize(1);
         sc = 1.87;


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
                    Draw::color_of_hash( (iface+1545)*151 );

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
	Vec3d p = camMat.c*-3.0;
	Draw3D::drawPointCross( p, 0.1 );

	uint8_t iface; double a,b;
	sampleOctahedron( p, iface, a, b );

	//printf( "iface %i | (%f,%f,%f) \n", iface, p.x, p.y, p.z );

	Vec3i ivs = Solids::Octahedron_tris[iface];
	Vec3d c;
	sampleTri( p, ivs, Solids::Octahedron_verts, c );

    if ( SDL_GetMouseState(NULL, NULL) & SDL_BUTTON(SDL_BUTTON_LEFT) ) {
        //SDL_Log("Mouse Button 1 (left) is pressed.");
        uint32_t * pixF = pix + nsamp*nsamp*iface;
        //float step = 1.0/nsamp;
        int ia = (int)(c.a*nsamp);
        int ib = (int)(c.b*nsamp);
        int i = ia*nsamp + ib;
        pixF[i] = 0xffffffff;
    }


	double fsc = 1.1;
	Vec3d p_;
	p_.set_lincomb( c.a,c.b,c.c,  Solids::Octahedron_verts[ivs.a], Solids::Octahedron_verts[ivs.b], Solids::Octahedron_verts[ivs.c] );
	p_.mul(fsc);
	Draw3D::drawLine( Solids::Octahedron_verts[ivs.a]*fsc, p_ );
	Draw3D::drawLine( Solids::Octahedron_verts[ivs.b]*fsc, p_ );
	Draw3D::drawLine( Solids::Octahedron_verts[ivs.c]*fsc, p_ );





	glColor3f( 1.0f,0.0f,1.0f );
	fsc = 1.1;
    //Draw3D::drawTriangle( Solids::Octahedron_verts[ivs.a]*fsc, Solids::Octahedron_verts[ivs.b]*fsc, Solids::Octahedron_verts[ivs.c]*fsc );

    glDisable ( GL_LIGHTING );

    glColor3f( 1.0f,0.0f,0.0f );


    //exit(0);

    glColor3f( 0.0f,1.0f,0.0f );
    for( int i=0; i<Solids::Octahedron_nverts; i++ ){
        sprintf( str, "%i", i );
        Draw3D::drawText(str, Solids::Octahedron_verts[i], fontTex, 0.03, 0);
    }


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
















