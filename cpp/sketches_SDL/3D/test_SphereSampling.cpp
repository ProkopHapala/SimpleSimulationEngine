
/// @file @brief This demo procedurally generates a planet or asteroid using `SphereSampling.h` and `Noise.h`. It applies one or more layers of a noise function (like Simplex noise) to a sphere to create a heightmap, resulting in features like continents, mountains, and craters. The result is rendered using `DrawSphereMap.h`, and the user can rotate the generated celestial body.
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

#include <cstdio>
#include <cmath>
#include <cassert>
#include <algorithm>

#include "SDL_utils.h"
#include "GLUtils.h"
#include "argparse.h"

using namespace SphereSampling;



// ======================  TestApp
int nsamp = 32;
int npix  = 10*nsamp*nsamp;

uint32_t  *pix;
float     *heights;
Vec3f     *samplePs;

float fontScale = 0.03;

const int nCrater = 50;
Vec3d  craterPos[nCrater];
double craterSz[nCrater];

bool debugPlot = false;
Vec3d curPos = (Vec3d){5.0,0.0,0.0};


char str[2048];
//int  fontTex;

// --- CLI mode selection
enum ModeSphereSampling { MODE_FULL=0, MODE_TINYRAND=1 };
static ModeSphereSampling g_mode = MODE_FULL; // default keeps existing behavior

// Evaluate height at an arbitrary unit vector by sampling the icosahedral per-face grid.
// Grid layout per face is rectangular: n.a x n.b, stored row-major with index ia*n.b + ib.
// We use a diagonal split inside each cell consistent with diTri2cartes():
// if (u<=v): use triangle (v00, v01, v11); else: (v00, v10, v11).
static inline float evalIcosaGrid(const Vec3d& p_in, Vec2i n, const float* heights){
    Vec3d p = p_in; p.normalize();
    uint8_t iface; double fa, fb; sampleIcosa2quads(p, iface, fa, fb);
    uint8_t face10 = iface >> 1; // sampleIcosa2quads returns 20 triangles; storage is 10 rhombi

    const int na = n.a, nb = n.b; const int nab = na*nb; const int total = 10*nab;

    // Basic sanity checks
    if(!(na>=2 && nb>=2)){
        std::printf("[evalIcosaGrid] ERROR: na/nb too small: na=%d nb=%d\n", na, nb);
        std::printf(" p=(%.6f,%.6f,%.6f)\n", p.x, p.y, p.z);
        std::fflush(stdout);
        assert(false);
    }
    if(!(std::isfinite(fa) && std::isfinite(fb))){
        std::printf("[evalIcosaGrid] ERROR: non-finite fa/fb: fa=%g fb=%g iface=%u\n", fa, fb, (unsigned)iface);
        std::printf(" p=(%.6f,%.6f,%.6f)\n", p.x, p.y, p.z);
        std::fflush(stdout);
        assert(false);
    }
    if(!(face10 < 10)){
        std::printf("[evalIcosaGrid] ERROR: face10 out of range: iface=%u face10=%u (expected 0..9)\n", (unsigned)iface, (unsigned)face10);
        std::printf(" fa=%g fb=%g p=(%.6f,%.6f,%.6f)\n", fa, fb, p.x, p.y, p.z);
        std::fflush(stdout);
        assert(false);
    }

    // map [0,1] to [0,na-1] and [0,nb-1]
    double sa = fa * (na-1);
    double sb = fb * (nb-1);
    int ia = (int)floor(sa);
    int ib = (int)floor(sb);
    if(ia<0) ia=0; if(ia>na-2) ia=na-2;
    if(ib<0) ib=0; if(ib>nb-2) ib=nb-2;
    double u = sa - ia;  // in [0,1)
    double v = sb - ib;  // in [0,1)

    // Bounds for the face block
    const int baseOff = face10*nab;
    const int i00 = ia*nb + ib;
    const int i10 = (ia+1)*nb + ib;
    const int i01 = ia*nb + (ib+1);
    const int i11 = (ia+1)*nb + (ib+1);
    if(!(baseOff>=0 && baseOff < total)){
        std::printf("[evalIcosaGrid] ERROR: baseOff out of total: baseOff=%d total=%d iface=%u face10=%u\n", baseOff, total, (unsigned)iface, (unsigned)face10);
        std::fflush(stdout);
        assert(false);
    }
    {
        int imin = std::min(std::min(i00,i10), std::min(i01,i11));
        int imax = std::max(std::max(i00,i10), std::max(i01,i11));
        if(!(imin >= 0 && imax < nab)){
            std::printf("[evalIcosaGrid] ERROR: local indices out of face range: i00=%d i10=%d i01=%d i11=%d nab=%d (ia=%d ib=%d na=%d nb=%d)\n", i00,i10,i01,i11,nab,ia,ib,na,nb);
            std::fflush(stdout);
            assert(false);
        }
    }
    if(!(baseOff + i11 < total)){
        std::printf("[evalIcosaGrid] ERROR: global index overflow: baseOff=%d i11=%d total=%d (iface=%u face10=%u)\n", baseOff, i11, total, (unsigned)iface, (unsigned)face10);
        std::fflush(stdout);
        assert(false);
    }

    const float* base = heights + baseOff;
    float h00 = base[i00];
    float h10 = base[i10];
    float h01 = base[i01];
    float h11 = base[i11];

    // Seam-aware fix for the top/bottom edge (fb direction):
    // If we are sampling the last row of cells (ib == nb-2), replace the samples that
    // lie on the outer row (i01,i11) with values from the corresponding neighbor face.
    // This mirrors the cap stitching used in drawDiTri_seam with views
    //   current:  {-1, n.a*n.b-1}
    //   neighbor: {-n.a, n.a*(n.a-1)}
    // For square grids (na==nb), this maps (ia, nb-1) -> neighbor linear index ia*na.
    if( ib == (nb-2) ){
        int band   = face10 / 5;         // 0: upper belt, 1: lower belt
        int i_ring = face10 % 5;         // position within the belt
        int i2     = (i_ring + 1); if(i2>=5) i2=0;
        int neigh_face10;
        if(band==0){
            // upper belt connects to lower belt face (i2+5) across +b
            neigh_face10 = i2 + 5;
        }else{
            // lower belt connects back to upper belt face i across +b
            neigh_face10 = i_ring; // i in upper belt
        }
        int neigh_off = neigh_face10 * nab;
        if(!(neigh_face10>=0 && neigh_face10<10)){
            std::printf("[evalIcosaGrid] ERROR: neigh_face10 out of range: %d (face10=%u)\n", neigh_face10, (unsigned)face10);
            std::fflush(stdout);
            assert(false);
        }
        if(!(neigh_off>=0 && neigh_off+nab <= total)){
            std::printf("[evalIcosaGrid] ERROR: neighbor block out of bounds: neigh_off=%d nab=%d total=%d\n", neigh_off, nab, total);
            std::fflush(stdout);
            assert(false);
        }
        const float* baseN = heights + neigh_off;

        // Using the same orientation as seam stitching: neighbor linear index = ia*na and (ia+1)*na
        // Note: this assumes na==nb as used elsewhere in this demo.
        int in01 = ia*na;
        int in11 = (ia+1)*na;
        if(!(in01>=0 && in11<nab)){
            std::printf("[evalIcosaGrid] ERROR: neighbor indices out of range: in01=%d in11=%d nab=%d (ia=%d na=%d)\n", in01, in11, nab, ia, na);
            std::fflush(stdout);
            assert(false);
        }
        h01 = baseN[in01];
        h11 = baseN[in11];
    }
    if(u <= v){
        // weights: w00=1-v, w01=v-u, w11=u
        return (float)((1.0-v)*h00 + (v-u)*h01 + u*h11);
    }else{
        // weights: w00=1-u, w10=u-v, w11=v
        return (float)((1.0-u)*h00 + (u-v)*h10 + v*h11);
    }
}

class TestAppSphereSampling : public AppSDL2OGL_3D {
	public:

    int ogl_base=-1;
	int ogl_points=-1;
	int ogl_asteroide=-1;


	// ---- function declarations

	virtual void draw   ();

	TestAppSphereSampling( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppSphereSampling::TestAppSphereSampling( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    //glEnable(GL_AUTO_NORMALS);

    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );


    //const int nsamp = 8;

    if(g_mode == MODE_FULL){
        // Original asteroid/crater example with nsamp=32

        nsamp = 32;
        npix  = 10*nsamp*nsamp;
        pix      = new uint32_t[npix];
        heights  = new float[npix];
        samplePs = new Vec3f[npix];

        srand(34646);
        for(int i=0; i<nCrater; i++){
            craterPos[i].fromRandomSphereSample();
            craterSz[i] = randf()+0.4;
        }

        ogl_base=glGenLists(1);
        glNewList( ogl_base, GL_COMPILE );
            for(int i=0; i<nCrater; i++) Draw3D::drawPointCross( craterPos[i], craterSz[i] );
            glPushMatrix();
            glDisable ( GL_LIGHTING );
            glColor3f( 1.0f, 0.0f, 1.0f );
            glPopMatrix();
            glPushMatrix();
            glEnable    ( GL_LIGHTING );
            glShadeModel( GL_FLAT     );
            glColor3f( 0.8f, 0.8f, 0.8f );
            glScalef(0.5,0.5,0.5);
            glRotatef( 31.7174744, 0,0,1 );
            glDisable( GL_LIGHTING ); glColor3f(1.0,1.0,1.0); Draw3D::drawLines( Solids::Icosahedron_nedges, (int*)Solids::Icosahedron_edges, Solids::Icosahedron_verts );
            glRotatef( -31.7174744, 0,0,1 );
        glEndList();
    

        ogl_asteroide=glGenLists(1);
        glNewList( ogl_asteroide, GL_COMPILE );
            double sc = 1.90211303259;
            glColor3f(0.0,0.0,0.0);
            glDisable( GL_LIGHTING );
            glScalef(sc,sc,sc);
            glDisable(GL_LIGHTING);
            glShadeModel( GL_SMOOTH     );
            glColor3f(1.0,1.0,1.0);
            drawIcosaMap( (Vec2i){nsamp,nsamp}, heights, 0.0 );
            glPopMatrix();
        glEndList();

        Vec2d dab = (Vec2d){ 1.0/nsamp, 1.0/nsamp };
        for(int iface=0; iface<10; iface++ ){
            for(int ia=0; ia<nsamp; ia++ ){
                for(int ib=0; ib<nsamp; ib++ ){
                    Vec3d p;
                    icosa2cartes( (Vec2i){nsamp,nsamp}, iface, ia*dab.a, ib*dab.b, p );
                    int i = iface*nsamp*nsamp + ia*nsamp + ib;
                    p.normalize();
                    double h = 0.0;
                    h += Noise::getCraterHeight( p, nCrater, 1.0, craterPos, craterSz )*0.3;
                    heights[i] = h;
                }
            }
        }
    }else{
        nsamp = 3;
        npix  = 10*nsamp*nsamp;
        pix      = new uint32_t[npix];
        heights  = new float   [npix];
        samplePs = new Vec3f   [npix];
        // MODE_TINYRAND: nsample=3 random heights to verify interpolation

        //srand(1337); for(int i=0; i<npix; i++ ){ heights[i] = (randf()*2.0f - 1.0f); }

        for(int i=0; i<npix; i++ ){ heights[i] = -1.0f; }; heights[0] = 1.0f;heights[1] = 2.0f;heights[2] = 3.0f;
        


        ogl_points=glGenLists(1);
        glNewList( ogl_points, GL_COMPILE );
            glDisable(GL_LIGHTING);
            glPointSize(3.0f);
            const int nPoints = 50000;
            const double sc_pts = 1.0;
            glBegin(GL_POINTS);
            for(int i=0; i<nPoints; i++){
                Vec3d p; p.fromRandomSphereSample();
                float h = evalIcosaGrid(p, (Vec2i){nsamp,nsamp}, heights);
                heightColor(h);
                //glVertex3f( (float)(sc_pts*p.x), (float)(sc_pts*p.y), (float)(sc_pts*p.z) );
                glVertex3f( (float)(sc_pts*-p.x), (float)(sc_pts*-p.y), (float)(sc_pts*-p.z) );
            }
            glEnd();

            glEnable    ( GL_LIGHTING );
            glShadeModel( GL_FLAT     );
            glColor3f( 0.8f, 0.8f, 0.8f );
            glScalef(0.5,0.5,0.5);
            glRotatef( 31.7174744, 0,0,1 );
            glDisable( GL_LIGHTING ); glColor3f(1.0,1.0,1.0); Draw3D::drawLines( Solids::Icosahedron_nedges, (int*)Solids::Icosahedron_edges, Solids::Icosahedron_verts );
            glRotatef( -31.7174744, 0,0,1 );
        glEndList();
    }



    zoom = 3.0;

}

void TestAppSphereSampling::draw   (){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glEnable(GL_DEPTH_TEST);
	glCallList( ogl_base );
	glCallList( ogl_points );
    glCallList( ogl_asteroide );


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

    double fsc = 1.1;

    // debugPlot = true;
    // sampleIcosa2quads( curPos, iface, a, b );
    // sprintf( str, "(%.3f,%.3f)", a,b );   Draw3D::drawText(str,curPos, fontTex, fontScale, 0);
    // debugPlot = false;
	// Vec3d p_;
	// p_.set_lincomb( c.a,c.b,c.c,  Solids::Octahedron_verts[ivs.a], Solids::Octahedron_verts[ivs.b], Solids::Octahedron_verts[ivs.c] );
	// p_.mul(fsc);
	// Draw3D::drawLine( Solids::Octahedron_verts[ivs.a]*fsc, p_ );
	// Draw3D::drawLine( Solids::Octahedron_verts[ivs.b]*fsc, p_ );
	// Draw3D::drawLine( Solids::Octahedron_verts[ivs.c]*fsc, p_ );


    // ---- view icosahedron polar vertices
    // for( int i=0; i<5; i++ ){
    //     Vec3d p;
    //     glColor3f(1.0,0.0,0.0);
    //     p=  ((Vec3d*)icosa_polar_verts)[i]*2.0;
    //     //printf("%i p: %f %f %f\n", i, p.x, p.y, p.z );
    //     sprintf( str, "%i", i );   Draw3D::drawText(str, p, fontTex, fontScale, 0);
    //     p =  ((Vec3d*)icosa_polar_verts)[i+5]*2.0;
    //     sprintf( str, "%i", i+5 ); Draw3D::drawText(str, p, fontTex, fontScale, 0);
    // }

    for( int i=0; i<12; i++ ){
        Vec3d p;
        glColor3f(1.0,0.0,0.0);
        p =  ((Vec3d*)icosa_polar_verts)[i]*2.0;
        sprintf( str, "%i", i ); 
        Draw3D::drawText(str, p, fontTex, fontScale, 0);
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


TestAppSphereSampling * testApp;

int main(int argc, char *argv[]){
    // no buffering
    setbuf(stdout, NULL);
    setbuf(stderr, NULL);
    // Parse CLI flags to choose mode
    LambdaDict funcs;
    funcs["-full"]     = {0, [&](const char**){ g_mode = MODE_FULL;     printf("test_SphereSampling.cpp: arg: -full\n"); }};
    funcs["-testRand"] = {0, [&](const char**){ g_mode = MODE_TINYRAND; printf("test_SphereSampling.cpp: arg: -testRand \n"); }};
    process_args(argc, argv, funcs, false);

    SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestAppSphereSampling( junk , 1600, 1200 );
	testApp->loop( 1000000 );
	return 0;
}
