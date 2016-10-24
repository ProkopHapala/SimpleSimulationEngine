
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
#include "Besier.h"
#include "Solids.h"

#include "Draw.h"
#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"

// ======================  TestApp

inline void subEdge( Vec3d& result, const Vec3d& up, const Vec3d& down, const Vec3d& left, const Vec3d& right ){
    constexpr double c_para = 0.125d;
    constexpr double c_perp = 0.375d;
    result.set_lincomb( c_para, up, c_para, down );
    result.add_mul(left, c_perp);
    result.add_mul(right, c_perp);
};


// Loop subdivision : not finished
// TODO : point order should be modified to be consistent with the rim
// later this should be optimized especially for case of uniform grid
// see here  : https://graphics.stanford.edu/~mdfisher/subdivision.html
void subdivTrinagle( int level,
    const Vec3d& p0, const Vec3d& p1,  const Vec3d& p2,
    const Vec3d& p3, const Vec3d& p4,  const Vec3d& p5,
    const Vec3d& p6, const Vec3d& p7,  const Vec3d& p8,
    const Vec3d& p9, const Vec3d& p10, const Vec3d& p11
    ){
    level--;
    Vec3d h0,h1,h2;
    subEdge( h0,  p0,p1,  p2,p4  );
    subEdge( h1,  p1,p2,  p0,p7  );
    subEdge( h2,  p0,p2,  p1,p10 );
    if( level>0 ){
        Vec3d h3,h4,h5,h6, h7,h8,h9,h10, h11,h12,h13,h14;
        subEdge( h3 ,  p0,p10,  p2,p11);
        subEdge( h4 ,  p0,p11,  p3,p10);
        subEdge( h5 ,  p0,p3,   p4,p11);
        subEdge( h6 ,  p0,p4,   p1,p3 );
        subEdge( h7 ,  p1,p4,   p0,p5 );
        subEdge( h8 ,  p1,p5,   p4,p6 );
        subEdge( h9 ,  p1,p6,   p5,p7 );
        subEdge( h10,  p1,p7,   p2,p6 );
        subEdge( h11,  p2,p7,   p1,p8 );
        subEdge( h12,  p2,p8,   p7,p9 );
        subEdge( h13,  p2,p9,   p8,p10);
        subEdge( h14,  p2,p10,  p0,p9 );
        subdivTrinagle( level, p0,h1,h2,  h5,h6,h7 ,p1 ,h1 ,p2 ,h14,h3 ,h4 );
        subdivTrinagle( level, h0,p1,h1,  h6,h7,h8 ,h9 ,h10,h11,p2 ,h2 ,p0 );
        subdivTrinagle( level, h2,h1,p2,  p0,h0,p1 ,h10,h11,h12,h13,h14,h3 );
        subdivTrinagle( level, h0,h1,h2,  h7,p1,h10,h11,p2 ,h14,h3 ,p0 ,h6 );
    }else{
        Draw3D::drawTriangle( p0, h0, h2 );
        Draw3D::drawTriangle( p1, h1, h0 );
        Draw3D::drawTriangle( p2, h2, h1 );
        Draw3D::drawTriangle( h0, h1, h2 );
    }
}


void drawBesierHull( const BesierTriangle& btri, uint32_t inner, uint32_t outer, uint32_t normals ){
    glDisable(GL_LIGHTING);
    if( normals != 0 ){
        Draw::setRGBA( normals );

        //Draw3D::drawVecInPos( btri.ns[0], btri.ps[0] );
        //Draw3D::drawVecInPos( btri.ns[1], btri.ps[3] );
        //Draw3D::drawVecInPos( btri.ns[2], btri.ps[6] );

        Draw3D::drawVecInPos( btri.ns[3], (btri.ps[1]+btri.ps[4])*0.5 );   // ab
        Draw3D::drawVecInPos( btri.ns[4], (btri.ps[2]+btri.ps[7])*0.5 );  // ac
        Draw3D::drawVecInPos( btri.ns[5], (btri.ps[5]+btri.ps[8])*0.5 );    // bc


    }
    if( outer != 0 ){
        glBegin(GL_LINE_LOOP);
        Draw::setRGBA( outer );
        for( int i=0; i<9; i++ ){
            Vec3f p; convert(btri.ps[ BesierTriangle_borderHull[i] ],p);
            //p.add(randf(-nR,nR),randf(-nR,nR),randf(-nR,nR));
            glVertex3f(p.x,p.y,p.z);
            //printf(" %i %i (%3.3f,%3.3f,%3.3f)\n", i, btriEdge[i], p.x,p.y,p.z  );
        }
        glEnd();
    }
    if( inner != 0 ){
        glBegin(GL_LINE_LOOP);
        Draw::setRGBA( inner  );
        for( int i=0; i<9; i++ ){
            Vec3f p;
            convert(btri.ps[ BesierTriangle_innerHull[i] ],p);
            //p.add(randf(-nR,nR),randf(-nR,nR),randf(-nR,nR));
            glVertex3f(p.x,p.y,p.z);
            //printf(" %i %i (%3.3f,%3.3f,%3.3f)\n", i, btriIn[i], p.x,p.y,p.z  );
        }
        glEnd();
    }
}

void drawBesierTriangle( int nsub, const BesierTriangle& btri, uint32_t faces, uint32_t edges, uint32_t normals ){

    double da=1.0d/nsub;
    double db=1.0d/nsub;

    if( edges != 0 ){
        Draw::setRGBA( edges );
        glDisable(GL_LIGHTING);
        for(int ia=0; ia<nsub; ia++){
            Vec3d p;
            glBegin(GL_LINE_STRIP); for(int ib=0; ib<=(nsub-ia); ib++){ p = btri.getPoint(  ia*da, ib*db ); glVertex3d(p.x,p.y,p.z);    } glEnd();
        }
        for(int ib=0; ib<nsub; ib++){
            Vec3d p;
            glBegin(GL_LINE_STRIP); for(int ia=0; ia<=(nsub-ib); ia++){ p = btri.getPoint( ia*da, ib*db); glVertex3d(p.x,p.y,p.z); }     glEnd();
        }
        for(int ic=0; ic<nsub; ic++){
            Vec3d p;
            glBegin(GL_LINE_STRIP); for(int ia=0; ia<=(nsub-ic); ia++){ p = btri.getPoint( ia*da, (nsub-ic-ia)*db ); glVertex3d(p.x,p.y,p.z); }     glEnd();
        }
    }

    if( normals != 0 ){
        Draw::setRGBA( normals );
        glBegin(GL_LINES);
        glDisable(GL_LIGHTING);
        for(int ia=0; ia<nsub; ia++){
            Vec3d p,n;
            for(int ib=0; ib<=(nsub-ia); ib++){
                p = btri.getPoint(  ia*da, ib*db );
                glVertex3d(p.x,p.y,p.z);
                n = btri.getNormal(  ia*da, ib*db );
                p.add_mul(n,0.2);
                glVertex3d(p.x,p.y,p.z);
            }
        }
        glEnd();
    }

    if( faces != 0 ){
        Draw::setRGBA( faces );
        glEnable(GL_LIGHTING);
        for(int ia=0; ia<nsub; ia++){
            glBegin(GL_TRIANGLE_STRIP);
            //glNormal3f(nA.x,nA.y,nA.z);
            Vec3d p,n;
            for(int ib=0; ib<(nsub-ia); ib++){
                double a,b;
                a=ia*da; b=ib*db;
                n = btri.getNormal(  a,b ); glNormal3d(n.x,n.y,n.z);
                p = btri.getPoint (  a,b ); glVertex3d(p.x,p.y,p.z);
                a+=da;
                n = btri.getNormal( a,b );  glNormal3d(n.x,n.y,n.z);
                p = btri.getPoint ( a,b );  glVertex3d(p.x,p.y,p.z);
                //printf( " %i %i (%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) \n",  ia, ib, ia*da, ib*db, p1.x,p1.y,p1.z, p2.x,p2.y,p2.z  );
            }
            p = btri.getNormal( ia*da, (nsub-ia)*db );  glNormal3d(n.x,n.y,n.z);
            p = btri.getPoint ( ia*da, (nsub-ia)*db );  glVertex3d(p.x,p.y,p.z);
            glEnd();
        }
    }

}


class TestAppPatches : public AppSDL2OGL_3D {
	public:

	int point_cloud;
	int shape;

	Vec3d A,B,C,nA,nB,nC,center;
	BesierTriangle btri;

	bool refresh=true;

	// ---- function declarations

	virtual void draw   ();
	virtual void eventHandling ( const SDL_Event& event  );

	TestAppPatches( int& id, int WIDTH_, int HEIGHT_ );


};


TestAppPatches::TestAppPatches( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {


    // initialization
    center.set(0.0);
    srand(123);
    double pR=1.0;
    double nR=0.5;

    A.set(randf(-pR,pR),randf(-pR,pR),randf(-pR,pR));
    B.set(randf(-pR,pR),randf(-pR,pR),randf(-pR,pR));
    C.set(randf(-pR,pR),randf(-pR,pR),randf(-pR,pR));
    Vec3d cog; cog.set( (A+B+C)*0.33333333 ); A.sub(cog); B.sub(cog); C.sub(cog);

    nA.set_cross(B-A,C-A); nA.normalize(); nB.set(nA); nC.set(nB);
    nA.add(randf(-nR,nR),randf(-nR,nR),randf(-nR,nR)); nA.normalize();
    nB.add(randf(-nR,nR),randf(-nR,nR),randf(-nR,nR)); nB.normalize();
    nC.add(randf(-nR,nR),randf(-nR,nR),randf(-nR,nR)); nC.normalize();


    btri.fromPNs( A, B, C,  nA, nB, nC );

    for(int i=0; i<10; i++){
        printf(" p %i (%3.3f,%3.3f,%3.3f)\n", i, btri.ps[i].x, btri.ps[i].y, btri.ps[i].z );
    };
    for(int i=0; i<6; i++){
        printf(" n %i (%3.3f,%3.3f,%3.3f)\n", i, btri.ns[i].x, btri.ns[i].y, btri.ns[i].z );
    };

     /*
     p 0 (-0.594,0.577,-0.880)
     p 1 (-0.671,0.338,-0.650)
     p 2 (-0.735,0.276,-0.627)
     p 3 (-0.731,-0.277,-0.303)
     p 4 (-0.728,0.037,-0.380)
     p 5 (-0.784,-0.350,-0.304)
     p 6 (-0.911,-0.481,-0.248)
     p 7 (-0.903,-0.099,-0.375)
     p 8 (-0.881,-0.404,-0.241)
     p 9 (-0.803,-0.020,-0.406)
     n 0 (0.464,-0.688,-0.558)
     n 1 (0.337,-0.230,-0.913)
     n 2 (0.741,-0.227,-0.632)
     n 3 (0.403,-0.591,-0.699)
     n 4 (0.601,-0.561,-0.569)
     n 5 (0.231,-0.645,-0.728)

    p 0 (0.151,0.637,-0.403)
     p 1 (0.075,0.399,-0.173)
     p 2 (0.010,0.337,-0.150)
     p 3 (0.015,-0.216,0.174)
     p 4 (0.018,0.097,0.096)
     p 5 (-0.038,-0.289,0.173)
     p 6 (-0.166,-0.421,0.229)
     p 7 (-0.158,-0.038,0.101)
     p 8 (-0.136,-0.343,0.236)
     p 9 (-0.057,0.040,0.071)
     n 0 (0.464,-0.688,-0.558)
     n 1 (0.337,-0.230,-0.913)
     n 2 (0.741,-0.227,-0.632)
     n 3 (0.403,-0.591,-0.699)
     n 4 (0.601,-0.561,-0.569)
     n 5 (0.231,-0.645,-0.728)


     */

    zoom=2.0;
}

void TestAppPatches::draw   (){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    if( refresh ){
        printf("Refreshing!!!\n");

        glDeleteLists(shape,1);
        shape=glGenLists(1);
        glNewList( shape, GL_COMPILE );
        /*
            btri.fromPNs( A, B, C,  nA, nB, nC );
            drawBesierHull    (     btri,  0xFF008000, 0xFF00FF00, 0xFF0000FF );
            drawBesierTriangle( 10, btri,  0xFFFFFFFF, 0xFF000000, 0xFF000080 );
        */

            for ( int i=0; i<Solids::Tetrahedron_ntris; i++ ){
                Vec3d A,B,C, nA,nB,nC;


                int i3=3*i;
                A.set(Solids::Tetrahedron_verts[ Solids::Tetrahedron_tris[i3  ] ]);  nA.set(A); nA.normalize(); A.add(center);
                B.set(Solids::Tetrahedron_verts[ Solids::Tetrahedron_tris[i3+1] ]);  nB.set(B); nB.normalize(); B.add(center);
                C.set(Solids::Tetrahedron_verts[ Solids::Tetrahedron_tris[i3+2] ]);  nC.set(C); nC.normalize(); C.add(center);
                //printf("",A.x);
                btri.fromPNs( A, B, C,  nA, nB, nC );
                drawBesierHull    (     btri,  0xFF008000, 0xFF00FF00, 0xFF0000FF );
                drawBesierTriangle( 20, btri,  0xFFFFFFFF, 0xFF000000, 0xFF000080 );
            };


        /*
        TO DO : problem that control points are the same
        glColor3f(1.0f,1.0f,1.0f);
        subdivTrinagle( 1,
            Solids::Tetrahedron_verts[0],Solids::Tetrahedron_verts[1],Solids::Tetrahedron_verts[2],
            Solids::Tetrahedron_verts[3],Solids::Tetrahedron_verts[3],Solids::Tetrahedron_verts[3],
            Solids::Tetrahedron_verts[3],Solids::Tetrahedron_verts[3],Solids::Tetrahedron_verts[3],
            Solids::Tetrahedron_verts[3],Solids::Tetrahedron_verts[3],Solids::Tetrahedron_verts[3]
        );
        */
        glEndList();
        refresh=false;
	}

	glCallList( shape );

	glDisable ( GL_LIGHTING );

	//Draw3D::drawAxis ( 3.0f );

};


void TestAppPatches::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_a: center.add(camMat.a*0.1); refresh=true; break;
                case SDLK_d: center.sub(camMat.a*0.1); refresh=true; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            printf("keyPressed\n");
            break;

    };
    AppSDL2OGL::eventHandling( event );
}

// ===================== MAIN

TestAppPatches * testApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestAppPatches( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}
















