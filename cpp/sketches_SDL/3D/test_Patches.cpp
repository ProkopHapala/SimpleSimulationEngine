
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

#include "grids2D.h"

#include "TerrainSimplex.h"
#include "Draw.h"
#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"

#include "spline_triC1.h"

// ======================  TestApp

inline void subEdge( Vec3d& result, const Vec3d& up, const Vec3d& down, const Vec3d& left, const Vec3d& right ){
    constexpr double c_para = 0.375d;
    constexpr double c_perp = 0.125d;
    result.set_lincomb( c_para, up, c_para, down );
    result.add_mul(left,  c_perp);
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
    constexpr double c_on  = 0.625d;
    constexpr double c_off = 0.0625d;
    Vec3d v0,v1,v2;
    v0.set_lincomb( c_on, p0, c_off, p1+p2 +p3+p4+p10+p11 );
    v1.set_lincomb( c_on, p1, c_off, p0+p2 +p4+p5+p6+p7   );
    v2.set_lincomb( c_on, p2, c_off, p0+p1 +p7+p8+p9+p10  );
    //printf("h0 (%3.3f,%3.3f,%3.3f) h1 (%3.3f,%3.3f,%3.3f) h2 (%3.3f,%3.3f,%3.3f) \n", h0.x, h0.y, h0.z, h1.x, h1.y, h1.z, h2.x, h2.y, h2.z );
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

        subdivTrinagle( level, v0,h0,h2,  h5,h6,h7 ,v1 ,h1 ,v2 ,h14,h3 ,h4 );
        subdivTrinagle( level, h0,v1,h1,  h6,h7,h8 ,h9 ,h10,h11,v2 ,h2 ,v0 );
        subdivTrinagle( level, h2,h1,v2,  v0,h0,v1 ,h10,h11,h12,h13,h14,h3 );
        subdivTrinagle( level, h0,h1,h2,  h7,v1,h10,h11,v2 ,h14,h3 ,v0 ,h6 );
        /*
        subdivTrinagle( level, p0,h0,h2,  h5,h6,h7 ,p1 ,h1 ,p2 ,h14,h3 ,h4 );
        subdivTrinagle( level, h0,p1,h1,  h6,h7,h8 ,h9 ,h10,h11,p2 ,h2 ,p0 );
        subdivTrinagle( level, h2,h1,p2,  p0,h0,p1 ,h10,h11,h12,h13,h14,h3 );
        subdivTrinagle( level, h0,h1,h2,  h7,p1,h10,h11,p2 ,h14,h3 ,p0 ,h6 );
        */
        /*
        glDisable(GL_LIGHTING);
        Draw3D::drawPointCross(h0,0.05);
        Draw3D::drawPointCross(h1,0.05);
        Draw3D::drawPointCross(h2,0.05);

        Draw3D::drawPointCross(h3,0.05);
        Draw3D::drawPointCross(h4,0.05);
        Draw3D::drawPointCross(h5,0.05);
        Draw3D::drawPointCross(h6,0.05);

        Draw3D::drawPointCross(h7,0.05);
        Draw3D::drawPointCross(h8,0.05);
        Draw3D::drawPointCross(h9,0.05);
        Draw3D::drawPointCross(h10,0.05);

        Draw3D::drawPointCross(h11,0.05);
        Draw3D::drawPointCross(h12,0.05);
        Draw3D::drawPointCross(h13,0.05);
        Draw3D::drawPointCross(h14,0.05);
        */
    }else{
        glEnable(GL_LIGHTING);
        Draw3D::drawTriangle( v0, h0, h2 );
        Draw3D::drawTriangle( v1, h1, h0 );
        Draw3D::drawTriangle( v2, h2, h1 );
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

void interpolateTriC1( int na, int nb, int nsub, double* CPs, double* out ){
    double d = 1.0/nsub;
    int na_ = na*nsub;
    for(int ib=0; ib<nb; ib++){
        for(int ia=0; ia<na; ia++){
            int i0=ib*na+ia;
            double A=0,B=0,C=0,D=0,E=0,F=0,E_=0,F_=0;
            if((ia<(na-1))&&(ib<(nb-1)) ) A = CPs[i0+na+1];
            if((ia<(na-1))              ) B = CPs[i0   +1];
            if(             (ib<(nb-1)) ) C = CPs[i0+na  ];
                                          D = CPs[i0     ];
            if(             (ib<(nb-2)) ) E = CPs[i0+na*2];
            if((ia<(na-2))              ) F = CPs[i0   +2];
            if((ia>1     )&&(ib<(nb-1)) ) E_= CPs[i0+na-1];
            if((ia<(na-1))&&(ib>1)      ) F_= CPs[i0-na+1];
            for(int isb=0;isb<nsub;isb++){
                for(int isa=0;isa<nsub;isa++){
                    double v = d*isb;
                    double u = d*isa;
                    double w = 1-u-v;
                    if (w>0){
                        out[ (ib*nsub+isb)*na_ + ia*nsub+isa ] = Spline_triC1::val<double>( u,v,w,  D, B, C, A, E_, F_ );
                        //out[ (ib*nsub+isb)*na_ + ia*nsub+isa ] = 1;
                    }else{
                        //u=1-u;
                        //v=1-v;
                        //w=1-u-v;
                        //out[ (ib*nsub+isb)*na_ + ia*nsub+isa ] = Spline_triC1::val<double>( u,v,w,  A, C, B, D, F, E  );
                        out[ (ib*nsub+isb)*na_ + ia*nsub+isa ] = Spline_triC1::val<double>( 1-u, 1-v, u+v-1,  A, C, B, D, F, E  );
                        //out[ (ib*nsub+isb)*na_ + ia*nsub+isa ] = 0;
                    }
                }
            }
            //exit(0);
        }
    }
}

/////////////////////////////////////////////////////////
//          TestAppPatches
/////////////////////////////////////////////////////////


class TestAppPatches : public AppSDL2OGL_3D {
	public:

	int point_cloud;
	int shape;

	Vec3d A,B,C,nA,nB,nC,center;
	BesierTriangle btri;

	int nsub = 10;
	double* val_interp = NULL;

	bool refresh=true;

	double test_tris[3*12] = {
        -0.5, 0.0, 0.0,
         0.5, 0.0, 0.0,
         0.0, 1.0, 0.0,
        -1.0,-1.0, 0.0,
         0.0,-1.0, 0.0,
         1.0,-1.0, 0.0,
         1.5, 0.0, 0.0,
         1.0, 1.0, 0.0,
         0.5, 2.0, 0.0,
        -0.5, 2.0, 0.0,
        -1.0, 1.0, 0.0,
        -1.5, 0.0, 0.0,
	};

	int      na,nb;
	double * vals;

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

    na=8; nb=8;
    vals = new double[ na*nb];
    for(int ia=0;ia<na; ia++){
        for(int ib=0;ib<nb; ib++){
            vals[ia*nb+ib] = randf(0.0,1.0)*randf(0.0,1.0);
            //vals[ia*nb+ib] = 0;
        }
    }
    //vals[2*nb+3] = 1;

    val_interp = new double[nsub*nsub*na*nb];
    interpolateTriC1( na, nb, nsub, vals, val_interp);

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
        /*

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
        */


        //TO DO : problem that control points are the same
        /*
        glColor3f(1.0f,1.0f,1.0f);
        subdivTrinagle( 3,
            Solids::Tetrahedron_verts[0],Solids::Tetrahedron_verts[1],Solids::Tetrahedron_verts[2],
            Solids::Tetrahedron_verts[3],Solids::Tetrahedron_verts[3],Solids::Tetrahedron_verts[3],
            Solids::Tetrahedron_verts[3],Solids::Tetrahedron_verts[3],Solids::Tetrahedron_verts[3],
            Solids::Tetrahedron_verts[3],Solids::Tetrahedron_verts[3],Solids::Tetrahedron_verts[3]
        );
        */




/*
        //glEnable(GL_LIGHTING);
        glDisable(GL_LIGHTING);
        glColor3f(1.0f,1.0f,1.0f);
        Vec3d * tris = (Vec3d*)&test_tris[0];
        double dmax = 0.5;
        for(int i=0; i<12; i++){
            //tris[i].z += randf(-dmax,dmax);
            //tris[i].y += randf(-dmax,dmax);
            //tris[i].x += randf(-dmax,dmax);
            if(i>2) tris[i].z +=-1;
            Draw3D::drawPointCross(tris[i],0.1);
        }
        subdivTrinagle( 3,
            tris[0],tris[1],tris[2],
            tris[3],tris[4],tris[5],
            tris[6],tris[7],tris[8],
            tris[9],tris[10],tris[11]
        );
        */


        Vec2d da,db;
        da.set( 0.5, 0.86602540378 );
        db.set( 1.0,           0.0 );

        /*
        //Draw3D::drawSimplexGrid( na, nb, da, db, vals, vals, 0, NULL );
        glColor3f(0.1f,0.1f,0.8f); Draw3D::drawSimplexGridLines( na, nb, da, db,  vals );

        int na_=na<<1,nb_=nb<<1;
        double * vals_ = new double[na_*nb_];
        subdivideLoopGrid( nb, na, vals, vals_ );
        //Draw3D::drawSimplexGrid( na_, nb_, da*0.5, db*0.5, vals_, vals_, 0, NULL );
        glColor3f(0.1f,0.8f,0.1f); Draw3D::drawSimplexGridLines( na_, nb_, da*0.5, db*0.5,  vals_ );

        int na__=na_<<1,nb__=nb_<<1;
        double * vals__= new double[na__*nb__];
        subdivideLoopGrid( nb_, na_, vals_, vals__ );
        //Draw3D::drawSimplexGrid( na__, nb__, da*0.25, db*0.25, vals__, vals__, 0, NULL );
        glColor3f(0.8f,0.1f,0.1f); Draw3D::drawSimplexGridLines( na__, nb__, da*0.25, db*0.25,  vals__ );
        */

        //glColor3f(0.8f,0.1f,0.1f); Draw3D::drawSimplexGridLinesToned( na, nb, da, db,  vals );
        glColor3f(0.8f,0.1f,0.1f); Draw3D::drawSimplexGridLinesToned( na*nsub, nb*nsub, da*0.1, db*0.1,  val_interp );


        //delete vals_,vals__;

        //Draw3D::drawSimplexGrid( na, nb, da, db, vals, vals );
        //Draw3D::drawSimplexGrid( na, nb, da, db );
        //Draw3D::drawSimplexGrid( );


        glEndList();
        refresh=false;
	}

	glCallList( shape );

	glDisable ( GL_LIGHTING );

	Draw3D::drawAxis ( 3.0f );

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
















