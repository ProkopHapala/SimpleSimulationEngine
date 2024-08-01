
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw3D.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "Mat4.h"
//#include "VecN.h"


#include "AppSDL2OGL_3D.h"
#include "testUtils.h"
#include "SDL_utils.h"
#include "Plot2D.h"

//#include "MMFF.h"
//#define R2SAFE  1.0e-8f

//#include "eFF.h"
//#include "e2FF.h"


class TestAppSp3Space: public AppSDL2OGL_3D { public:

    //RigidAtom     atom1;
    //RigidAtomType type1,type2;

    bool bRun = false;


    std::vector<Vec3d> ps{ 1,Vec3dZero };


    Mat4d orbs;

    int ipicked  = -1, ibpicked = -1;

    //Plot2D plot1;

    //double Emin,Emax;
    //int     npoints;
    //Vec3d*  points  =0;
    //double* Energies=0;
    //Vec3d * Forces  =0;

    int      fontTex;

    virtual void draw   ();
    virtual void drawHUD();
    //virtual void mouseHandling( );
    virtual void eventHandling   ( const SDL_Event& event  );
    //virtual void keyStateHandling( const Uint8 *keys );

    TestAppSp3Space( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppSp3Space::TestAppSp3Space( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );

    Mat4d UU;

    /*
    orbs.setOne();
    double d = 1.0;
    for(int i=0;i<16;i++){orbs.array[i]+=randf(-d,d);}
    UU.set_mmul_NT( orbs, orbs );
    printf("orbs before\n"); orbs.print();
    printf("UU before\n"); UU.print();
    */

    // tetrahedron
    double s = 0.81649658092;
    double c = 0.57735026919;

    //orbs.vecs[0] = (Quat4d){.0,-s,+c,1.0d};
    //orbs.vecs[1] = (Quat4d){.0,+s,+c,1.0d};
    //orbs.vecs[2] = (Quat4d){-s,.0,-c,1.0d};
    //orbs.vecs[3] = (Quat4d){+s,.0,-c,1.0d};


    double sqrt3 = sqrt(3.0); // normalization of spherical harmonsi
    double invsqrt3 = 1/sqrt3;
    orbs.vecs[0] = (Quat4d){-1,-1,+1,1};
    orbs.vecs[1] = (Quat4d){+1,+1,+1,1};
    orbs.vecs[2] = (Quat4d){+1,-1,-1,1};
    orbs.vecs[3] = (Quat4d){-1,+1,-1,1};

    orbs.normalize();

    /*
    double errMax = 1e-6;
    int niter = orbs.orthoRun( sq(errMax), 100 );
    printf("orthoRun n=%i err<%g \n", niter, errMax );
    */

    UU.set_mmul_NT( orbs, orbs );
    printf("UU after\n"); UU.print();
    printf("orbs after\n"); orbs.print();


    Vec3d p=Vec3dZero;
    for(Vec3d& pi : ps){
        p.addRandomCube(0.2);
        pi=p.normalized();
    }
    for(Vec3d& pi : ps){ print(pi); printf("\n"); }

    /*

    // ============ TO DO:

    1) ==============
    Orthogonalization can be done more efficietly.
    We can orthogonalize just two vectors vi, vj than we don't have to recalculate the overlap matrix Cij=<vi|vj>

        algorithm ( similar to Jacobi Rotation ):
            1) O(n^3) : calculate initial matrix Cij=<vi|vj>
            2) O(n^2) : find the maximum value  i,j = argmax( Cij ), make list cimax[j]
            3) O(n^1) : orthogonalize these two vectors
                e.g.   vi_ = vi - 0.5*cij*vj;   vj_ = vj - 0.5*cij*vi
            4) O(n^1) : update
                cii <= 1 - cij^2/4     = <(vi - 0.5*cij*vj)^2>                 = <vi|vi> - 0.5cij<vi|vj> + 0.25cij^2<vj|vj>   = 1 - 0.5*cij^2 + 0.25cij^2
                cij <=     cij^3/4     = < vi - 0.5*cij*vj| vj - 0.5*cij*vi >  = <vi|vj> - 0.5*cij*(<vi|vi>+<vj|vj>) + 0.25*cij^2*<vi|vj> = cij-cij  + 0.25*cij^3
                cik <=     cik - cjk*cij/2
            5) O(n^1) : update cimax[j] and find new i,j = argmin(cij)
            6) goto 3)

    2) ==============

    For symmetric orthogonalization of two vectors {vi,vj} = {a,b} we can use better update than vi_ = vi - 0.5*cij*vj;   vj_ = vj - 0.5*cij*vi
    we can make symetrized (half-angle) coordinate system
        u=(a+b)/(2c);   v=(a-b)/(2s)
        where 2c = |a+b| = 2*cos(a/2) = 2*sqrt((1+cij)/2)
              2s = |a-b| = 2*sin(a/2) = 2*sqrt((1-cij)/2)
        then a_ = (u+v)/sqrt2 = sqrt2*( (a+b)/(2c) + (a-b)/(2s) )  =  sqrt2*( s*(a+b) + c*(a-b) )/(2s*c) = sqrt2*( a*(s+c) + b*(s-c) )/(2s*c)
             b_ = (u-v)/sqrt2 = sqrt2*( (a+b)/(2c) - (a-b)/(2s) )  =  sqrt2*( s*(a+b) - c*(a-b) )/(2s*c) = sqrt2*( a*(s-c) + b*(s+c) )/(2s*c)

    3) ==============
    We can also make sure forces keeps orthogonality == purely rotational forces
        ui,uj are orthogonal vectors
        fi = dui/dt   fj = duj/dt
        if ui,uj should be kept orthogonal than any rotation applied to ui must be also applied to uj and wise versa
        rotation is charactrerized by sine of angle sij = <fi|uj>*dt;   sji = <fj|ui>*dt;
        fi_ = fi - Sum_j{ sji * uj } = fi - Sum_j{ <fj|ui>uj }




    */



}

void TestAppSp3Space::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);


    glColor3f(0.0,0.0,0.0);




    Mat4d f; f.set(.0d);

    double K = 1.1;

    for(int i=0; i<4; i++){
        for(const Vec3d& p : ps ){
            double cik = orbs.vecs[i].s + orbs.vecs[i].p.dot( p )*M_SQRT2;
            double E   = K*cik*cik;
            f.vecs[i].add_mul( {p.x,p.y,p.z,1}, K*cik );
            //f.vecs[i].add_mul( {p,1}, K*cik );
        }
    }

    f.makeOrthoU( orbs );
    //for(int i=0; i<4; i++){  // remove force component which breaks normality
    //    f.vecs[i].add_mul( orbs.vecs[i], -orbs.vecs[i].dot( f.vecs[i] ) );
    //}

    if(bRun){
        orbs.add_mul( f, 0.1 );

        orbs.orthoRun( 1e-6, 10 );
    }
    //orbs.orthoRun( 1e-6, 10 );


    for(int i=0; i<4; i++){
        // 0.5 0.86
        double s  = orbs.vecs[i].s;
        double c  = sqrt( 1-s*s );
        //printf( "%i : sc %g %g \n", i, s, c );
        double l1 =  c*(1+s);   // s+p
        double l2 = -c*(1-s);   // s-p
        glColor3f(1.0,1.0,1.0);
        Draw3D::drawArrow( orbs.vecs[i].p*l2, orbs.vecs[i].p*l1, 0.1 );
        glColor3f(1.0,0.0,0.0);
        Draw3D::drawVecInPos( f.vecs[i].p, orbs.vecs[i].p*l1 );

    }

    glColor3f(0.0,0.0,0.0);
    Draw3D::drawPoints(ps.size(),&ps[0],0.1);

    //Draw3D::drawAxis(1.5);


};


void TestAppSp3Space::drawHUD(){


	glTranslatef( 100.0,100.0,0.0 );
	glScalef    ( 10.0,10.00,1.0  );
	//plot1.view();

}

/*
void TestAppSp3Space::mouseHandling( ){
    int mx,my; Uint32 buttons = SDL_GetRelativeMouseState( &mx, &my);
    if ( buttons & SDL_BUTTON(SDL_BUTTON_LEFT)) {

    }
    AppSDL2OGL_3D::mouseHandling( );
};
*/

void TestAppSp3Space::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
                case SDLK_SPACE: bRun = !bRun;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //ipicked = pickParticle( ff.natoms, ff.apos, ray0, (Vec3d)cam.rot.c , 0.5 );
                break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //ibpicked = pickParticle( ff.natoms, ff.apos, ray0, (Vec3d)cam.rot.c , 0.5 );
                    //printf( "dist %i %i = ", ipicked, ibpicked );
                    break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}

// ===================== MAIN

TestAppSp3Space* thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppSp3Space( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















