
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw3D.h"
//#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
//#include "Mat3.h"
//#include "quaternion.h"
//#include "raytrace.h"
//#include "Body.h"

//#include "geom3D.h"

#include "kBoxes.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

// ============= Application

/*

ToDo: kBoxes.h seems to be good choice for dipole distribution
      (we should not reinvetn the wheel unless we can do much better efficieny)
      but it is probably not hierarchical !

*/


class MultipoleAccel{
    int nps=0;
    Vec3d*  ps=0; // points
    double* qs=0; // point charges
    int*    p2c =0; // maping points to centers
    int*    p2cT=0;
    std::vector<Vec3d> centers;
    std::vector<int>   levels;
    double* charges=0;
    Vec3d * dipoles=0;

    int getLevelN( int i ){
        int n1;
        if(i<levels.size()-1){ n1=levels[i+1]; }else{ n1=centers.size(); };
        return n1-levels[i];
    }

    void swapp(int i,int j){
        //Vec3d tp; double td; int ti;
        //ps[i]=
        _swap(ps[i],ps[j]);
        _swap(qs[i],qs[j]);
    }

    void KmeansRecur( int k, int ip0, int np, int nlevel ){
        int ic0=centers.size();
        // initialze new centers
        centers.resize(ic0+k);
        //Vec3d cogs [k];
        int   npbox[k];
        int   i0box[k];
        for(int ic=0; ic<k; ic++){
            int ip = rand()&(np-ic)+ic;
            swapp(ip0+ip,ip0+ic);
            Vec3d p = ps[ip0+ic];
            centers[ic0+ic] = p;
            //cogs [ic]=p;
            npbox [ic]=0;
            i0box [ic]=0;
        }
        // assign points to centers
        int ipmax=ip0+np;
        for(int ip=ip0; ip<ipmax; ip++){
            int    icbest =0;
            double r2best =1e+300;
            Vec3d  p  = ps[ip];
            double q2 = qs[ip]; q2*=q2;
            for(int ic=0; ic<k; ic++){
                Vec3d d = centers[ic0+ic]; d.sub(p);
                double r2 = d.norm2();
                if(r2<r2best){ r2best=r2; icbest=ic; }
            }
            npbox[icbest]++;
            //cogs [icbest].add( p );
            p2c  [ip]=icbest+ic0;
            //centers[icbest];
        }
        // prepare boxes size and starting index
        int nsum=0;
        for(int ic=0; ic<k; ic++){
            int ni    = npbox[ic];
            i0box[ic] = nsum;
            nsum     += ni;
            npbox[ic] = 0;
        }
        // find permut index
        for(int ip=ip0; ip<ipmax; ip++){
            int ic   = p2c  [ip];
            int ip0c = i0box[ic];
            int ni   = npbox[ic];
            //int npc  = npbox [ic];
            //swapp(ip,ip0c+ni); // ToDo : Warrning - this can delete some unprocessed cells
            //p2cT[ip0c+ni] = ip;
            // ToDO : perhaps we need allocate new array for points
            npbox[ic]++;
        }
        //
    }

};

class TestAppMultipoleAccel : public AppSDL2OGL_3D { public:


	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppMultipoleAccel( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppMultipoleAccel::TestAppMultipoleAccel( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

}

void TestAppMultipoleAccel::draw(){
    //rintf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

};


void TestAppMultipoleAccel::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}


void TestAppMultipoleAccel::drawHUD(){
    glDisable ( GL_LIGHTING );
}

// ===================== MAIN

TestAppMultipoleAccel * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppMultipoleAccel( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















