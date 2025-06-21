
/// @file @brief  This demo illustrates an advanced algorithm for accelerating N-body simulations. It uses the Fast Multipole Method (FMM) via `kBoxes.h` and `HierarchicalKpivot` to approximate the influence of distant particles, reducing the computational complexity from O(N^2) to O(N). The visualization would likely show a system of particles and the hierarchical grid used to group them for the multipole expansion.
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

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

#include "kBoxes.h"

// ============= Application

/*

ToDo: kBoxes.h seems to be good choice for dipole distribution
      (we should not reinvetn the wheel unless we can do much better efficieny)
      but it is probably not hierarchical !

*/

void drawCells( int nc, const Vec2i* cellNI, const Vec3d* ps, float sz, int ik, const Vec3d* pivots=0, bool bLines=true ){
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    for(int i=0; i<nc; i++){
        Vec3f col;
        //Draw::color_of_hash( i*5454 + 646*ik + 454 );
        Draw::color_of_hash( i*5454 + 646*ik + 454, col );
        glColor4f(col.x,col.y,col.z,1.0);
        const Vec2i& c = cellNI[i];
        for(int j=0; j<c.a; j++){
            glColor4f(col.x,col.y,col.z,1.0);
            Draw3D::drawPointCross( ps[c.b+j], sz );
            glColor4f(col.x,col.y,col.z,0.3);
            if(pivots&&bLines){ Draw3D::drawLine(pivots[i],ps[c.b+j]); }
        }
        if(pivots){
            Draw3D::drawPointCross( pivots[i], sz*3 );
        }
    }
}

/*
void KpivotsRecur( int np, Vec3d* ps, Vec3d* ps_, int* o2c, int* c2o, int perCellMax, int kmax, int lev ){
    if(lev>3)return;
    int K = _min( 1+np/perCellMax, kmax );
    if(K<=1)return;
    printf( "%0*cK %i np %i \n", lev+1, '+', K, np );
    //if(np<100)for(int i=0; i<np; i++){ printf( "ps[%i] (%g,%g,%g)\n", i, ps[i].x,ps[i].y,ps[i].z ); }
    KPivots kpiv(K);
    kpiv.build( np, ps, o2c, c2o, true, false, 1.0, 2.0 );
    applyPermutTmp( np, c2o, ps, ps_ );
    //if(np<100)for(int i=0; i<np; i++){ printf( "ps_[%i] c2o[] %i (%g,%g,%g) \n", i, c2o[i], ps[i].x,ps[i].y,ps[i].z ); }
    //return;
    int nsum=0;
    for(int k=0; k<kpiv.K; k++){
        int ni = kpiv.cellNIs[k].a;
        nsum+=ni;
        printf( "%0*c[%i] ni %i \n", lev+2, '.', k, ni );
        if( ni>perCellMax ){
            int i0=kpiv.cellNIs[k].b;
            KpivotsRecur( ni, ps+i0, ps_+i0, o2c+i0, c2o+i0, perCellMax, kmax, lev+1 );
        }
    }
    //printf( "%0*cnsum %i np %i \n", lev, '-', nsum, np );
}
*/

class MultipoleAccel{ public:
    int     np=0;
    Vec3d*  ps=0; // points
    double* qs=0; // point charges
    int*    p2c =0; // maping points to centers
    int*    p2cT=0;
    std::vector<Vec3d> centers;
    std::vector<int>   levels;
    double* charges=0;
    Vec3d * dipoles=0;
};


class TestAppMultipoleAccel : public AppSDL2OGL_3D { public:

    MultipoleAccel mpol;
    KPivots* kpiv=0;
    HierarchicalKpivot* hkpiv=0;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppMultipoleAccel( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppMultipoleAccel::TestAppMultipoleAccel( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {
    mpol.np  = 40;
    mpol.ps  = new Vec3d[mpol.np];

    Box span; span.setSymetric(0.5);
    Vec3d pos = Vec3dZero;
    for(int i=0;i<mpol.np;i++){
        pos.add( span.genRandomSample() );
        mpol.ps[i] = pos;
    }

    /*
    //kpiv = new KPivots( (int)(sqrt(mpol.np)) );
    kpiv = new KPivots( 10 );
    //kpiv->build( mpol.np, mpol.ps );
    //kpiv->build( mpol.np, mpol.ps, 0, 0, true, false, -5.0, 5.0 );
    kpiv->build( mpol.np, mpol.ps, 0, 0, true, true, -5.0, 5.0 );
    applyPermut_relloc( mpol.np, mpol.ps, kpiv->c2o );
    */

    /*
    // ========= Test Hierarchical   .... Currently not working  ... WIP
    Vec3d* ps_ = new Vec3d[mpol.np];
    int*   o2c = new int  [mpol.np];
    int*   c2o = new int  [mpol.np];
    KpivotsRecur( mpol.np, mpol.ps, ps_, o2c, c2o, 4, 4, 0 );
    */

    hkpiv = new HierarchicalKpivot(mpol.np, mpol.ps);
    hkpiv->run( 4, 4, 3 );

}

void TestAppMultipoleAccel::draw(){
    //rintf( " ==== frame %i \n", frameCount );
    glClearColor( 0.9f, 0.9f, 0.9f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	//drawCells( kpiv->K, kpiv->cellNIs, kpiv->c2o , mpol.ps, 0.1, 0*frameCount/100, kpiv->pivots, true );

	drawCells( hkpiv->cells.size(), &hkpiv->cells[0], mpol.ps, 0.05, 0*frameCount/100, &hkpiv->pivots[0], true );

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
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppMultipoleAccel( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}





/*
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
*/
