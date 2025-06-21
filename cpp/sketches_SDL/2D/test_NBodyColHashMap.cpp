/// @file @brief This demo simulates 2D N-body collisions, leveraging `HashMap2D.h` for optimized broad-phase collision detection. It visualizes the dynamic interaction of numerous circular bodies, where the hash map efficiently identifies potential collision pairs, and the simulation resolves their physical responses. The simulation can be toggled with the `SPACE` key, and bodies can be dragged with the mouse.

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

//#include <unordered_set>

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"
#include "HashMap2D.h"
#include "Body2D.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "AppSDL2OGL.h"
#include "testUtils.h"

// ===========================
//          NBodyWorld
//============================

inline void stringForce( const Vec2d& pa, const Vec2d& pb, double k, Vec2d& Fout ){
    Vec2d d;
    d.set_sub( pb, pa );
    Fout.set( d.x*k, d.y*k );
};

inline bool interactPairwise( PointBody2D& pa, PointBody2D& pb ){
    const double r2max = 1.0d;
    Vec2d d;
    d.set_sub( pa.pos, pb.pos );
    double r2 = d.norm2();
    if( r2 < r2max ){
        d.normalize();
        double ca = d.dot( pa.vel );
        double cb = d.dot( pb.vel );
        if( ca < cb ){ // closing in
            pa.vel.add_mul( d, cb-ca );
            pb.vel.add_mul( d, ca-cb );
        }
    }
    return false;
};

class NBodyWorld{
	public:
    double dt_frame  = 0.1;
    int    per_frame = 10;
    double damping   = 0.2;

    double anchorStiffness = 0.5;

    double damp;
    double dt;

    double v2max;
    double f2max;

    Rect2d boundBox;
    UINT    nBuckets;
    ULONG * buckets;


	int nParticles;
	PointBody2D* particles;
	HashMap2D<PointBody2D> map;

    //std::unordered_set<ULONG> activeCells;
    //std::unordered_set<ULONG> activeCellsNeighbors;
    //int nActiveParticles;
    //PointBody2D** activeParticles;

    Vec2d anchor;
    PointBody2D* picked = NULL;

    int n_moves, n_interactions, n_cellReads;
    int n_moves_sum, n_interactions_sum, n_cellReads_sum;

    void init();
    void update();
    void simulationStep( double dt );
    void assembleForces( ULONG i );
    void assembleForces_offside( ULONG i, ULONG j, UINT ni, PointBody2D** buf_i );
    bool moveParticle     ( PointBody2D* pi );

    inline void setSimParams( double dt_frame_, double per_frame_, double damping_ ){
        dt_frame  = dt_frame_;
        per_frame = per_frame_;
        damping   = damping_;
        evalAuxSimParams();
        //printf( " dt_frame, per_frame,  dt, damp " );
    }

    inline void evalAuxSimParams(){
        const double dampMin = 0.5;
        dt   = dt_frame / per_frame;
        damp = ( 1 - damping * dt );
        if( damp < dampMin ){ damp = dampMin; }
    };

};

void NBodyWorld::update( ){
    for( int i=0; i<per_frame; i++  ){
        n_moves=0; n_interactions=0;   n_cellReads=0;
        simulationStep( dt );
        n_moves_sum        += n_moves;        n_interactions_sum += n_interactions;        n_cellReads_sum    += n_cellReads;
        printf( " ==== DONE sub_step %i moves %i pairs %i cells %i \n", i, n_moves, n_interactions, n_cellReads );
    }
}

void NBodyWorld::simulationStep( double dt ){
    /*
    for( int i=0; i<nParticles; i++ ){
        particles[i].force.set( 0.0, 0.0 );
    }
    */

    for( int i=0; i<nBuckets; i++ ){ assembleForces( buckets[i] ); }

    //exit(0);
/*
    if( picked != NULL ){
        Vec2d fstring;
        stringForce( picked->pos, anchor, anchorStiffness, fstring );
        picked->force.add( fstring );
    }
*/
    v2max=0; f2max=0;
    for( int i=0; i<nParticles; i++ ){
        PointBody2D* pi = particles+i;
        bool active     = moveParticle( pi );
    }
};



void NBodyWorld::assembleForces( ULONG i ){
    // BE WARE : particle->force should be cleaned before we start
    // onside step
    PointBody2D*  buf_i_[256];
    PointBody2D** buf_i = &buf_i_[0];
    UINT ni = map.HashMap<PointBody2D>::getBucketObjects( i, buf_i );
    n_cellReads++;
    for(int ii=0; ii<ni; ii++ ){
        PointBody2D* pi = buf_i[ii];
        // confine box
        if( ( pi->pos.x < boundBox.x0 ) && ( pi->vel.x < 0 ) ){  pi->vel.x *= -1; };
        if( ( pi->pos.x > boundBox.x1 ) && ( pi->vel.x > 0 ) ){  pi->vel.x *= -1; };
        if( ( pi->pos.y < boundBox.y0 ) && ( pi->vel.y < 0 ) ){  pi->vel.y *= -1; };
        if( ( pi->pos.y > boundBox.y1 ) && ( pi->vel.y > 0 ) ){  pi->vel.y *= -1; };
        for(int jj=0; jj<ii; jj++ ){
            PointBody2D* pj = buf_i[jj];
            bool interacts = interactPairwise( *pi, *pj );
            n_interactions++;
            //DEBUG_PLOT_INTERACTION( pi, pj, 0.1f, 0.9f, 0.1f )
        }
    }
    // offside part
    UHALF ix,iy;
    map.unfoldBucketInt( i, ix, iy );
    assembleForces_offside( i, map.getBucketInt( ix-1, iy-1 ), ni, buf_i );
    assembleForces_offside( i, map.getBucketInt( ix  , iy-1 ), ni, buf_i );
    assembleForces_offside( i, map.getBucketInt( ix+1, iy-1 ), ni, buf_i );
    assembleForces_offside( i, map.getBucketInt( ix-1, iy   ), ni, buf_i );
    //         onside                          ix   iy
    assembleForces_offside( i, map.getBucketInt( ix+1, iy   ), ni, buf_i );
    assembleForces_offside( i, map.getBucketInt( ix-1, iy+1 ), ni, buf_i );
    assembleForces_offside( i, map.getBucketInt( ix  , iy+1 ), ni, buf_i );
    assembleForces_offside( i, map.getBucketInt( ix+1, iy+1 ), ni, buf_i );
};

void NBodyWorld::assembleForces_offside( ULONG i, ULONG j, UINT ni, PointBody2D** buf_i ){
    if( i < j ){
        PointBody2D* buf_j[256];
        UINT nj = map.HashMap<PointBody2D>::getBucketObjects( j, buf_j );
        n_cellReads++;
        for(int ii=0; ii<ni; ii++ ){
            PointBody2D* pi = buf_i[ii];
            for(int jj=0; jj<nj; jj++ ){
                PointBody2D* pj = buf_j[jj];
                bool interacts = interactPairwise( *pi, *pj );
                //printf( " %i %i   %i %i  (%3.3f,%3.3f)(%3.3f,%3.3f)\n",   i, j,  ii, jj,  pi->pos.x,pi->pos.y,  pj->pos.x,pj->pos.y );
                n_interactions++;
                //DEBUG_PLOT_INTERACTION( pi, pj, 0.9f, 0.1f, 0.9f )
            }
        }
    }
};

bool NBodyWorld::moveParticle( PointBody2D* pi ){

    //pi->vel.mul( damp );
    ULONG old_index = map.getBucket( pi->pos.x, pi->pos.y );
    pi->move_PointBody2D( dt );
    ULONG new_index = map.getBucket( pi->pos.x, pi->pos.y );
    n_moves++;
    if( old_index != new_index ){
        bool removed = map.HashMap<PointBody2D>::tryRemove  ( pi, old_index );
        if( removed ){
            int iinsert = map.HashMap<PointBody2D>::insertIfNew( pi, new_index );
            if( iinsert < 0 ){
                printf( " cannot insert! inconsistent HashMap! \n" );
                exit(0);
            }
        }else{
            printf( " cannot remove! inconsistent HashMap! \n" );
            exit(0);
        }
    }
    return true;
}

void NBodyWorld::init(){
    evalAuxSimParams();

    //int power = 8; int nside = 1;
    int power = 16; int nside = 60;
    //int power = 16; int nside = 300;
    nParticles = (2*nside+1)*(2*nside+1);
    //nParticles = 4*nside*nside;
	particles  = new PointBody2D[nParticles];

    double cellsize = 2.5d;
    map.init( cellsize, power );
	printf( "map: %i %i %i %i \n", map.power, map.mask, map.capacity, map.filled );

    boundBox.set( -nside*cellsize, -nside*cellsize,  nside*cellsize, nside*cellsize );

    buckets  = new ULONG[ 65536 ];
    nBuckets = map.getBucketsInRect( boundBox.x0-cellsize*0.5f, boundBox.y0-cellsize*0.5f, boundBox.x1+cellsize*0.5f, boundBox.y1+cellsize*0.5f, buckets );
    printf( "nBuckets: %i boundBox %f %f %f %f \n", boundBox.x0, boundBox.y0, boundBox.x1, boundBox.y1, nBuckets );

	int i = 0;
	for( int iy=-nside; iy<nside; iy++ ){
		for( int ix=-nside; ix<nside; ix++ ){
			particles[i].vel.set( 0.0, 0.0 );
			particles[i].setMass( 1.0 );
			particles[i].pos.set(
                ( ix + randf(0.2,0.8) ) * map.step*0.4,
                ( iy + randf(0.2,0.8) ) * map.step*0.4 );
			particles[i].vel.set( randf(-1.0,1.0), randf(-1.0,1.0) );
			//map.insertNoTest( &(points[i]), points[i].x, points[i].y  );
			map.insertIfNew( &(particles[i]), particles[i].pos.x, particles[i].pos.y  );
			//printf( " insering %03i-th particle  (%3.3f,%3.3f) to (%i,%i) %i \n", i,  particles[i].pos.x,  particles[i].pos.y, map.getIx( particles[i].pos.x), map.getIy( particles[i].pos.y), map.getBucket( particles[i].pos.x, particles[i].pos.y) );
            i++;
		};
	};
	nParticles = i;
	printf( "map: %i %i %i %i \n", map.power, map.mask, map.capacity, map.filled );

	//exit(1);

};










// ===========================
//          NBodyWorldApp
//============================


class NBodyWorldApp : public AppSDL2OGL {
	public:
    NBodyWorld world;

	virtual void draw   ();
	virtual void drawHUD();
	virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );

	void pickParticle( PointBody2D*& picked );

	NBodyWorldApp( int& id, int WIDTH_, int HEIGHT_ );

};

NBodyWorldApp::NBodyWorldApp( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {
    world.init();
}

void NBodyWorldApp::draw(){
    long tstart = getCPUticks();
    glClearColor( 0.9f, 0.9f, 0.9f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    world.n_moves_sum = 0; world.n_interactions_sum = 0; world.n_cellReads_sum = 0;
    long t1 = getCPUticks();
    world.update();
    long t2 = getCPUticks();

    glColor3f( 0.6f, 0.6f, 0.6f );
    for( int i=0; i<world.nBuckets; i++ ){
        ULONG icell = world.buckets[ i ];
        double x,y;
        world.map.unfoldBucket( icell, x, y );
        Draw2D::drawRectangle( (float)x, (float)y, (float)(x+world.map.step), (float)(y+world.map.step), false );
    }

    PointBody2D* screenObjects[65536];
    //float camXmin_ =-1; float camXmax_ =+1;
    //float camYmin_ =-1; float camYmax_ =+1;
    float camXmin_ =camXmin; float camXmax_ =camXmax;
    float camYmin_ =camYmin; float camYmax_ =camYmax;

    UHALF ix0 = world.map.getIx( camXmin_ );  UHALF iy0 = world.map.getIy( camYmin_ );
    UHALF ix1 = world.map.getIx( camXmax_ );  UHALF iy1 = world.map.getIy( camYmax_ );
    UINT nfound_tot = 0;
    int ncells  = 0;
    for( UHALF iy = iy0; iy<=iy1; iy++ ){
        for( UHALF ix = ix0; ix<=ix1; ix++ ){
            UINT nfoundi = world.map.getBucketObjectsInt( ix, iy, screenObjects );
            nfound_tot += nfoundi;
            for( int i=0; i<nfoundi; i++ ){
                PointBody2D* p = screenObjects[i];
                glColor3f( 0.2f, 0.2f, 0.2f );
                Draw2D::drawCircle_d( p->pos, 0.5, 8, true );
                glColor3f( 0.8f, 0.2f, 0.2f );
                Vec2d pvel; pvel.set_mul( p->vel, 5.0 ); pvel.add( p->pos );
                Draw2D::drawLine_d( p->pos, pvel );
            }
        }
    }

    Draw2D::drawPointCross_d( world.anchor, 0.5 );
    if( world.picked != NULL ) Draw2D::drawLine_d( world.anchor, world.picked->pos );

    long tend = getCPUticks();

    //printf( " ======== frame %i DONE ( map.filled %i nfound_tot %i )\n", frameCount, world.map.filled, nfound_tot );
    printf( " ======== frame %i DONE T=%3.3f Mticks/frame( %3.3f sim. ) moves %i pairs %i cells %i \n",
           frameCount, (tend-tstart)*1.0e-6, (t2-t1)*1.0e-6, world.n_moves_sum, world.n_interactions_sum, world.n_cellReads_sum );

    //n_moves_sum        += n_moves;        n_interactions_sum += n_interactions;        n_cellReads_sum    += n_cellReads;
    //printf( " ==== DONE sub_step %i moves %i pairs %i cells \n", i, v2max, f2max, n_moves, n_interactions, n_cellReads );

	//STOP = true;

};

void NBodyWorldApp::mouseHandling( ){
    Uint32 buttons = SDL_GetMouseState( &mouseX, &mouseY );
    defaultMouseHandling( mouseX, mouseY );
    world.anchor.set( mouse_begin_x, mouse_begin_y );
}


void NBodyWorldApp::pickParticle( PointBody2D*& picked ){
    PointBody2D* screenObjects[256];
    // mouse picking
    double rmax  = 2.0d;
    double r2max = rmax*rmax;
    Vec2d vmouse;
    vmouse.set( mouse_begin_x, mouse_begin_y );
    UINT mfound = world.map.getObjectsInRect( (float)(vmouse.x - rmax ), (float)(vmouse.y - rmax ), (float)(vmouse.x + rmax ), (float)(vmouse.y + rmax ), &(screenObjects[0]) );
    //printf( "mfound  %i \n", mfound );
    int imin = -1;
    double r2min = 1e+300;
    for( int i=0; i<mfound; i++ ){
        PointBody2D* p = screenObjects[i];
        Vec2d d;
        d.set_sub( p->pos, vmouse );
        double r2 = d.norm2();
        if( r2 < r2max ){
            if( r2 < r2min ){  r2min = r2; imin = i; }
        }
        //printf( " r2 %3.3f r2max %3.3f r2min %3.3f imin %i \n", r2, r2max, r2min, imin );
    }
    if( imin >= 0 ){
        picked = screenObjects[imin];
    }else{
        picked = NULL;
    }
    //if( world.picked != NULL ) printf( " MOUSE picked (%f,%f) \n", world.picked->pos.x, world.picked->pos.y );
}

void NBodyWorldApp::eventHandling ( const SDL_Event& event  ){
    //printf( "NBodyWorldApp::eventHandling() \n" );
    switch( event.type ){
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    pickParticle( world.picked );
                break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    world.picked = NULL;
                    break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}

void NBodyWorldApp::drawHUD(){

}

// ===================== MAIN

NBodyWorldApp * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	thisApp = new NBodyWorldApp( junk , 800, 600 );
	thisApp->camStep = 10;
	thisApp->delay = 0;
	thisApp->zoom = 30;
	thisApp->loop( 1000000 );
	return 0;
}
