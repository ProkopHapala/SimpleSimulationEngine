
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <unordered_map>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"
#include "Solids.h"

#include "DynamicOpt.h"
#include "SoftBody.h"
#include "CubicRuler.h"


#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"


// ======================  TestApp

class TestAppSphereTree : public AppSDL2OGL_3D {
	public:

    bool running = false;
    int perFrame = 10;

	int nSpheres = 0;
	const static int nMaxSpheres = 1024;
	Vec3d sphere_pos[nMaxSpheres];
	//std::vector<Vec3d>

	std::unordered_multimap<int_fast64_t,int>  grid;

	CubicRuler ruler;
    SoftBody   truss;

    int sphereShape;

    const double rStartOff = 2.0;
    double rStart          = rStartOff;
    double driftSpeed      = 0.001;
    double collisionDist2  = 1.0;
    double maxStep         = 0.25;
    Vec3d maxPos,minPos;

    bool DLArunning = true;
    Vec3d DLA_pos;

	// ---- function declarations

	bool insertSphere( const Vec3d& pos, int_fast64_t ind );
    bool insertSphere( const Vec3d& pos );

	bool DLAstep ( Vec3d& pos );
	void DLAstart( );

	virtual void draw   ();
    //virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );

	TestAppSphereTree ( int& id, int WIDTH_, int HEIGHT_ );

};

bool TestAppSphereTree::insertSphere( const Vec3d& pos, int_fast64_t ind ){
    if( nSpheres >= (nMaxSpheres-2) ) return false;
    sphere_pos[nSpheres] = pos;
    grid.insert({ind,nSpheres});
    nSpheres++;
    return true;
};

bool TestAppSphereTree::insertSphere( const  Vec3d& pos ){
    int_fast64_t ind = ruler.pos2index( pos );
    return insertSphere( pos, ind );
};

void TestAppSphereTree::DLAstart( ){
    double z   = randf( -1.0,1.0 );
    double xy  = rStart * sqrt (  1 - z*z );
    z         *= rStart;
    double phi = randf( 0.0,M_PI*2.0);
    DLA_pos.set( xy*cos(phi), xy*sin(phi), z );
}



bool TestAppSphereTree::DLAstep( Vec3d& pos ){
    Vec3d dpos;

    dpos.set( randf( -maxStep, maxStep ), randf( -maxStep, maxStep ), randf( -maxStep, maxStep ) );
    pos.add( dpos );
    double r = pos.norm();
    if( r > rStart ){ pos.mul( rStart/r ); }

/*

    if ( !pos.isBetween( minPos, maxPos ) ){
        //printf( " out of box \n" );
        return false;
    }
*/

    for( int i=0; i<nSpheres; i++ ){
        Vec3d posi = sphere_pos[i];
        double dr2 = pos.dist2( posi );
        if ( dr2 < collisionDist2 ){
            if( (r + rStartOff) > rStart ) rStart = r + rStartOff;
            insertSphere( pos );
            printf( " insert (%3.3f,%3.3f,%3.3f) %i \n", pos.x,pos.y,pos.z, nSpheres );
            return false;
        }
    }


/*
    // with HashMap acceleration

    Vec3d dabc;
    Vec3i iabc;
    ruler.pos2index( pos, dabc, iabc );
    for( int ineigh=0; ineigh<CubicRuler_nNeighs; ineigh++ ){
        int ix = iabc.x + CubicRuler_neighs[ineigh][0];
        int iy = iabc.y + CubicRuler_neighs[ineigh][1];
        int iz = iabc.z + CubicRuler_neighs[ineigh][2];
        int_fast64_t ind = xyz2id( ix , iy, iz );
        //int_fast64_t ind = ruler.pos2index( pos );
        auto range = grid.equal_range( ind );
        //printf( " %i  (%i,%i,%i)  %i \n", ineigh, ix, iy, iz, ind );
        //printf( " %f  (%3.3f,%3.3f,%3.3f)  %i \n", r, pos.x, pos.y, pos.z, ind );
        for ( auto it = range.first; it != range.second; ++it ){
            Vec3d posi = sphere_pos[it->second];
            double dr2 = pos.dist2( posi );
            //printf( " > %f  (%3.3f,%3.3f,%3.3f) \n", dr2, posi.x, posi.y, posi.z );
            if ( dr2 < collisionDist2 ){
                if( r > rStart ) rStart = r;
                insertSphere( pos, ind );
                printf( " insert (%i,%i,%i) %i %i \n", ix,iy,iz, ind, nSpheres );
                return false;
            }
        }
    }
*/

    return true;
}

TestAppSphereTree::TestAppSphereTree( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    maxPos.set(  100.0, 100.0, 100.0 );
    minPos.set( -100.0,-100.0,-100.0 );

    ruler.setup( minPos*-2.0, { 1.0,0.0,0.0,  0.0,1.0,0.0,  0.0,0.0,1.0 } );

    //insertSphere( {0.0,0.0,0.0} );
    insertSphere( {1.0,-1.0,1.0} );

    sphereShape = glGenLists(1);
    glNewList( sphereShape , GL_COMPILE );
        //glPushMatrix();
        //glDisable ( GL_LIGHTING );
        //Draw3D::drawAxis ( 3.0f );
        //glColor3f( 1.0f, 0.0f, 1.0f ); Draw3D::drawLines   ( Solids::Icosahedron_nedges, Solids::Icosahedron_edges, Solids::Icosahedron_verts                             );
        //glEnable( GL_LIGHTING );
        //glColor3f( 0.8f, 0.8f, 0.8f ); Draw3D::drawPolygons( Solids::Icosahedron_nfaces, Solids::Icosahedron_ngons, Solids::Icosahedron_faces,  Solids::Icosahedron_verts );
        glEnable( GL_LIGHTING );
        //glColor3f( 0.8f, 0.8f, 0.8f );
        Draw3D::drawSphere_oct( 4, 0.6, {0.0,0.0,0.0} );
        //glPopMatrix();
    glEndList();

    DLAstart( );

};

void TestAppSphereTree::draw   (){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glEnable(GL_DEPTH_TEST);

    float ambient  [] = { 0.8f, 0.8f, 0.8f, 1.0f };
	float diffuse  [] = { 0.2f, 0.2f,  0.2f,  1.0f };
	float specular [] = { 0.5f, 0.5f,  0.5f,  1.0f };
	//float shininess[] = { 80.0f                    };
	//glMaterialfv ( GL_FRONT_AND_BACK, GL_AMBIENT,   ambient);
	//glMaterialfv ( GL_FRONT_AND_BACK, GL_DIFFUSE,   diffuse);
	//glMaterialfv ( GL_FRONT_AND_BACK, GL_SPECULAR,  specular);
	glEnable     ( GL_COLOR_MATERIAL   );
	glLightfv    ( GL_LIGHT0, GL_AMBIENT,  ambient  );
	glLightfv    ( GL_LIGHT0, GL_DIFFUSE,  diffuse  );
	glLightfv    ( GL_LIGHT0, GL_SPECULAR, specular  );
	//glMaterialfv ( GL_FRONT_AND_BACK, GL_SHININESS, shininess);

    Mat3d lrot; lrot.setOne();

    int perFrame = 10000;
    glColor3f(0.8f,0.2f,0.2f);
    for( int i=0; i<perFrame; i++ ){
        if( nSpheres < (nMaxSpheres-2) ){
            if( DLArunning ){
                DLArunning = DLAstep( DLA_pos );
                //Draw3D::drawShape  ( DLA_pos, lrot, sphereShape );
            }else{
                DLAstart( );
                DLArunning = true;
            }
        }
    }

/*
    int ip = 0;
    for(auto p : grid ) {
        if( p.first != 0 ) printf( " %i %i \n", p.first, p.second );
        ip++;
    }
*/

    //for( auto o : world.objects ) {
    //glColor3f(0.8f,0.8f,0.8f);
    for( int i=0; i<nSpheres; i++ ){
        Vec3d lpos = sphere_pos[i];
        double r = lpos.norm();
        float c  = (float)(r/rStart);
        glColor3f( c, 2*c*(1-c), 0 );
        Draw3D::drawShape    ( lpos, lrot, sphereShape );
    }

	glDisable(GL_DEPTH_TEST);

};

void TestAppSphereTree::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    int rot;
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_d:  ix ++; if( ix >= builder.nMax ) ix = 0;                break;
                //case SDLK_a:  ix --; if( ix <  0            ) ix = builder.nMax-1;   break;
            }
            break;
        /*
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
        */
    };
    AppSDL2OGL_3D::eventHandling( event );
}



// ===================== MAIN

TestAppSphereTree * testApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestAppSphereTree( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}
