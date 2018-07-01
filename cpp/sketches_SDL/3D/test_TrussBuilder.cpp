
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
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
#include "TrussBuilder.h"

#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"


// ======================  TestApp

class TestAppTrussBuilder : public AppSDL2OGL_3D {
	public:

	//const static int nMaxTypes = 256;
	//int nTypes = 0;
	//int iType  = 0;
    //BondType bondTypes[nMaxTypes];

    bool running = false;
    int perFrame = 10;

    int iType  = 0;

    bool startSet = false;
    int_fast16_t ix=TrussBuilder::ioff,iy=TrussBuilder::ioff,iz=TrussBuilder::ioff;
    int_fast16_t oix,oiy,oiz;
    int cursorShape;

    TrussBuilder builder;
    SoftBody     truss;

	// ---- function declarations
	virtual void draw   ();
    //virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );

	void op_enter     ();
	void op_backspace ();
	void op_fixNode   ();
	void drawTruss  ( bool DEBUG );
	void drawBuilder( bool DEBUG );

	TestAppTrussBuilder( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppTrussBuilder::TestAppTrussBuilder( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    builder.init( 256, 256*8, 16 );

    builder.insertBox( ix, iy, iz, 0b111, &builder.bondTypes[iType] );
    printf( " nnodes %i  nbonds %i \n", builder.nodes.size(), builder.bonds.size() );

    cursorShape = glGenLists(1);
    glNewList( cursorShape , GL_COMPILE );
        glPushMatrix();
        glScalef( 0.5f, 0.5f, 0.5f );
        glColor3f(0.01f,0.01f,0.01f); Draw3D::drawLines    ( Solids::Cube_nedges, (int*)Solids::Cube_edges, Solids::Cube_verts );
        //Draw3D::drawAxis( 0.5f );
        glPopMatrix();
    glEndList();

};

void TestAppTrussBuilder::drawTruss( bool DEBUG ){
    if( ( truss.bonds != NULL )&&( truss.points != NULL )  ){
        glBegin( GL_LINES );
        for( int i=0; i<truss.nbonds; i++ ){
            Bond& bond = truss.bonds[i];
            Vec3d& pi  = truss.points[bond.i];
            Vec3d& pj  = truss.points[bond.j];
            glVertex3f( (float)pi.x, (float)pi.y, (float)pi.z );
            glVertex3f( (float)pj.x, (float)pj.y, (float)pj.z );
            if( DEBUG ){
                printf( " %i  %i %i   (%3.3f,%3.3f,%3.3f)  (%3.3f,%3.3f,%3.3f)\n",  i,   bond.i, bond.j,  pi.x, pi.y, pi.z,  pj.x, pj.y, pj.z  );
            }
        }
        glEnd();
    }
}

void TestAppTrussBuilder::drawBuilder( bool DEBUG ){
    glColor3f( 0.0f, 0.0f, 0.0f );
    for( int i=0; i<builder.nodes.size(); i++ ){
        //truss.points[ node.id ] = node.pos; // do we need this at all ?
        GridNode& node = builder.nodes[i];
        Vec3d pos;
        builder.index2pos( {node.ix,node.iy,node.iz}, pos );
        //printf(  " (%i,%i,%i)  (%3.3f,%3.3f,%3.3f) \n" , node.ix,node.iy,node.iz, pos );
        Draw3D::drawPointCross( pos, 0.2 );

        if( node.fixed ){
            Draw3D::drawSphereOctLines  ( 16, 0.2, pos );
        }
    }
    glBegin(GL_LINES);
    for( auto it : builder.bonds ){
        Bond& bond = it.second;
        //printf( " %i %i %i \n", bond.i, bond.j, builder.nodes.size()  );
        GridNode& node1 = builder.nodes[bond.i];
        GridNode& node2 = builder.nodes[bond.j];
        Vec3d p1,p2;
        builder.index2pos( {node1.ix,node1.iy,node1.iz}, p1 );
        builder.index2pos( {node2.ix,node2.iy,node2.iz}, p2 );
        glVertex3f( (float)p1.x, (float)p1.y, (float)p1.z );
        glVertex3f( (float)p2.x, (float)p2.y, (float)p2.z );
        //exit(0);
    }
    glEnd();
}


void TestAppTrussBuilder::draw   (){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glEnable(GL_DEPTH_TEST);

    //printf( " ==== frame %i \n", frameCount );
    //printf( " perspective %i first_person %i \n", perspective, first_person );
    //printf( " %i %i %i   %i \n", ix, iy, iz, builder.nMax );

    if(running){
        truss.dt   = 0.005;
        truss.damp = 0.5;
        for( int i=0; i<perFrame; i++ ){
             truss.step( );
        }
        drawTruss  ( false );
        //exit(0);
    }else{
        drawBuilder( false );
    }

    if( startSet ){

        Vec3d p1,p2;
        builder.index2pos( {oix,oiy,oiz}, p1 );
        builder.index2pos( {ix,iy,iz}, p2 );
        glColor3f( 1.0f, 1.0f, 1.0f ); Draw3D::drawLine( p1, p2 );
    }

    int ioff = (TrussBuilder::ioff);
    glDisable(GL_LIGHTING);
    glColor3f(0.5f,0.0f,0.0f); Draw3D::drawScale( {    0.0,    0.0,0.0}, {ix-ioff,0.0    ,    0.0}, {0.0,1.0,0.0}, 1.0, 0.1, 0.1 );
    glColor3f(0.0f,0.5f,0.0f); Draw3D::drawScale( {ix-ioff,    0.0,0.0}, {ix-ioff,iy-ioff,    0.0}, {1.0,0.0,0.0}, 1.0, 0.1, 0.1 );
    glColor3f(0.0f,0.0f,0.5f); Draw3D::drawScale( {ix-ioff,iy-ioff,0.0}, {ix-ioff,iy-ioff,iz-ioff}, {1.0,0.0,0.0}, 1.0, 0.1, 0.1 );
	Draw3D::drawAxis ( 3.0f );

    //printf( " %e   (%3.3f,%3.3f,%3.3f) \n", t, hitPos.x, hitPos.y, hitPos.z );
    //glColor3f(0.0f,1.0f,1.0f); Draw3D::drawPointCross   ( hitPos, 0.1 );
    //glColor3f(1.0f,0.0f,1.0f); Draw3D::drawVecInPos( normal,   hitPos );

	glDisable(GL_DEPTH_TEST);

};


void TestAppTrussBuilder::op_enter(){
    if( startSet ){
        //printf( " insert \n" );
        if( (ix!=oix) || (iy!=oiy) || (iz!=oiz) ){
            builder.insertBond( oix, oiy, oiz, ix, iy, iz, &builder.bondTypes[iType] );
        }
        startSet=false;
    }else{
        oix=ix; oiy=iy; oiz=iz;
        startSet=true;
    }
}

void TestAppTrussBuilder::op_backspace(){
    if( startSet ){
        if( (ix!=oix) || (iy!=oiy) || (iz!=oiz) ){
            builder.removeBond( oix, oiy, oiz, ix, iy, iz );
        }
        startSet=false;
    }else{
        oix=ix; oiy=iy; oiz=iz;
        startSet=true;
    }
}

void TestAppTrussBuilder::op_fixNode(){
    int_fast32_t inode = builder.getNodeIndex( ix, iy, iz );
    if( inode >= 0 ){
        if( builder.nodes[inode].fixed ){
            builder.nodes[inode].fixed = false;
            builder.nfixed--;
        }else{
            builder.nodes[inode].fixed = true;
            builder.nfixed++;
        }
    }
}

void TestAppTrussBuilder::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    int rot;
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_d:  ix ++; if( ix >= builder.nMax ) ix = 0;                break;
                case SDLK_a:  ix --; if( ix <  0            ) ix = builder.nMax-1;   break;
                case SDLK_w:  iy ++; if( iy >= builder.nMax ) iy = 0;                break;
                case SDLK_s:  iy --; if( iy <  0            ) iy = builder.nMax-1;   break;
                case SDLK_q:  iz ++; if( iz >= builder.nMax ) iz = 0;                break;
                case SDLK_e:  iz --; if( iz <  0            ) iz = builder.nMax-1;   break;
                case SDLK_x:  op_fixNode(); break;
                case SDLK_LEFTBRACKET:  iType ++; if( iType>=builder.bondTypes.size() ) iType = 0;                            break;
                case SDLK_RIGHTBRACKET: iType --; if( iType<0                         ) iType = builder.bondTypes.size()-1;   break;
                case SDLK_RETURN:    op_enter    (); break;
                case SDLK_BACKSPACE: op_backspace(); break;
                case SDLK_u:  qCamera.setOne();  break;
                case SDLK_k: printf("to   file \n"); builder.toFile  ( "truss.txt" ); break;
                case SDLK_l: printf("from file \n"); builder.fromFile( "truss.txt" ); break;
                case SDLK_SPACE: running = !running; if(running){ builder.toSoftBody( truss ); }; break;
                case SDLK_KP_1: builder.insertBox( ix, iy, iz, 1, &builder.bondTypes[iType] ); break;
                case SDLK_KP_2: builder.insertBox( ix, iy, iz, 2, &builder.bondTypes[iType] ); break;
                case SDLK_KP_3: builder.insertBox( ix, iy, iz, 3, &builder.bondTypes[iType] ); break;
                case SDLK_KP_4: builder.insertBox( ix, iy, iz, 4, &builder.bondTypes[iType] ); break;
                case SDLK_KP_5: builder.insertBox( ix, iy, iz, 5, &builder.bondTypes[iType] ); break;
                case SDLK_KP_6: builder.insertBox( ix, iy, iz, 6, &builder.bondTypes[iType] ); break;
                case SDLK_KP_7: builder.insertBox( ix, iy, iz, 7, &builder.bondTypes[iType] ); break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
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

TestAppTrussBuilder * testApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestAppTrussBuilder( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}
