
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "Draw3D.h"
#include "SDL_utils.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"

#include "Truss.h"
#include "SpaceCraft.h"
#include "SoftBody.h"

#include "AppSDL2OGL_3D.h"
#include "GUI.h"
#include "testUtils.h"

#include "EditSpaceCraft.h"

using namespace SpaceCrafting;

enum class EDIT_MODE : int { vertex=0, edge=1, size }; // http://www.cprogramming.com/c++11/c++11-nullptr-strongly-typed-enum-class.html

void drawTruss( Truss& truss ){
    //printf("=============\n");
    for( int i=0; i<truss.blocks.size(); i++ ){
        Draw::color_of_hash(i+15454);
        Vec2i bj,bi = truss.blocks[i];
        if( i==(truss.blocks.size()-1) ){ bj = {truss.points.size(),truss.edges.size()}; }else{ bj = truss.blocks[i+1]; };
        //Draw3D::drawPoints( bj.a-bi.a, &truss.points[bi.a], 0.1 );
        glBegin( GL_LINES );
        for( int j=bi.b; j<bj.b; j++ ){
            Vec3f a,b;
            TrussEdge& ed = truss.edges[j];
            convert( truss.points[ed.a], a ); glVertex3f( a.x, a.y, a.z );
            convert( truss.points[ed.b], b ); glVertex3f( b.x, b.y, b.z );

            //if(i==(truss.blocks.size()-1)){ printf( "%i %i (%i,%i) (%g,%g,%g) (%g,%g,%g) \n", i, j, ed.a, ed.b, a.x, a.y, a.z,  b.x, b.y, b.z ); }
        }
        glEnd();
    }
    //printf( "points.size() %i \n", truss.points.size() );
};

class SpaceCraftEditGUI : public AppSDL2OGL_3D {
	public:

	SpaceCraft craft;
    Truss      truss;

    EDIT_MODE edit_mode = EDIT_MODE::vertex;
    //EDIT_MODE edit_mode = EDIT_MODE::edge;
    int picked = -1;
    Vec3d mouse_ray0;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );
    //virtual void mouseHandling( );

	SpaceCraftEditGUI( int& id, int WIDTH_, int HEIGHT_ );

};

SpaceCraftEditGUI::SpaceCraftEditGUI( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {
    //truss.loadXYZ(  "data/octShip.xyz" );

    //Lua1.init();

    theSpaceCraft = new SpaceCraft();

    initSpaceCraftingLua();
    //luaL_dofile(theLua, "data/spaceshil1.lua");
    Lua::dofile(theLua,"data/spaceshil1.lua");
    //luaL_dostring(theLua, "print('LuaDEBUG 1'); n1 = Node( {-100.0,0.0,0.0} ); print('LuaDEBUG 2'); print(n1)");

    //exit(0);

    Truss trussPlan;
    trussPlan.loadXYZ(  "data/octShip.xyz" );

    //trussPlan.affineTransform( (Mat3d){5.5,0.0,0.0, 0.0,5.5,0.0, 0.0,0.0,5.5}, false );

    trussPlan.affineTransform( (Mat3d){6,0.0,0.0, 0.0,6,0.0, 0.0,0.0,6}, false );

    GirderParams * gpar = new GirderParams [trussPlan.edges.size()];
    //Vec3d ups  = new Vec3d[trussPlan.edges.size()];
    Vec3d ups[] = {
        (Vec3d){0.0,-1.0,0.0},
        (Vec3d){0.0,+1.0,0.0},
        (Vec3d){0.0,0.0,-1.0},
        (Vec3d){0.0,0.0,+1.0},
        (Vec3d){-1.0,0.0,0.0},
        (Vec3d){+1.0,0.0,0.0}
    };

    //truss.makeGriders( 6, &trussPlan.edges[0], &trussPlan.points[0], gpar, ups );
    //truss.makeGriders( trussPlan, gpar, ups, NULL );
    std::vector<Vec2i> ends;
    truss.makeGriders( trussPlan, gpar, ups, &ends );
    truss.autoBridge(ends.size(), &ends[0], 0.8, 0 );

    delete [] gpar;
    //delete [] ups;

    //truss.girder1( (Vec3d){-5.0,0.0,0.0}, (Vec3d){5.0,0.0,0.0}, (Vec3d){0.0,1.0,0.0}, 5, 1.0 );
}

void SpaceCraftEditGUI::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glDisable(GL_DEPTH_TEST);

    drawTruss( truss );

    //printf( "%i\n", EDIT_MODE::vertex );
    if(picked>=0){
        switch(edit_mode){
            case EDIT_MODE::vertex:
                glColor3f(1.0,1.0,1.0); Draw3D::drawPointCross( truss.points[picked], 0.3 );
                if ( SDL_GetMouseState(NULL, NULL) & SDL_BUTTON(SDL_BUTTON_LEFT) ){ Draw3D::drawLine(truss.points[picked], mouse_ray0); }
                break;
            case EDIT_MODE::edge  : glColor3f(1.0,1.0,1.0); auto ed = truss.edges[picked]; Draw3D::drawLine( truss.points[ed.a], truss.points[ed.b] ); break;
        }

    }




    mouse_ray0 = camMat.a*mouse_begin_x + camMat.b*mouse_begin_y;

    //glColor3f(0.0f,0.0f,0.0f); drawTruss( truss.edges.size(), &truss.edges[0], &truss.points[0] );
    //glColor3f(1.0f,1.0f,1.0f); Draw3D::drawPoints( truss.points.size(), &truss.points[0], 0.1 );

};

void SpaceCraftEditGUI::drawHUD(){
    glDisable ( GL_LIGHTING );
    //glColor3f(1.0f,1.0f,1.0f);   txtStatic.view3D( {5,5}, fontTex, 8 );
    glPopMatrix();
}

//void SpaceCraftEditGUI::keyStateHandling( const Uint8 *keys ){ };
/*
void SpaceCraftEditGUI::mouseHandling( ){
    SDL_GetMouseState( &mouseX, &mouseY ); mouseY=HEIGHT-mouseY;
    mouse_t   = (mouseX-20 )/timeScale + tstart;
    mouse_val = (mouseY-200)/valScale;
    sprintf( curCaption, "%f %f\0", mouse_t, mouse_val );
    int ipoint_ = binSearchFrom<double>(mouse_t,splines.n,splines.ts);
    if( (splines.ts[ipoint_+1]-mouse_t)<(mouse_t-splines.ts[ipoint_]) ) ipoint_++;
    if(ipoint_!=ipoint){
        ipoint=ipoint_;
        char buff[100];
        Vec3d r,v;
        r.set( splines.CPs[0][ipoint],splines.CPs[1][ipoint],splines.CPs[2][ipoint] );
        v.set( splines.getPointDeriv(ipoint,0), splines.getPointDeriv(ipoint,1), splines.getPointDeriv(ipoint,2) );
        sprintf(buff, "%i %f r(%3.3f,%3.3f,%3.3f) v(%3.3f,%3.3f,%3.3f)", ipoint, splines.ts[ipoint], r.x,r.y,r.z,   v.x,v.y,v.z );
        txtStatic.inputText = buff;
    }
    //printf("%i %i %3.3f %3.3f \n", mouseX,mouseY, mouse_t, mouse_val );
}
*/

void SpaceCraftEditGUI::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_m:  edit_mode = (EDIT_MODE)((((int)edit_mode)+1)%((int)EDIT_MODE::size)); printf("edit_mode %i\n", (int)edit_mode); break;
                //case SDLK_h:  warrior1->tryJump(); break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    switch(edit_mode){
                        case EDIT_MODE::vertex: picked = truss.pickVertex( mouse_ray0, camMat.c, 0.5  ); printf("picked %i\n", picked); break;
                        case EDIT_MODE::edge  : picked = truss.pickEdge  ( mouse_ray0, camMat.c, 0.25 ); printf("picked %i\n", picked); break;
                    }; break;
                case SDL_BUTTON_RIGHT: break;
            }
            break;

        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    switch(edit_mode){
                        case EDIT_MODE::vertex: int ip2 = truss.pickVertex( mouse_ray0, camMat.c, 0.5  ); if((picked>=0)&(ip2!=picked)); truss.edges.push_back((TrussEdge){picked,ip2,0}); break;
                        //case EDIT_MODE::edge  : picked = truss.pickEdge  ( mouse_ray0, camMat.c, 0.25 ); printf("picked %i\n", picked); break;
                    }; break;
                case SDL_BUTTON_RIGHT:break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}

// ===================== MAIN

SpaceCraftEditGUI * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new SpaceCraftEditGUI( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















