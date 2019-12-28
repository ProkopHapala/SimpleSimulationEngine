
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
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
//#include "raytrace.h"


#include "geom2D.h"


//#include "SphereSampling.h"
//#include "DrawSphereMap.h"
#include "Draw3D_Surf.h"

#include "AppSDL2OGL_3D.h"
#include "GUI.h"
#include "IO_utils.h"
#include "testUtils.h"

//#include "Table.h"
//#include "Tree.h"

#include "CAD2D.h"

using namespace CAD;


int fontTex = 0 ;

class CAD2DGUI : public AppSDL2OGL_3D { public:

	//DropDownList lstLuaFiles;
    GUI gui;

    //DropDownList* lstLuaFiles=0;
    //OnSelectLuaShipScript onSelectLuaShipScript;

    double gridStep = 1.0;

    int picked = -1;

    void selectCompGui();

    // https://stackoverflow.com/questions/29145476/requiring-virtual-function-overrides-to-use-override-keyword

	virtual void draw   () override;
	virtual void drawHUD() override;
	//virtual void mouseHandling( );
	//virtual void camera();
	virtual void eventHandling   ( const SDL_Event& event  ) override;
	virtual void keyStateHandling( const Uint8 *keys ) override;
    //virtual void mouseHandling( );

	CAD2DGUI( int& id, int WIDTH_, int HEIGHT_ );

};

CAD2DGUI::CAD2DGUI( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex     = makeTexture    ( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    GUI_fontTex = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );

}

void drawPointGrid(const Vec2i& ns,const Vec2d& p0,const Vec2d& da,const Vec2d& db ){
    glBegin(GL_POINTS);
    for(int ia=0;ia<ns.a;ia++){
        for(int ib=0;ib<ns.b;ib++){
            Vec2d p = p0 + da*ia + db*ib;
            glVertex2f(p.x,p.y);
        }
    }
    glEnd();
}


void CAD2DGUI::draw(){
    //printf( " ==== frame %i \n", frameCount );
    //glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    //glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
    //glClearColor( 0.0f, 0.0f, 0.0f, 1.0f );
    glClearColor( 0.8f, 0.8f, 0.8f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	gridStep = 0.01*pow( 10., round( log10( zoom ) ) );

	double mx = gridStep * round( mouse_begin_x/gridStep );
	double my = gridStep * round( mouse_begin_y/gridStep );

    Vec2d p0 = { -ASPECT_RATIO*zoom, -zoom };
    Vec2d p1 = { +ASPECT_RATIO*zoom,  zoom };

	glColor3f(0.75,0.75,0.75);
	Draw2D::drawLine( {mx,my-1000}, {mx,my+1000} );
	Draw2D::drawLine( {mx-1000,my}, {mx+1000,my} );

	//glColor3f(0.25,0.25,0.25);

    int nx= (int)(p1.x-p0.x)/gridStep;
    int ny= (int)(p1.y-p0.y)/gridStep;
    //printf( "%i %i \n", nx, ny );
	//drawPointGrid( {nx, ny}, {mx,my}, {gridStep,0.}, {0.,gridStep} );
	drawPointGrid( {nx, ny}, p0+((Vec3d)cam.pos).xy(), {gridStep,0.}, {0.,gridStep} );

	/*
    {
        Vec2d p1 = { 0.2, .8};
        Vec2d p2 = {-0.9, .5};
        Vec2d p3 = { 1.2,-.8};

        Circle2d c1;
        c1.from3points( p1,p2,p3 );
        //c1.from2points( p1,p2 );
        //c1.fromCenterAndPoint( p1,p2 );

        glColor3f(0.,0.,0.);
        Draw2D::drawPointCross_d( p1, 0.1 );
        Draw2D::drawPointCross_d( p2, 0.1 );
        Draw2D::drawPointCross_d( p3, 0.1 );
        Draw2D::drawCircle_d( c1.p0, c1.r, 64, false );
    }
    */

    /*
    {
        //Vec2d p1 = { 0.2, .8};
        //Vec2d p3 = { mouse_begin_x, mouse_begin_y };
        Vec2d p1 = { 0.2, .8};
        Vec2d p2 = {-0.9, .5};
        Vec2d p3 = { mouse_begin_x, mouse_begin_y };

        Circle2d c1;
        Arc2d    arc;

        c1.from3points( p1,p2,p3 );
        //arc.fromCenter2points( c1.p0, p1, p2 );
        arc.fromCenterPoints( c1.p0, p1, p2, p3 );

        Vec2d v;
        glColor3f( 1.0,0.0,0.0 ); v.fromAngle( arc.angs[0] ); Draw2D::drawVecInPos_d( v*c1.r, c1.p0 );
        glColor3f( 0.0,0.0,1.0 ); v.fromAngle( arc.angs[1] ); Draw2D::drawVecInPos_d( v*c1.r, c1.p0 );

        glColor3f( 0.7,0.7,0.7 ); Draw2D::drawPointCross_d( p1, 0.1 );
        glColor3f( 0.7,0.7,0.7 ); Draw2D::drawPointCross_d( p2, 0.1 );
        glColor3f( 0.7,0.0,0.7 ); Draw2D::drawPointCross_d( p3, 0.1 );

        glColor3f( 0.7,0.7,0.7 ); Draw2D::drawCircle_d( c1.p0, c1.r,  64, false );
        glColor3f( 0.0,0.0,0.0 ); Draw2D::drawArc_d   ( c1.p0, c1.r, arc.angs[0], arc.angs[1], 0.1, false );

    }
    */

    {
        Vec2d pc = { 0.2, .8};
        Vec2d d1 = {-0.9, .5};
        Vec2d d2 = { mouse_begin_x, mouse_begin_y };
        d2.sub(pc);
        d1.normalize();
        d2.normalize();
        double r = 0.4;

        Circle2d c1;
        Arc2d    arc;

        c1.fromCorner ( pc, d1, d2, r );
        arc.fromCorner( d1, d2 );

        glColor3f( 0.7,0.7,0.7 ); Draw2D::drawPointCross_d( pc, 0.1 );
        glColor3f( 1.0,0.0,0.0 ); Draw2D::drawVecInPos_d  ( d1*5, pc );
        glColor3f( 0.0,0.0,1.0 ); Draw2D::drawVecInPos_d  ( d2*5, pc );


        Vec2d d = (d1+d2);
        glColor3f( 0.0,1.0,0.0 ); Draw2D::drawVecInPos_d( d*( r/d.norm() ), pc );

        glColor3f( 0.7,0.7,0.7 ); Draw2D::drawCircle_d( c1.p0, c1.r,  64, false );
        glColor3f( 0.0,0.0,0.0 ); Draw2D::drawArc_d   ( c1.p0, c1.r, arc.angs[0], arc.angs[1], 0.1, false );

    }


};

void CAD2DGUI::drawHUD(){
    glDisable( GL_LIGHTING );
    glDisable(GL_DEPTH_TEST);

    gui.draw();
    //glColor3f(1.0f,1.0f,1.0f);   txtStatic.view3D( {5,5}, fontTex, 8 );
    //glPopMatrix();
}

//void CAD2DGUI::keyStateHandling( const Uint8 *keys ){ };
/*
void CAD2DGUI::mouseHandling( ){
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



void CAD2DGUI::keyStateHandling( const Uint8 *keys ){

    //Mat3d camMat;
    //qCamera.toMatrix_T(camMat);
    //Draw3D::drawMatInPos(  (Mat3f)camMat, (Vec3f){0.0,0.0,0.0} );
    //qCamera.toMatrix(camMat);

    if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.dyaw  (  keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.dyaw  ( -keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.dpitch(  keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.dpitch( -keyRotSpeed ); }

    if( keys[ SDL_SCANCODE_W ] ){ cam.pos.add_mul( cam.rot.b, +0.05*zoom ); }
	if( keys[ SDL_SCANCODE_S ] ){ cam.pos.add_mul( cam.rot.b, -0.05*zoom );  }
	if( keys[ SDL_SCANCODE_A ] ){ cam.pos.add_mul( cam.rot.a, -0.05*zoom );  }
	if( keys[ SDL_SCANCODE_D ] ){ cam.pos.add_mul( cam.rot.a, +0.05*zoom );  }

};

void CAD2DGUI::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    //if(event.type == SDL_MOUSEWHEEL){
    //    if     (event.wheel.y > 0){ zoom*=VIEW_ZOOM_STEP; }
    //    else if(event.wheel.y < 0){ zoom/=VIEW_ZOOM_STEP; }
    //}

    if(event.type == SDL_MOUSEWHEEL){
        if     (event.wheel.y > 0){ zoom/=VIEW_ZOOM_STEP; }
        else if(event.wheel.y < 0){ zoom*=VIEW_ZOOM_STEP; }
    }
    gui.onEvent(mouseX,mouseY,event);
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_m:  edit_mode = (EDIT_MODE)((((int)edit_mode)+1)%((int)EDIT_MODE::size)); printf("edit_mode %i\n", (int)edit_mode); break;
                //case SDLK_h:  warrior1->tryJump(); break;
                //case SDLK_l:    onSelectLuaShipScript.GUIcallback(lstLuaFiles); break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //switch(edit_mode){
                    //    case EDIT_MODE::vertex    : picked = truss.pickVertex( mouse_ray0, (Vec3d)cam.rot.c, 0.5  ); printf("picked %i\n", picked); break;
                    //    case EDIT_MODE::edge      : picked = truss.pickEdge  ( mouse_ray0, (Vec3d)cam.rot.c, 0.25 ); printf("picked %i\n", picked); break;
                    //    case EDIT_MODE::component :
                    //        picked = theSpaceCraft->pick( mouse_ray0, (Vec3d)cam.rot.c, 3.0 );
                    //        //printf("picked %i %i\n", picked, theSpaceCraft->pickedTyp );
                    //        selectCompGui();
                    //        //printf("onEvent %li %li \n", girderGui, compGui );
                    //        if( compGui ){
                    //            if(picked>=0){ compGui->bindLoad( theSpaceCraft->getPicked(picked) ); }
                    //            //else         { compGui->unbind( ); };
                    //        }
                    //        break;
                    //}; break;
                case SDL_BUTTON_RIGHT: break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //switch(edit_mode){
                        //case EDIT_MODE::vertex: int ip2 = truss.pickVertex( mouse_ray0, (Vec3d)cam.rot.c, 0.5  ); if((picked>=0)&(ip2!=picked)); truss.edges.push_back((TrussEdge){picked,ip2,0}); break;
                        //case EDIT_MODE::edge  : picked = truss.pickEdge  ( mouse_ray0, camMat.c, 0.25 ); printf("picked %i\n", picked); break;
                    //}; break;
                case SDL_BUTTON_RIGHT:break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );

    //printf( "compGui %li \n", compGui );
    //if(compGui) if( compGui->check() ){ renderShip(); }

}

// ===================== MAIN

CAD2DGUI * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);


	// https://www.opengl.org/discussion_boards/showthread.php/163904-MultiSampling-in-SDL
	//https://wiki.libsdl.org/SDL_GLattr
	//SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 8);
    //SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 8);
    //SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 8);
    //SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE, 8);
    //SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    //SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    //SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 1);
    //SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 2);
    //SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 4);
    //SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 8);
    //SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 16);
    //glEnable(GL_MULTISAMPLE);


	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
    SDL_DisplayMode dm;
    SDL_GetDesktopDisplayMode(0, &dm);
	thisApp = new CAD2DGUI( junk , dm.w-150, dm.h-100 );
	//thisApp = new CAD2DGUI( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















