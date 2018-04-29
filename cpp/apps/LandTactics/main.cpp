
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include <SDL2/SDL_image.h>
//#include <SDL2/SDL_ttf.h>
//#include "Texture.h"
#include "Draw.h"
#include "Draw2D.h"

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"
#include "AppSDL2OGL.h"

#include "testUtils.h"
#include "SDL_utils.h"

#include "TerrainCubic.h"
#include "TiledView.h"

#include "LTcommon.h"
#include "LTUnitType.h"
#include "LTUnit.h"
#include "LTShelter.h"
#include "LTFaction.h"
#include "LTWorld.h"

#include "GUI.h"
#include "Plot2D.h"

// font rendering:
//  http://www.willusher.io/sdl2%20tutorials/2013/12/18/lesson-6-true-type-fonts-with-sdl_ttf
//  http://stackoverflow.com/questions/28880562/rendering-text-with-sdl2-and-opengl

/*
//   Units should have 2 scales of resolution
// - Company scale (individual units are abstracted out- not rendered, not evaluated)
// - Detailed scale ( individual units are rendered and evaluated )
*/

int   default_font_texture;
char strBuf[0x10000];

void cmapHeight(double g){

    //glColor3f(0.2+0.8*g*g,0.2+0.3*(1-g)*g,0.2);
    //double snow = g*g*g*g;
    //glColor3f(g*0.8+0.2,0.5+0.5*snow,0.2+0.8*snow);

    glColor3f(g,g,g);

    //return g;
}

/*
void drawStaticObject( LTStaticObject& o ){
    glBegin(GL_LINE_LOOP);
    glVertex3f( );
    glVertex3f( );
    glVertex3f( );
    glVertex3f( );
    glEnd();
    Draw2D::drawShape( o.pos, o.rot, o.type->glo );
}
*/

void plotSiteFittness( const LTsurrounding& sur, const LTUnit& u, const Vec2d& pos ){
    double E = sur.unitPosFittness( &u, pos );
    sprintf( strBuf, "%3.3f", E );
    Draw2D::drawText(strBuf, pos, {100.0,20.0}, default_font_texture, 0.5 );
}

void plotSurrounding( const LTsurrounding& sur, const Vec2d& pos ){
    glColor3f( 0.0,0.0,1.0 ); for( LTUnit* u : sur.coleagues     ){ Draw2D::drawLine_d( pos, u->pos ); }
    glColor3f( 1.0,0.0,0.0 ); for( LTUnit* u : sur.enemies       ){ Draw2D::drawLine_d( pos, u->pos ); }
    glColor3f( 0.5,1.0,0.0 ); for( LTLinearObject* l : sur.lobjs ){ Draw2D::drawLine_d( pos, (l->p1+l->p2)*0.5 ); }
    glColor3f( 0.0,1.0,0.5 ); for( LTStaticObject* o : sur.objs  ){ Draw2D::drawLine_d( pos, o->pos ); }
}


class FormationTacticsApp : public AppSDL2OGL {
	public:
    LTWorld world;

    int formation_view_mode = 0;
    LTSquad   * currentSquad      = NULL;
    LTFaction * currentFaction   = NULL;
    int       ifaction = 0;

    bool bDrawing = false;

    double xsc,ysc;

    //GLuint       itex;

    // ==== function declaration
    void printASCItable( int imin, int imax  );
    //GLuint makeTexture( char * fname );
    //GLuint renderImage( GLuint itex, const Rect2d& rec );
    //void drawString( char * str, int imin, int imax, float x, float y, float sz, int itex );
    //void drawString( char * str, float x, float y, float sz, int itex );

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling( const SDL_Event& event  );
	void debug_buffinsert( );
	//void pickParticle( Particle2D*& picked );
	//virtual int tileToList( float x0, float y0, float x1, float y1 );

	double polyLinIntercept( int n, double* xs, double* ys, double y0, double a );
	void   drawVisibilityIsolines( Vec2d ray0, int ndhs, int nDirs, double phiMin, double phiMax, double dhmin, double dhmax, double tmax );

	FormationTacticsApp( int& id, int WIDTH_, int HEIGHT_ );

};

void FormationTacticsApp::printASCItable( int imin, int imax  ){
    int len = imax-imin;
    char str[len];
    for ( int i=0; i<len; i++ ){
        str[i] = (char)(i+imin);
    }
    printf("%s\n", str );
};

FormationTacticsApp::FormationTacticsApp( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    //SDL_MaximizeWindow( window);  // does not work if not resizable windows // https://stackoverflow.com/questions/311818/maximize-sdl-window

    /*
    for(int i=0; i<100; i++){
        //Vec3d hs = {0.0,1.0,1.0};
        Vec3d hs = {randf(-1.0,1.0),randf(-1.0,1.0),randf(-1.0,1.0)};
        Vec2d p  = {0.2, 0.2};
        double h   = trinagleInterp( toBaricentric(p), hs );
        double hdx = trinagleInterp( toBaricentric({p.x+0.1,p.y }), hs );
        double hdy = trinagleInterp( toBaricentric({p.x,p.y+0.1 }), hs );
        Vec2d  dh  = trinagleDeriv ( hs );
        printf( " %i (%g,%g) (%g,%g) \n", i, (hdx-h)/0.1, (hdy-h)/0.1, dh.x, dh.y );
    }
    exit(0);
    */

    printASCItable( 33, 127  );

    world.init();

    camX0 = world.map_center.x;
    camY0 = world.map_center.y;

    currentFaction = world.factions[0];  printf( "currentFaction: %s\n", currentFaction->name );
    currentSquad   = currentFaction->squads[0];

    //TiledView::init( 6, 6 );
    //tiles    = new int[ nxy ];
    //TiledView::renderAll( -10, -10, 10, 10 );

    //default_font_texture = makeTexture(  "common_resources/dejvu_sans_mono.bmp" );
    default_font_texture = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    //default_font_texture = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    //itex = makeTexture(  "data/tank.bmp" );
    //itex = makeTexture(  "data/nehe.bmp" );
    printf( "default_font_texture :  %i \n", default_font_texture );

}

double FormationTacticsApp::polyLinIntercept( int n, double* xs, double* ys, double y0, double a ){
    double ox = xs[0];
    double oy = ys[0];
    //for( int i=0; i<n; i++ ){ Draw2D::drawPointCross_d( { camXmin + xs[i]*xsc, camYmin + ysc*ys[i] }, 5.0 );  };
    for(int i=1; i<n; i++){
        double x  = xs[i];
        double y  = ys[i];
        double h  = (y0 + a*x ) -  y;
        glColor3f(0.0f,0.0f,1.0f); Draw2D::drawPointCross_d( { camXmin+x*xsc, camYmin + ysc*y     }, 1.0 );
        glColor3f(1.0f,0.0f,0.0f); Draw2D::drawPointCross_d( { camXmin+x*xsc, camYmin + ysc*(y+h) }, 1.0 );
        glColor3f(1.0f,0.0f,1.0f); Draw2D::drawLine( { camXmin+x*xsc, camYmin + ysc*y     },{ camXmin+x*xsc, camYmin + ysc*(y+h) } );
        //printf( "... %i x=%f %f \n", i, x, oy );
        if(h<0){
            double oh = (y0 + a*ox) - oy;
            double f  = oh/(oh-h);
            //printf( "OOO %i x=%f dx=%f %f %f oh=%f \n", i, x, x-ox, a, oy, oh );
            x  = f*(x-ox) + ox;
            //printf( "<<< %i %f %f %f dx=%f   %f %f %f \n", i, y0, h, x, (x-ox),   oh, f, x );

            glColor3f(0.0f,1.0f,0.0f); Draw2D::drawPointCross_d( { camXmin+x*xsc, camYmin + ysc*(y0 + a*x ) }, 1.0 );
            return x;
        }
        ox=x; oy=y;
    }
    return -1e-8;
}


void FormationTacticsApp::drawVisibilityIsolines( Vec2d ray0, int ndhs, int nDirs, double phiMin, double phiMax, double dhMin, double dhMax, double tmax ){
    //const int  ndhs = 5;
    //const int nDirs = 50;
    double dhs[ndhs];
    double ts [ndhs];
    Vec2d  hRays    [nDirs];
    double horizonts[nDirs*ndhs];
    // renerate slopes
    double ddhs = (dhMax-dhMin)/(ndhs-1);
    for(int i=0; i<ndhs ; i++ ){ dhs[i]=dhMin+ddhs*i; };
    // generate directions
    double dphi = (phiMax-phiMin)/(nDirs-1);
    for(int i=0; i<nDirs; i++ ){ double a = phiMin+dphi*i; hRays[i]=(Vec2d){cos(a),sin(a)}; };
    // raytrace terrain
    //double tmax = 200.0;
    for(int i=0; i<nDirs; i++ ){
        int nhit = world.ruler.rayHorizonts( ray0, hRays[i], world.ground, 10.0, ndhs, dhs, ts, tmax );
        for(int j=0; j<ndhs ; j++ ){
            int ij = j*nDirs + i;
            if( j<nhit ){ horizonts[ij]=ts[j];    }
            else        { horizonts[ij]=tmax+1.0; }
        }
    };
    // plot lines
    float cstep = 1.0f/(ndhs-1);
    for(int j=0; j<ndhs; j++ ){
        //glBegin(GL_LINE_LOOP);
        float c = cstep*j;
        glBegin(GL_LINE_STRIP);
            for(int i=0; i<nDirs; i++ ){
                int ij = j*nDirs + i;
                double t = horizonts[ij];
                Vec2d p = ray0 + hRays[i]*t;
                if( t>tmax ){ glColor3f(0.0f,0.0f,0.0f); }else{ glColor3f(1.0f-c,0.0f,c); };
                glVertex3f( (float)p.x, (float)p.y, 100.0);
            }
        glEnd();
    }
}

void FormationTacticsApp::draw(){
    //long tTot = getCPUticks();
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

	if(bDrawing){
        Vec2i ind;Vec2d dind;
        //int i = world.ruler.simplexIndex({mouse_begin_x,mouse_begin_y},ind,dind)>>1;
        world.ruler.simplexIndexBare({mouse_begin_x,mouse_begin_y},ind);
        ind = world.ruler.wrap_index(ind);
        int i = world.ruler.ip2i(ind);
        //printf("i %i (%i,%i) \n",i, ind.x, ind.y );
        world.ground[i] = randf(0.0,1.0);
	}

    //glDisable    ( GL_LIGHTING   );
    //glDisable    ( GL_DEPTH_TEST );
    glShadeModel ( GL_SMOOTH     );

	glPushMatrix();
	glScalef(world.ruler.step,world.ruler.step,1.0);
	Draw2D::drawTriaglePatch<cmapHeight>( {0,0}, {128,128}, world.ruler.na, world.ground, 0.0, world.maxHeight );
	glPopMatrix();

	//printf( "world.objects.size %i \n", world.objects.size() );
	glColor3f(0.5,0.5,0.5);
	for( LTStaticObject& o : world.objects ){
        //Draw2D::drawShape( o.pos, o.dir, o.type->glo );
        o.view();
        //o.type->render( o.pos, o.dir );
    }
    //exit(0);

    for( LTLinearObject& o : world.linObjects ){
        //printf( " (%f,%f) (%f,%f) \n", o.p1.x, o.p1.y, o.p2.x, o.p2.y );
        Draw2D::drawLine_d( o.p1, o.p2 );
    }

    /*
    LTLinearObject l1,l2;
    l1.p1 = (Vec2d){0.0,-1.0};   l1.p1.add(world.map_center);
    l1.p2 = (Vec2d){0.0, 1.0};   l1.p2.add(world.map_center);
    l2.p1 = (Vec2d){-1.0,-1.0};  l2.p1.add(world.map_center);
    l2.p2 = (Vec2d){ 1.0, 1.0};  l2.p2.add(world.map_center);

    Draw2D::drawLine_d(l1.p1, l1.p2);
    Draw2D::drawLine_d(l2.p1, l2.p2);
    Vec2d X;
    glColor3f(0.5,0.9,0.5);
    char c = l1.intersection( l2.p1, l2.p2, X );
    printf( "intersection mask %i \n", c );
    if( c==0 ){
        Draw2D::drawPointCross_d(X, 5.0 );
    }
    */

    //float camMargin = ( camXmax - camXmin )*0.1;
    //float camMargin = 0;
    //TiledView::draw(  camXmin-camMargin, camYmin-camMargin, camXmax+camMargin, camYmax+camMargin  );
    //printf( " camRect  %f %f %f %f \n", camXmin-camMargin, camYmin-camMargin, camXmax+camMargin, camYmax+camMargin );
    //long tComp = getCPUticks();
    world.update( );
    //tComp = getCPUticks() - tComp;

    long tDraw = getCPUticks();
    int i=0;
    for( LTSquad* u : world.squads ){
        if( (u!= NULL) ){
            // TODO : check if on screen
            //printf( "squad %i \n", i ); i++;
            if  ( u == currentSquad ){ u->render( u->faction->color, 1 );   }
            else                     { u->render( u->faction->color, 1 );   }
        }
    }
    tDraw = getCPUticks() - tDraw;

    if( currentSquad != 0 ){
        //glColor3f(1.0,0.0,1.0);
        glColor3f(0.0,1.0,0.0);
        Draw2D::drawCircle_d( currentSquad->pos, 0.5, 16, false );
        currentSquad->renderJob( currentSquad->faction->color );
        drawVisibilityIsolines( currentSquad->pos, 5, 50, 0, 2*M_PI, -0.1, +0.1, 500.0 );
    }

    //world.tmpSur.bConstr=false;
    world.tmpSur.bConstr=true;
    world.tmpSur.ConstrPos = currentSquad->goal;
    world.tmpSur.ConstrRad = currentSquad->goalRadius;
    world.tmpSur.ConstrE   = -1.0;

    world.tmpSur.clear();
    world.getSurroundings( world.tmpSur, currentSquad->faction, currentSquad->pos, 50.0 );
    plotSurrounding( world.tmpSur, {mouse_begin_x,mouse_begin_y} );
    plotSiteFittness( world.tmpSur, currentSquad->units[0], {mouse_begin_x,mouse_begin_y} );

    world.optimizeDeployment( currentSquad, 50.0, 5, 5, true );
    //for( LTUnit& u: currentSquad->units ){ u.pos = u.goal_pos; }


    //Vec2d ray0 = (Vec2d){mouse_begin_x,mouse_begin_y};
    //Vec2d hray = (Vec2d){0.0,1.0};   hray.normalize();
    //drawVisibilityIsolines( ray0, 5, 50, 0, 2*M_PI, -0.1, +0.1, 500.0 );

};

/*
int FormationTacticsApp::tileToList( float x0, float y0, float x1, float y1 ){
	int ilist=glGenLists(1);
	glNewList( ilist, GL_COMPILE );
		world.terrain.renderRect( x0, y0, x1, y1, 31 );
		//glColor3f(0.9f,0.2f,0.2f); Draw2D::drawRectangle( x0+0.1, y0+0.1, x1-0.1, y1-0.1, false );
	glEndList();
	return ilist;
}
*/

void FormationTacticsApp::drawHUD(){
    if(currentSquad){
        glPushMatrix();
        //printf( "currentSquad.type %i %s \n", currentSquad->type, currentSquad->type->name.c_str() );
        //Draw::drawText( "abcdefghijklmnopqrstuvwxyz \n0123456789 \nABCDEFGHIJKLMNOPQRTSTUVWXYZ \nxvfgfgdfgdfgdfgdfgdfg", fontTex, 8, {10,5} );
        glTranslatef( 10.0,HEIGHT-20,0.0  ); glColor3f(1.0,0.0,1.0); currentSquad->type->toStrCaptioned(strBuf,true); Draw::drawText( strBuf, default_font_texture, 8, {80,50} );
        //glTranslatef( 300.0,      0.0,0.0 ); glColor3f(0.0,0.5,0.0); currentFormation->reportStatus(strBuf); Draw::drawText( strBuf, fontTex, 8, {80,15} );
        //glTranslatef( 300.0,      0.0,0.0 ); glColor3f(0.0,0.5,0.8); currentFormation->soldiers[0].type->toStrCaptioned(strBuf); Draw::drawText( strBuf, fontTex, 8, {80,15} );
        glPopMatrix();
    }
};

void FormationTacticsApp::eventHandling ( const SDL_Event& event  ){
    //printf( "NBodyWorldApp::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_0:  formation_view_mode = 0;            printf( "view : default\n" ); break;
                //case SDLK_1:  formation_view_mode = VIEW_INJURY;  printf( "view : injury\n"  ); break;
                //case SDLK_2:  formation_view_mode = VIEW_STAMINA; printf( "view : stamina\n" ); break;
                //case SDLK_3:  formation_view_mode = VIEW_CHARGE;  printf( "view : charge\n"  ); break;
                //case SDLK_4:  formation_view_mode = VIEW_MORAL;   printf( "view : moral\n"   ); break;
                case   SDLK_n: ifaction++; if(ifaction>=world.factions.size()) ifaction=0; currentFaction = world.factions[ifaction]; printf("ifaction %i\n",ifaction); break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    if( currentFaction != NULL ) currentSquad = currentFaction->getUnitAt( { mouse_begin_x, mouse_begin_y } );
                    //bDrawing=true;
                break;
                case SDL_BUTTON_RIGHT:
                    //printf( "left button pressed !!!! " );
                    if( currentSquad != NULL ){
                        int imin = world.getUnitAt( { mouse_begin_x, mouse_begin_y }, currentFaction );
                        if( imin > -1 ) {
                            printf( "target selected %i %i\n", imin, world.squads[imin] );
                            currentSquad->setOpponent( world.squads[imin] );
                        }else{
                            printf( "goal selected (%3.3f,%3.3f)\n", mouse_begin_x, mouse_begin_y );
                            currentSquad->setGoal  ( { mouse_begin_x, mouse_begin_y } );
                        }
                    }
                break;
            }
            break;
            /*
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
    AppSDL2OGL::eventHandling( event );
    camStep = zoom*0.05;
}

// ===================== MAIN

FormationTacticsApp * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);

    SDL_DisplayMode dm;
    SDL_GetDesktopDisplayMode(0, &dm);
    //printf( "Display size W: %i H:%i \n", dm.w, dm.h );

	int junk;
	thisApp = new FormationTacticsApp( junk , dm.w-150, dm.h-100 );
	thisApp->zoom = 30;
	thisApp->loop( 1000000 );
	return 0;
}
















