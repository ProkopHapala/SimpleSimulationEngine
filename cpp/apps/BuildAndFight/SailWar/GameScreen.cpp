
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>



#include "drawMath2D.h"
//#include "drawDrawUtils.h"

#include "GameScreen.h" // THE HEADER

void GameScreen::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    //printf( " GameScreen::draw \n" );

	glEnable (GL_LIGHTING);
	glShadeModel(GL_FLAT);

	world->update();

    //printf( " world->update DONE \n" );

	//for( std::vector<Projectile*>::iterator it = world->projectiles.begin(); it != world->projectiles.end(); ++it ) {
	for( auto ship : world->ships ) {
        //printf( "projectile draw \n " );
        glColor3f( 0.8f, 0.8f, 0.8f ); 	ship->drawHitBox( );
        glColor3f( 0.8f, 0.8f, 0.8f ); 	ship->draw_shape( );
        glColor3f( 0.2f, 0.2f, 0.2f );  ship->draw( );
	}

	for( auto p : world->projectiles ) {
        //printf( "projectile draw \n " );
		p->draw();
		//p->update_old_pos();
	}


/*

    //printf( " mouse x,y : %i %i %f %f \n", mouseX, mouseY, mouseRight(mouseX), mouseUp(mouseY) );
    Vec2d mvec;
    mvec.set( mouseRight(mouseX), mouseUp(mouseY) );
    //mvec.set( 0.0, 0.0 );
    Draw2D::drawPointCross_d( mvec, 1.0 );

	for( auto isle : world->isles ) {
        //printf( "draw isle: %i %f %f \n ", isle->n, isle->corners[0].x, isle->corners[0].y );
        if( isle->pointIn( mvec ) ){
            glColor3f( 1.0f, 0.2f, 0.2f );
            //printf( " mouse in isle \n " );
        }else{
            glColor3f( 0.2f, 0.2f, 0.2f );
        }
		Draw2D::drawConvexPolygon( isle->n, isle->corners, true );

		//exit(0);
	}

*/

    Convex2d base( 5 );
    /*
	base.corners[0].set( -1.0, +1.0 );
    base.corners[1].set( +1.0, -1.5 );
    base.corners[2].set( +2.0,  0.5 );
    */
	base.corners[0].set( -1.0, -1.0 );
    base.corners[1].set( +1.0, -1.0 );
    base.corners[2].set( +2.0,  0.0 );
	base.corners[3].set( +1.0, +1.0 );
	base.corners[4].set( -1.0, +1.0 );

    glColor3f( 0.9f, 0.9f, 0.9f ); Draw2D::drawConvexPolygon( base.n, base.corners, false );


    Vec2d Acut,Bcut;
    Acut.set( -1.0, -5.0 );
    Bcut.set( +3.0, +5.0 );
    Line2d cutline; cutline.set( Acut, Bcut );
    glColor3f( 0.2f, 0.9f, 0.2f ); Draw2D::drawLine_d( Acut, Bcut );

/*
    Vec2d p1,p2;
    cutline.intersectionPoint(      base.corners[0], base.corners[1], p1 );
    //intersection_point( Acut, Bcut, base.corners[0], base.corners[1], p2 );
    glColor3f( 0.9f, 0.2f, 0.2f ); Draw2D::drawLine_d( base.corners[0], base.corners[1] );
    glColor3f( 0.9f, 0.9f, 0.2f ); Draw2D::drawPointCross_d( p1, 0.3 );
    //glColor3f( 0.8f, 0.8f, 0.2f ); Draw2D::drawPointCross_d( p2, 0.3 );
*/


    Convex2d left;
    Convex2d right;

    //Convex2d * pLeft  = new Convex2d();
    //Convex2d * pRight = new Convex2d();

    //base.funcking_empty1( cutline, right );
    //base.funcking_empty2( cutline, left, right );
    //base.funcking_empty3( cutline, &left, &right );
    //base.funcking_empty3( cutline, pLeft, pRight );
    base.cut( cutline, left, right );
    glColor3f( 0.9f, 0.2f, 0.9f ); Draw2D::drawPointCross_d( left.corners[0], 0.3 );
    glColor3f( 0.2f, 0.9f, 0.9f ); Draw2D::drawPointCross_d( left.corners[1], 0.3 );
    glColor3f( 0.9f, 0.2f, 0.2f ); Draw2D::drawConvexPolygon( left.n,  left.corners,  true );
    glColor3f( 0.2f, 0.2f, 0.9f ); Draw2D::drawConvexPolygon( right.n, right.corners, true );

	Vec2d compass_pos; compass_pos.set( 0.8*ASPECT_RATIO*zoom, 0.8*zoom );

	glColor3f( 0.2f, 0.2f, 0.2f );  Draw2D::drawPointCross_d( compass_pos, zoom*0.1 );
	glColor3f( 0.2f, 0.5f, 0.2f );  Draw2D::drawVecInPos_d( ( *(Vec2d*)&(world->wind_speed)) *zoom*0.1,   compass_pos );
	glColor3f( 0.2f, 0.2f, 0.8f );  Draw2D::drawVecInPos_d( world->watter_speed*zoom*0.1, compass_pos );

	glDisable  (GL_LIGHTING);

	//exit(0);
	//drawAxis( 10 );
};


void GameScreen::drawHUD(){
    float bar_y =  10.0f;
    float bar_x = 100.0f;
    if( thisShip != NULL ){
        Draw2D::drawRectangle( { 0.0f,  0.0f }, { bar_x * thisShip->gunload_left,    bar_y } );
        Draw2D::drawRectangle( { 0.0f, bar_y }, { bar_x * thisShip->gunload_right, 2*bar_y } );
        Draw2D::drawLine     ( { bar_x, 0.0f }, { bar_x, 2*bar_y } );
    }
};




GameScreen::GameScreen( int& id, int WIDTH_, int HEIGHT_ ) : ScreenSDL2OGL( id, WIDTH_, HEIGHT_ ) {
//	init( id, WIDTH_, HEIGHT_ );
}

