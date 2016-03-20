
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "Draw2D.h"
#include "Voronoi.h"

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

    // ============== plot Voronoi

/*
    glColor3f( 0.9f, 0.9f, 0.9f );
    int iii=0;
    for ( Vec2d * p : *(world->voronoi->places) ){
        iii++;
        Draw2D::drawPointCross_d( *p, 0.1 );
	}
	//printf( " Vertices N %i \n", iii );

    glColor3f( 0.9f, 0.9f, 0.9f );
    iii=0;
    for ( VoronoiNamespace::VEdge * pe : *(world->voronoi->edges) ){
        iii++;
        Vec2d a,b;
        a.set( *(pe->start) );
        b.set( *(pe->end)   );
        Draw2D::drawLine_d( a, b );
        //printf( " %f %f   %f %f \n", a.x, a.y,    b.x, b.y  );
	}
	//exit(0);
	//printf( " Edges N %i \n", iii );

*/

	// ============== rest

	glDisable  (GL_LIGHTING);

	//exit(0);
	//drawAxis( 10 );
};


void GameScreen::drawHUD(){
    float bar_y =  10.0f;
    float bar_x = 100.0f;
    if( thisShip != NULL ){
        Draw2D::drawRectangle( { 0.0f,  0.0f }, { bar_x * thisShip->gunload_left,    bar_y }, true );
        Draw2D::drawRectangle( { 0.0f, bar_y }, { bar_x * thisShip->gunload_right, 2*bar_y }, true );
        Draw2D::drawLine     ( { bar_x, 0.0f }, { bar_x, 2*bar_y } );
    }

    Vec2d compass_pos; compass_pos.set( WIDTH*0.9, HEIGHT*0.9 );
	glColor3f( 0.2f, 0.2f, 0.2f );  Draw2D::drawPointCross_d( compass_pos, 10 );
	glColor3f( 0.2f, 0.5f, 0.2f );  Draw2D::drawVecInPos_d  ( ( *(Vec2d*)&(world->wind_speed))*10.1,   compass_pos );
	glColor3f( 0.2f, 0.2f, 0.8f );  Draw2D::drawVecInPos_d  ( world->watter_speed*10.1, compass_pos );
};




GameScreen::GameScreen( int& id, int WIDTH_, int HEIGHT_ ) : ScreenSDL2OGL( id, WIDTH_, HEIGHT_ ) {
//	init( id, WIDTH_, HEIGHT_ );
}

