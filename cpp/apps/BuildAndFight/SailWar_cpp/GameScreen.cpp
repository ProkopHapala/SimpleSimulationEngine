
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "drawMath2D.h"
//#include "drawDrawUtils.h"

#include "GameScreen.h" // THE HEADER

void GameScreen::draw(){

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

	Vec2d compass_pos; compass_pos.set( 0.8*ASPECT_RATIO*zoom, 0.8*zoom );

	glColor3f( 0.2f, 0.2f, 0.2f );  drawPointCross( compass_pos, zoom*0.1 );
	glColor3f( 0.2f, 0.5f, 0.2f );  drawVecInPos( ( *(Vec2d*)&(world->wind_speed)) *zoom*0.1,   compass_pos );
	glColor3f( 0.2f, 0.2f, 0.8f );  drawVecInPos( world->watter_speed*zoom*0.1, compass_pos );

	glDisable  (GL_LIGHTING);

	//exit(0);
	//drawAxis( 10 );
};

GameScreen::GameScreen( int& id, int WIDTH_, int HEIGHT_ ) : Screen2D( id, WIDTH_, HEIGHT_ ) {
//	init( id, WIDTH_, HEIGHT_ );
}

