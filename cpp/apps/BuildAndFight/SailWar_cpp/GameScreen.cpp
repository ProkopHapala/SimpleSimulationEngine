
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "drawMath2D.h"
//#include "drawDrawUtils.h"

#include "GameScreen.h" // THE HEADER

void GameScreen::draw(){
	glEnable (GL_LIGHTING);
	glShadeModel(GL_FLAT);

	world->update();

	glColor3f( 0.8f, 0.8f, 0.8f ); 	world->ship1.draw_shape( );  
	glColor3f( 0.2f, 0.2f, 0.2f );  world->ship1.draw( ); 

	glColor3f( 0.8f, 0.8f, 0.8f );  world->ship2.draw_shape( );  
	glColor3f( 0.2f, 0.2f, 0.2f );  world->ship2.draw( ); 

	for( std::vector<Projectile*>::iterator it = world->projectiles.begin(); it != world->projectiles.end(); ++it ) {
		(*it) -> draw();
	}

	Vec2d compass_pos; compass_pos.set( 0.8*ASPECT_RATIO*zoom, 0.8*zoom );

	glColor3f( 0.2f, 0.2f, 0.2f );  drawPointCross( compass_pos, zoom*0.1 );
	glColor3f( 0.2f, 0.5f, 0.2f );  drawVecInPos( ( *(Vec2d*)&(world->wind_speed)) *zoom*0.1,   compass_pos );
	glColor3f( 0.2f, 0.2f, 0.8f );  drawVecInPos( world->watter_speed*zoom*0.1, compass_pos );

	glDisable  (GL_LIGHTING);
	//drawAxis( 10 );
};

GameScreen::GameScreen( int& id, int WIDTH_, int HEIGHT_ ) : Screen2D( id, WIDTH_, HEIGHT_ ) {
//	init( id, WIDTH_, HEIGHT_ );
}

