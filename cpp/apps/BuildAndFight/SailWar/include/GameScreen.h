
class GameScreen : public Screen2D {
	public:

	virtual void draw(){
		glEnable (GL_LIGHTING);
		glShadeModel(GL_FLAT);

		for( int i=0; i<perFrame; i++ ){
			ship1.clean_temp( );
			ship1.applySailForces(  windSpeed,  watterSpeed );
			ship1.move( dt );

			ship2.clean_temp( );
			ship2.applySailForces(  windSpeed,  watterSpeed );
			ship2.move( dt );

			std::vector<Projectile*>::iterator it = projectiles.begin();
			while( it != projectiles.end() ) {
				Projectile * p = *it; 
				p -> evalForce(    );
				p -> move     ( dt );
				if( p -> check_hit( ) ){ it = projectiles.erase( it ); }
				else                   { ++it;                  }
			}

		}


		glColor3f( 0.8f, 0.8f, 0.8f ); 	ship1.draw_shape( );  
		glColor3f( 0.2f, 0.2f, 0.2f );  ship1.draw( ); 

		glColor3f( 0.8f, 0.8f, 0.8f ); 	ship2.draw_shape( );  
		glColor3f( 0.2f, 0.2f, 0.2f );  ship2.draw( ); 


		for( std::vector<Projectile*>::iterator it = projectiles.begin(); it != projectiles.end(); ++it ) {
			(*it) -> draw();
		}


		Vec2d compass_pos; compass_pos.set( 0.8*ASPECT_RATIO*zoom, 0.8*zoom );

		glColor3f( 0.2f, 0.2f, 0.2f );  drawPointCross( compass_pos, zoom*0.1 );
		glColor3f( 0.2f, 0.5f, 0.2f );  drawVecInPos( windSpeed*zoom*0.1,   compass_pos );
		glColor3f( 0.2f, 0.2f, 0.8f );  drawVecInPos( watterSpeed*zoom*0.1, compass_pos );

		glDisable  (GL_LIGHTING);
		drawAxis( 10 );
	}

/*
	virtual void update( ){
		camera();
		glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
		draw();
		glPopMatrix();
		SDL_RenderPresent(renderer);
	};
*/

	GameScreen( int& id, int WIDTH_, int HEIGHT_ ) : Screen2D( id, WIDTH_, HEIGHT_ ) {
	//	init( id, WIDTH_, HEIGHT_ );
	}

/*
	virtual void projectMouse       (){ Screen2D::projectMouse();        };
	virtual void inputHanding       (){ Screen2D::inputHanding();        };
	virtual void getCameraDirections(){ Screen2D::getCameraDirections(); };
	virtual void update             (){	Screen2D::update();              };
*/





};





