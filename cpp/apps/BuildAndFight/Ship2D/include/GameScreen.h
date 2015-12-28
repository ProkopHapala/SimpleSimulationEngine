class GameScreen : public Screen2D {
	public:

	virtual void draw(){
		glEnable (GL_LIGHTING);
		glShadeModel(GL_FLAT);

		for( int i=0; i<perFrame; i++ ){
			yacht1.clean_temp( );
			yacht1.applySailForces(  windSpeed,  watterSpeed );
			yacht1.move( dt );
		}

		chain1->move( yacht1.pos );
		glColor3f( 0.3f, 0.3f, 0.3f ); chain1->draw( yacht1.pos );

		glColor3f( 0.8f, 0.8f, 0.8f ); 	yacht1.draw_shape( );  
		glColor3f( 0.2f, 0.2f, 0.2f );  yacht1.draw( ); 

		Vec2d compass_pos; compass_pos.set( 0.8*ASPECT_RATIO*zoom, 0.8*zoom );

		glColor3f( 0.2f, 0.2f, 0.2f );  drawPointCross( compass_pos, zoom*0.1 );
		glColor3f( 0.2f, 0.5f, 0.2f );  drawVecInPos( windSpeed*zoom*0.1,   compass_pos );
		glColor3f( 0.2f, 0.2f, 0.8f );  drawVecInPos( watterSpeed*zoom*0.1, compass_pos );

		glDisable  (GL_LIGHTING);
		drawAxis( 10 );
	}

	GameScreen( int& id, int WIDTH_, int HEIGHT_ ) : Screen2D( id, WIDTH_, HEIGHT_ )  {	};

};
