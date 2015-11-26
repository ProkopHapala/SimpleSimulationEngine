
class Screen{
	public:
	int   WIDTH;
	int   HEIGHT;
	float VIEW_DEPTH;
	float ASPECT_RATIO;
	float zoom;

	float qCamera   [4];
	float qCameraOld[4]; 
	Vec3d camDir, camUp, camRight;

	int   mouseX, mouseY;
	bool  mouse_spinning;
	float mouse_begin_x ;
	float mouse_begin_y ;

	bool hasFocus;
	SDL_Window*      window;
	SDL_Renderer*    renderer;

// ============ Setup functions

	void initRenderer();
	void setupRenderer();

	void setDefaults(){
		VIEW_DEPTH   = VIEW_DEPTH_DEFAULT;
		ASPECT_RATIO = WIDTH/(float)HEIGHT;
		zoom = VIEW_ZOOM_DEFAULT;
		printf(" %f %f %f \n", zoom, ASPECT_RATIO, VIEW_DEPTH  );
		mouse_spinning = false;
		mouse_begin_x  = 0;
		mouse_begin_y  = 0;
	}

	Screen( int& id, int WIDTH_, int HEIGHT_ ){
		WIDTH  = WIDTH_;
		HEIGHT = HEIGHT_;
		setDefaults();
		SDL_CreateWindowAndRenderer(WIDTH, HEIGHT, SDL_WINDOW_OPENGL, &window, &renderer);
		id = SDL_GetWindowID(window); printf( " win id %i \n", id );
		char str[40];  sprintf(str, " Window id = %d", id );
		SDL_SetWindowTitle( window, str );

		setupRenderer();
		trackball (qCamera, 0, 0, 0, 0);
		printf( " ASPECT_RATIO %f \n", ASPECT_RATIO );
	}

// =========== control Handling

	void mouse_camera ( float x, float y );
	void startSpining ( float x, float y ){ mouse_spinning = true; mouse_begin_x  = x; mouse_begin_y  = y;	}
	void endSpining   ()                  { mouse_spinning = false;	}
	float mouseUp     ( float mY ){ return 2*zoom*( 0.5 -mY/float(HEIGHT)                    ); };
	float mouseRight  ( float mX ){ return 2*zoom*(      mX/float(HEIGHT) - 0.5*ASPECT_RATIO ); };
	void  projectMouse( float mX, float mY, Vec3d& mp ){	mp.set_lincomb( mouseRight(mX), camRight,  mouseUp(mY), camUp ); };

// ============ per frame functions

	void camera();

	void draw();
	void projectMouse();
	void inputHanding();
	void getCameraDirections();

	void update( ){
		//SDL_RenderPresent(renderer);
		camera();
		getCameraDirections();
		glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
		draw();
		SDL_RenderPresent(renderer);
	};

};


// ============ Setup ops

void Screen::setupRenderer(){
	float ambient  [] = { 0.1f, 0.15f, 0.25f, 1.0f };
	float diffuse  [] = { 0.9f, 0.8f,  0.7f,  1.0f };
	float specular [] = { 1.0f, 1.0f,  1.0f,  1.0f };
	float shininess[] = { 80.0f                    };
	float lightPos [] = { 1.0f, 1.0f, +1.0f, 0.0f  };
	glMaterialfv ( GL_FRONT_AND_BACK, GL_AMBIENT,   ambient);
	glMaterialfv ( GL_FRONT_AND_BACK, GL_DIFFUSE,   diffuse);
	glMaterialfv ( GL_FRONT_AND_BACK, GL_SPECULAR,  specular);
	glMaterialfv ( GL_FRONT_AND_BACK, GL_SHININESS, shininess);
	glEnable     ( GL_COLOR_MATERIAL    );
	glLightfv    ( GL_LIGHT0, GL_POSITION, lightPos);
	glLightfv    ( GL_LIGHT0, GL_AMBIENT, ambient);
	glEnable     ( GL_LIGHTING         );
	glEnable     ( GL_LIGHT0           );
	glEnable     ( GL_NORMALIZE        );
	glEnable     ( GL_DEPTH_TEST       );
	glHint       ( GL_LINE_SMOOTH_HINT, GL_NICEST );	
	glShadeModel ( GL_SMOOTH           );
	glPolygonMode( GL_FRONT_AND_BACK,GL_FILL );
}

// =================== Control Handling

void Screen::mouse_camera (float x, float y){
	if (mouse_spinning){
		trackball ( qCameraOld,
			(2.0f*mouse_begin_x-WIDTH)/WIDTH,  (HEIGHT-2.0f*mouse_begin_y)/HEIGHT, 
			(2.0f*x-WIDTH            )/WIDTH,  (HEIGHT-2.0f*y            )/HEIGHT 
		);
		add_quats ( qCameraOld, qCamera, qCamera );
		mouse_begin_x = x; mouse_begin_y = y; 
	}
}

void Screen::getCameraDirections( ){
	float mat[4][4];
	glGetFloatv (GL_MODELVIEW_MATRIX, &mat[0][0]);
	camRight.set( mat[0][0], mat[1][0], mat[2][0] );
	camUp   .set( mat[0][1], mat[1][1], mat[2][1] );
	camDir  .set( mat[0][2], mat[1][2], mat[2][2] );
	camDir.mul( -1 ); // for some reason it is inverted
}

// ================ PER FRAME OPS

void Screen::draw(){
	glEnable (GL_LIGHTING);
	glShadeModel(GL_FLAT);
	glColor3f( 0.2f, 0.8f, 0.2f );
	for( int i=0; i<glObjects->n; i++ ){
		//printf( " screen %i object %i %i \n", id, i, glObjects->data[i] );
		int obj = glObjects->data[i];  
		if( obj > 0 ){ 
			glCallList(obj); 
			//printf( " render screen %i object %i %i \n", id, i, glObjects->data[i] ); 
		};
	}
	world->draw();

	glDisable  (GL_LIGHTING);
	drawAxis( 10 );
}


void Screen::camera(){
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
	glOrtho ( -zoom*ASPECT_RATIO, zoom*ASPECT_RATIO, -zoom, zoom, -VIEW_DEPTH, +VIEW_DEPTH );
	//glOrtho ( -zoom, zoom, -zoom, zoom, -VIEW_DEPTH, +VIEW_DEPTH );
	//glOrtho ( -zoom, zoom, -zoom*ASPECT_RATIO, zoom*ASPECT_RATIO, -VIEW_DEPTH, +VIEW_DEPTH );
	glMatrixMode (GL_MODELVIEW);
	float camMatrix[4][4];
	build_rotmatrix (camMatrix, qCamera );
	glLoadMatrixf(&camMatrix[0][0]);
}





