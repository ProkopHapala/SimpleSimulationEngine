
#include "ScreenSDL2OGL.h" // THE HEADER

// ============== per frame

void ScreenSDL2OGL::update( ){
	//SDL_RenderPresent(renderer);
	//glPushMatrix();
	camera();
	draw();
	cameraHUD();
	drawHUD();
	//glPopMatrix();
	SDL_RenderPresent(renderer);
};

void ScreenSDL2OGL::draw   (){
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
};

void ScreenSDL2OGL::drawHUD(){ };

void ScreenSDL2OGL::inputHanding(){};

void ScreenSDL2OGL::camera(){
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
	glOrtho ( -zoom*ASPECT_RATIO, zoom*ASPECT_RATIO, -zoom, zoom, -VIEW_DEPTH, +VIEW_DEPTH );
	glMatrixMode (GL_MODELVIEW);
}

void ScreenSDL2OGL::cameraHUD(){
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
	glOrtho ( -0, WIDTH, 0, HEIGHT, -VIEW_DEPTH, +VIEW_DEPTH );
	glMatrixMode (GL_MODELVIEW);
}

// ============== initialization

void ScreenSDL2OGL::setupRenderer(){
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
	glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
}

void ScreenSDL2OGL::setDefaults(){
	VIEW_DEPTH   = VIEW_DEPTH_DEFAULT;
	ASPECT_RATIO = WIDTH/(float)HEIGHT;
	zoom         = VIEW_ZOOM_DEFAULT;
	//printf(" %f %f %f \n", zoom, ASPECT_RATIO, VIEW_DEPTH  );
	mouse_begin_x  = 0;
	mouse_begin_y  = 0;
}

void ScreenSDL2OGL::init( int& id, int WIDTH_, int HEIGHT_ ){
	WIDTH  = WIDTH_;
	HEIGHT = HEIGHT_;
	setDefaults();
	SDL_CreateWindowAndRenderer(WIDTH, HEIGHT, SDL_WINDOW_OPENGL, &window, &renderer);
	id = SDL_GetWindowID(window); printf( " win id %i \n", id );
	char str[40];  sprintf(str, " Window id = %d", id );
	SDL_SetWindowTitle( window, str );
	setupRenderer();
	//printf( " ASPECT_RATIO %f \n", ASPECT_RATIO );
}

ScreenSDL2OGL::ScreenSDL2OGL( int& id, int WIDTH_, int HEIGHT_ ){
	init( id, WIDTH_, HEIGHT_ );
};

