
#include "ScreenSDL2OGL3.h" // THE HEADER

// ============== per frame

/*
void ScreenSDL2OGL3::update( ){
	draw( );
	if( deffered ) drawDeffered( );
	SDL_GL_SwapWindow(window);
	//SDL_RenderPresent(renderer);
};
*/

void ScreenSDL2OGL3::draw(){
    //printf( "DEBUG ScreenSDL2OGL3::draw\n" );
    SDL_GL_MakeCurrent(window,context);
    //glClearColor(1.0, 1.0, 1.0, 1.0);
    //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );
    for(SceneOGL3* scene: scenes){
        if( scene){
            //printf( "DEBUG ScreenSDL2OGL3[%i]::draw\n", id );
            scene->draw( cam );
        }
    }
    //if( deffered ){ glBindFramebuffer(GL_FRAMEBUFFER, FramebufferName); }
    //else          { glBindFramebuffer(GL_FRAMEBUFFER, 0);               }
    //glBindFramebuffer(GL_FRAMEBUFFER, 0);
    //glEnable   ( GL_DEPTH_TEST );
    //glDepthFunc( GL_LESS );
    //scene->render( camPos, camRot, tgFrustrum, 1.0 );
    SDL_GL_SwapWindow(window);
}

void ScreenSDL2OGL3::drawDeffered( ){

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST_MIPMAP_LINEAR);
    // Generate mipmaps, by the way.
    glGenerateMipmap(GL_TEXTURE_2D);

    GLuint uloc;

    glUseProgram(canvasShader->shaderprogram);
    uloc = glGetUniformLocation( canvasShader->shaderprogram, "canvas_depth_tex");   glUniform1i(uloc, 0);
    uloc = glGetUniformLocation( canvasShader->shaderprogram, "canvas_color_tex"); glUniform1i(uloc, 1);

    glActiveTexture(GL_TEXTURE0 );    glBindTexture  (GL_TEXTURE_2D, canvas_depth_tex );
    glActiveTexture(GL_TEXTURE1 );    glBindTexture  (GL_TEXTURE_2D, canvas_color_tex );

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    //glViewport(0,0,WIDTH,HEIGHT);                    // instead of drawing canvasQuad ?
    glEnableVertexAttribArray(0); canvasQuad->draw();  // To Do - canvasQuad can render just aread which changed

}

void ScreenSDL2OGL3::setupDefferedRender( ){

    GLuint uloc;

    //GLfloat vertexes[4][2] = {
    float * quadVerts = new float [4*2] {
	  -0.9f,  -0.9f,
	  -0.9f,   0.9f,
	   0.9f,  -0.9f,
	   0.9f,   0.9f
    };

    canvasQuad = new GLObject( );
    /*
	canvasQuad->nVert    = 4;
	canvasQuad->vertexes = quadVerts;
	canvasQuad->init();
	*/

    canvasShader=new Shader();
	canvasShader->init( "shaders/plain_vert.c", "shaders/SSAO_frag.c" );

    GLfloat resolution[2]; resolution[0]=WIDTH; resolution[1]=HEIGHT;

	glUseProgram(canvasShader->shaderprogram);

    uloc = glGetUniformLocation( canvasShader->shaderprogram, "resolution" );	glUniform2fv(uloc, 1, resolution  );

    // ------------- texture

    glGenTextures(1, &canvas_color_tex );
    glBindTexture(GL_TEXTURE_2D, canvas_color_tex );
    glTexImage2D (GL_TEXTURE_2D, 0,GL_RGB,             WIDTH, HEIGHT, 0, GL_RGB,               GL_UNSIGNED_BYTE, 0);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    glGenTextures(1, &canvas_depth_tex);
    glBindTexture(GL_TEXTURE_2D, canvas_depth_tex);
    glTexImage2D(GL_TEXTURE_2D, 0,GL_DEPTH_COMPONENT, WIDTH, HEIGHT, 0, GL_DEPTH_COMPONENT, GL_FLOAT,         0);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    // ------------- frameBuffer

    glGenFramebuffers(1, &FramebufferName);
    glBindFramebuffer(GL_FRAMEBUFFER, FramebufferName);
    // The depth buffer

    glGenRenderbuffers(1, &depthrenderbuffer);
    glBindRenderbuffer(GL_RENDERBUFFER, depthrenderbuffer);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, WIDTH, HEIGHT );
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthrenderbuffer);

    // Set "renderedTexture" as our colour attachement #0

    glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,  canvas_depth_tex,   0 );
    glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, canvas_color_tex, 0 );
    GLenum DrawBuffers[2] = {GL_DEPTH_ATTACHMENT, GL_COLOR_ATTACHMENT0};
    glDrawBuffers(2, DrawBuffers);

/*
    glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,  canvas_depth_tex,   0 );
    GLenum DrawBuffers[1] = {GL_DEPTH_ATTACHMENT};
    glDrawBuffers(1, DrawBuffers);
*/

    if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE){
        printf(" problem in FBO ! \n ");
        checkFramebufferStatus();
    }
}


// ============== initialization

void ScreenSDL2OGL3::setDefaults(){
	VIEW_DEPTH   = VIEW_DEPTH_DEFAULT;
	ASPECT_RATIO = WIDTH/(float)HEIGHT;
	zoom         = VIEW_ZOOM_DEFAULT;
	//printf(" %f %f %f \n", zoom, ASPECT_RATIO, VIEW_DEPTH  );
	mouse_begin_x  = 0;
	mouse_begin_y  = 0;
}

/*
void ScreenSDL2OGL3::init( int& id, int WIDTH_, int HEIGHT_ ){
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
*/

int ScreenSDL2OGL3::init( int WIDTH_, int HEIGHT_ ){

    WIDTH  = WIDTH_;
	HEIGHT = HEIGHT_;
	setDefaults();

    //printf("S DEBUG 1\n");
    window  = SDL_CreateWindow("Tutorial2", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WIDTH, HEIGHT, SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN);
    //printf("S DEBUG 2\n");
    context = SDL_GL_CreateContext( window );
    //printf("S DEBUG 3\n");
    SDL_GL_SetSwapInterval(1);
    //printf("S DEBUG 4\n");

    id = SDL_GetWindowID(window); printf( " win id %i \n", id );
	char str[40];  sprintf(str, " Window id = %d", id );
	//printf("S DEBUG 4\n");
	SDL_SetWindowTitle( window, str );
    if ( !window ){  printf( "Unable to initialize SDL\n" ); return -1; };
    //printf("S DEBUG 6\n");
	return id;
}

ScreenSDL2OGL3::ScreenSDL2OGL3( int WIDTH_, int HEIGHT_ ){
	init( WIDTH_, HEIGHT_ );
	qCamera.setOne();
};


// ========= MISC

bool ScreenSDL2OGL3::checkFramebufferStatus(){
    // check FBO status
    GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
    switch(status)    {
    case GL_FRAMEBUFFER_COMPLETE:
        printf( "Framebuffer complete.\n" );
        return true;
    case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT:
        printf( "[ERROR] Framebuffer incomplete: Attachment is NOT complete.\n" );
        return false;
    case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT:
        printf( "[ERROR] Framebuffer incomplete: No image is attached to FBO.\n" );
        return false;
    case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER:
        printf( "[ERROR] Framebuffer incomplete: Draw buffer.\n" );
        return false;
    case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER:
        printf( "[ERROR] Framebuffer incomplete: Read buffer.\n" );
        return false;
    case GL_FRAMEBUFFER_UNSUPPORTED:
        printf( "[ERROR] Framebuffer incomplete: Unsupported by FBO implementation.\n" );
        return false;
    default:
        printf( "[ERROR] Framebuffer incomplete: Unknown error.\n" );
        return false;
    }
}



