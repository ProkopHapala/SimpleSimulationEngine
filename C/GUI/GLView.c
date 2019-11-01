#ifndef  GLView_c
#define  GLView_c

#include "GLView.h"

#include "Vec3math.h"
#include "Vec4math.h"
#include "Mat3math.h"

//#include "macros.h"
//#include "Camera"

//#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "Camera.c"

#include "Draw3D.c"

/*
#include "Vec3.h"
#include "Vec2.h"
#include "Mat3.h"
#include "quaternion.h"

#include "Camera.h"
#include "Draw3D.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "GLView.h"  // THE HEADER
*/

//namespace Cam{

void setupRenderer(){
    //float white    [] = { 1.0f, 1.0f,  1.0f,  1.0f };
    float ambient  [] = { 0.1f, 0.15f, 0.25f, 1.0f };
    float diffuse  [] = { 0.9f, 0.8f,  0.7f,  1.0f };
    float specular [] = { 1.0f, 1.0f,  1.0f,  1.0f };
    float shininess[] = { 80.0f                    };
    float lightPos [] = { 1.0f, -1.0f, 1.0f, 0.0f  };

    //glMaterialfv ( GL_FRONT_AND_BACK, GL_AMBIENT,   ambient);
    //glMaterialfv ( GL_FRONT_AND_BACK, GL_DIFFUSE,   diffuse);
    glMaterialfv ( GL_FRONT_AND_BACK, GL_SPECULAR,  specular);
    glMaterialfv ( GL_FRONT_AND_BACK, GL_SHININESS, shininess);

    glEnable     ( GL_COLOR_MATERIAL    );
    glLightfv    ( GL_LIGHT0, GL_POSITION,  lightPos );
    glLightfv    ( GL_LIGHT0, GL_DIFFUSE,   diffuse  );
    glLightfv    ( GL_LIGHT0, GL_AMBIENT,   ambient  );
    glLightfv    ( GL_LIGHT0, GL_SPECULAR,  specular );
    //glLightfv    ( GL_LIGHT0, GL_AMBIENT,  ambient  );
    glEnable     ( GL_LIGHTING         );
    glEnable     ( GL_LIGHT0           );
    glEnable     ( GL_NORMALIZE        );

    //glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, 1 );

    glEnable     ( GL_DEPTH_TEST       );
    glHint       ( GL_LINE_SMOOTH_HINT, GL_NICEST );
    glShadeModel ( GL_SMOOTH           );
    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
}

void drawCrosshair( float whalf, float hhalf,  float sz ){
    glBegin( GL_LINES );
    //float whalf = this->WIDTH *0.5;
    //float hhalf = this->HEIGHT*0.5;
    glVertex3f( whalf-10,hhalf, 0 ); glVertex3f( whalf+10,hhalf, 0 );
    glVertex3f( whalf,hhalf-10, 0 ); glVertex3f( whalf,hhalf+10, 0 );
    glEnd();
}

void Cam_ortho( const Camera cam, bool zsym ){
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    float zmin = cam.zmin; if(zsym) zmin=-cam.zmax;
    glOrtho( -cam.zoom*cam.aspect, cam.zoom*cam.aspect, -cam.zoom, cam.zoom, zmin, cam.zmax );
    float glMat[16];
    Draw3f_toGLMatCam( f3_Zero, cam.rot, glMat );
    glMultMatrixf( glMat );

    glMatrixMode ( GL_MODELVIEW );
    glLoadIdentity();
    glTranslatef(-cam.pos.x,-cam.pos.y,-cam.pos.z);
}

void Cam_perspective( const Camera cam ){
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    //glFrustum( -ASPECT_RATIO, ASPECT_RATIO, -1, 1, camDist/zoom, VIEW_DEPTH );
    glFrustum( -cam.aspect*cam.zoom, cam.aspect*cam.zoom, -cam.zoom, cam.zoom, cam.zmin, cam.zmax );
    //glFrustum( -cam.zoom*cam.aspect, cam.zoom*cam.aspect, -cam.zoom, cam.zoom, cam.zmin, cam.zmax );
    float glMat[16];
    Draw3f_toGLMatCam( f3_Zero, cam.rot, glMat );
    glMultMatrixf( glMat );

    glMatrixMode ( GL_MODELVIEW );
    glLoadIdentity();
    glTranslatef(-cam.pos.x,-cam.pos.y,-cam.pos.z);
    //glTranslatef ( -camPos.x+camMat.cx*camDist, -camPos.y+camMat.cy*camDist, -camPos.z+camMat.cz*camDist );
    //glTranslatef ( -cam.pos.x+camMat.cx*camDist, -camPos.y+camMat.cy*camDist, -camPos.z+camMat.cz*camDist );
}

// ===================== INPUTS

void GLView_quit( GLView* this ){ SDL_Quit(); this->loopEnd=true; /*exit(1);*/ };

void GLView_wait( GLView* this, int ms){
    for( int i=0; i<ms; i+=this->timeSlice ){
        uint32_t tnow=SDL_GetTicks();
        //printf( "GLView_wait %i %i \n", tnow, this->upTime+ms );
        if(tnow>=(this->upTime+ms)){this->upTime=tnow; break; }
        SDL_Delay(this->timeSlice);
    }
}

_inline float GLView_mouseUp     ( GLView* this, float mY ){ return 2*this->zoom*( 0.5 -mY/((float)this->HEIGHT)                          ); };
_inline float GLView_mouseUp_    ( GLView* this, float mY ){ return 2*this->zoom*(      mY/((float)this->HEIGHT) - 0.5                    ); };
_inline float GLView_mouseRight  ( GLView* this, float mX ){ return 2*this->zoom*(      mX/((float)this->HEIGHT) - 0.5*this->ASPECT_RATIO ); };


void GLView_updateMousePos ( GLView* this,  int x, int y ){
    this->mouse_begin_x = GLView_mouseRight( this, x );
    this->mouse_begin_y = GLView_mouseUp   ( this, y );
}

void GLView_eventHandling ( GLView* this, const SDL_Event* event  ){
    switch( event->type ){
        case SDL_KEYDOWN :
            switch( event->key.keysym.sym ){
                case SDLK_ESCAPE: GLView_quit(this); break;
                //case SDLK_SPACE:    STOP = !STOP; printf( STOP ? " STOPED\n" : " UNSTOPED\n"); break;
                case SDLK_KP_MINUS: this->zoom*=VIEW_ZOOM_STEP; break;
                case SDLK_KP_PLUS:  this->zoom/=VIEW_ZOOM_STEP; break;
                case SDLK_o:  this->perspective   = !this->perspective; break;
                case SDLK_p:  this->first_person  = !this->first_person ;   break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event->button.button ){
                case SDL_BUTTON_LEFT:  this->LMB = true;  break;
                case SDL_BUTTON_RIGHT: this->RMB = true;  break;
            };  break;
        case SDL_MOUSEBUTTONUP:
            switch( event->button.button ){
                case SDL_BUTTON_LEFT:  this->LMB = false; break;
                case SDL_BUTTON_RIGHT: this->RMB = false; break;
            }; break;
        case SDL_QUIT: GLView_quit(this); break;
    };
}

void GLView_keyStateHandling( GLView* this, const Uint8 *keys ){

    //if( keys[ SDL_SCANCODE_LEFT  ] ){ this->qCamera.dyaw  (  keyRotSpeed ); }
    //if( keys[ SDL_SCANCODE_RIGHT ] ){ this->qCamera.dyaw  ( -keyRotSpeed ); }
    //if( keys[ SDL_SCANCODE_UP    ] ){ this->qCamera.dpitch(  keyRotSpeed ); }
    //if( keys[ SDL_SCANCODE_DOWN  ] ){ this->qCamera.dpitch( -keyRotSpeed ); }

    if( keys[ SDL_SCANCODE_A ] ){ f3_fma_p( & this->cam.pos, this->cam.rot.a, -this->cameraMoveSpeed ); }
    if( keys[ SDL_SCANCODE_D ] ){ f3_fma_p( & this->cam.pos, this->cam.rot.a,  this->cameraMoveSpeed ); }
    if( keys[ SDL_SCANCODE_W ] ){ f3_fma_p( & this->cam.pos, this->cam.rot.b,  this->cameraMoveSpeed ); }
    if( keys[ SDL_SCANCODE_S ] ){ f3_fma_p( & this->cam.pos, this->cam.rot.b, -this->cameraMoveSpeed ); }
    if( keys[ SDL_SCANCODE_Q ] ){ f3_fma_p( & this->cam.pos, this->cam.rot.c, -this->cameraMoveSpeed ); }
    if( keys[ SDL_SCANCODE_E ] ){ f3_fma_p( & this->cam.pos, this->cam.rot.c,  this->cameraMoveSpeed ); }

    //printf( "frame %i keyStateHandling cam.pos (%g,%g,%g) \n", frameCount, cam.pos.x, cam.pos.y, cam.pos.z );

}

void GLView_mouseHandling( GLView* this ){
    SDL_GetMouseState( &this->mouseX, &this->mouseY );   this->mouseY=this->HEIGHT-this->mouseY;
    this->mouse_begin_x = (2*this->mouseX-this->WIDTH )*this->zoom/this->HEIGHT;
    this->mouse_begin_y = (2*this->mouseY-this->HEIGHT)*this->zoom/this->HEIGHT;
    int mx,my; Uint32 buttons = SDL_GetRelativeMouseState( &mx, &my);
    //printf( " %i %i \n", mx,my );
    if ( buttons & SDL_BUTTON(SDL_BUTTON_RIGHT)) {
        float4 q = f4_fromTrackball( 0, 0, -mx*this->mouseRotSpeed, my*this->mouseRotSpeed );
        this->qCamera = f4_qmul( this->qCamera, q );
    }
    //qCamera.qmul( q );
}

void GLView_inputHanding( GLView* this ){
    //printf( "GLView_inputHanding \n" );
    const Uint8 *keys = SDL_GetKeyboardState(NULL);
    GLView_keyStateHandling( this, keys );
    GLView_mouseHandling( this );
    SDL_Event		 event;
    while(SDL_PollEvent(&event)){
        GLView_eventHandling( this, &event );
    }
}





// ======== Update

void GLView_camera( GLView* this ){
    this->cam.rot = f4_toMat3( this->qCamera );
    this->cam.zoom   = this->zoom;
    this->cam.aspect = this->ASPECT_RATIO;
    if (this->perspective){ Cam_perspective( this->cam ); }
    else                 { Cam_ortho( this->cam, true ); }
}

void GLView_cameraHUD( GLView* this ){
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    glOrtho ( 0, this->WIDTH, 0, this->HEIGHT, -this->VIEW_DEPTH, +this->VIEW_DEPTH );
    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity();
}

void GLView_draw(){
    //printf( "default_draw \n" );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    glEnable    ( GL_LIGHTING );
    glShadeModel( GL_FLAT     );

    Draw3f_box  ( -1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 1.0f, 0.8f, 0.8f, 0.8f );
    //Draw3D_box       ( -1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 1.0f, 0.8f, 0.8f, 0.8f );

    glShadeModel( GL_SMOOTH     );
    //Draw3D_sphere_oct( 5, 1.0f, (float3){3.0,3.0,3.0} );

    glDisable ( GL_LIGHTING );
    //Draw3D_axis ( 3.0f );
    Draw3f_axis( 3.0f );
}

void GLView_drawHUD(){

}

bool GLView_pre_draw( GLView* this ){
    //thisApp->update();
    if( this->GL_LOCK ) return true;
    this->inputHanding(this);
    this->GL_LOCK = true;
    this->camera(this);

    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glEnable    ( GL_LIGHTING );
	glShadeModel( GL_SMOOTH );

    //printf( "DEBU.add_mul(meCount );
    return this->GL_LOCK;
}

bool GLView_post_draw( GLView* this  ){
    //printf( "DEBUG post_draw[fame=%i] \n", frameCount );
    SDL_RenderPresent(this->renderer);
    this->frameCount++;
    this->GL_LOCK = false;
    return this->loopEnd;
}

bool GLView_update( GLView* this ){
    //printf( "GLView::update %i \n", this->GL_LOCK );
    if( this->GL_LOCK ) return true;
    this->GL_LOCK = true;
    //printf( "DEBUG 1 %li %li %li \n", (long)this->inputHanding, (long)this->camera, (long)this->draw );
    this->inputHanding(this);
    this->camera(this);
    //this->pre_draw();
    this->draw(this);
    this->cameraHUD(this);
    this->drawHUD(this);
    GLView_post_draw( this );
    return this->GL_LOCK;
}

int GLView_loop( GLView* this, int n ){
    this->loopEnd = false;
    int iframe=0;
    for( ; iframe<n; iframe++ ){
        //printf( "GLView::loop[%i] \n", iframe );
        //GLView_inputHanding();
        //if(!STOP){update();} // DEPRECATED: usually we want to stop physics, not drawing
        GLView_update(this);
        //printf(" %i \n", iframe );
        GLView_wait(this, this->delay);
        //if( delay>0 ) SDL_Delay( delay );
        if(this->loopEnd) break;
    }
    return iframe;
}



// ================ Constructors

void GLView_setDefaults( GLView* this ){
    this->VIEW_DEPTH   = VIEW_DEPTH_DEFAULT;
    this->ASPECT_RATIO = this->WIDTH/(float)this->HEIGHT;
    this->zoom         = VIEW_ZOOM_DEFAULT;
    //printf(" %f %f %f \n", zoom, ASPECT_RATIO, VIEW_DEPTH  );
    this->mouse_begin_x  = 0;
    this->mouse_begin_y  = 0;
}

void GLView_init( GLView* this, int id, int WIDTH_, int HEIGHT_ ){
    this->WIDTH  = WIDTH_;
    this->HEIGHT = HEIGHT_;
    GLView_setDefaults(this);
    SDL_CreateWindowAndRenderer(this->WIDTH, this->HEIGHT, SDL_WINDOW_OPENGL, &this->window, &this->renderer);
    this->id = SDL_GetWindowID(this->window); 
    printf( " win id %i \n", id );
    char str[40];  
    sprintf(str, " Window id = %d", this->id );
    SDL_SetWindowTitle( this->window, str );
    setupRenderer();
    //printf( " ASPECT_RATIO %f \n", ASPECT_RATIO );
}




GLView GLView_default( ){
    GLView this;
    this.WIDTH=800; this.HEIGHT=600;
    this.VIEW_DEPTH=VIEW_DEPTH_DEFAULT;
    this.ASPECT_RATIO=8.f/6.f;
    this.zoom        =10.f;

    this.mouseX=0,this.mouseY=0;
    this.mouse_begin_x=0;
    this.mouse_begin_y=0;

    this.GL_LOCK = false;

    this.hasFocus=true;
    this.window=0;
    this.renderer=0;

    this.mouseRotSpeed   = 0.001;
    this.keyRotSpeed     = 0.01;
    this.cameraMoveSpeed = 0.2f;

    this.qCamera    = f4_Eye;
    this.qCameraOld = f4_Eye;
    //Mat3f  camMat     = Mat3fIdentity;

    this.mouse_spinning    = false;

    this.cam = Camera_default();

    this.camDist = 50.0;
    //this.spinning_start;

    this.perspective  = false;
    this.first_person = false;

    //this.LMB=false,this.false;
    this.upTime=0;
    this.delay=20;
    this.timeSlice=5;
    this.frameCount=0;
    this.loopEnd = false; 
    this.STOP    = false;
    this.camStep = VIEW_MOVE_STEP;

    
    this.camera       = GLView_camera;
    this.cameraHUD    = GLView_cameraHUD;

    this.draw         = GLView_draw;
    this.drawHUD      = GLView_drawHUD;
    this.inputHanding = GLView_inputHanding;
    
    /*
    this.camera       = 0;
    this.cameraHUD    = 0;

    this.draw         = 0;
    this.drawHUD      = 0;
    this.inputHanding = 0;
    */
    //printf( "DEBUG 1 %li %li %li \n", (long)this.inputHanding, (long)this.camera, (long)this.draw );
    return this;
}

GLView GLView_( int id, int WIDTH_, int HEIGHT_ ){
    GLView this = GLView_default( );
    //printf( "DEBUG 1 %li %li %li \n", (long)this.inputHanding, (long)this.camera, (long)this.draw );
    GLView_init(&this, id,WIDTH_,HEIGHT_);
    this.qCamera = f4_Eye;
    this.cam.rot = f4_toMat3( this.qCamera );
    this.cam.pos = f3_newf( 0 );
    //GLbyte* s;
    // http://stackoverflow.com/questions/40444046/c-how-to-detect-graphics-card-model
    printf( "GL_VENDOR  : %s \n", glGetString(GL_VENDOR)  );
    printf( "GL_VERSION : %s \n", glGetString(GL_VERSION) );
    return this;
}


//GLView * thisApp = 0;

//extern "C"{

/*
//void init( int w, int h, void* world_, char* work_dir_ ){
void init( int w, int h ){
    //world = (FlightWorld*)world_;
    //work_dir = work_dir_;
    SDL_Init(SDL_INIT_VIDEO);
    SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
    //SDL_SetRelativeMouseMode( SDL_TRUE );
    //SDL_ShowCursor(SDL_DISABLE);
    int junk;
    thisApp = new GLView( junk , w, h );
    SDL_SetWindowPosition(thisApp->window, 100, 0 );
    thisApp->loopEnd = false;
}
*/

/*
bool draw(){ return thisApp->update(); }

bool pre_draw(){ return thisApp->pre_draw(); }

bool post_draw(){ return thisApp->post_draw(); }

int run_Nframes(int nframes){ thisApp->loop( nframes ); }

void set_draw_function( ProcedurePointer draw_func ){ thisApp->draw_func_ptr = draw_func; }
*/

#endif