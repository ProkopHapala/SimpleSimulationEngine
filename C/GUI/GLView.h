#ifndef  GLView_h
#define  GLView_h

#include <stdbool.h>
#include "VectorTypes.h" 

#include "Camera.h"

#include <SDL2/SDL.h>
//#include <SDL2/SDL_opengl.h>

typedef struct GLView GLView;

typedef void (*ProcedurePointer)();
typedef void (*GLViewMethod)(GLView*);
//typedef void (*Method)(void*);

const float	VIEW_ZOOM_STEP       = 1.2f;
const float	VIEW_ZOOM_DEFAULT    = 10.0f;
const float	VIEW_DEPTH_DEFAULT   = 1000.0f;
const float	VIEW_MOVE_STEP       = 0.2f;

typedef struct GLView{
    int id;

    int   WIDTH,HEIGHT;
    float VIEW_DEPTH;
    float ASPECT_RATIO;
    float zoom;

    int   mouseX,mouseY;
    float mouse_begin_x;
    float mouse_begin_y;

    bool GL_LOCK;

    bool hasFocus;
    SDL_Window*      window;
    SDL_Renderer*    renderer;

    float  mouseRotSpeed  ;
    float  keyRotSpeed    ;
    float  cameraMoveSpeed;

    float4 qCamera    ;
    float4 qCameraOld ;
    //Mat3f  camMat     = Mat3fIdentity;

    bool  mouse_spinning;

    Camera cam;

    float camDist;
    //Vec2i spinning_start;

    bool perspective  ;
    bool first_person ;

    bool LMB,RMB;
    int  upTime,delay,timeSlice;
    int frameCount;
    bool loopEnd,STOP;
    float camStep;

    // Function Pointers

    GLViewMethod camera;
    GLViewMethod cameraHUD;
    GLViewMethod draw;
    GLViewMethod drawHUD;
    GLViewMethod inputHanding;

} GLView;


void init( int w, int h );
bool draw();
bool pre_draw();
bool post_draw();
int  run_Nframes(int nframes);
void set_draw_function( ProcedurePointer draw_func );

#endif
