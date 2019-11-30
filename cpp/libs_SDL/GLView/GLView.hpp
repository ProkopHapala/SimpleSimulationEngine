#include "Vec3.h"
#include "Vec2.h"
#include "Mat3.h"
#include "quaternion.h"

#include "Camera.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "GLView.h" 

class GLView{ public:

    constexpr static const float	VIEW_ZOOM_STEP       = 1.2f;
    constexpr static const float	VIEW_ZOOM_DEFAULT    = 10.0f;
    constexpr static const float	VIEW_DEPTH_DEFAULT   = 1000.0f;
    constexpr static const float	VIEW_MOVE_STEP       = 0.2f;

    int   WIDTH=800,HEIGHT=600;
    float VIEW_DEPTH=VIEW_DEPTH_DEFAULT;
    float ASPECT_RATIO=8.f/6.f;
    float zoom=10.f;

    //float camX0=0.0f,camY0=0.0f;
    //float fWIDTH, fHEIGHT, camXmin, camYmin, camXmax, camYmax;

    int   mouseX=0,mouseY=0;
    float mouse_begin_x=0;
    float mouse_begin_y=0;

    bool GL_LOCK = false;

    bool hasFocus;
    SDL_Window*      window;
    SDL_Renderer*    renderer;

    float  mouseRotSpeed   = 0.001;
    float  keyRotSpeed     = 0.01;
    float  cameraMoveSpeed = 0.2f;

    Quat4f qCamera    = Quat4fIdentity;
    Quat4f qCameraOld = Quat4fIdentity;
    //Mat3f  camMat     = Mat3fIdentity;

    bool  mouse_spinning    = false;

    Camera cam;

    float camDist = 50.0;
    Vec2i spinning_start;

    bool perspective  = false;
    bool first_person = false;


	public:
	bool LMB=false,RMB=false;
	int  upTime=0,delay=20,timeSlice=5,frameCount=0;
	bool loopEnd    = false, STOP = false;
	float camStep   = VIEW_MOVE_STEP;


// ============ function declarations

    void wait        (float ms);
    virtual void quit(       );
    void         wait(int ms);
    virtual int loop( int n );
    void defaultMouseHandling    ( const int& mouseX, const int& mouseY );

    // ==== Functions Virtual

    // ---- setup
    virtual void setupRenderer ();
    virtual void setDefaults   ();

    // ---- Inputs

    virtual void keyStateHandling( const Uint8 *keys       );
    virtual void eventHandling   ( const SDL_Event& event  );
    virtual void mouseHandling   (                         );
    virtual void inputHanding();
    virtual void updateMousePos ( int x, int y );

    // ---- Draw & Update
    virtual void camera      ();
    virtual void cameraHUD   ();
    virtual void draw        ();
    virtual void drawHUD     ();
    
    //ProcedurePointer draw_func_ptr = default_draw;
    ProcedurePointer draw_func_ptr = 0;

    // ==== Functions

    // ---- Init
    void init( int& id, int WIDTH_, int HEIGHT_ );
    GLView ( int& id, int WIDTH_, int HEIGHT_ );

    // ---- Draw & Update
    void  startSpining ( float x, float y              ){ mouse_spinning = true; mouse_begin_x  = x; mouse_begin_y  = y;	}
    void  endSpining   (                               ){ mouse_spinning = false;	                                    }
    //void  projectMouse ( float mX, float mY, Vec3d& mp ){ mp.set_lincomb( mouseRight(mX), camRight,  mouseUp(mY), camUp ); };
    void  projectMouse ( float mX, float mY, Vec3f& mp ){ mp.set_lincomb( mouseRight(mX), cam.rot.a,  mouseUp(mY), cam.rot.b ); };
    void drawCrosshair( float sz );

    // ---- Camera
    //void camera_FPS       ( const Vec3d& pos, const Mat3d& rotMat );
    //void camera_FwUp      ( const Vec3d& pos, const Vec3d& fw, const Vec3d& up, bool upDominant );
    //void camera_FreeLook  ( const Vec3d& pos );
    //void camera_OrthoInset( const Vec2d& p1, const Vec2d& p2, const Vec2d& zrange, const Vec3d& fw, const Vec3d& up, bool upDominant );

    bool update( );
    bool pre_draw ( );
    bool post_draw( );

    inline float mouseUp     ( float mY ){ return 2*zoom*( 0.5 -mY/float(HEIGHT)                    ); };
    inline float mouseUp_    ( float mY ){ return 2*zoom*(      mY/float(HEIGHT) - 0.5              ); };
    inline float mouseRight  ( float mX ){ return 2*zoom*(      mX/float(HEIGHT) - 0.5*ASPECT_RATIO ); };

    inline bool wait_LOCK( int n, int ms ){ if(!GL_LOCK) return true; for(int i=0; i<n; i++){ SDL_Delay(ms); if(!GL_LOCK) return true; } return false; }

};





GLView* getGLViewInstance();