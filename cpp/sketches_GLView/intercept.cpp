#if 0
target="$0"
binout=$target.x
PREDIR="../"
PREI="-I"$PREDIR
IFLAGS=$PREI"common/math "$PREI"common/CombatModels "$PREI"common_SDL/SDL2OGL "$PREI"common_SDL/ "$PREI"common/utils "$PREI"common/dataStructures "$PREI"common_SDL/SDL2OGL "$PREI"libs_SDL/GLView -I/usr/include" 
echo "PREI " $PREI
echo "IFLAGS " $IFLAGS
ERRFLAGS="-Wno-missing-braces -Werror=return-type"
CFLAGS="-std=c++17 -Ofast -march=native -mtune=native"
LFLAGS="-L"$PREDIR"Build/libs_SDL/GLView -lGLView -lGL -lSDL2"
g++ -o $binout $target $CFLAGS $ERRFLAGS $IFLAGS $LFLAGS 
#clang++ -o $binout $target $CFLAGS $ERRFLAGS $IFLAGS $LFLAGS 
ls
bin_path=$PREDIR/Build/libs_SDL/GLView
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$bin_path
echo $LD_LIBRARY_PATH
$binout
exit
#endif

#include <cmath>
#include <cstdio>
#include "GLView.h"

#include "fastmath.h"

#include "Draw.h"
#include "Draw2D.h"
#include "Draw3D.h"
//#include "Plot2D.h"
//#include "PlotScreen2D.h"
#include "Solids.h"
#include "testUtils.h"

#include "SDL_utils.h"

//#include "gonioApprox.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

// ==================== Clases

struct Operation{
    char in1,in2;
    char op;
    char out1;
};

class Interpeter{
    double scalars[128];
    Vec3d  vecs   [128];

    inline void exec( char i1, char i2, char ret, char op  ){
        //Vec3d v1,v2;
        double f1,f2,fo;
        f1=scalars[i1];
        f2=scalars[i2];
        switch(op){
            case '+': fo=f1+f2; break;
            case '-': fo=f1-f2; break;
            case '*': fo=f1*f2; break;
            case '/': fo=f1/f2; break;
            case '^': fo=pow(f1,f2); break;
        }
        scalars[ret]=fo;
    }
    
};


class Pray{ public:
    Vec3d pos;
    Vec3d vel;
    void move(double dt){
        pos = pos + vel*dt;
    }
};


class Interceptor{ public:
    Vec3d pos;
    Vec3d vel;
    Vec3d thrust;
    double thrustMax2;
    void limitTrhust(){
        double f2 = thrust.norm2();
        if(f2>thrustMax2){ thrust.mul(sqrt(thrustMax2/f2)); }
    }
    void control( const Pray& pray, double Dt ){
        Vec3d dpos   = pray.pos - pos;
        Vec3d dvel   = pray.vel - vel;
        Vec3d dpos_t = dpos + dvel*Dt;
        thrust = (dpos_t*2.).mul( 1/(Dt*Dt) ); 
    }
    void move(double dt){
        limitTrhust();
        vel = vel + thrust*dt;
        pos = pos + vel*dt;
    }
};



// ==================== Variables

int ogl=0;
//Plot2D plot1;
//int fontTex;

double dt=0.01;
double Dt=1.25; 
Interceptor hunter;
Pray        pray;


// ==================== Functions


void my_draw(){

    //printf(" my_draw ! %i %i %g \n", plot1.lines.size(), plot1.lines.back()->n, plot1.lines.back()->ys[2]  );
    glClearColor( 1.0f, 1.0f, 1.0f, 0.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glDisable( GL_DEPTH_TEST );

    pray  .move(dt);
    hunter.control( pray, Dt );
    hunter.move(dt);

    glColor3f(0.0,0.0,1.0);
    Draw3D::drawPointCross(pray.pos, 0.25 );
    Draw3D::drawVecInPos  (pray  .vel*Dt, pray  .pos );

    glColor3f(1.0,0.0,0.0);
    Draw3D::drawPointCross(hunter.pos, 0.25 );
    Draw3D::drawVecInPos  (hunter.vel*Dt, hunter.pos );
    glColor3f(0.0,0.7,0.0);
    Draw3D::drawVecInPos  (hunter.thrust, hunter.pos );

}

void setup(){
    pray.pos      = (Vec3d){ -5.0, 2.0, 0.0 };
    pray.vel      = (Vec3d){ 1.0, 0.0, 0.0 };
    hunter.pos    = (Vec3d){ 0.0, -2.0, 0.0 };
    hunter.vel    = (Vec3d){ 0.0, 0.0, 0.0 };
    hunter.thrust = Vec3dZero;
    hunter.thrustMax2 = sq(0.4);
}

int main(){
     printf( "HERE !!!! \n"  );
    init( 800, 600 );
    set_draw_function( my_draw );
    setup();
    run_Nframes(1000000);
}
