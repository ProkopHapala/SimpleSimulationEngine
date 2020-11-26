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
#g++ -o $binout $target $CFLAGS $ERRFLAGS $IFLAGS $LFLAGS 
clang++ -o $binout $target $CFLAGS $ERRFLAGS $IFLAGS $LFLAGS 
./$binout
rm -f $binout
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


int ogl=0;
//Plot2D plot1;
//int fontTex;

// ==================== MAIN

void my_draw(){

    //printf(" my_draw ! %i %i %g \n", plot1.lines.size(), plot1.lines.back()->n, plot1.lines.back()->ys[2]  );
    glClearColor( 1.0f, 1.0f, 1.0f, 0.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    //glDisable( GL_DEPTH_TEST );

    Draw3D::drawLine((Vec3d){0,0,0},(Vec3d){1,1,1});
    glCallList(ogl);
}

void setup(){
    uint32_t cpoly = 0xFFFFFFFF;
    uint32_t cwire = 0xFF000000;

    double ang = atan2( sqrt(2), 1 )/M_PI;
    printf( "angle %g \n", ang );

    double sz = 3.2;
    Vec3d a = {    0.00000000000,     sz     , 0};
    Vec3d b = { sz*0.86602540378,     sz*0.5 , 0};
    Vec3d c = { sz*0.86602540378*0.5*2./3., sz*0.75*2./3., sz*0.86602540378};

    ogl=Draw::list(ogl);
    //glDisable(GL_LIGHTING);
    /*
    float sc=2.0;
    glScalef(sc,sc,sc);
    Draw3D::drawMesh( Solids::Octahedron, 0xFF0000FF, 0xFFFFFFFF );
    glScalef(1/sc,1/sc,1/sc);
    Draw3D::drawMesh( Solids::Cube      , 0xFFFF0000, 0xFFFFFFFF );
    Draw3D::drawMesh( Solids::RhombicDodecahedron, 0, 0xFF00FF00 );
    */

   int n=2;

    Draw3D::drawAxis(5.0);
    Draw3D::drawVec( (Vec3d){10,10,1} );
    for(int ic=0; ic<=1; ic++){
        for(int ib=-n; ib<=n-ic; ib++){
            for(int ia=-n; ia<=n-ic; ia++){
                float x = a.x*ia + b.x*ib + c.x*ic;
                float y = a.y*ia + b.y*ib + c.y*ic;    
                float z = a.z*ia + b.z*ib + c.z*ic;            
                glPushMatrix();
                glTranslatef(x,y,z);
                glRotatef( 45,      0,0,1 );
                glRotatef( 180*ang, 1,1,0 );
                Draw3D::drawMesh( Solids::RhombicDodecahedron, 0xFF808080, 0xFF00FF00 );
                glPopMatrix();
            }
        }
    }
    
    glEndList();

}

int main(){
    init( 800, 600 );
    set_draw_function( my_draw );
    setup();
    run_Nframes(1000000);
}
