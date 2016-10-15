
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "Draw3D.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"

#include "spline_hermite.h"
#include "SplineManager.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

class OrbitEditor : public AppSDL2OGL_3D {
	public:
    SplineManager splines;

    int     nsamp=0;
    double tstart, tend, dt;
    Vec3d * positions=NULL;
    Vec3d * velocities=NULL;
    Vec3d * accelerations=NULL;
    int which[3] = {0,1,2};

    int iedit=0, idraged=0;
    double mouse_t,mouse_val;
    double timeScale=100.0,valScale=100.0f;

    Vec3d *  dpos;
    Vec3d * ddpos;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );
    virtual void mouseHandling( );

	OrbitEditor( int& id, int WIDTH_, int HEIGHT_ );

};

OrbitEditor::OrbitEditor( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    double dt_ = 1.0;
    bool derivs = true;
    splines.allocate( 3, 10, derivs );
    for(int i=0; i<splines.n; i++){
        splines.ts [i]     = i*dt_;
        splines.CPs[0][i]  = randf( -1.0f, 1.0f );
        splines.CPs[1][i]  = randf( -1.0f, 1.0f );
        splines.CPs[2][i]  = randf( -1.0f, 1.0f );
        if(derivs){
            splines.dCPs[0][i]  = randf( -1.0f, 1.0f )/dt_;
            splines.dCPs[1][i]  = randf( -1.0f, 1.0f )/dt_;
            splines.dCPs[2][i]  = randf( -1.0f, 1.0f )/dt_;

            //splines.dCPs[0][i]  = 0;
            //splines.dCPs[1][i]  = 0;
            //splines.dCPs[2][i]  = 0;
            printf("%i %3.3f (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) \n", i, splines.ts[i], splines.CPs[0][i], splines.CPs[1][i], splines.CPs[2][i],  splines.dCPs[0][i], splines.dCPs[1][i], splines.dCPs[2][i]    );
        }else{
            printf("%i %3.3f (%3.3f,%3.3f,%3.3f)\n", i, splines.ts[i], splines.CPs[0][i], splines.CPs[1][i], splines.CPs[2][i]    );
        }

    }
    splines.ts[2]     -= 0.3;

    int nsamp_max = 1000;
    positions      = new Vec3d[nsamp_max ];
    velocities     = new Vec3d[nsamp_max ];
    accelerations  = new Vec3d[nsamp_max ];

    tstart = 1.15;
    tend   = 3.9;
    dt     = 0.01;
    nsamp = splines.evalUniform( tstart, tend, dt, 3, which, (double*)positions, (double*)velocities, (double*)accelerations );

    dpos   = new Vec3d[nsamp_max ];
    ddpos  = new Vec3d[nsamp_max ];

    double invdt = 1/dt;
    for(int i=0;i<nsamp-1;i++){ dpos [i]=( positions[i+1]-positions[i] )*invdt;  };
    for(int i=0;i<nsamp-2;i++){ ddpos[i]=( dpos[i+1]-dpos[i] )*invdt;             };

}

void OrbitEditor::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    //glEnable( GL_LIGHTING );
    //glCallList( world.defaultObjectShape );
	glDisable ( GL_LIGHTING );
	Draw3D::drawAxis ( 3.0f );

    glColor3f(1.0f,1.0f,1.0f);
    if(nsamp>0) Draw3D::drawPolyLine( nsamp, positions );

    //printf( "camMat.a (%3.3f,%3.3f,%3.3f) \n", camMat.a.x, camMat.a.y, camMat.a.z );

};

void OrbitEditor::drawHUD(){
    glDisable ( GL_LIGHTING );


    glPushMatrix();

        glTranslatef(20, 200, 0);
        glScalef    (timeScale,valScale,1);
        glTranslatef(-tstart, 0, 0);

        glColor3f( 0.0f, 1.0f, 0.0f );
        Draw2D::drawPointCross( {mouse_t,mouse_val} , 0.1 );
        if(idraged!=0){ Draw2D::drawLine( {mouse_t,mouse_val}, {splines.ts[idraged],splines.CPs[iedit][idraged]} ); };

        glColor3f( 1.0f, 1.0f, 1.0f );
        Draw2D::drawLine( { 0, 0.0}, {WIDTH, 0.0} );

        glColor3f( 1.0f, 1.0f, 1.0f ); Draw2D::plot_cross( splines.n, splines.ts, splines.CPs[0], 0.03 );

        glColor3f( 1.0f, 1.0f, 1.0f ); glBegin(GL_LINE_STRIP); for(int i=0; i<nsamp; i++){ glVertex3d( (tstart+i*dt), (float)positions[i].x, 0.0 ); }; glEnd();

        /*
        glColor3f( 1.0f, 0.0f, 0.0f ); glBegin(GL_LINE_STRIP); for(int i=0; i<nsamp; i++){ glVertex3d( (tstart+i*dt), 0.5*(float)velocities[i].x, 0.0 ); }; glEnd();
        glColor3f( 0.0f, 1.0f, 0.0f ); glBegin(GL_LINE_STRIP); for(int i=0; i<nsamp; i++){ glVertex3d( (tstart+i*dt), 0.25*(float)accelerations[i].x, 0.0 ); }; glEnd();

        glColor3f( 1.0f, 0.0f, 1.0f ); glBegin(GL_LINE_STRIP); for(int i=0; i<nsamp; i++){ glVertex3d( (tstart+i*dt), 0.1+0.5*(float)dpos[i].x, 0.0 ); }; glEnd();
        glColor3f( 0.0f, 1.0f, 1.0f ); glBegin(GL_LINE_STRIP); for(int i=0; i<nsamp; i++){ glVertex3d( (tstart+i*dt), 0.1+0.25*(float)ddpos[i].x, 0.0 ); }; glEnd();
        */

        glColor3f( 1.0f, 1.0f, 1.0f );
        glBegin(GL_LINES);
        double * dCPi=splines.dCPs[iedit];
        for(int i=1; i<splines.n-1; i++){
            double t,d,p;
            t=splines.ts[i];
            p=splines.CPs[iedit][i];
            if(dCPi){d=dCPi[i];}else{  d=splines.inferDeriv(i,iedit); }
            glVertex3d( t-0.5, p-0.5*d, 0.0 );
            glVertex3d( t+0.5, p+0.5*d, 0.0 );
        }
        glEnd();

    glPopMatrix();
}

//void OrbitEditor::keyStateHandling( const Uint8 *keys ){ };

void OrbitEditor::mouseHandling( ){
    SDL_GetMouseState( &mouseX, &mouseY ); mouseY=HEIGHT-mouseY;
    mouse_t   = (mouseX-20 )/timeScale + tstart;
    mouse_val = (mouseY-200)/valScale;
    printf("%i %i %3.3f %3.3f \n", mouseX,mouseY, mouse_t, mouse_val );
}

void OrbitEditor::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_f:  warrior1->tryJump(); break;
                //case SDLK_h:  warrior1->tryJump(); break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
                case SDLK_ESCAPE:   quit(); break;
                //case SDLK_SPACE:    STOP = !STOP; printf( STOP ? " STOPED\n" : " UNSTOPED\n"); break;
                case SDLK_KP_MINUS: zoom*=VIEW_ZOOM_STEP; break;
                case SDLK_KP_PLUS:  zoom/=VIEW_ZOOM_STEP; break;
            }
            break;

        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                case SDL_BUTTON_RIGHT:
                    float t   = (mouseX+20 )/timeScale + tstart;
                    //float val = (mouseY-200)/valScale;
                    idraged = binSearchFrom<double>(t,splines.n,splines.ts);
                    break;
            }
            break;

        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    if(idraged!=0){
                        splines.CPs[iedit][idraged]=mouse_val;
                        splines.ts        [idraged]=mouse_t;
                        nsamp = splines.evalUniform( tstart, tend, dt, 3, which, (double*)positions, (double*)velocities, (double*)accelerations );
                    }
                    idraged=0;
                    break;
                case SDL_BUTTON_RIGHT:
                    if(idraged!=0){
                        double * dCPi = splines.dCPs[iedit];
                        if(dCPi==0) break;
                        dCPi       [idraged]=(mouse_val-splines.CPs[iedit][idraged])/(mouse_t-splines.ts[idraged]);
                        //splines.ts [idraged]=mouse_t;
                        nsamp = splines.evalUniform( tstart, tend, dt, 3, which, (double*)positions, (double*)velocities, (double*)accelerations );
                    }
                    idraged=0;
                    break;
            }
            break;
    };
    //AppSDL2OGL::eventHandling( event );
}




// ===================== MAIN

OrbitEditor * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new OrbitEditor( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















