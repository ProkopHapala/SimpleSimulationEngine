
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "Draw3D.h"
#include "SDL_utils.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"

#include "spline_hermite.h"
#include "SplineManager.h"

#include "AppSDL2OGL_3D.h"
#include "GUI.h"
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

    bool dragging = false;
    int iedit=0, ipoint=0;
    double mouse_t,mouse_val;
    double timeScale=100.0,valScale=100.0f;

    int          fontTex;
    GUITextInput txtStatic;
    //GUITextInput txtPop;
    char curCaption[100];

    //DEGUB
    //Vec3d *  dpos;
    //Vec3d * ddpos;

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
    tend   = 8.9;
    dt     = 0.01;
    nsamp = splines.evalUniform( tstart, tend, dt, 3, which, (double*)positions, (double*)velocities, (double*)accelerations );

    //DEGUB
    //dpos   = new Vec3d[nsamp_max ];
    //ddpos  = new Vec3d[nsamp_max ];
    //double invdt = 1/dt;
    //for(int i=0;i<nsamp-1;i++){ dpos [i]=( positions[i+1]-positions[i] )*invdt;  };
    //for(int i=0;i<nsamp-2;i++){ ddpos[i]=( dpos[i+1]-dpos[i] )*invdt;             };

    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    //txtPop   .inputText = "txtPop";
    txtStatic.inputText = "txtStatic";

    zoom = 2;

}

void OrbitEditor::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    //glEnable( GL_LIGHTING );
    //glCallList( world.defaultObjectShape );

    //glMatrixMode(GL_MODELVIEW);
    glMatrixMode(GL_PROJECTION);
    glTranslatef(0,1.0,0);

	glDisable ( GL_LIGHTING );
	Draw3D::drawAxis ( 3.0f );

    glColor3f(0.0f,0.0f,0.0f);
    if(nsamp>0) Draw3D::drawPolyLine( nsamp, positions );
    for(int i=0; i<splines.n; i++){
        Draw3D::drawPointCross( {splines.CPs[0][i],splines.CPs[1][i],splines.CPs[2][i]}, 0.02 );
    }

    Vec3d p,v,a;
    splines.eval(mouse_t,0, 3, which, (double*)&p, (double*)&v, (double*)&a );
    glColor3f(1.0f,1.0f,1.0f); Draw3D::drawPointCross( p, 0.1 );
    glColor3f(1.0f,1.0f,0.0f); Draw3D::drawVecInPos  ( v, p );
    glColor3f(0.0f,1.0f,1.0f); Draw3D::drawVecInPos  ( a, p );

    //printf( "camMat.a (%3.3f,%3.3f,%3.3f) \n", camMat.a.x, camMat.a.y, camMat.a.z );

};

void OrbitEditor::drawHUD(){
    glDisable ( GL_LIGHTING );

    glColor3f(1.0f,1.0f,1.0f);   txtStatic.view3D( {5,5}, fontTex, 8 );

    glPushMatrix();

        glTranslatef(20, 200, 0);
        glScalef    (timeScale,valScale,1);
        glTranslatef(-tstart, 0, 0);

        // axis-zero
        glColor3f( 0.0f, 0.0f, 0.0f ); Draw2D::drawLine( { 0, 0.0}, {WIDTH, 0.0} );

        // curves
        glColor3f( 0.0f, 0.0f, 0.0f ); Draw2D::plot_cross( splines.n, splines.ts, splines.CPs[iedit], 0.03 );

        // curves
        //glColor3f( 1.0f, 1.0f, 1.0f ); glBegin(GL_LINE_STRIP); for(int i=0; i<nsamp; i++){ glVertex3d( (tstart+i*dt), (float)positions[i].array[iedit], 0.0 ); }; glEnd();
        glColor3f( 1.0f, 0.0f, 0.0f ); glBegin(GL_LINE_STRIP); for(int i=0; i<nsamp; i++){ glVertex3d( (tstart+i*dt), (float)positions[i].x, 0.0 ); }; glEnd();
        glColor3f( 0.0f, 1.0f, 0.0f ); glBegin(GL_LINE_STRIP); for(int i=0; i<nsamp; i++){ glVertex3d( (tstart+i*dt), (float)positions[i].y, 0.0 ); }; glEnd();
        glColor3f( 0.0f, 0.0f, 1.0f ); glBegin(GL_LINE_STRIP); for(int i=0; i<nsamp; i++){ glVertex3d( (tstart+i*dt), (float)positions[i].z, 0.0 ); }; glEnd();

        //DEGUB
        //glColor3f( 1.0f, 0.0f, 0.0f ); glBegin(GL_LINE_STRIP); for(int i=0; i<nsamp; i++){ glVertex3d( (tstart+i*dt), 0.5*(float)velocities[i].x, 0.0 ); }; glEnd();
        //glColor3f( 0.0f, 1.0f, 0.0f ); glBegin(GL_LINE_STRIP); for(int i=0; i<nsamp; i++){ glVertex3d( (tstart+i*dt), 0.25*(float)accelerations[i].x, 0.0 ); }; glEnd();
        //glColor3f( 1.0f, 0.0f, 1.0f ); glBegin(GL_LINE_STRIP); for(int i=0; i<nsamp; i++){ glVertex3d( (tstart+i*dt), 0.1+0.5*(float)dpos[i].x, 0.0 ); }; glEnd();
        //glColor3f( 0.0f, 1.0f, 1.0f ); glBegin(GL_LINE_STRIP); for(int i=0; i<nsamp; i++){ glVertex3d( (tstart+i*dt), 0.1+0.25*(float)ddpos[i].x, 0.0 ); }; glEnd();

        // tangents
        glColor3f( 0.0f, 0.0f, 0.0f );
        glBegin(GL_LINES);
        double * dCPi=splines.dCPs[iedit];
        for(int i=1; i<splines.n-1; i++){
            double t,d,p;
            t =splines.ts[i];
            p =splines.CPs[iedit][i];
            d =splines.getPointDeriv(i,iedit);
            //if(dCPi){d=dCPi[i];}else{  d=splines.inferDeriv(i,iedit); }
            glVertex3d( t-0.5, p-0.5*d, 0.0 );
            glVertex3d( t+0.5, p+0.5*d, 0.0 );
        }
        glEnd();

        // cursor
        //Vec2d pcur;    pcur.set( mouse_t,mouse_val );
        Vec2f ppoint;  ppoint.set( splines.ts[ipoint],splines.CPs[iedit][ipoint] );
        glColor3f( 1.0f, 1.0f, 1.0f );
        //Draw2D::drawPointCross( {mouse_t,mouse_val} , 0.1 );
        Draw2D::drawLine( {mouse_t-0.1,mouse_val}, {mouse_t+0.1,mouse_val} );
        Draw2D::drawLine( {mouse_t,valScale}, {mouse_t,-valScale} );
        Draw2D::drawCircle( ppoint, 0.1, 16, false );
        if(dragging){ Draw2D::drawLine( {mouse_t,mouse_val}, ppoint ); };

        //txtPop.view3D( {splines.ts[ipoint],splines.CPs[iedit][ipoint],1.0}, fontTex, 10 );
        glTranslatef(mouse_t,0,0);
        glRotatef(90,0,0,1);
        Draw::drawText(curCaption, fontTex, 0.1, 0, 0 );

    glPopMatrix();
}

//void OrbitEditor::keyStateHandling( const Uint8 *keys ){ };

void OrbitEditor::mouseHandling( ){
    SDL_GetMouseState( &mouseX, &mouseY ); mouseY=HEIGHT-mouseY;
    mouse_t   = (mouseX-20 )/timeScale + tstart;
    mouse_val = (mouseY-200)/valScale;
    sprintf( curCaption, "%f %f\0", mouse_t, mouse_val );
    int ipoint_ = binSearchFrom<double>(mouse_t,splines.n,splines.ts);
    if( (splines.ts[ipoint_+1]-mouse_t)<(mouse_t-splines.ts[ipoint_]) ) ipoint_++;
    if(ipoint_!=ipoint){
        ipoint=ipoint_;
        char buff[100];
        Vec3d r,v;
        r.set( splines.CPs[0][ipoint],splines.CPs[1][ipoint],splines.CPs[2][ipoint] );
        v.set( splines.getPointDeriv(ipoint,0), splines.getPointDeriv(ipoint,1), splines.getPointDeriv(ipoint,2) );
        sprintf(buff, "%i %f r(%3.3f,%3.3f,%3.3f) v(%3.3f,%3.3f,%3.3f)", ipoint, splines.ts[ipoint], r.x,r.y,r.z,   v.x,v.y,v.z );
        txtStatic.inputText = buff;
    }
    //printf("%i %i %3.3f %3.3f \n", mouseX,mouseY, mouse_t, mouse_val );
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

                case SDLK_x:  iedit=0; break;
                case SDLK_y:  iedit=1; break;
                case SDLK_z:  iedit=2; break;
            }
            break;

        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                case SDL_BUTTON_RIGHT:
                    //float t   = (mouseX+20 )/timeScale + tstart;
                    //float val = (mouseY-200)/valScale;
                    //ipoint = binSearchFrom<double>(mouse_t,splines.n,splines.ts);
                    dragging=true;
                    break;
            }
            break;

        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    if(dragging){
                        splines.CPs[iedit][ipoint]=mouse_val;
                        splines.ts        [ipoint]=mouse_t;
                        nsamp = splines.evalUniform( tstart, tend, dt, 3, which, (double*)positions, (double*)velocities, (double*)accelerations );
                    }
                    dragging=false;
                    break;
                case SDL_BUTTON_RIGHT:
                    if(dragging){
                        double * dCPi = splines.dCPs[iedit];
                        if(dCPi==0) break;
                        dCPi       [ipoint]=(mouse_val-splines.CPs[iedit][ipoint])/(mouse_t-splines.ts[ipoint]);
                        //splines.ts [idraged]=mouse_t;
                        nsamp = splines.evalUniform( tstart, tend, dt, 3, which, (double*)positions, (double*)velocities, (double*)accelerations );
                    }
                    dragging=false;
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
















