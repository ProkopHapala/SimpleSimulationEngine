
// see
// https://www.opengl.org/discussion_boards/showthread.php/171319-glFlush-or-glFinish-with-mulithreading


#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <thread>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"

#include "Draw.h"
#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"



double  glob_var               = 0.0;
int     default_sphere         = 0;

class Vis3DApp : public AppSDL2OGL_3D {
	public:

    Mesh mesh;

    bool  dragging;
    Vec2f mouse0;
    int   ipicked;

    std::vector<int>  displayLists;

	virtual void draw   ();
	//virtual void drawHUD();
	//virtual void mouseHandling( );
	//virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	Vis3DApp( int& id, int WIDTH_, int HEIGHT_ );

};

Vis3DApp::Vis3DApp( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {};

void Vis3DApp::draw(  ) {
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    //printf("displayLists.size() %i \n", displayLists.size() );
    for(int i=0; i<displayLists.size(); i++){
        //printf( " %i \n", i, displayLists[i] );
        int ilist = displayLists[i];
        if(ilist>0){
            //printf(" %i %i \n", i, ilist );
            glCallList( ilist );
        }
    }
    //printf( "glob var = %f\n", glob_var );
};


Vis3DApp * thisApp;

extern "C"{

    void setGobVar( double f ){
        glob_var = f;
    }

    void printHello(){
        printf("Hello!\n");
    }

    void initWindow(){
        SDL_Init(SDL_INIT_VIDEO);
        SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
        int junk;
        thisApp = new Vis3DApp( junk , 800, 600 );
    }

    void loop( int n ){
        thisApp->loop( n );
    }

    void asyncLoop( int n){
        std::thread ( loop, n ).detach();
    }

    int erase( int iobj ){
        if( !thisApp->wait_LOCK(10,5) ) return -1;  thisApp->GL_LOCK = true;
        printf( "erase iobj %i ilist %i \n", iobj, thisApp->displayLists[iobj]  );
        glDeleteLists( thisApp->displayLists[iobj], 1 );
        thisApp->displayLists[iobj] = 0;
        thisApp->GL_LOCK = false;
        return 0;
    }

    int defaultSphere( int n ){
        //thisApp->GL_LOCK = true; glFinish();
        //if( !thisApp->wait_LOCK(10,5) ) return -1;  thisApp->GL_LOCK = true;
        if(default_sphere) glDeleteLists(default_sphere,1);
        default_sphere = glGenLists(1);
        glNewList( default_sphere, GL_COMPILE);
            Draw3D::drawSphere_oct( n, 1.0, {0.0,0.0,0.0} );
        glEndList();
        //thisApp->GL_LOCK = false;
        //glFinish();  thisApp->GL_LOCK = false;
        return 0;
    }

    int spheres( int n, double * poss_, double * colors_, double * radius ){
        if( !thisApp->wait_LOCK(10,5) ) return -1; thisApp->GL_LOCK = true;
        Vec3d * poss   = (Vec3d*)poss_;
        Vec3d * colors = (Vec3d*)colors_;
        if( !default_sphere ) defaultSphere( 4 );
        int ilist = glGenLists(1);
        glNewList(ilist, GL_COMPILE);
        glEnable (GL_LIGHTING);
        for( int i=0; i<n; i++ ){
            Vec3f clr; convert(colors[i],clr);
            Vec3f pos; convert(poss[i],pos);
            glColor3f   (clr.x,clr.y,clr.z);
            glPushMatrix();
            glTranslatef(pos.x,pos.y,pos.z);
            float r = radius[i];
            glScalef(r,r,r);
            glCallList(default_sphere);
            glPopMatrix();
            //printf( " %i (%3.3f,%3.3f,%3.3f) \n", i, pos.x, pos.y, pos.z );
        }
        glEndList();
        thisApp->displayLists.push_back( ilist );
        thisApp->GL_LOCK = false;
        return thisApp->displayLists.size()-1;
    }

    int polyline( int n, double * points_, int closed, uint32_t icolor ){
        if( !thisApp->wait_LOCK(10,5) ) return -1; thisApp->GL_LOCK = true;
        int ilist = glGenLists(1);
        glNewList(ilist, GL_COMPILE);
            glDisable (GL_LIGHTING);
            Draw::setRGBA(icolor);
            Draw3D::drawPolyLine( n, (Vec3d*)points_, closed );
        glEndList();
        thisApp->displayLists.push_back( ilist );
        thisApp->GL_LOCK = false;
        return thisApp->displayLists.size()-1;
    }

    int lines( int nedges, int * edges, double * points_, uint32_t icolor ){
        if( !thisApp->wait_LOCK(10,5) ) return -1; thisApp->GL_LOCK = true;
        int ilist = glGenLists(1);
        glNewList(ilist, GL_COMPILE);
            glDisable (GL_LIGHTING);
            Draw::setRGBA(icolor);
            Draw3D::drawLines( nedges, edges, (Vec3d *)points_ );
        glEndList();
        thisApp->displayLists.push_back( ilist );
        thisApp->GL_LOCK = false;
        return thisApp->displayLists.size()-1;
    }
    
    int vectors( int n, double * vecs, double * poss, uint32_t icolor ){
        if( !thisApp->wait_LOCK(10,5) ) return -1; thisApp->GL_LOCK = true;
        int ilist = glGenLists(1);
        glNewList(ilist, GL_COMPILE);
            glDisable (GL_LIGHTING);
            Draw::setRGBA(icolor);
            for(int i=0; i<n; i++){
                int i3 = i*3;
                glBegin(GL_LINES);
                glVertex3f((float) poss[i3],           (float) poss[i3+1],             (float) poss[i3+2]             );
                glVertex3f((float)(poss[i3]+vecs[i3]),(float)(poss[i3+1]+vecs[i3+1]),(float)(poss[i3+2]+vecs[i3+2]));
                glEnd();
           };
        glEndList();
        thisApp->displayLists.push_back( ilist );
        thisApp->GL_LOCK = false;
        return thisApp->displayLists.size()-1;
    }

    int triangles( int ntris, int * tris, double * points_, uint32_t icolor ){
        if( !thisApp->wait_LOCK(10,5) ) return -1; thisApp->GL_LOCK = true;
        int ilist = glGenLists(1);
        glNewList(ilist, GL_COMPILE);
            glEnable (GL_LIGHTING);
            Draw::setRGBA(icolor);
            Draw3D::drawTriangles( ntris, tris, (Vec3d *)points_ );
        glEndList();
        thisApp->displayLists.push_back( ilist );
        thisApp->GL_LOCK = false;
        return thisApp->displayLists.size()-1;
    }

}
