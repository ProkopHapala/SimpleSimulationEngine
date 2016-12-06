
#include <stdlib.h>
#include <stdio.h>
#include <vector>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"

#include "Draw.h"
#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"

#include <thread>  

double  glob_var               = 0.0;
bool    modifications_finished = true;
int     default_sphere         = 0;

class TestAppMesh : public AppSDL2OGL_3D {
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

	TestAppMesh( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppMesh::TestAppMesh( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {};

void TestAppMesh::draw(  ) {
    if( modifications_finished ){
        glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	    //printf("displayLists.size() %i \n", displayLists.size() );
	    for(int i=0; i<displayLists.size(); i++){
            //printf( " %i \n", i, displayLists[i] );
            int ilist = displayLists[i];
            if(ilist>0){  printf(" %i %i \n", i, ilist ); glCallList( ilist ); }
	    }
	
	    //printf( "glob var = %f\n", glob_var );
	}
};


TestAppMesh * thisApp;

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
        thisApp = new TestAppMesh( junk , 800, 600 );
    }

    void loop( int n ){
        thisApp->loop( n );
    }
    
    void asyncLoop( int n){
        std::thread ( loop, n ).detach();
    }
    
    void erase( int iobj ){
        modifications_finished = false;
        printf( "%i %i \n", iobj, thisApp->displayLists[iobj]  );
        glDeleteLists( thisApp->displayLists[iobj], 1 );
        thisApp->displayLists[iobj] = 0;
        modifications_finished = true;
    }

    int defaultSphere( int n ){
        modifications_finished = false;
        if(default_sphere) glDeleteLists(default_sphere,1);
        default_sphere = glGenLists(1);
        glNewList( default_sphere, GL_COMPILE);
            Draw3D::drawSphere_oct( n, 1.0, {0.0,0.0,0.0} );
        glEndList();
        modifications_finished = true;
    }

    int spheres( int n, double * poss_, double * colors_, double * radius ){
        modifications_finished = false;
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
        modifications_finished = true;
        return thisApp->displayLists.size()-1;
    }

    int polyline( int n, double * points_, int closed, uint32_t icolor ){
        modifications_finished = false;
        int ilist = glGenLists(1);
        glNewList(ilist, GL_COMPILE);
            glDisable (GL_LIGHTING);
            Draw::setRGBA(icolor);
            Draw3D::drawPolyLine( n, (Vec3d*)points_, closed );
        glEndList();
        thisApp->displayLists.push_back( ilist );
        modifications_finished = true;
        return thisApp->displayLists.size()-1;
    }

    int lines( int nedges, int * edges, double * points_, uint32_t icolor ){
        modifications_finished = false;
        int ilist = glGenLists(1);
        glNewList(ilist, GL_COMPILE);
            glDisable (GL_LIGHTING);
            Draw::setRGBA(icolor);
            Draw3D::drawLines( nedges, edges, (Vec3d *)points_ );
        glEndList();
        thisApp->displayLists.push_back( ilist );
        modifications_finished = true;
        return thisApp->displayLists.size()-1;
    }

    int triangles( int ntris, int * tris, double * points_, uint32_t icolor ){
        modifications_finished = false;
        int ilist = glGenLists(1);
        glNewList(ilist, GL_COMPILE);
            glEnable (GL_LIGHTING);
            Draw::setRGBA(icolor);
            Draw3D::drawTriangles( ntris, tris, (Vec3d *)points_ );
        glEndList();
        thisApp->displayLists.push_back( ilist );
        modifications_finished = true;
        return thisApp->displayLists.size()-1;
    }

}
