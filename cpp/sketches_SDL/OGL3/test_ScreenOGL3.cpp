#include <stdlib.h>
#include <stdio.h>

#include <vector>
#include <GL/glew.h>
//#define GL_GLEXT_PROTOTYPES
//#include <GL/gl.h>
//#include <SDL2/SDL.h>

#include <fastmath.h>
#include <Vec2.h>
#include <Vec3.h>
#include <Mat3.h>
#include <quaternion.h>
#include <raytrace.h>
//#include <Body.h>

#include "Shader.h"
#include "GLObject.h"
#include "SceneOGL3.h"
#include "ScreenSDL2OGL3.h"
#include "AppSDL2OGL3.h"


#include "Mesh.h"
#include "Solids.h"
#include "GLfunctions.h"
#include "GLobjects.h"
#include "GLObject.h"
#include "Shader.h"

// ========== functions

class TestAppScreenOGL3: public AppSDL2OGL3, public SceneOGL3 { public:
    //virtual void draw(){}

    Shader *sh1;
    GLMesh *glmesh,*gledges;

    //int npoints=0;
    //Vec3f* points=NULL;

    TestAppScreenOGL3():AppSDL2OGL3(),SceneOGL3(){

        // FIXME: second window does not work :((((
        //SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
        //screens.push_back( new ScreenSDL2OGL3( 800, 600) );

        for( ScreenSDL2OGL3* screen: screens ) screen->scenes.push_back( this );

        sh1=new Shader();
        sh1->init( "common_resources/shaders/const3D.glslv",   "common_resources/shaders/const3D.glslf"   );
        //sh1->init( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/color3D.glslf"   );
        //sh1->init( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/cut3DTexture.glslf"   );
        sh1->getDefaultUniformLocation();

        /*
        Mesh mesh;
        mesh.fromFileOBJ( "common_resources/turret.obj" );
        mesh.polygonsToTriangles(false);
        mesh.tris2normals(true);
        mesh.findEdges( );
        */

        /*
        CMesh mesh = Solids::Octahedron;
        nVerts = countVerts( mesh.nfaces, mesh.ngons );
        Vec3f * model_vpos = new Vec3f[nVerts];
        Vec3f * model_vnor = new Vec3f[nVerts];
        hardFace( mesh.nfaces, mesh.ngons, mesh.faces, mesh.verts, (GLfloat*)model_vpos, (GLfloat*)model_vnor );
        glmesh = new GLMesh();
        glmesh->init_d( mesh.points.size(), mesh.triangles.size()*3, ((int*)&mesh.triangles[0]), (double*)&(mesh.points [0]), (double*)&(mesh.normals[0]), NULL, NULL );
        */

        /*
        const CMesh& cmsh = Solids::Icosahedron;
        glmesh->init_d( cmsh.nvert, cmsh.ntri*3, cmsh., (double*)&(mesh.points [0]), (double*)&(mesh.normals[0]), NULL, NULL );
        */
        glmesh = new GLMesh();
        //glmesh.draw_mode = GL_LINES;
        //glmesh->init_wireframe( Solids::Icosahedron );
        glmesh->init_wireframe( Solids::Cube );

        /*
        npoints = 100;
        points  = new Vec3f[npoints];
        float span=50.0;
        for(int i=0; i<npoints; i++){ points[i]={randf(-span,span),randf(-span,span),randf(-span,span)}; }
        */

        Camera& cam = screens[0]->cam;
        //cam.zmin   = -1000.0;
        //cam.zmin = -1000.0; cam.zmax = 1000.0; cam.zoom = 20.00f;
        cam.zmin = 10.0; cam.zmax = 1000.0; cam.zoom = 20.00f;
        cam.aspect = screens[0]->HEIGHT/(float)screens[0]->WIDTH;
        //cam.aspect = (float)screens[0]->WIDTH/(float)screens[0]->HEIGHT;
    }

    virtual void draw( Camera& cam ){

        glClearColor(1.0, 1.0, 1.0, 1.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT  );

        sh1->use();

        Mat3f mrot; mrot.setOne();
        sh1->set_modelMat( (GLfloat*)&mrot );
        sh1->set_modelPos( (const GLfloat[]){0.0f,0.0f,0.0f} );
        setCamera( *sh1, cam );

        /*
        setCameraOrtho   (*sh1, cam);
        sh1->set_modelMat( (GLfloat*)&cam.rot );
        sh1->set_modelPos( (const GLfloat[]){0.0f,0.0f,0.0f} );
        //glmesh->draw();
        */

        float span = 100.0;
        srand(15454);
        GLuint ucolor = sh1->getUloc("baseColor");
        glmesh->preDraw ();
        for( int i=0; i<100; i++ ){
            //sh1->set_modelPos( (GLfloat*)(points+i) );
            Vec3f pos = (Vec3f){ randf(-span,span),randf(-span,span),randf(-span,span) };
            sh1->set_modelPos( (GLfloat*)&pos );
            glUniform4f( ucolor, randf(0,1), randf(0,1), randf(0,1), 1.0 );
            //glUniform4fv( sh1->getUloc("baseColor"), 1,  (const float[]){1.0, 0.0, 0.0, 1.0} );
            //glmesh->draw();
            glmesh->drawRaw();
        }

    };



};


// ================== main

TestAppScreenOGL3 * app;

int main(int argc, char *argv[]){
    app = new TestAppScreenOGL3( );
    app->loop( 1000000 );
    app->quit();
    return 0;
}

