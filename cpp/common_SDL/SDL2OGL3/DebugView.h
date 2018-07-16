
#ifndef  DebugView_h
#define  DebugView_h

#include <GL/glew.h>
#include "Vec3.h"
#include "DrawOGL3.h"

extern GLMeshBuilder* DEBUG_mesh;
extern Shader*        DEBUG_shader;

//#define DEBUG_VIEW_INIT() GLMeshBuilder* DEBUG_mesh = 0; Shader* DEBUG_shader = 0;
#define DEBUG_VIEW_DEFINE()  GLMeshBuilder* DEBUG_mesh; Shader* DEBUG_shader;
//#define DEBUG_VIEW_INIT() DEBUG_mesh = new GLMeshBuilder(); DEBUG_shader = new Shader( DEFAULT_vertex_shader_code, DEFAULT_fragment_shader_code, false);
#define DEBUG_VIEW_INIT() DEBUG_mesh = new GLMeshBuilder(); DEBUG_shader = new Shader( "common_resources/shaders/color3D.glslv",   "common_resources/shaders/color3D.glslf", true);

void DEBUG_drawCamera( const Camera& cam,  Vec3f c ){
    GLMesh* msh = new GLMesh();
    msh->draw_mode = GL_LINES;

    printf( "zoom %g aspect %g zmin %g zmax %g\n", cam.zoom, cam.aspect, cam.zmin, cam.zmax  );

    float tgx  = cam.getTgX();
    float tgy  = cam.getTgY();
    float zmin = cam.zmin;
    float zmax = cam.zmax;
    float xzmin = tgx*zmin;
    float xzmax = tgx*zmax;
    float yzmin = tgy*zmin;
    float yzmax = tgy*zmax;
    Vec3f verts[] = {
        {-xzmin,-yzmin,zmin},
        {-xzmax,-yzmax,zmax},
        {-xzmin,+yzmin,zmin},
        {-xzmax,+yzmax,zmax},
        {+xzmin,-yzmin,zmin},
        {+xzmax,-yzmax,zmax},
        {+xzmin,+yzmin,zmin},
        {+xzmax,+yzmax,zmax}
    };
    int edges[] = { 0,1, 0,2, 0,4,  1,3,1,5, 2,3, 2,6,    7,5, 7,6, 7,3,  5,4, 6,4    };
    //msh->init( 8, 12*2, edges, verts, NULL, NULL, NULL );
    //return msh;
    GLMeshBuilder* mesh = DEBUG_mesh;
    int i0   = mesh->vpos.size();
    mesh->addLines( 12, edges, verts, c );
    int i1   = DEBUG_mesh->vpos.size();
    mesh->applyMatrixT( {i0,i1}, cam.rot );
    mesh->move( {i0,i1}, cam.pos );
};

void DEBUG_draw( const Camera& cam, const Vec3d& pos, const Mat3d& rotMat ){
    GLMesh * msh = DEBUG_mesh->makeLineMesh();
    DEBUG_shader->use();
    setCamera(*DEBUG_shader, cam);
    DEBUG_shader->setModelPoseT( pos, rotMat );
    msh->draw(GL_LINES);
    //msh.drawRaw(GL_LINES);
    delete msh;
    DEBUG_mesh->clear();
};

//#define DEBUG_line( , )  ;

#endif
