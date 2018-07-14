
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
