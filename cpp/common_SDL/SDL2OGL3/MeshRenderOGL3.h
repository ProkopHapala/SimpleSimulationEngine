#ifndef MeshRenderOGL3_h
#define MeshRenderOGL3_h

#include "Vec3.h"
#include "Mat3.h"
#include "Mat4.h"

#include "Camera.h"

#include "GLobjects.h"
#include "Shader.h"

class MeshRenderOGL3{
public:
    GLMesh* mesh_tri   = nullptr;
    GLMesh* mesh_lines = nullptr;

    Shader sh_solid;
    Shader sh_const;

    Vec3f modelPos;
    Mat3f modelMat;

    void initDefaultShaders();

    void setModelPos(const Vec3f& pos);
    void setModelMat(const Mat3f& mat);

    void uploadMesh_d( int nVerts, int nTris,  const int* tris,  const double* verts, const double* nors );
    void uploadLines_d( int nVerts, int nEdges, const int* edges, const double* verts );

    void draw(const Camera& cam);
};

#endif
